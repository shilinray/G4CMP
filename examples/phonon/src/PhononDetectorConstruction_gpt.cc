/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc \brief
/// Implementation of the PhononDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.
// 20221006  Remove unused features; add phonon sensor pad with use of
//		G4CMPPhononElectrode to demonstrate KaplanQP.

#include "PhononDetectorConstruction.hh"
#include "PhononSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::PhononDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), topSurfProp(0), botSurfProp(0), wallSurfProp(0),
    electrodeSensitivity(0), fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::~PhononDetectorConstruction() {
  delete topSurfProp;
  delete botSurfProp;
  delete wallSurfProp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* PhononDetectorConstruction::Construct()
{
  if (fConstructed) {
    if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    G4LatticeManager::GetLatticeManager()->Reset();
    // Clear all LogicalSurfaces
    // NOTE: No need to redefine the G4CMPSurfaceProperties
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm); // half (-16,16)
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0); // physical placement
  
  //                               
  // Germanium cylinder - this is the volume in which we will propagate phonons
  //  
  const G4double geHalfZ = 0.05*mm;
  G4VSolid* fGermaniumSolid = new G4Box("fGermaniumSolid", 0.5*cm, 0.5*cm,
                                         geHalfZ);
  G4LogicalVolume* fGermaniumLogical =
    new G4LogicalVolume(fGermaniumSolid,fGermanium,"fGermaniumLogical");
  G4VPhysicalVolume* GePhys = 
    new G4PVPlacement(0,G4ThreeVector(),fGermaniumLogical,"fGermaniumPhysical",
                      worldLogical,false,0); // placing physical volume at center of world logical

  //
  //Germanium lattice information
  //

  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* GeLogical = LM->LoadLattice(fGermanium, "Ge");

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* GePhysical = new G4LatticePhysical(GeLogical);
  GePhysical->SetMillerOrientation(1,0,0); // how crystal is oriented, also 4 coord orient. (online calculator)
  LM->RegisterLattice(GePhys, GePhysical); // connects physical lattice to volume

  // NOTE:  Above registration can also be done in single step:
  // G4LatticlePhysical* GePhysical = LM->LoadLattice(GePhys, "Ge");

  //
  // Aluminum. This is where phonon hits are registered

  // Aluminum meander: multiple long horizontal traces, plus left short pad comb,
  // a top bus, and a bottom semicircular coupling arc; all lie in XY plane on top of Ge.
  const G4double alHalfZ      = 0.05*mm;          // same thickness as before
  const G4double zPos         = geHalfZ + alHalfZ;
  // Meander long horizontal traces
  const G4int    nLongTraces  = 24;
  const G4double traceWidth   = 0.05*mm;          // full width = 0.10 mm
  const G4double traceGap     = 0.05*mm;          // full gap between traces
  const G4double traceHalfY   = traceWidth/2.;
  const G4double longTraceHalfX = 3.0*mm;         // full length 6 mm
  // Left comb short pads
  const G4int    nShortPads   = 12;
  const G4double shortPadHalfX = 0.6*mm;          // full length 1.2 mm
  const G4double padPitch     = 0.30*mm;          // vertical spacing (center-to-center)
  // Top bus
  const G4double busHalfX     = 3.2*mm;
  const G4double busHalfY     = 0.06*mm;
  // Bottom coupling semicircle
  const G4double arcInnerR    = 0.40*mm;
  const G4double arcWidth     = 0.10*mm;
  const G4double arcOuterR    = arcInnerR + arcWidth;
  const G4double arcHalfZ     = alHalfZ;
  const G4double arcPhiStart  = 180.*deg;
  const G4double arcPhiSpan   = 180.*deg;

  // Center offsets
  const G4double patternCenterX = 0.*mm;
  const G4double patternCenterY = 0.*mm;

  // Solids/logicals reused
  G4VSolid* longTraceSolid =
    new G4Box("longTraceSolid", longTraceHalfX, traceHalfY, alHalfZ);
  G4LogicalVolume* longTraceLogical =
    new G4LogicalVolume(longTraceSolid, fAluminum, "longTraceLogical");

  G4VSolid* shortPadSolid =
    new G4Box("shortPadSolid", shortPadHalfX, traceHalfY, alHalfZ);
  G4LogicalVolume* shortPadLogical =
    new G4LogicalVolume(shortPadSolid, fAluminum, "shortPadLogical");

  G4VSolid* busSolid =
    new G4Box("busSolid", busHalfX, busHalfY, alHalfZ);
  G4LogicalVolume* busLogical =
    new G4LogicalVolume(busSolid, fAluminum, "busLogical");

  G4VSolid* arcSolid =
    new G4Tubs("couplingArcSolid", arcInnerR, arcOuterR, arcHalfZ,
               arcPhiStart, arcPhiSpan);
  G4LogicalVolume* arcLogical =
    new G4LogicalVolume(arcSolid, fAluminum, "couplingArcLogical");

  std::vector<G4VPhysicalVolume*> alPhysParts;

  // Place long traces (stacked along +Y)
  const G4double totalHeight =
      nLongTraces*(2.*traceHalfY) + (nLongTraces-1)*traceGap;
  const G4double firstTraceY =
      patternCenterY - 0.5*totalHeight + traceHalfY;
  for (G4int i = 0; i < nLongTraces; ++i) {
    G4double y = firstTraceY + i*(2.*traceHalfY + traceGap);
    G4String name = "longTracePhys_" + std::to_string(i);
    alPhysParts.push_back(new G4PVPlacement(
        0, G4ThreeVector(patternCenterX, y, zPos),
        longTraceLogical, name, worldLogical, false, i));
  }

  // Place short pads (left comb), aligned roughly to upper section
  const G4double combStartY = firstTraceY + 2.*mm; // shift start
  for (G4int i = 0; i < nShortPads; ++i) {
    G4double y = combStartY + i*padPitch;
    G4String name = "shortPadPhys_" + std::to_string(i);
    alPhysParts.push_back(new G4PVPlacement(
        0, G4ThreeVector(patternCenterX - longTraceHalfX + shortPadHalfX - 0.4*mm,
                         y, zPos),
        shortPadLogical, name, worldLogical, false, i));
  }

  // Place top bus
  G4VPhysicalVolume* busPhys = new G4PVPlacement(
      0, G4ThreeVector(patternCenterX,
                       firstTraceY + nLongTraces*(2.*traceHalfY + traceGap)/2. + busHalfY + 0.05*mm,
                       zPos),
      busLogical, "busPhys", worldLogical, false, 0);
  alPhysParts.push_back(busPhys);

  // Place bottom coupling semicircle
  G4VPhysicalVolume* arcPhys = new G4PVPlacement(
      0, G4ThreeVector(patternCenterX,
                       firstTraceY - arcOuterR - 0.2*mm,
                       zPos),
      arcLogical, "couplingArcPhys", worldLogical, false, 0);
  alPhysParts.push_back(arcPhys);

  //
  // detector -- Note : "sensitive detector" is attached to Germanium crystal
  // want a phonon sensitive detector, attached to Ge crystal
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!electrodeSensitivity)
    electrodeSensitivity = new PhononSensitivity("PhononElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fGermaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  // Assign SD to all aluminum logicals
  longTraceLogical->SetSensitiveDetector(electrodeSensitivity);
  shortPadLogical->SetSensitiveDetector(electrodeSensitivity);
  busLogical->SetSensitiveDetector(electrodeSensitivity);
  arcLogical->SetSensitiveDetector(electrodeSensitivity);

  //
  // surface between Al and Ge determines phonon reflection/absorption
  //
  // 
  if (!fConstructed) {
    const G4double GHz = 1e9 * hertz; 

    //the following coefficients and cutoff values are not well-motivated
    //the code below is used only to demonstrate how to set these values.
    const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
    const std::vector<G4double> diffCoeffs =
      {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
    const std::vector<G4double> specCoeffs =
      {0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};

    const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

    topSurfProp = new G4CMPSurfaceProperty("TopAlSurf", 1.0, 0.0, 0.0, 0.0,  
					  	        0.3, 1.0, 0.0, 0.0);   // absorption and reflection are the first two, opposite for the wall
    topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);
    AttachPhononSensor(topSurfProp);

    botSurfProp = 0;

    wallSurfProp = new G4CMPSurfaceProperty("WallSurf", 0.0, 1.0, 0.0, 0.0,
					    	          0.0, 1.0, 0.0, 0.0);
    wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);

  }

  // Connects the inner volume, outer volume, and physics that applies at the surface
  // Logical border surface applies the specified physics for ANYWHERE the two volumes touch
  //
  for (size_t i = 0; i < alPhysParts.size(); ++i) {
    new G4CMPLogicalBorderSurface("AlCurveSurface_" + std::to_string(i),
                                  GePhys, alPhysParts[i], topSurfProp);
  }

  // Remove old: new G4CMPLogicalBorderSurface("feedlineTop", ...), ("sensorTop", ...)

  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  fGermaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  // fAluminumFeedlineLogical->SetVisAttributes(simpleBoxVisAtt);
  // fAluminumSensorLogical->SetVisAttributes(simpleBoxVisAtt);

  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0.7,0.7,0.2));
  alVis->SetForceSolid(true);
  longTraceLogical->SetVisAttributes(alVis);
  shortPadLogical->SetVisAttributes(alVis);
  busLogical->SetVisAttributes(alVis);
  arcLogical->SetVisAttributes(alVis);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void PhononDetectorConstruction::
AttachPhononSensor(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  // Specify properties of aluminum sensor, same on both detector faces
  // See G4CMPPhononElectrode.hh or README.md for property keys

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", 0.20);    // True sensor area
  sensorProp->AddConstProperty("filmThickness", 600.*nm);
  sensorProp->AddConstProperty("gapEnergy", 173.715e-6*eV);
  sensorProp->AddConstProperty("lowQPLimit", 3.);
  sensorProp->AddConstProperty("phononLifetime", 242.*ps);
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
  sensorProp->AddConstProperty("vSound", 3.26*km/s);
  sensorProp->AddConstProperty("subgapAbsorption", 0.1);

  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

