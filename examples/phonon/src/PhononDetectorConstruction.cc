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
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
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
  G4VSolid* worldSolid = new G4Box("World",8.*cm,8.*cm,8.*cm); // half (-16,16)
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

  // Aluminum feedline
  const G4double alFeedlineHalfX = 0.5*cm;
  const G4double alFeedlineHalfY = 1.5*um;
  const G4double alFeedlineHalfZ = 0.1*um;

  // Outer feedline dimensions (100um width)
  const G4double alOuterFeedlineHalfY = 50.0*um;
  const G4double feedlineGap = 2.0*um;

  G4Box* solidCenter = new G4Box("feedlineCenter", alFeedlineHalfX, alFeedlineHalfY, alFeedlineHalfZ);
  G4Box* solidOuter  = new G4Box("feedlineOuter",  alFeedlineHalfX, alOuterFeedlineHalfY, alFeedlineHalfZ);

  G4MultiUnion* fAluminumFeedlineSolid = new G4MultiUnion("feedlineMultiUnion");

  G4RotationMatrix rot;
  G4ThreeVector posCenter(0, 0, 0);
  G4Transform3D trCenter(rot, posCenter);
  fAluminumFeedlineSolid->AddNode(*solidCenter, trCenter);

  G4double yOffset = alFeedlineHalfY + feedlineGap + alOuterFeedlineHalfY;

  G4ThreeVector posTop(0, yOffset, 0);
  G4Transform3D trTop(rot, posTop);
  fAluminumFeedlineSolid->AddNode(*solidOuter, trTop);

  G4ThreeVector posBot(0, -yOffset, 0);
  G4Transform3D trBot(rot, posBot);
  fAluminumFeedlineSolid->AddNode(*solidOuter, trBot);

  fAluminumFeedlineSolid->Voxelize();

  G4LogicalVolume* fAluminumFeedlineLogical =
    new G4LogicalVolume(fAluminumFeedlineSolid,fAluminum,"feedlineLogical");

  G4VPhysicalVolume* aluminumFeedlinePhysical = new G4PVPlacement(
    0, G4ThreeVector(0.,0., geHalfZ + alFeedlineHalfZ), fAluminumFeedlineLogical, "feedlinePhysical",
    worldLogical, false, 0);

  // Aluminum sensor
  const G4double alSensorHalfXY = 0.5*mm;     // full side = 1.0 mm
  const G4double alSensorHalfZ  = 0.1*um;  
  const G4double gaptosensor = .1*mm;

  G4VSolid* fAluminumSensorSolid = new G4Box("sensor", alSensorHalfXY, alSensorHalfXY, alSensorHalfZ);
  G4LogicalVolume* fAluminumSensorLogical =
    new G4LogicalVolume(fAluminumSensorSolid,fAluminum,"sensorLogical");
  
  G4double ycentertosensor = yOffset + alOuterFeedlineHalfY + gaptosensor + alSensorHalfXY;

  G4VPhysicalVolume* aluminumSensorPhysical = new G4PVPlacement(
    0, G4ThreeVector(0., ycentertosensor, geHalfZ + alSensorHalfZ), fAluminumSensorLogical, "sensorPhysical",
    worldLogical, false, 0);

  // QPD
  
  //
  // detector -- Note : "sensitive detector" is attached to Germanium crystal
  // want a phonon sensitive detector, attached to Ge crystal
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!electrodeSensitivity)
    electrodeSensitivity = new PhononSensitivity("PhononElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fGermaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  //
  // surface between Al and Ge determines phonon reflection/absorption
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
  new G4CMPLogicalBorderSurface("feedlineTop", GePhys, aluminumFeedlinePhysical,
				topSurfProp);
  new G4CMPLogicalBorderSurface("sensorTop", GePhys, aluminumSensorPhysical,
				topSurfProp);
  new G4CMPLogicalBorderSurface("detWall", GePhys, fWorldPhys,
				wallSurfProp);

  //                                        
  // Visualization attributes
  //
  // World remains invisible
  G4VisAttributes* wrldVis = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  wrldVis->SetForceSolid(true);
  worldLogical->SetVisAttributes(wrldVis);

  // Germanium crystal: light gray, solid
  G4VisAttributes* geVis = new G4VisAttributes(G4Colour(0.85,0.85,0.85));
  geVis->SetForceSolid(true);
  fGermaniumLogical->SetVisAttributes(geVis);

  // Aluminum patterned parts
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  alVis->SetForceSolid(true);
  fAluminumSensorLogical->SetVisAttributes(alVis);
  fAluminumFeedlineLogical->SetVisAttributes(alVis);
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

