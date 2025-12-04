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
  const G4double thickness = 0.1*um;

  // Outer feedline dimensions (100um width)
  const G4double alUPGhalfy = 50.0*um;
  const G4double alLPGhalfy = 2.0*um;
  const G4double feedlineGap = 2.0*um;

  G4Box* feedline = new G4Box("feedlineCenter", alFeedlineHalfX, alFeedlineHalfY, thickness);
  G4Box* uppergroundplane  = new G4Box("uppergroundplane",  alFeedlineHalfX, alUPGhalfy, thickness);
  G4Box* lowergroundplane  = new G4Box("lowergroundplane",  alFeedlineHalfX, alLPGhalfy, thickness);
  
  G4LogicalVolume* alFLlogical = new G4LogicalVolume(feedline,fAluminum,"alFLlogical"); // logical feedline
  G4LogicalVolume* alUGPlogical = new G4LogicalVolume(uppergroundplane,fAluminum,"alUGPlogical"); // logical feedline
  G4LogicalVolume* alLGPlogical = new G4LogicalVolume(lowergroundplane,fAluminum,"alLGPlogical"); // logical feedline

  G4VPhysicalVolume* alFLphysical = new G4PVPlacement(
    0, G4ThreeVector(0.,0., geHalfZ + thickness), alFLlogical, "alFLphysical",
    worldLogical, false, 0); // physical feedline
  
  G4double yUGPoffset = alFeedlineHalfY + feedlineGap + alUPGhalfy;
  G4double yLGPoffset = alFeedlineHalfY + feedlineGap + alLPGhalfy;

  G4VPhysicalVolume* alUGPphysical = new G4PVPlacement(
    0, G4ThreeVector(0., yUGPoffset, geHalfZ + thickness), alUGPlogical, "alUGPphysical",
    worldLogical, false, 0); // physical UPG
  G4VPhysicalVolume* alLGPphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -yLGPoffset, geHalfZ + thickness), alLGPlogical, "alLGPphysical",
    worldLogical, false, 0); // physical LGP


  // QPD

  // coupling capacitor
  const G4double couplCapy = 3.0*um;
  const G4double couplCapx = 300*um;
  const G4double couplCapgap = 2.0*um;

  G4Box* couplCap = new G4Box("couplCap", couplCapx, couplCapy, thickness);
  G4LogicalVolume* couplCaplogical = new G4LogicalVolume(couplCap,fAluminum,"couplCaplogical");

  G4double ycouplcapoffset = yLGPoffset + alLPGhalfy + couplCapgap+ couplCapy;

  G4VPhysicalVolume* couplCapphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -ycouplcapoffset, geHalfZ + thickness), couplCaplogical, "couplCapphysical",
    worldLogical, false, 0);

  // connector
  const G4double connectory = 5.0*um;
  const G4double connectorx = 3.0*um;

  G4Box* connector = new G4Box("connector", connectorx, connectory, thickness);
  G4LogicalVolume* connectorlogical = new G4LogicalVolume(connector,fAluminum,"connectorlogical");

  G4double yconnectoroffset = ycouplcapoffset + couplCapy + connectory;

  G4VPhysicalVolume* connectorphysical = new G4PVPlacement(
    0, G4ThreeVector(0., -yconnectoroffset, geHalfZ + thickness), connectorlogical, "connectorphysical",
    worldLogical, false, 0);

  // cap ind connector
  const G4double capindconnectory = 3.0*um;
  const G4double capindconnectorx = 200.5*um;

  G4Box* capindconnector = new G4Box("capindconnector", capindconnectorx, capindconnectory, thickness);
  G4LogicalVolume* capindconnectorlogical = new G4LogicalVolume(capindconnector,fAluminum,"capindconnectorlogical");

  G4double ycapindconnectoroffset = yconnectoroffset + connectory + capindconnectory;
  G4double xcapindconnectoroffset = -32.5*um;

  G4VPhysicalVolume* capindconnectorphysical = new G4PVPlacement(
    0, G4ThreeVector(xcapindconnectoroffset, -ycapindconnectoroffset, geHalfZ + thickness), capindconnectorlogical, "capindconnectorphysical",
    worldLogical, false, 0); 
    
  // indvert1 
  const G4double indcenter = 0 - 208*um;

  const G4double indhorzx = 23*um;
  const G4double indhorzy = 1.0*um;

  const G4double indvertx = 1.0*um;
  const G4double indverty = 8.0*um;

  G4Box* indvert1 = new G4Box("indvert1", indvertx, indverty, thickness);
  G4LogicalVolume* indvert1logical = new G4LogicalVolume(indvert1, fAluminum, "indvert1logical");
  
  G4double xindvert1offset = indcenter - indhorzx - indvertx;
  G4double yindvert1offset = ycapindconnectoroffset + capindconnectory + indverty;

  G4VPhysicalVolume* indvert1physical = new G4PVPlacement(
    0, G4ThreeVector(xindvert1offset, -yindvert1offset, geHalfZ + thickness), indvert1logical, "indvert1physical",
    worldLogical, false, 0);

  // Generate meandering inductor pattern from indhorz1 to indhorz40
  std::vector<G4VPhysicalVolume*> inductorPhys;
  std::vector<G4LogicalVolume*> inductorLog;

  G4double currentYVertOffset = yindvert1offset;
  G4double yIndHorzLastPos = 0.0;

  // Total Inductor Height (Top of indhorz1 to Bottom of indhorz40):
  // Pitch = 2 * (indverty - indhorzy) = 2 * (8.0 - 1.0) = 14.0 um
  // Distance = 39 * Pitch + 2 * indhorzy
  //          = 39 * 14.0 + 2 * 1.0 = 548.0 um

  for (int i = 1; i <= 40; ++i) {
    // indhorz_i
    G4String hName = "indhorz" + std::to_string(i);
    G4Box* solidH = new G4Box(hName, indhorzx, indhorzy, thickness);
    G4LogicalVolume* logicH = new G4LogicalVolume(solidH, fAluminum, hName + "logical");
    inductorLog.push_back(logicH);

    G4double yHorzOffset = currentYVertOffset + indverty - indhorzy;
    if (i == 40) yIndHorzLastPos = yHorzOffset;
    
    G4VPhysicalVolume* physH = new G4PVPlacement(
      0, G4ThreeVector(indcenter, -yHorzOffset, geHalfZ + thickness), 
      logicH, hName + "physical", worldLogical, false, 0);
    inductorPhys.push_back(physH);

    // indvert_{i+1}
    if (i < 40) {
      int nextVertIdx = i + 1;
      G4String vName = "indvert" + std::to_string(nextVertIdx);
      G4Box* solidV = new G4Box(vName, indvertx, indverty, thickness);
      G4LogicalVolume* logicV = new G4LogicalVolume(solidV, fAluminum, vName + "logical");
      inductorLog.push_back(logicV);

      G4double xVertOffset = (nextVertIdx % 2 == 0) ? 
                             (indcenter + indhorzx + indvertx) : // Even -> Right
                             (indcenter - indhorzx - indvertx);  // Odd -> Left
      
      G4double yVertOffset = yHorzOffset - indhorzy + indverty;
      currentYVertOffset = yVertOffset; 

      G4VPhysicalVolume* physV = new G4PVPlacement(
        0, G4ThreeVector(xVertOffset, -yVertOffset, geHalfZ + thickness), 
        logicV, vName + "physical", worldLogical, false, 0);
      inductorPhys.push_back(physV);
    }
  }

  // leftcapwall
  const G4double leftcapwally = 284*um;
  const G4double leftcapwallx = 1*um;

  G4Box* leftcapwall = new G4Box("leftcapwall", leftcapwallx, leftcapwally, thickness);
  G4LogicalVolume* leftcapwalllogical = new G4LogicalVolume(leftcapwall, fAluminum, "leftcapwalllogical");
  
  G4double xleftcapwalloffset = 0 - 161*um;
  G4double yleftcapwalloffset = ycapindconnectoroffset + 305*um;

  G4VPhysicalVolume* leftcapwallphysical = new G4PVPlacement(
    0, G4ThreeVector(xleftcapwalloffset, -yleftcapwalloffset, geHalfZ + thickness), leftcapwalllogical, "leftcapwallphysical",
    worldLogical, false, 0);

  // rightcapwall
  const G4double rightcapwally = 291*um;
  const G4double rightcapwallx = 1*um;

  G4Box* rightcapwall = new G4Box("rightcapwall", rightcapwallx, rightcapwally, thickness);
  G4LogicalVolume* rightcapwalllogical = new G4LogicalVolume(rightcapwall, fAluminum, "rightcapwalllogical");
  
  G4double xrightcapwalloffset = 0 + 167*um;
  G4double yrightcapwalloffset = ycapindconnectoroffset + 294*um;

  G4VPhysicalVolume* rightcapwallphysical = new G4PVPlacement(
    0, G4ThreeVector(xrightcapwalloffset, -yrightcapwalloffset, geHalfZ + thickness), rightcapwalllogical, "rightcapwallphysical",
    worldLogical, false, 0);

  // Interdigitated Capacitor (IDC) Fingers: caphorz1 to caphorz24
  const G4double xcapcenter = 0 + 0;
  const G4double xcapoffset = 6*um;
  const G4double caphorzy = 1*um;
  const G4double caphorzx = 160*um;
  const G4double captocapgap = 8*um;
  const G4double ycapstart = ycapindconnectoroffset + 22*um;
  // Total IDC height (top of 1 to bottom of 24) = 2*caphorzy + 23*captocapgap = 186.0 um

  std::vector<G4VPhysicalVolume*> idcPhys;
  std::vector<G4LogicalVolume*> idcLog;

  for (int i = 1; i <= 57; ++i) {
    G4String name = "caphorz" + std::to_string(i);
    G4Box* solid = new G4Box(name, caphorzx, caphorzy, thickness);
    G4LogicalVolume* logic = new G4LogicalVolume(solid, fAluminum, name + "logical");
    idcLog.push_back(logic);

    // Odd indices shift left, Even indices shift right
    G4double xpos = (i % 2 != 0) ? (xcapcenter) : (xcapcenter + xcapoffset);
    G4double ypos = ycapstart + (i - 1) * captocapgap;

    G4VPhysicalVolume* phys = new G4PVPlacement(
      0, G4ThreeVector(xpos, -ypos, geHalfZ + thickness), 
      logic, name + "physical", worldLogical, false, 0);
    idcPhys.push_back(phys);
  }

  // botcapindconnector1
  const G4double botcapindconnector1x = 1*um;
  const G4double botcapindconnector1y = 8*um;

  G4Box* botcapindconnector1 = new G4Box("botcapindconnector1", botcapindconnector1x, botcapindconnector1y, thickness);
  G4LogicalVolume* botcapindconnector1logical = new G4LogicalVolume(botcapindconnector1, fAluminum, "botcapindconnector1logical");
  
  G4double xbotcapindconnector1offset = xindvert1offset;
  G4double ybotcapindconnector1offset = yIndHorzLastPos + 7*um;

  G4VPhysicalVolume* botcapindconnector1physical = new G4PVPlacement(
    0, G4ThreeVector(xbotcapindconnector1offset, -ybotcapindconnector1offset, geHalfZ + thickness), 
    botcapindconnector1logical, "botcapindconnector1physical",
    worldLogical, false, 0);

  // botcapindconnector2
  // Connects botcapindconnector1 and leftcapwall
  const G4double botcapindconnector2x = 34.5*um; // 69um / 2
  const G4double botcapindconnector2y = 1.0*um;  // 2um / 2

  G4Box* botcapindconnector2 = new G4Box("botcapindconnector2", botcapindconnector2x, botcapindconnector2y, thickness);
  G4LogicalVolume* botcapindconnector2logical = new G4LogicalVolume(botcapindconnector2, fAluminum, "botcapindconnector2logical");

  G4double xbotcapindconnector2offset = (xleftcapwalloffset + xbotcapindconnector1offset) / 2.0;
  // Align bottom edge with botcapindconnector1 bottom edge to ensure connection
  G4double ybotcapindconnector2offset = (ybotcapindconnector1offset + botcapindconnector1y) - botcapindconnector2y;

  G4VPhysicalVolume* botcapindconnector2physical = new G4PVPlacement(
    0, G4ThreeVector(xbotcapindconnector2offset, -ybotcapindconnector2offset, geHalfZ + thickness),
    botcapindconnector2logical, "botcapindconnector2physical",
    worldLogical, false, 0);

  // capjunctconnect
  const G4double capjunctconnectx = 150*um; 
  const G4double capjunctconnecty = 1.0*um; 

  G4Box* capjunctconnect = new G4Box("capjunctconnect", capjunctconnectx, capjunctconnecty, thickness);
  G4LogicalVolume* capjunctconnectlogical = new G4LogicalVolume(capjunctconnect, fAluminum, "capjunctconnectlogical");

  // Right align with rightcapwall: (xrightcapwalloffset + rightcapwallx) - capjunctconnectx
  // Top align with rightcapwall bottom: (yrightcapwalloffset + rightcapwally) + capjunctconnecty
  G4double xcapjunctconnectoffset = (xrightcapwalloffset + rightcapwallx) - capjunctconnectx;
  G4double ycapjunctconnectoffset = (yrightcapwalloffset + rightcapwally) + capjunctconnecty;

  G4VPhysicalVolume* capjunctconnectphysical = new G4PVPlacement(
    0, G4ThreeVector(xcapjunctconnectoffset, -ycapjunctconnectoffset, geHalfZ + thickness),
    capjunctconnectlogical, "capjunctconnectphysical",
    worldLogical, false, 0);


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

    wallSurfProp = new G4CMPSurfaceProperty("WallSurf", 0.0, 1.0, 0.0, 0.0,
					    	          0.0, 1.0, 0.0, 0.0);
    wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);

  }

  // Connects the inner volume, outer volume, and physics that applies at the surface
  // Logical border surface applies the specified physics for ANYWHERE the two volumes touch
  //
  new G4CMPLogicalBorderSurface("Al", GePhys, alFLphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, alUGPphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, alLGPphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, capindconnectorphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, connectorphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, couplCapphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, indvert1physical, topSurfProp);
  for (auto phys : inductorPhys) {
    new G4CMPLogicalBorderSurface("Al", GePhys, phys, topSurfProp);
  }
  new G4CMPLogicalBorderSurface("detWall", GePhys, fWorldPhys, wallSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, leftcapwallphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, rightcapwallphysical, topSurfProp);
  for (auto phys : idcPhys) {
    new G4CMPLogicalBorderSurface("Al", GePhys, phys, topSurfProp);
  }
  new G4CMPLogicalBorderSurface("Al", GePhys, botcapindconnector1physical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, botcapindconnector2physical, topSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, capjunctconnectphysical, topSurfProp);


  //                                        
  // Visualization attributes
  //
  // World remains invisible
  G4VisAttributes* wrldVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.1, 0.1));
  wrldVis->SetVisibility(false);
  worldLogical->SetVisAttributes(wrldVis);

  // Germanium crystal: light gray, solid
  G4VisAttributes* geVis = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85, 0.4));
  geVis->SetVisibility(true);
  fGermaniumLogical->SetVisAttributes(geVis);
  

  // Aluminum patterned parts
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  alVis->SetVisibility(true);
  alUGPlogical->SetVisAttributes(alVis);
  alLGPlogical->SetVisAttributes(alVis);
  alFLlogical->SetVisAttributes(alVis);
  connectorlogical->SetVisAttributes(alVis);
  couplCaplogical->SetVisAttributes(alVis);
  capindconnectorlogical->SetVisAttributes(alVis);
  indvert1logical->SetVisAttributes(alVis);
  for (auto log : inductorLog) {
    log->SetVisAttributes(alVis);
  }
  leftcapwalllogical->SetVisAttributes(alVis);
  rightcapwalllogical->SetVisAttributes(alVis);
  for (auto log : idcLog) {
    log->SetVisAttributes(alVis);
  }
  botcapindconnector1logical->SetVisAttributes(alVis);
  botcapindconnector2logical->SetVisAttributes(alVis);
  capjunctconnectlogical->SetVisAttributes(alVis);

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
  sensorProp->AddConstProperty("filmAbsorption", 1.0);    // True sensor area
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

