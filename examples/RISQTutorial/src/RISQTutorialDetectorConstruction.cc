// QPD

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
#include "RISQTutorialConfigManager.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "CADMesh.hh"

using namespace RISQTutorialDetectorParameters;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::RISQTutorialDetectorConstruction()
  : fAir(0), fVacuum(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), AlSurfProp(0), polishedwallSurfProp(0), sidewallSurfProp(0), 
    fSuperconductorSensitivity(0), fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {
  delete AlSurfProp;
  delete polishedwallSurfProp;
  delete sidewallSurfProp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* RISQTutorialDetectorConstruction::Construct()
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

void RISQTutorialDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fAir = nistManager->FindOrBuildMaterial("G4_AIR"); 
  fVacuum = new G4Material("VACUUM", 
        1.,
		1.008*CLHEP::g/CLHEP::mole,
		1.0e-25*CLHEP::g/CLHEP::cm3,
		kStateGas,
		0.01*CLHEP::kelvin,
	    3.0e-18*pascal);
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",8.*cm,8.*cm,8.*cm); // half (-16,16)
  G4LogicalVolume* worldLogical =
    new G4LogicalVolume(worldSolid,fVacuum,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0); // physical placement
  
  //                               
  // Germanium cylinder - this is the volume in which we will propagate phonons
  //  
  const G4double geHalfX = 0.5*cm;
  const G4double geHalfY = 0.5*cm;
  const G4double geHalfZ = 0.05*mm;
  G4VSolid* fGermaniumSolid = new G4Box("fGermaniumSolid", geHalfX, geHalfY, geHalfZ);
  G4LogicalVolume* fGermaniumLogical = new G4LogicalVolume(fGermaniumSolid,fGermanium,"fGermaniumLogical");
  G4VPhysicalVolume* GePhys = new G4PVPlacement(0,G4ThreeVector(),fGermaniumLogical,"fGermaniumPhysical", worldLogical,false,0); 
  // placing physical volume at center of world logical

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
  // Air boxes touching the 4 side faces of the Ge crystal (+X, -X, +Y, -Y)
  //
  const G4double airSideThickness = 0.5*cm; // Distance the air extends away from the chip
  
  // X boxes cover the +X and -X faces of the Ge crystal
  G4VSolid* airSideSolidX = new G4Box("airSideSolidX", airSideThickness/2, geHalfY, geHalfZ);
  G4LogicalVolume* airSideLogicalX = new G4LogicalVolume(airSideSolidX, fAir, "airSideLogicalX");
  G4VPhysicalVolume* airSideRightXPhys = new G4PVPlacement(0, G4ThreeVector(geHalfX + airSideThickness/2, 0, 0), airSideLogicalX, "airSideRightX", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftXPhys = new G4PVPlacement(0, G4ThreeVector(-geHalfX - airSideThickness/2, 0, 0), airSideLogicalX, "airSideLeftX", worldLogical, false, 0);

  // Y boxes cover the +Y and -Y faces of the Ge crystal
  G4VSolid* airSideSolidY = new G4Box("airSideSolidY", geHalfX, airSideThickness/2, geHalfZ);
  G4LogicalVolume* airSideLogicalY = new G4LogicalVolume(airSideSolidY, fAir, "airSideLogicalY");
  G4VPhysicalVolume* airSideRightYPhys = new G4PVPlacement(0, G4ThreeVector(0, geHalfY + airSideThickness/2, 0), airSideLogicalY, "airSideRightY", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftYPhys = new G4PVPlacement(0, G4ThreeVector(0, -geHalfY - airSideThickness/2, 0), airSideLogicalY, "airSideLeftY", worldLogical, false, 0);

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

  // Placeholder removed for Side Ground Planes

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

  // junct1
  const G4double junct1x = 5.0*um; // 10um / 2
  const G4double junct1y = 0.05*um; // 0.1um / 2

  const G4double capjunctconnectjunct1gap = 0.25*um;

  G4Box* junct1 = new G4Box("junct1", junct1x, junct1y, thickness);
  G4LogicalVolume* junct1logical = new G4LogicalVolume(junct1, fAluminum, "junct1logical");

  G4double xjunct1offset = 2.0*um;
  // 0.25um below capjunctconnect
  G4double yjunct1offset = ycapjunctconnectoffset + capjunctconnecty + capjunctconnectjunct1gap + junct1y;

  G4VPhysicalVolume* junct1physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct1offset, -yjunct1offset, geHalfZ + thickness),
    junct1logical, "junct1physicalshunt",
    worldLogical, false, 0);

  // junct2
  const G4double junct2x = 0.05*um; // 0.1um / 2
  const G4double junct2y = 0.5*um;  // 1.0um / 2

  G4Box* junct2 = new G4Box("junct2", junct2x, junct2y, thickness);
  G4LogicalVolume* junct2logical = new G4LogicalVolume(junct2, fAluminum, "junct2logical");

  G4double xjunct2offset = xjunct1offset;
  // Connects to center of junct1 and extends downward
  G4double yjunct2offset = yjunct1offset + junct2y;

  G4VPhysicalVolume* junct2physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct2offset, -yjunct2offset, geHalfZ + thickness),
    junct2logical, "junct2physicalshunt",
    worldLogical, false, 0);

  // junct3
  const G4double junct3x = 1.0*um; // 2.0um / 2
  const G4double junct3y = 0.05*um; // 0.1um / 2

  G4Box* junct3 = new G4Box("junct3", junct3x, junct3y, thickness);
  G4LogicalVolume* junct3logical = new G4LogicalVolume(junct3, fAluminum, "junct3logical");

  // Left side at x=0 -> Center at +1.0um
  G4double xjunct3offset = 1.0*um;
  // Intersect junct2 0.2um above bottom of junct2
  // Bottom of junct2 is at -(yjunct2offset + junct2y)
  // Target Y is -(yjunct2offset + junct2y - 0.2*um)
  G4double yjunct3offset = yjunct2offset + junct2y - 0.2*um;

  G4VPhysicalVolume* junct3physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct3offset, -yjunct3offset, geHalfZ + thickness),
    junct3logical, "junct3physicalshunt",
    worldLogical, false, 0);

  // junct4
  const G4double junct4x = 0.05*um; // 0.1um / 2
  const G4double junct4y = 1.0*um;  // 2.0um / 2

  G4Box* junct4 = new G4Box("junct4", junct4x, junct4y, thickness);
  G4LogicalVolume* junct4logical = new G4LogicalVolume(junct4, fAluminum, "junct4logical");

  // Left align with junct3: (xjunct3offset - junct3x) + junct4x
  G4double xjunct4offset = (xjunct3offset - junct3x) + junct4x;
  // Top connects to bottom of junct3: (yjunct3offset + junct3y) + junct4y
  G4double yjunct4offset = yjunct3offset + junct3y + junct4y;

  G4VPhysicalVolume* junct4physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct4offset, -yjunct4offset, geHalfZ + thickness),
    junct4logical, "junct4physicalshunt",
    worldLogical, false, 0);

  // junct5
  const G4double junct5x = 0.5*um; // 1.0um / 2
  const G4double junct5y = 1.0*um; // 2.0um / 2

  G4Box* junct5 = new G4Box("junct5", junct5x, junct5y, thickness);
  G4LogicalVolume* junct5logical = new G4LogicalVolume(junct5, fAluminum, "junct5logical");

  // Centered with junct4
  G4double xjunct5offset = xjunct4offset;
  // Top connects to bottom of junct4
  G4double yjunct5offset = yjunct4offset + junct4y + junct5y;

  G4VPhysicalVolume* junct5physical = new G4PVPlacement(
    0, G4ThreeVector(xjunct5offset, -yjunct5offset, geHalfZ + thickness),
    junct5logical, "junct5physicalshunt",
    worldLogical, false, 0);

  // absorber
  const G4double absorberRadius = 350.0*um;
  // Semicircle extending downwards (180 to 360 degrees)
  G4Tubs* absorber = new G4Tubs("absorber", 0.0, absorberRadius, thickness, 180.0*deg, 180.0*deg);
  G4LogicalVolume* absorberlogical = new G4LogicalVolume(absorber, fAluminum, "absorberlogical");

  G4double xabsorberoffset = xjunct5offset;
  // Flat face connects to bottom of junct5
  G4double yabsorberoffset = yjunct5offset + junct5y;

  G4VPhysicalVolume* absorberphysical = new G4PVPlacement(
    0, G4ThreeVector(xabsorberoffset, -yabsorberoffset, geHalfZ + thickness),
    absorberlogical, "absorberphysicalshunt",
    worldLogical, false, 0);

  // Side Ground Planes
  const G4double sideGroundx = 10.0*um; // 20um width / 2
  
  // They extend from the bottom of the lowergroundplane to the top flat edge of the bridge
  G4double yLGPBottom = yLGPoffset + alLPGhalfy;
  G4double ySideGroundTop = yLGPBottom;
  // The bridge top flat edge is at yabsorberoffset (same y as the absorber flat edge)
  G4double ySideGroundBottom = yabsorberoffset;
  
  G4double sideGroundy = (ySideGroundBottom - ySideGroundTop) / 2.0;
  // Note: Y coordinates are negative in placement
  G4double ySideGroundCenter = ySideGroundTop + sideGroundy;

  G4Box* sideGround = new G4Box("sideGround", sideGroundx, sideGroundy, thickness);
  G4LogicalVolume* sideGroundlogical = new G4LogicalVolume(sideGround, fAluminum, "sideGroundlogical");

  // Centers of the side ground planes
  // They wrap around the absorber, so their inner edge aligns with the absorber outer edge.
  // Left: center = xabsorberoffset - absorberRadius - sideGroundx
  // Right: center = xabsorberoffset + absorberRadius + sideGroundx
  G4double xLeftGround = xabsorberoffset - absorberRadius - sideGroundx;
  G4double xRightGround = xabsorberoffset + absorberRadius + sideGroundx;

  G4VPhysicalVolume* leftSideGroundPhysical = new G4PVPlacement(
    0, G4ThreeVector(xLeftGround, -ySideGroundCenter, geHalfZ + thickness), 
    sideGroundlogical, "leftSideGroundPhysical",
    worldLogical, false, 0);

  G4VPhysicalVolume* rightSideGroundPhysical = new G4PVPlacement(
    0, G4ThreeVector(xRightGround, -ySideGroundCenter, geHalfZ + thickness), 
    sideGroundlogical, "rightSideGroundPhysical",
    worldLogical, false, 0);

  // Bottom bridge wrapping around the absorber
  const G4double bridgeInnerRadius = absorberRadius; 
  const G4double bridgeOuterRadius = absorberRadius + 2.0 * sideGroundx; // 20.0um wide
  
  G4Tubs* sideGroundBridge = new G4Tubs("sideGroundBridge", 
                                        bridgeInnerRadius, 
                                        bridgeOuterRadius, 
                                        thickness, 
                                        180.0*deg, 180.0*deg);
  G4LogicalVolume* sideGroundBridgeLogical =
    new G4LogicalVolume(sideGroundBridge, fAluminum, "sideGroundBridgeLogical");

  G4VPhysicalVolume* sideGroundBridgePhysical = new G4PVPlacement(
    0, G4ThreeVector(xabsorberoffset, -yabsorberoffset, geHalfZ + thickness),
    sideGroundBridgeLogical, "sideGroundBridgePhysical",
    worldLogical, false, 0);

  // ---- End of Geometry definition ----

  // 
  // detector -- Note : "sensitive detector" is attached to Germanium crystal
  // want a phonon sensitive detector, attached to Ge crystal
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  fGermaniumLogical->SetSensitiveDetector(fSuperconductorSensitivity);

  //
  // surface between Al and Ge determines phonon reflection/absorption
  //
  if (!fConstructed) {
    const G4double GHz = 1e9 * hertz; 

    //the following coefficients and cutoff values are not well-motivated
    //the code below is used only to demonstrate how to set these values.
    const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 0};
    const std::vector<G4double> diffCoeffs = {1, 0, 0, 0, 0, 0};
    const std::vector<G4double> specCoeffs = {0, 0, 0, 0, 0, 0};

    const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

    double pAbsProbAlSi = 0; //.488
    double pAbsProbSideWallSi = 0.01;//.01
    double pAbsProbPolishedWallSi = 0.0; //.0025

    AlSurfProp = new G4CMPSurfaceProperty("AlSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbAlSi, 1.0, 0.0, 0.0);
    AlSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);
    AttachPhononSensor_Al(AlSurfProp);

    NbSurfProp = new G4CMPSurfaceProperty("NbSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbAlSi, 1.0, 0.0, 0.0);
    NbSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);
    AttachPhononSensor_Nb(NbSurfProp);

    sidewallSurfProp = new G4CMPSurfaceProperty("SideWallSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbSideWallSi, 1.0, 0.0, 0.0 );
    sidewallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);
    
    
    polishedwallSurfProp = new G4CMPSurfaceProperty("polishedWallSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbPolishedWallSi, 1.0, 0.0, 0.0 );
    polishedwallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);

  }

  // Connects the inner volume, outer volume, and physics that applies at the surface
  // Logical border surface applies the specified physics for ANYWHERE the two volumes touch
  
  bool Al = true; 

  if (Al) {
    new G4CMPLogicalBorderSurface("Al", GePhys, alFLphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, alUGPphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, alLGPphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, capindconnectorphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, connectorphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, couplCapphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, leftSideGroundPhysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, rightSideGroundPhysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, sideGroundBridgePhysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, indvert1physical, AlSurfProp);
    for (auto phys : inductorPhys) {
      new G4CMPLogicalBorderSurface("Al", GePhys, phys, AlSurfProp);
    }
    
    new G4CMPLogicalBorderSurface("Al", GePhys, leftcapwallphysical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, rightcapwallphysical, AlSurfProp);
    for (auto phys : idcPhys) {
      new G4CMPLogicalBorderSurface("Al", GePhys, phys, AlSurfProp);
    }
    new G4CMPLogicalBorderSurface("Al", GePhys, botcapindconnector1physical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, botcapindconnector2physical, AlSurfProp);
    new G4CMPLogicalBorderSurface("Al", GePhys, capjunctconnectphysical, AlSurfProp);
  }

  if (!Al) {
    new G4CMPLogicalBorderSurface("Nb", GePhys, alFLphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, alUGPphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, alLGPphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, capindconnectorphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, connectorphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, couplCapphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, leftSideGroundPhysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, rightSideGroundPhysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, sideGroundBridgePhysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, indvert1physical, NbSurfProp);
    for (auto phys : inductorPhys) {
      new G4CMPLogicalBorderSurface("Nb", GePhys, phys, NbSurfProp);
    }
    
    new G4CMPLogicalBorderSurface("Nb", GePhys, leftcapwallphysical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, rightcapwallphysical, NbSurfProp);
    for (auto phys : idcPhys) {
      new G4CMPLogicalBorderSurface("Nb", GePhys, phys, NbSurfProp);
    }
    new G4CMPLogicalBorderSurface("Nb", GePhys, botcapindconnector1physical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, botcapindconnector2physical, NbSurfProp);
    new G4CMPLogicalBorderSurface("Nb", GePhys, capjunctconnectphysical, NbSurfProp);
  }

  // walls
  new G4CMPLogicalBorderSurface("GeToWorld", GePhys, fWorldPhys, polishedwallSurfProp);
  new G4CMPLogicalBorderSurface("GeToSideWall", GePhys, airSideRightXPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("GeToSideWall", GePhys, airSideLeftXPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("GeToSideWall", GePhys, airSideRightYPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("GeToSideWall", GePhys, airSideLeftYPhys, sidewallSurfProp);
  
  // sensor
  new G4CMPLogicalBorderSurface("Al", GePhys, junct1physical, AlSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, junct2physical, AlSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, junct3physical, AlSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, junct4physical, AlSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, junct5physical, AlSurfProp);
  new G4CMPLogicalBorderSurface("Al", GePhys, absorberphysical, AlSurfProp);


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
  
  // Air boxes: light blue, semi-transparent
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.1));
  airVis->SetVisibility(true);
  airSideLogicalX->SetVisAttributes(airVis);
  airSideLogicalY->SetVisAttributes(airVis);

  // Aluminum patterned parts
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  alVis->SetVisibility(true);
  alUGPlogical->SetVisAttributes(alVis);
  alLGPlogical->SetVisAttributes(alVis);
  alFLlogical->SetVisAttributes(alVis);
  connectorlogical->SetVisAttributes(alVis);
  couplCaplogical->SetVisAttributes(alVis);
  sideGroundlogical->SetVisAttributes(alVis);
  sideGroundBridgeLogical->SetVisAttributes(alVis);
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
  junct1logical->SetVisAttributes(alVis);
  junct2logical->SetVisAttributes(alVis);
  junct3logical->SetVisAttributes(alVis);
  junct4logical->SetVisAttributes(alVis);
  junct5logical->SetVisAttributes(alVis);
  absorberlogical->SetVisAttributes(alVis);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void RISQTutorialDetectorConstruction::
AttachPhononSensor_Al(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  double filmThicknessAl = 0;
  if (RISQTutorialConfigManager::GetfilmThicknessAl() != -1.0) filmThicknessAl = RISQTutorialConfigManager::GetfilmThicknessAl();

  std::cout<<"SR--Al filmThickness set to "<< filmThicknessAl <<std::endl;

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", .795);  //.795 
  sensorProp->AddConstProperty("filmThickness", filmThicknessAl*nm);
  sensorProp->AddConstProperty("gapEnergy", 173.715e-6*eV);
  sensorProp->AddConstProperty("lowQPLimit", 3.);
  sensorProp->AddConstProperty("phononLifetime", 242.*ps);
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
  sensorProp->AddConstProperty("vSound", 3.26*km/s);
  sensorProp->AddConstProperty("subgapAbsorption", 0.1);

  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

void RISQTutorialDetectorConstruction::
AttachPhononSensor_Nb(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  double filmThicknessNb = 0;
  if (RISQTutorialConfigManager::GetfilmThicknessNb() != -1.0) filmThicknessNb = RISQTutorialConfigManager::GetfilmThicknessNb();

  std::cout<<"SR--Nb filmThickness set to "<< filmThicknessNb <<std::endl;

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", .745);  //.745 
  sensorProp->AddConstProperty("filmThickness", filmThicknessNb*nm);
  sensorProp->AddConstProperty("gapEnergy",1538e-6*eV );            //SQD: From Eric, Dylan
  sensorProp->AddConstProperty("lowQPLimit",3.);                   //SQD: Taken from G4CMP phonon example (also aluminum).
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //SQD: From G4CMP phonon example (also aluminum), validated by Eric.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //REL: Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",2.44*CLHEP::km/CLHEP::s); //SQD: From Eric
  sensorProp->AddConstProperty("subgapAbsorption", 0.1);


  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

