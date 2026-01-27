// SQUAT

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
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
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), topSurfProp(0), botSurfProp(0), wallSurfProp(0), 
    fSuperconductorSensitivity(0), fConstructed(false) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {
  delete topSurfProp;
  delete botSurfProp;
  delete wallSurfProp;
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

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
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
  const G4double alFeedlineHalfY = 4.5*um;
  const G4double abs_thickness = 0.2*um;

  G4Box* feedline = new G4Box("feedlineCenter", alFeedlineHalfX, alFeedlineHalfY, abs_thickness);
 
  G4LogicalVolume* alFLlogical = new G4LogicalVolume(feedline,fAluminum,"alFLlogical"); // logical feedline

  G4VPhysicalVolume* alFLphysical = new G4PVPlacement(
    0, G4ThreeVector(0.,0., geHalfZ + abs_thickness), alFLlogical, "alFLphysical",
    worldLogical, false, 0); // physical feedline

  // Helper to print bounding box + volume for CADMesh tessellated solids
  const auto printSolidInfo = [](const char* label, G4VSolid* solid) {
    G4ThreeVector pMin, pMax;
    solid->BoundingLimits(pMin, pMax);
    const G4ThreeVector size = pMax - pMin;

    G4cout << label << " bbox min=" << (pMin/mm) << " mm"
           << " max=" << (pMax/mm) << " mm"
           << " size=" << (size/mm) << " mm" << G4endl;
    G4cout << label << " volume=" << (solid->GetCubicVolume()/mm3) << " mm^3" << G4endl;
  };

  // SQUAT - Load all STL parts first
  const G4double trap_thickness = .02*um;

  auto leftabs = CADMesh::TessellatedMesh::FromSTL("../../single_squat/single_squat_BE1.STL");
  leftabs->SetScale(1e-3);
  G4VSolid* leftabs_solid = leftabs->GetSolid();

  auto righttrap = CADMesh::TessellatedMesh::FromSTL("../../single_squat/single_squat_BE2.STL");
  righttrap->SetScale(1e-3);
  G4VSolid* righttrap_solid = righttrap->GetSolid();

  auto lefttrap = CADMesh::TessellatedMesh::FromSTL("../../single_squat/single_squat_BE3.STL");
  lefttrap->SetScale(1e-3);
  G4VSolid* lefttrap_solid = lefttrap->GetSolid();

  auto junction = CADMesh::TessellatedMesh::FromSTL("../../single_squat/single_squat_BE4.STL");
  junction->SetScale(1e-3);
  G4VSolid* junction_solid = junction->GetSolid();

  auto rightabs = CADMesh::TessellatedMesh::FromSTL("../../single_squat/single_squat_BE5.STL");
  rightabs->SetScale(1e-3);
  G4VSolid* rightabs_solid = rightabs->GetSolid();

  // Compute overall bounding box of all SQUAT parts to find the center
  std::vector<G4VSolid*> squatSolids = {leftabs_solid, righttrap_solid, lefttrap_solid, 
                                         junction_solid, rightabs_solid};
  G4ThreeVector overallMin(DBL_MAX, DBL_MAX, DBL_MAX);
  G4ThreeVector overallMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (auto* solid : squatSolids) {
    G4ThreeVector pMin, pMax;
    solid->BoundingLimits(pMin, pMax);
    overallMin.setX(std::min(overallMin.x(), pMin.x()));
    overallMin.setY(std::min(overallMin.y(), pMin.y()));
    overallMin.setZ(std::min(overallMin.z(), pMin.z()));
    overallMax.setX(std::max(overallMax.x(), pMax.x()));
    overallMax.setY(std::max(overallMax.y(), pMax.y()));
    overallMax.setZ(std::max(overallMax.z(), pMax.z()));
  }
  G4ThreeVector stlCenter = 0.5 * (overallMin + overallMax);
  
  G4cout << "SQUAT overall bbox min=" << (overallMin/um) << " um"
         << " max=" << (overallMax/um) << " um" << G4endl;
  G4cout << "SQUAT STL center=" << (stlCenter/um) << " um" << G4endl;

  // Compute offset to place SQUAT center at desired position
  // Target Z: bottom of SQUAT sits on top of Ge surface (at z = geHalfZ)
  const G4double targetCenterX = 0.0*um;    // Desired X position of SQUAT center
  const G4double targetCenterY = -200.0*um; // Desired Y position: 200um below feedline
  const G4double targetCenterZ = geHalfZ + (overallMax.z() - overallMin.z()) / 2.0; // Bottom of SQUAT at Ge surface
  const G4ThreeVector squatOffset(targetCenterX - stlCenter.x(), 
                                   targetCenterY - stlCenter.y(), 
                                   targetCenterZ - stlCenter.z());

  G4cout << "SQUAT placement offset=" << (squatOffset/um) << " um" << G4endl;

  // Print individual solid info
  printSolidInfo("leftabs_solid", leftabs_solid);
  printSolidInfo("righttrap_solid", righttrap_solid);
  printSolidInfo("lefttrap_solid", lefttrap_solid);
  printSolidInfo("junction_solid", junction_solid);
  printSolidInfo("rightabs_solid", rightabs_solid);

  // Create logical volumes and place physical volumes
  G4LogicalVolume* leftabslogical = new G4LogicalVolume(leftabs_solid,fAluminum,"leftabslogical"); 
  G4VPhysicalVolume* leftabsphysical = new G4PVPlacement(0, squatOffset, leftabslogical, "leftabsphysicalshunt", worldLogical, false, 0);

  G4LogicalVolume* righttraplogical = new G4LogicalVolume(righttrap_solid,fAluminum,"righttraplogical"); 
  G4VPhysicalVolume* righttrapphysical = new G4PVPlacement(0, squatOffset, righttraplogical, "righttrapphysicalshunt", worldLogical, false, 0);

  G4LogicalVolume* lefttraplogical = new G4LogicalVolume(lefttrap_solid,fAluminum,"lefttraplogical"); 
  G4VPhysicalVolume* lefttrapphysical = new G4PVPlacement(0, squatOffset, lefttraplogical, "lefttrapphysicalshunt", worldLogical, false, 0);

  G4LogicalVolume* junctionlogical = new G4LogicalVolume(junction_solid,fAluminum,"junctionlogical"); 
  G4VPhysicalVolume* junctionphysical = new G4PVPlacement(0, squatOffset, junctionlogical, "junctionphysicalshunt", worldLogical, false, 0);

  G4LogicalVolume* rightabslogical = new G4LogicalVolume(rightabs_solid,fAluminum,"rightabslogical"); 
  G4VPhysicalVolume* rightabsphysical = new G4PVPlacement(0, squatOffset, rightabslogical, "rightabsphysicalshunt", worldLogical, false, 0);


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

    double pAbsProbAlSi = 0.488;
    double pAbsProbSideWallSi = 0.0;

    topSurfProp = new G4CMPSurfaceProperty("TopAlSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbAlSi, 1.0, 0.0, 0.0);
    topSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					 diffCoeffs, specCoeffs, GHz, GHz, GHz);
    AttachPhononSensor(topSurfProp);

    wallSurfProp = new G4CMPSurfaceProperty("WallSurf", 0.0, 1.0, 0.0, 0.0,
                                                      pAbsProbSideWallSi, 1.0, 0.0, 0.0 );
    wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					  diffCoeffs, specCoeffs, GHz, GHz,GHz);

  }

  // Connects the inner volume, outer volume, and physics that applies at the surface
  // Logical border surface applies the specified physics for ANYWHERE the two volumes touch
  //
  // Ge -> Al (phonons start in Ge and enter Al)
  new G4CMPLogicalBorderSurface("GeToAl_FL", GePhys, alFLphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("GeToAl_LA", GePhys, leftabsphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("GeToAl_RT", GePhys, righttrapphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("GeToAl_LT", GePhys, lefttrapphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("GeToAl_JN", GePhys, junctionphysical, topSurfProp);
  new G4CMPLogicalBorderSurface("GeToAl_RA", GePhys, rightabsphysical, topSurfProp);

  // Ge -> World (bare Ge where there is no Al coverage)
  new G4CMPLogicalBorderSurface("GeToWorld", GePhys, fWorldPhys, wallSurfProp);

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
  alFLlogical->SetVisAttributes(alVis);
  leftabslogical->SetVisAttributes(alVis);
  righttraplogical->SetVisAttributes(alVis);
  lefttraplogical->SetVisAttributes(alVis);
  junctionlogical->SetVisAttributes(alVis);
  rightabslogical->SetVisAttributes(alVis);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void RISQTutorialDetectorConstruction::
AttachPhononSensor(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  // Specify properties of aluminum sensor, same on both detector faces
  // See G4CMPPhononElectrode.hh or README.md for property keys

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", 0); // 0.745
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


