// OW200127 11KID device

#include "RISQTutorialDetectorConstruction.hh"
#include "RISQTutorialSensitivity.hh"
#include "RISQTutorialQubitHousing.hh"
#include "RISQTutorialPad.hh"
#include "RISQTutorialTransmissionLine.hh"
#include "RISQTutorialStraightFluxLine.hh"
#include "RISQTutorialCornerFluxLine.hh"
#include "RISQTutorialResonatorAssembly.hh"
#include "RISQTutorialConfigManager.hh" //added for batch
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
  : fAir(0), fVacuum(0), fSilicon(0), fAluminum(0), fTungsten(0), fNiobium(0),
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
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",8.*cm,8.*cm,8.*cm); // half (-16,16)
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,fVacuum,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,
                                 false,0); // physical placement
  
  //                               
  // Silicon crystal - this is the volume in which we will propagate phonons
  //  
  G4double numsensors = 0;
  if (RISQTutorialConfigManager::Getnumsensors() != -1.0) numsensors = RISQTutorialConfigManager::Getnumsensors();
    
  const G4double siHalfX = 1.1*cm / std::sqrt(numsensors);
  const G4double siHalfY = 1.1*cm / std::sqrt(numsensors);
  const G4double siHalfZ = 0.05*mm;
  G4VSolid* fSiliconSolid = new G4Box("fSiliconSolid", siHalfX, siHalfY, siHalfZ);
  G4LogicalVolume* fSiliconLogical = new G4LogicalVolume(fSiliconSolid,fSilicon,"fSiliconLogical");
  G4VPhysicalVolume* SiPhys = new G4PVPlacement(0,G4ThreeVector(),fSiliconLogical,"fSiliconPhysical", worldLogical,false,0); 
  // placing physical volume at center of world logical

  //
  //Silicon lattice information
  //

  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* SiLogical = LM->LoadLattice(fSilicon, "Si");

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
  SiPhysical->SetMillerOrientation(1,0,0); // how crystal is oriented, also 4 coord orient. (online calculator)
  LM->RegisterLattice(SiPhys, SiPhysical); // connects physical lattice to volume

  // NOTE:  Above registration can also be done in single step:
  // G4LatticlePhysical* SiPhysical = LM->LoadLattice(SiPhys, "Si");

  //
  // Air boxes touching the 4 side faces of the Si crystal (+X, -X, +Y, -Y)
  //
  const G4double airSideThickness = 0.5*cm; // Distance the air extends away from the chip
  
  // X boxes cover the +X and -X faces of the Si crystal
  G4VSolid* airSideSolidX = new G4Box("airSideSolidX", airSideThickness/2, siHalfY, siHalfZ);
  G4LogicalVolume* airSideLogicalX = new G4LogicalVolume(airSideSolidX, fAir, "airSideLogicalX");
  G4VPhysicalVolume* airSideRightXPhys = new G4PVPlacement(0, G4ThreeVector(siHalfX + airSideThickness/2, 0, 0), airSideLogicalX, "airSideRightX", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftXPhys = new G4PVPlacement(0, G4ThreeVector(-siHalfX - airSideThickness/2, 0, 0), airSideLogicalX, "airSideLeftX", worldLogical, false, 0);

  // Y boxes cover the +Y and -Y faces of the Si crystal
  G4VSolid* airSideSolidY = new G4Box("airSideSolidY", siHalfX, airSideThickness/2, siHalfZ);
  G4LogicalVolume* airSideLogicalY = new G4LogicalVolume(airSideSolidY, fAir, "airSideLogicalY");
  G4VPhysicalVolume* airSideRightYPhys = new G4PVPlacement(0, G4ThreeVector(0, siHalfY + airSideThickness/2, 0), airSideLogicalY, "airSideRightY", worldLogical, false, 0);
  G4VPhysicalVolume* airSideLeftYPhys = new G4PVPlacement(0, G4ThreeVector(0, -siHalfY - airSideThickness/2, 0), airSideLogicalY, "airSideLeftY", worldLogical, false, 0);

  // OW200127 - Load all STL parts first

  auto sensor = CADMesh::TessellatedMesh::FromSTL("../../OW200127/OW200127_1.STL");
  sensor->SetScale(1e-3);
  G4VSolid* sensor_solid = sensor->GetSolid();

  auto kidFeedline = CADMesh::TessellatedMesh::FromSTL("../../OW200127/OW200127_2.STL");
  kidFeedline->SetScale(1e-3 / 25.4); // corrects accidental um-to-inches (not mm) export
  G4VSolid* kidFeedline_solid = kidFeedline->GetSolid();

  auto otherKIDs = CADMesh::TessellatedMesh::FromSTL("../../OW200127/OW200127_3.STL");
  otherKIDs->SetScale(1e-3);
  G4VSolid* otherKIDs_solid = otherKIDs->GetSolid();

  // Compute overall max and min of all OW200127 parts to find the center
  std::vector<G4VSolid*> ow200127Solids = {sensor_solid, kidFeedline_solid, otherKIDs_solid};
  G4ThreeVector overallMin(DBL_MAX, DBL_MAX, DBL_MAX);
  G4ThreeVector overallMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (auto* solid : ow200127Solids) {
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

  // Compute offset to place OW200127 center at desired position
  // Target Z: bottom of OW200127 sits on top of Si surface (at z = siHalfZ)
  const G4double targetCenterX = 0.0*um;    // Desired X position of OW200127 center
  const G4double targetCenterY = 0.0*um; // Desired Y position: 200um below feedline
  const G4double targetCenterZ = siHalfZ + (overallMax.z() - overallMin.z()) / 2.0; // Bottom of OW200127 at Si surface
  
  G4ThreeVector targetCenter(targetCenterX, targetCenterY, targetCenterZ);
  const G4ThreeVector ow200127Offset = targetCenter - stlCenter - G4ThreeVector(0, 0, 0.1 * nm); // Force 0.1 nm geometric overlap;

  // Create logical volumes and place physical volumes
  G4LogicalVolume* sensorlogical = new G4LogicalVolume(sensor_solid,fAluminum,"sensorlogical"); 
  G4VPhysicalVolume* sensorphysical = new G4PVPlacement(0, ow200127Offset, sensorlogical, "sensorphysicalshunt", worldLogical, false, 0);

  G4LogicalVolume* kidFeedlinelogical = new G4LogicalVolume(kidFeedline_solid,fNiobium,"kidFeedlinelogical"); 
  G4VPhysicalVolume* kidFeedlinephysical = new G4PVPlacement(0, ow200127Offset, kidFeedlinelogical, "kidFeedlinephysical", worldLogical, false, 0);

  G4LogicalVolume* otherKIDslogical = new G4LogicalVolume(otherKIDs_solid,fNiobium,"otherKIDslogical"); 
  G4VPhysicalVolume* otherKIDsphysical = new G4PVPlacement(0, ow200127Offset, otherKIDslogical, "otherKIDsphysical", worldLogical, false, 0);


  // 
  // detector -- Note : "sensitive detector" is attached to Silicon crystal
  // want a phonon sensitive detector, attached to Si crystal
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  fSiliconLogical->SetSensitiveDetector(fSuperconductorSensitivity);

  //
  // surface between Al and Si determines phonon reflection/absorption
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
    double pAbsProbSideWallSi = 0;//.01
    if (RISQTutorialConfigManager::GetpAbsProbSideWallSi() != -1.0) pAbsProbSideWallSi = RISQTutorialConfigManager::GetpAbsProbSideWallSi();
    double pAbsProbPolishedWallSi = 0; //.0025
    if (RISQTutorialConfigManager::GetpAbsProbPolishedWallSi() != -1.0) pAbsProbPolishedWallSi = RISQTutorialConfigManager::GetpAbsProbPolishedWallSi();

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
  //
  // Si -> Al/Nb (phonons start in Si and enter the superconductor)
  new G4CMPLogicalBorderSurface("SiToAl_Sensor", SiPhys, sensorphysical, AlSurfProp);
  new G4CMPLogicalBorderSurface("SiToNb_Feedline", SiPhys, kidFeedlinephysical, NbSurfProp);
  new G4CMPLogicalBorderSurface("SiToNb_OtherKIDs", SiPhys, otherKIDsphysical, NbSurfProp);
  


  // Si -> World (bare Si where there is no Al coverage)
  new G4CMPLogicalBorderSurface("SiToWorld", SiPhys, fWorldPhys, polishedwallSurfProp);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideRightXPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideLeftXPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideRightYPhys, sidewallSurfProp);
  new G4CMPLogicalBorderSurface("SiToSideWall", SiPhys, airSideLeftYPhys, sidewallSurfProp);



//                                        
// Visualization attributes
//


  // World remains invisible
  G4VisAttributes* wrldVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.1, 0.1));
  wrldVis->SetVisibility(false);
  worldLogical->SetVisAttributes(wrldVis);

  // Silicon crystal: light gray, solid
  G4VisAttributes* siVis = new G4VisAttributes(G4Colour(0.85, 0.85, 0.85, 0.4));
  siVis->SetVisibility(true);
  fSiliconLogical->SetVisAttributes(siVis);
  
  // Air boxes: light blue, semi-transparent
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.1));
  airVis->SetVisibility(true);
  airSideLogicalX->SetVisAttributes(airVis);
  airSideLogicalY->SetVisAttributes(airVis);

  // Aluminum/Niobium patterned parts
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  alVis->SetVisibility(true);
  sensorlogical->SetVisAttributes(alVis);
  kidFeedlinelogical->SetVisAttributes(alVis);
  otherKIDslogical->SetVisAttributes(alVis);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void RISQTutorialDetectorConstruction::
AttachPhononSensor_Al(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  double filmThicknessAl = 30*nm;
  if (RISQTutorialConfigManager::GetfilmThicknessAl() != -1.0) filmThicknessAl = RISQTutorialConfigManager::GetfilmThicknessAl();

  std::cout<<"SR--Al filmThickness set to "<< filmThicknessAl <<std::endl;

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", .795);  //.795 
  sensorProp->AddConstProperty("filmThickness", filmThicknessAl*nm);
  sensorProp->AddConstProperty("gapEnergy", 173.715e-6*eV);
  sensorProp->AddConstProperty("lowQPLimit", 1.1);
  sensorProp->AddConstProperty("phononLifetime", 242.*ps);
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
  sensorProp->AddConstProperty("vSound", 3.26*km/s);
  sensorProp->AddConstProperty("subgapAbsorption", 0.0);
  sensorProp->AddConstProperty("absorberGap", 0.0);

  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

void RISQTutorialDetectorConstruction::
AttachPhononSensor_Nb(G4CMPSurfaceProperty *surfProp) {
  if (!surfProp) return;		// No surface, nothing to do

  double filmThicknessNb = 30*nm;
  if (RISQTutorialConfigManager::GetfilmThicknessNb() != -1.0) filmThicknessNb = RISQTutorialConfigManager::GetfilmThicknessNb();

  std::cout<<"SR--Nb filmThickness set to "<< filmThicknessNb <<std::endl;

  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", .745);  //.745 
  sensorProp->AddConstProperty("filmThickness", filmThicknessNb*nm);
  sensorProp->AddConstProperty("gapEnergy",1538e-6*eV );            //SQD: From Eric, Dylan
  sensorProp->AddConstProperty("lowQPLimit",1.1);                   //SQD: Taken from G4CMP phonon example (also aluminum).
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //SQD: From G4CMP phonon example (also aluminum), validated by Eric.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //REL: Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",2.44*CLHEP::km/CLHEP::s); //SQD: From Eric
  sensorProp->AddConstProperty("subgapAbsorption", 0.0);
  sensorProp->AddConstProperty("absorberGap", 0.0);


  // Attach electrode object to handle KaplanQP interface
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}
