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

using namespace RISQTutorialDetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::RISQTutorialDetectorConstruction()
  : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
    fWorldPhys(0), 
    fSuperconductorSensitivity(0), fConstructed(false) {;}//, fIfField(true) {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RISQTutorialDetectorConstruction::~RISQTutorialDetectorConstruction() {;}

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
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
  fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RISQTutorialDetectorConstruction::SetupGeometry()
{



  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // First, define border surface properties that can be referenced later
  const G4double GHz = 1e9 * hertz; 

  //the following coefficients and cutoff values are not well-motivated
  //the code below is used only to demonstrate how to set these values.
  const std::vector<G4double> anhCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const std::vector<G4double> diffCoeffs = {1,0,0,0,0,0};//Explicitly make this 1 for now
  const std::vector<G4double> specCoeffs = {0,0,0,0,0,0};//Turn this off temporarily
  const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external
    
  
  //These are just the definitions of the interface TYPES, not the interfaces themselves. These must be called in a set of loops
  //below, and invoke these surface definitions.
  if( !fConstructed ){
    fSiNbInterface = new G4CMPSurfaceProperty("SiNbInterface",
					      1.0, 0.0, 0.0, 0.0,
					      0.5, 1.0, 0.0, 0.0);
    fSiCopperInterface = new G4CMPSurfaceProperty("SiCopperInterface",
					      1.0, 0.0, 0.0, 0.0,
					      0.5, 1.0, 0.0, 0.0);

    double pAbsCircuit = 0.01;
    double pAbsMount = 0.001;
    double pAbsSilicon = 0.0001; 

    fSiFeedInterface = new G4CMPSurfaceProperty("SiFeedInterface",
					      1.0, 0.0, 0.0, 0.0,
					      pAbsCircuit, 1.0, 0.0, 0.0);
    fSiIndInterface = new G4CMPSurfaceProperty("SiIndInterface",
						  1.0, 0.0, 0.0, 0.0,
						  pAbsCircuit, 1.0, 0.0, 0.0 );
    fSiCapInterface = new G4CMPSurfaceProperty("SiCapInterface",
						  1.0, 0.0, 0.0, 0.0,
						  pAbsCircuit, 1.0, 0.0, 0.0 );
    fSiRightMountInterface = new G4CMPSurfaceProperty("SiRightMountInterface",
						  1.0, 0.0, 0.0, 0.0,
						  pAbsMount, 1.0, 0.0, 0.0 );
    fSiLeftMountInterface = new G4CMPSurfaceProperty("SiLeftMountInterface",
						  1.0, 0.0, 0.0, 0.0,
						  pAbsMount, 1.0, 0.0, 0.0 );          
    fSiVacuumInterface = new G4CMPSurfaceProperty("SiVacuumInterface",
						  0.0, 1.0, 0.0, 0.0,
						  pAbsSilicon, 1.0, 0.0, 0.0 );
    

    fSiFeedInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
					    diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiIndInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);  
    fSiCapInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);
    fSiRightMountInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);
    fSiLeftMountInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);
    fSiVacuumInterface->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
						diffCoeffs, specCoeffs, GHz, GHz, GHz);

    //Add a phonon sensor to the interface properties here.
    AttachPhononSensor(fSiFeedInterface);
    AttachPhononSensor(fSiIndInterface);
    AttachPhononSensor(fSiCapInterface);
    AttachPhononSensor(fSiRightMountInterface);
    AttachPhononSensor(fSiLeftMountInterface);
    AttachPhononSensor(fSiVacuumInterface);
  }


  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we start constructing the various components and their interfaces  
  //     
  // World
  //
  G4VSolid* solid_world = new G4Box("World",55.*cm,55.*cm,55.*cm);
  G4LogicalVolume* log_world = new G4LogicalVolume(solid_world,fLiquidHelium,"World");
  //  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0,
				 G4ThreeVector(),
				 log_world,
				 "World",
				 0,
                                 false,
				 0);
  
  
  bool checkOverlaps = true;


  //-------------------------------------------------------------------------------------------------------------------
  //First, set up the qubit chip substrate. By default, assume that we're using this. Otherwise, it's hard to establish
  //a sensitivity object for this.
  G4Box * solid_siliconChip = new G4Box("QubitChip_solid",
					0.5*dp_siliconChipDimX,
					0.5*dp_siliconChipDimY,
					0.5*dp_siliconChipDimZ);
  
  //Now attribute a physical material to the chip
  G4LogicalVolume * log_siliconChip = new G4LogicalVolume(solid_siliconChip,
							  fSilicon,
							  "SiliconChip_log");
    
  //Now, create a physical volume and G4PVPlacement for storing as the final output
  // G4ThreeVector siliconChipTranslate(0,0,0.5*(dp_housingDimZ - dp_siliconChipDimZ) + dp_eps); 

  G4ThreeVector siliconChipTranslate(0,0,0); 

  G4VPhysicalVolume * phys_siliconChip = new G4PVPlacement(0,
							   siliconChipTranslate,
							   log_siliconChip,
							   "SiliconChip", 
							   log_world,
							   false,
							   0,
							   checkOverlaps);

  G4VisAttributes* siliconChipVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  siliconChipVisAtt->SetVisibility(true);
  log_siliconChip->SetVisAttributes(siliconChipVisAtt);



  //Set up the G4CMP silicon lattice information using the G4LatticeManager
  // G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical* log_siliconLattice = LM->LoadLattice(fSilicon, "Si");
    
  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical* phys_siliconLattice = new G4LatticePhysical(log_siliconLattice);
  phys_siliconLattice->SetMillerOrientation(1,0,0); 
  LM->RegisterLattice(phys_siliconChip,phys_siliconLattice);

  //Set up border surfaces
  G4CMPLogicalBorderSurface * border_siliconChip_world = new G4CMPLogicalBorderSurface("border_siliconChip_world", phys_siliconChip, fWorldPhys, fSiVacuumInterface);


  G4double siliconThickness = 0.05 * cm / 2.0;
  G4double siliconWidth = 2.0 * cm / 2.0; 
  float indLength = 0.72 * mm / 2.0f;
  float capLength = 0.43 * mm / 2.0f;
  float kidWidth = 1.13 * mm / 2.0f;
  float kidThickness = 0.03 * um / 2.0f;
  float feedThickness = 0.07 * um;
  float feedWidth = 0.2 * mm / 2.0f;
  float feedIndGap = 30 * um;
  float capIndSpacing = 80 * um;
  float mountWidth = 1.0 * mm / 2.0f;

  //Inductor. 
  G4VSolid* fIndSolid = new G4Box("Ind", kidWidth, indLength, kidThickness);
  G4LogicalVolume* fIndLogical = new G4LogicalVolume(fIndSolid, fAluminum, "fIndLogical");
  G4VPhysicalVolume* IndPhys = new G4PVPlacement(0, G4ThreeVector(0, indLength + feedWidth/2.0 + feedIndGap, siliconThickness + kidThickness),
    fIndLogical, "fIndLogical", log_world, false, 0, checkOverlaps);

  //Capacitor. 
  //capLength + 2.0 * indLength + capIndSpacing
  G4VSolid* fCapSolid = new G4Box("Cap", kidWidth, capLength, kidThickness);
  G4LogicalVolume* fCapLogical = new G4LogicalVolume(fCapSolid, fAluminum, "fCapLogical");
  G4VPhysicalVolume* CapPhys = new G4PVPlacement(0,  G4ThreeVector(0, capLength + feedWidth/2.0 + feedIndGap + (2 * indLength) + capIndSpacing, siliconThickness + kidThickness),
    fCapLogical, "fCapLogical", log_world, false, 0, checkOverlaps);
  
  // Feedline.
  G4VSolid* fFeedSolid = new G4Box("Feed", siliconWidth, feedWidth, feedThickness);
  G4LogicalVolume* fFeedLogical = new G4LogicalVolume(fFeedSolid, fAluminum, "fFeedLogical");
  G4VPhysicalVolume* FeedPhys = new G4PVPlacement(0,  G4ThreeVector(0,-feedWidth/2.0, siliconThickness + feedThickness),
  fFeedLogical, "fFeedLogical", log_world, false, 0, checkOverlaps);

  // Right Mount
  G4VSolid* fRightMountSolid = new G4Box("RightMount", mountWidth, mountWidth, mountWidth/2);
  G4LogicalVolume* fRightMountLogical = new G4LogicalVolume(fRightMountSolid, fTeflon, "fRightMountSolid");
  G4VPhysicalVolume* RightMountPhys = new G4PVPlacement(0, G4ThreeVector(siliconWidth - mountWidth, siliconWidth - mountWidth, siliconThickness + mountWidth/2), fRightMountLogical, "fRightMountLogical", log_world, false, 0, checkOverlaps);

  // Left Mount
  G4VSolid* fLeftMountSolid = new G4Box("LeftMount", mountWidth, mountWidth, mountWidth/2);
  G4LogicalVolume* fLeftMountLogical = new G4LogicalVolume(fLeftMountSolid, fTeflon, "fLeftMountSolid");
  G4VPhysicalVolume* LeftMountPhys = new G4PVPlacement(0, G4ThreeVector(mountWidth - siliconWidth, mountWidth - siliconWidth, siliconThickness + mountWidth/2), fLeftMountLogical, "fLeftMountLogical", log_world, false, 0, checkOverlaps);

  G4CMPLogicalBorderSurface * detInd = new G4CMPLogicalBorderSurface("detInd", phys_siliconChip, IndPhys, fSiIndInterface);
  G4CMPLogicalBorderSurface * detCap = new G4CMPLogicalBorderSurface("detCap", phys_siliconChip, CapPhys, fSiCapInterface);
  G4CMPLogicalBorderSurface * detFeed = new G4CMPLogicalBorderSurface("detFeed", phys_siliconChip, FeedPhys, fSiFeedInterface);
  G4CMPLogicalBorderSurface * detRightMount = new G4CMPLogicalBorderSurface("detRightMount", phys_siliconChip, RightMountPhys, fSiRightMountInterface);
  G4CMPLogicalBorderSurface * detLeftMount = new G4CMPLogicalBorderSurface("detLeftMount", phys_siliconChip, LeftMountPhys, fSiLeftMountInterface);

  //---------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------
  // Now we establish a sensitivity object
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQTutorialSensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  log_siliconChip->SetSensitiveDetector(fSuperconductorSensitivity);



}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Set up a phonon sensor for this surface property object. I'm pretty sure that this
// phonon sensor doesn't get stapled to individual geometrical objects, but rather gets
// stapled to a surface property, but I'm not sure... have to ask mKelsey
void RISQTutorialDetectorConstruction::AttachPhononSensor(G4CMPSurfaceProperty * surfProp)
{
  //If no surface, don't do anything
  if(!surfProp) return;

  //Specify properties of the niobium sensors
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption",0.0);              //NOT WELL MOTIVATED - probably parametrize and put on slider?
  sensorProp->AddConstProperty("filmThickness",90.*CLHEP::nm);     //Accurate for our thin film.
  sensorProp->AddConstProperty("gapEnergy",1.6e-3*CLHEP::eV);       //Reasonably motivated. Actually, looks like Novotny and Meincke are quoting 2Delta, and this is delta. Nuss and Goossen mention that Nb has a delta value closer to this.
  sensorProp->AddConstProperty("lowQPLimit",3.);                   //NOT WELL MOTIVATED YET -- Dunno how to inform this...
  sensorProp->AddConstProperty("phononLifetime",4.17*CLHEP::ps);   //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  sensorProp->AddConstProperty("phononLifetimeSlope",0.29);        //Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound",3.480*CLHEP::km/CLHEP::s); //True for room temperature, probably good to 10%ish - should follow up
  sensorProp->AddConstProperty("subgapAbsorption",0.0);            //Assuming that since we're mostly sensitive to quasiparticle density, phonon "heat" here isn't something that we're sensitive to? Unsure how to select this.

  //  sensorProp->AddConstProperty("gapEnergy",3.0e-3*CLHEP::eV);      //Reasonably motivated. Novotny and Meincke, 1975 (2.8-3.14 meV)
  //  sensorProp->AddConstProperty("phononLifetime",242.*ps);      //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  
  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
  
}
