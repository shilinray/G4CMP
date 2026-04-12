/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File:  RISQTutorialConfigMessenger.cc
//
// Description:	Macro command defitions to set user configuration in
//		RISQTutorialConfigManager.
//
// 20170816  Michael Kelsey

#include "RISQTutorialConfigMessenger.hh"
#include "RISQTutorialConfigManager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"


// Constructor and destructor

RISQTutorialConfigMessenger::RISQTutorialConfigMessenger(RISQTutorialConfigManager* mgr)
  : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
    theManager(mgr), primFileNameCmd(0), hitsFileNameCmd(0) {
  primFileNameCmd = CreateCommand<G4UIcmdWithAString>("primFileName",
			      "Set text file name for primary particle information");
  hitsFileNameCmd = CreateCommand<G4UIcmdWithAString>("hitsFileName",
			      "Set text file name for phonon hit information");

  filmThicknessAlCmd = CreateCommand<G4UIcmdWithADouble> ("filmThicknessAl", 
            "Set the filmThickness parameter for Al");
  filmThicknessNbCmd = CreateCommand<G4UIcmdWithADouble> ("filmThicknessNb", 
            "Set the filmThickness parameter for Nb");
  AlCmd = CreateCommand<G4UIcmdWithABool> ("Al", 
            "Set the boolean parameter for Al");
  pAbsProbSideWallSiCmd = CreateCommand<G4UIcmdWithADouble> ("pAbsProbSideWallSi", 
            "Set the parameter for pAbsProbSideWallSi");
  pAbsProbPolishedWallSiCmd = CreateCommand<G4UIcmdWithADouble> ("pAbsProbPolishedWallSi", 
            "Set the parameter for pAbsProbPolishedWallSi");
  numsensorsCmd = CreateCommand<G4UIcmdWithADouble> ("numsensors", 
            "Set the parameter for numsensors");
}


RISQTutorialConfigMessenger::~RISQTutorialConfigMessenger() {
  delete primFileNameCmd; primFileNameCmd=0;
  delete hitsFileNameCmd; hitsFileNameCmd=0;
  delete filmThicknessAlCmd; filmThicknessAlCmd=0;
  delete filmThicknessNbCmd; filmThicknessNbCmd=0;
  delete AlCmd; AlCmd=0;
  delete pAbsProbSideWallSiCmd; pAbsProbSideWallSiCmd=0;
  delete pAbsProbPolishedWallSiCmd; pAbsProbPolishedWallSiCmd=0;
  delete numsensorsCmd; numsensorsCmd=0;
}


// Parse user input and add to configuration

void RISQTutorialConfigMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == primFileNameCmd) theManager->SetPrimaryOutput(value);
  if (cmd == hitsFileNameCmd) theManager->SetHitOutput(value);
  if (cmd == filmThicknessAlCmd) theManager->SetfilmThicknessAl(filmThicknessAlCmd->GetNewDoubleValue(value));
  if (cmd == filmThicknessNbCmd) theManager->SetfilmThicknessNb(filmThicknessNbCmd->GetNewDoubleValue(value));
  if (cmd == AlCmd) theManager->SetAl(AlCmd->GetNewBoolValue(value));
  if (cmd == pAbsProbSideWallSiCmd) theManager->SetpAbsProbSideWallSi(pAbsProbSideWallSiCmd->GetNewDoubleValue(value));
  if (cmd == pAbsProbPolishedWallSiCmd) theManager->SetpAbsProbPolishedWallSi(pAbsProbPolishedWallSiCmd->GetNewDoubleValue(value));
  if (cmd == numsensorsCmd) theManager->Setnumsensors(numsensorsCmd->GetNewDoubleValue(value));
 
}
