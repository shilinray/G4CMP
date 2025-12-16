//---------------------------------------------------------
//
// RISQTutorialAnalysis.cc
// linehan3@fnal.gov
// 5/23/2024
//
// This is an analysis macro that is built to analyze
// the output of the G4CMP simulations that are run in
// the RISQ G4CMP Tutorial code.
//
//---------------------------------------------------------

//C++ includes
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>

//ROOT includes
#include "TH2F.h"
#include "TH1F.h"

//---------------------------------------------------------------------------------------
// Define a set of structs for use interpreting the output from G4CMP
struct Hit
{
  int runID;
  int eventID;
  int trackID;
  std::string particleName;
  double startEnergy_eV;
  double startX_mm;
  double startY_mm;
  double startZ_mm;
  double startT_ns;
  double eDep_eV;
  double trackWeight;
  double endX_mm;
  double endY_mm;
  double endZ_mm;
  double endT_ns;
};

struct PrimaryInfo
{
  int runID;
  int eventID;
  int trackID;
  std::string particleName;
  double energy_eV;
  double X_mm;
  double Y_mm;
  double Z_mm;
  double T_ns;  
};

struct Event
{
  int runID;
  int eventID;
  std::vector<Hit> hitVect;
  PrimaryInfo thePrim;
};

//Forward declarations
std::vector<Event> ReadInG4CMPPrimaryAndHitFiles(std::string hitTextFilename, std::string primaryTextFilename);
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename);
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename);


bool compareByEndTime(const Hit &a, Hit &b) {
  return a.endT_ns < b.endT_ns;
}

//---------------------------------------------------------------------------------------
// Analysis Script: Phonon Collection Efficiency
void PCEStudy(std::string primariesFilename, std::string hitsFilename)
{
  //Read in the text files
  std::vector<Event> eventList = ReadInG4CMPPrimaryAndHitFiles(hitsFilename,primariesFilename);
  std::cout << "Done reading in events." << std::endl;


  //Define an outfile
  TFile * fOut = new TFile("AnalysisOutput.root","RECREATE");


  std::vector<double> etaPhonons;

  std::vector<double> etaPhononSilicon;
  std::vector<double> etaPhononCapacitor;
  std::vector<double> etaPhononFeedline;
  std::vector<double> etaPhononMounts;

  //Loop over events
  for( int iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;
    
    double eDepInd = 0.0;
    double eDepSilicon = 0.0;
    double eDepCapacitor = 0.0;
    double eDepFeedline = 0.0;
    double eDepMounts = 0.0;

    //Get the event
    Event tE = eventList[iE];

    std::sort(tE.hitVect.begin(), tE.hitVect.end(), compareByEndTime);
    
    //Add to the primary vector   
    
    //Plot a number of hit-related things: hit multiplicity, hit locations in XYZ, hits in XYZ weighted by energy, etc.
    for ( int iH = 0; iH < tE.hitVect.size(); ++iH ){

      //Gather hit information
      double hitX = tE.hitVect[iH].endX_mm;
      double hitY = tE.hitVect[iH].endY_mm;
      double hitZ = tE.hitVect[iH].endZ_mm;
      double startX = tE.hitVect[iH].startX_mm;
      double startY = tE.hitVect[iH].startY_mm;
      double startZ = tE.hitVect[iH].startZ_mm;
      double logEnergy_eV = TMath::Log10(tE.hitVect[iH].eDep_eV);
      double eDep = tE.hitVect[iH].eDep_eV;
      double initKE = tE.hitVect[iH].startEnergy_eV;
      double endTime = tE.hitVect[iH].endT_ns;


      //Isolate above 2-Delta (aluminum)
      if (eDep > 0.00036){
        //Inequality to isolate Inductor 
        if (hitZ == 0.25 && hitY <= 0.8 && hitY >= 0.08 && hitX <= 0.565 && hitX >= -0.5655) {
          eDepInd += eDep;
        }
        // Inequality to isolate capacitor
        else if (hitZ == 0.25 && hitY >= 0.88 && hitY <= 1.31 && hitX <= 0.565 && hitX >= -0.5655) {
          eDepCapacitor += eDep;
        }
        // Inequality to isolate feedline
        else if (hitZ == 0.25 && hitY <= 0.05 && hitY >= -0.15) {
          eDepFeedline += eDep;
        }
        // Inequality to isolate mounts
        else if (hitZ == 0.25 && ((hitX <= 10 && hitX >= 9 && hitY <= 10 && hitY >= 9) || (hitX >= -10 && hitX <= -9 && hitY >= -10 && hitY <= -9))) {
          eDepMounts += eDep;
        }
        // Silicon
        else {
          eDepSilicon += eDep;
        }
      }

    }
    etaPhonons.push_back((100 * 0.6 * eDepInd)/0.5);
    etaPhononSilicon.push_back((100 * 0.6 * eDepSilicon)/0.5);
    etaPhononCapacitor.push_back((100 * 0.6 * eDepCapacitor)/0.5);
    etaPhononFeedline.push_back((100 * 0.6 * eDepFeedline)/0.5);
    etaPhononMounts.push_back((100 * 0.6 * eDepMounts)/0.5);

  }
  //calculate the average eta Phonon

  double etaPhononSum = 0.0;
  double etaPhononSiSum = 0.0;
  double etaPhononCapSum = 0.0;
  double etaPhononFeedSum = 0.0;
  double etaPhononMountSum = 0.0;

  for (int i = 0; i < etaPhonons.size(); i++) {
    etaPhononSum += etaPhonons.at(i);
    etaPhononSiSum += etaPhononSilicon.at(i);
    etaPhononCapSum += etaPhononCapacitor.at(i);
    etaPhononFeedSum += etaPhononFeedline.at(i);
    etaPhononMountSum += etaPhononMounts.at(i);
  }

  double etaPhononAvg = etaPhononSum/etaPhonons.size();
  double etaPhononSiAvg = etaPhononSiSum/etaPhononSilicon.size();
  double etaPhononCapAvg = etaPhononCapSum/etaPhononCapacitor.size();
  double etaPhononFeedAvg = etaPhononFeedSum/etaPhononFeedline.size();
  double etaPhononMountAvg = etaPhononMountSum/etaPhononMounts.size();

  double stdSum = 0.0;
  double stdSiSum = 0.0;
  double stdCapSum = 0.0;
  double stdFeedSum = 0.0;
  double stdMountSum = 0.0;

  for (int i = 0; i < etaPhonons.size(); i++) {
    stdSum += pow(etaPhonons.at(i) - etaPhononAvg, 2);
    stdSiSum += pow(etaPhononSilicon.at(i) - etaPhononSiAvg, 2);
    stdCapSum += pow(etaPhononCapacitor.at(i) - etaPhononCapAvg, 2);
    stdFeedSum += pow(etaPhononFeedline.at(i) - etaPhononFeedAvg, 2);
    stdMountSum += pow(etaPhononMounts.at(i) - etaPhononMountAvg, 2);
  }

  double stdEtaPhonon = sqrt(stdSum/etaPhonons.size());
  double stdEtaPhononSi = sqrt(stdSiSum/etaPhononSilicon.size());
  double stdEtaPhononCap = sqrt(stdCapSum/etaPhononCapacitor.size());
  double stdEtaPhononFeed = sqrt(stdFeedSum/etaPhononFeedline.size());
  double stdEtaPhononMounts = sqrt(stdMountSum/etaPhononMounts.size());

  std::string siProbStr = hitsFilename.substr(hitsFilename.find("_Si") + 3, hitsFilename.find("Mount") - (hitsFilename.find("_Si") + 3));
  std::string mountProbStr = hitsFilename.substr(hitsFilename.find("Mount") + 5, hitsFilename.find(".txt") - (hitsFilename.find("Mount") + 5));
  
  // Convert the extracted strings to doubles
  double siProb = std::stod(siProbStr);
  double mountProb = std::stod(mountProbStr);

  const std::string filename = "etaPhonons.csv";
  
  // Check if the file exists
  struct stat buffer;
  bool fileExists = (stat(filename.c_str(), &buffer) == 0);

  // Open file in append mode
  std::ofstream csvFile(filename, std::ios_base::app);

  // Write header only if the file does not exist
  if (!fileExists) {
    csvFile << "Detector Part,Silicon Absorption Probability,Mount Absorption Probability,Average eta Phonon,Standard Deviation\n";
  }

  // Write the data to the CSV file
  csvFile << "Inductor," << siProb << "," << mountProb << "," << etaPhononAvg << "," << stdEtaPhonon << "\n";
  csvFile << "Silicon," << siProb << "," << mountProb << "," << etaPhononSiAvg << "," << stdEtaPhononSi << "\n";
  csvFile << "Capacitor," << siProb << "," << mountProb << "," << etaPhononCapAvg << "," << stdEtaPhononCap << "\n";
  csvFile << "Feedline," << siProb << "," << mountProb << "," << etaPhononFeedAvg << "," << stdEtaPhononFeed << "\n";
  csvFile << "Mount," << siProb << "," << mountProb << "," << etaPhononMountAvg << "," << stdEtaPhononMounts << "\n";

  // Close the file
  csvFile.close();

  //Write stuff currently established
  fOut->Write();
}

  
  




  



//---------------------------------------------------------------------------------------
// Parsing function
std::vector<Event> ReadInG4CMPPrimaryAndHitFiles(std::string hitTextFilename, std::string primaryTextFilename)
{
  //Output: a vector of event objects nicely organizing our data.
  std::vector<Event> output;

  //First, let's open up our primary file and parse it. We'll get a map of int (event ID)
  //to primary info and a map of int (eventID) to a list of hits
  std::map<int,PrimaryInfo> primaryInfo = ParsePrimaryTextFileForPrimaries(primaryTextFilename);
  std::map<int,std::vector<Hit> > hitInfo = ParseHitTextFileForHits(hitTextFilename);

  //Now we do a final loop over event ID to merge these into actual events.
  for( std::map<int,PrimaryInfo>::iterator it = primaryInfo.begin(); it != primaryInfo.end(); ++it ){
    Event theEvent;
    theEvent.runID = it->second.runID;
    theEvent.eventID = it->first;
    theEvent.hitVect = hitInfo[it->first];
    theEvent.thePrim = it->second;
    output.push_back(theEvent);
  }
  return output;

}
  

//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename)
{
  std::map<int,std::vector<Hit> > output;
  std::vector<Hit> dummy;
  
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  int eventID = -1;
  int runID = -1;
  int counter = 0;
  while(1){
    if(!infile.good()) break;
    if(infile.is_open()){
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between spaces)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while(getline(check1,token,' ')){
	tokens.push_back(token);
      }
      if( tokens.size() == 0 ) break;
	  
      //If we're on the first line of the file
      if( tokens[0].find("Run") != std::string::npos ){
	continue;
      }

      //Check the runID and eventID, and if different than existing one,
      //push back a new event into the map
      if( std::atoi(tokens[0].c_str()) != runID || std::atoi(tokens[1].c_str()) != eventID ){
	output.emplace(std::atoi(tokens[1].c_str()),dummy);
	runID = std::atoi(tokens[0].c_str());
	eventID = std::atoi(tokens[1].c_str());
	counter++; 
	if( counter % 1000 == 0 ) std::cout << "Done reading " << counter << " events for hits." << std::endl;
      }

      //Log the hit information and push back into the most recently-created event in the vector
      Hit theHit;      
      theHit.runID = std::atoi(tokens[0].c_str());
      theHit.eventID = std::atoi(tokens[1].c_str());
      theHit.trackID = std::atoi(tokens[2].c_str());
      theHit.particleName = tokens[3];
      theHit.startEnergy_eV = std::atof(tokens[4].c_str());
      theHit.startX_mm = std::atof(tokens[5].c_str());
      theHit.startY_mm = std::atof(tokens[6].c_str());
      theHit.startZ_mm = std::atof(tokens[7].c_str());
      theHit.startT_ns = std::atof(tokens[8].c_str());
      theHit.eDep_eV = std::atof(tokens[9].c_str());
      theHit.trackWeight = std::atof(tokens[10].c_str());
      theHit.endX_mm = std::atof(tokens[11].c_str());
      theHit.endY_mm = std::atof(tokens[12].c_str());
      theHit.endZ_mm = std::atof(tokens[13].c_str());
      theHit.endT_ns = std::atof(tokens[14].c_str());
      output[eventID].push_back(theHit);
    }
  }
  return output;
}


//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename)
{
  std::map<int,PrimaryInfo> output;
  std::ifstream infile;
  infile.open(filename.c_str());
  std::string wholeLine;

  //Begin loop through file
  while(1){
    if(!infile.good()) break;
    if(infile.is_open()){
      std::getline(infile,wholeLine);
      
      //Tokenize the string (split between spaces)
      stringstream check1(wholeLine);
      string token;
      std::vector<std::string> tokens;
      while(getline(check1,token,' ')){
	tokens.push_back(token);
      }
      if( tokens.size() == 0 ) break;
	  
      //If we're on the first line of the file
      if( tokens[0].find("Run") != std::string::npos ){
	continue;
      }

      //Here, it's simpler since there's one line per event (assuming only a single run)
      //So we can use the event id as an index.
      PrimaryInfo thePrim;
      thePrim.runID = std::atoi(tokens[0].c_str());
      thePrim.eventID = std::atoi(tokens[1].c_str());
      thePrim.particleName = tokens[2];
      thePrim.energy_eV = std::atof(tokens[3].c_str());
      thePrim.X_mm = std::atof(tokens[4].c_str());
      thePrim.Y_mm = std::atof(tokens[5].c_str());
      thePrim.Z_mm = std::atof(tokens[6].c_str());
      thePrim.T_ns = std::atof(tokens[7].c_str());
      output.emplace(thePrim.eventID,thePrim);
    }
  }
  return output;
}



