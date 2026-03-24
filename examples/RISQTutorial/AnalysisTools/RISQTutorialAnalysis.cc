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
#include <map>
#include <sstream>
#include <string>
#include <vector>

//ROOT includes
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

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
/**
 * @brief Read and merge G4CMP primary-particle and hit information from two text files into a list of events.
 *
 * This routine is intended to parse two companion text outputs:
 *  - a "hit" text file containing detector/step hit records, and
 *  - a "primary" text file containing primary particle records,
 * and then associate the records by event (and any other identifiers used by the file formats) to produce
 * a unified collection of @c Event objects.
 *
 * @param hitTextFilename     Path to the text file containing hit records.
 * @param primaryTextFilename Path to the text file containing primary/vertex records.
 * @return std::vector<Event> A vector of reconstructed events containing merged primary + hit data.
 *
 * @throws std::runtime_error (or another exception) if either file cannot be opened, is malformed, or if the
 *         contents cannot be consistently associated into events (implementation-defined).
 *
 * @note From the snippet provided, this function is only *declared* in
 *       @c /Users/shilinray/Library/CloudStorage/OneDrive-CaliforniaInstituteofTechnology/Research/G4CMP/examples/RISQTutorial/AnalysisTools/RISQTutorialAnalysis.cc.
 *       Its *definition* is not shown; to locate it, search the codebase for
 *       @c "ReadInG4CMPPrimaryAndHitFiles(" (it may be defined later in the same .cc file or in another translation unit).
 */
std::vector<Event> ReadInG4CMPPrimaryAndHitFiles(std::string hitTextFilename, std::string primaryTextFilename);
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename);
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename);
int FindClosestQubitID(double hitX_mm, double hitY_mm);

  
void AnalyzeMuonEvent(std::string primariesFilename, std::string hitsFilename,double scaleFactorEHPairs)
{
  //Read in the text files
  std::vector<Event> eventList = ReadInG4CMPPrimaryAndHitFiles(hitsFilename,primariesFilename);
  std::cout << "Done reading in events." << std::endl;

  //Define an outfile
  TFile * fOut = new TFile("AnalysisOutput.root","RECREATE");
  
  //Define a number of histograms for the hits
  TH1F * h_eDep = new TH1F("h_eDep","Hit EDeps; log10(eDep[eV]); nEvents",200,-6,1);
  TH1F * h_qubitTotalHitEnergy_singleEvent = new TH1F("h_qubitTotalHitEnergy_singleEvent","In-Qubit Total Energy; Qubit ID; Energy [eV]",6,0,6);
  TH2F * h_hitXY = new TH2F("h_nHitsXY","XY Locations of Hits; X [mm]; Y [mm]; NHits/bin",200,-5,5,200,-5,5);
  TH2F * h_hitYZ = new TH2F("h_nHitsYZ","YZ Locations of Hits; Y [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_hitXZ = new TH2F("h_nHitsXZ","XZ Locations of Hits; X [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  
  //Loop over events
  for( int iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;
    
    //Get the event
    Event tE = eventList[iE];
    
    //Plot a number of hit-related things: 
    for ( int iH = 0; iH < tE.hitVect.size(); ++iH ){

      //Gather hit information
      double hitX_mm = tE.hitVect[iH].endX_mm;
      double hitY_mm = tE.hitVect[iH].endY_mm;
      double hitZ_mm = tE.hitVect[iH].endZ_mm;      
      double energy_eV = tE.hitVect[iH].eDep_eV;

      //From hit XY, figure out which qubit we're closest to, and add the
      //energy to that qubit id
      int closestQubitID = FindClosestQubitID(hitX_mm,hitY_mm);
      h_qubitTotalHitEnergy_singleEvent->Fill(closestQubitID,energy_eV);
      
      //Plot other hit information
      h_hitXY->Fill(hitX_mm,hitY_mm);
      h_hitYZ->Fill(hitY_mm,hitZ_mm);
      h_hitXZ->Fill(hitX_mm,hitZ_mm);
    }
  }

  //Last up, we need to remember that we did downsampling, so we need to scale our simulations back up to match
  h_qubitTotalHitEnergy_singleEvent->Scale(1.0/scaleFactorEHPairs);
  
  fOut->Write();
}

//---------------------------------------------------------------------------------------
// This is a bit of a kludge. What we really should do is pass the hit VOLUME out of
// G4CMP in the hit info. However, this involves some broader modifications to G4CMP
// beyond this example so we're kludging it together with this function.
int FindClosestQubitID(double hitX_mm, double hitY_mm)
{
  const int nQubits = 6;
  double qubitXY[nQubits][2] = {{-1.85,1.8},
				{0,1.8},
				{1.85,1.8},
				{-1.85,-1.8},
				{0,-1.8},
				{1.85,-1.8}};
  double lowestDeltaR2_mm = 1000000;
  int theQubitID = -1;
  for( int iQ = 0; iQ < nQubits; ++iQ ){
    double deltaR2_mm = TMath::Power((qubitXY[iQ][0] - hitX_mm),2) + TMath::Power((qubitXY[iQ][1] - hitY_mm),2);
    if( deltaR2_mm < lowestDeltaR2_mm ){
      lowestDeltaR2_mm = deltaR2_mm;
      theQubitID = iQ;
    }
  }
  return theQubitID;	 
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
  
  //Define a number of histograms for the hits
  TH1F * h_nHits = new TH1F("h_nHits","Number of Hits Per Event; nHits; nEvents",100,0,100);
  TH1F * h_eDep = new TH1F("h_eDep","Hit EDeps; log10(eDep[eV]); nEvents",200,-6,1);  
  TH2F * h_hitXY = new TH2F("h_nHitsXY","XY Locations of Hits; X [mm]; Y [mm]; NHits/bin",200,-5,5,200,-5,5);
  TH2F * h_hitYZ = new TH2F("h_nHitsYZ","YZ Locations of Hits; Y [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_hitXZ = new TH2F("h_nHitsXZ","XZ Locations of Hits; X [mm]; Z [mm]; NHits/bin",200,-5,5,200,4.6,5.02);
  
  //Define a number of histograms for the primaries
  TH2F * h_primXY = new TH2F("h_primXY","XY Location of Primary; X [mm]; Y [mm]; Primaries/bin",200,-5,5,200,-5,5);
  TH2F * h_primYZ = new TH2F("h_primYZ","YZ Location of Primary; Y [mm]; Z [mm]; Primaries/bin",200,-5,5,200,4.6,5.02);
  TH2F * h_primXZ = new TH2F("h_primXZ","XZ Location of Primary; X [mm]; Z [mm]; Primaries/bin",200,-5,5,200,4.6,5.02);

  //Define histograms useful for calculating PCE
  int nPCEBinsX = 50;
  int nPCEBinsY = 50;
  TH2F * h_totalHitEnergyAtPrimaryXY = new TH2F("h_totalHitEnergyAtPrimaryXY","Total Energy of Hits Generated by Primaries at this XY; X [mm]; Y [mm]; Energy/bin [eV]",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  TH2F * h_totalPrimaryEnergyAtPrimaryXY = new TH2F("h_totalPrimaryEnergyAtPrimaryXY","Total Energy of Primaries at this XY; X [mm]; Y [mm]; Energy/bin [eV]",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  TH2F * h_pceVsXY = new TH2F("h_pceVsXY","Phonon Collection Efficiency vs. Primary XY; X [mm]; Y [mm]; PCE",nPCEBinsX,-5,5,nPCEBinsY,-5,5);
  
  
  //Loop over events
  for( int iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;

    //Get the event
    Event tE = eventList[iE];
    
    //Add to the primary vector
    PrimaryInfo thePrim = tE.thePrim;
    h_primXY->Fill(thePrim.X_mm,thePrim.Y_mm);
    h_primXZ->Fill(thePrim.X_mm,thePrim.Z_mm);
    h_primYZ->Fill(thePrim.Y_mm,thePrim.Z_mm);
    h_totalPrimaryEnergyAtPrimaryXY->Fill(thePrim.X_mm,thePrim.Y_mm,thePrim.energy_eV);      
    
    //Plot a number of hit-related things: hit multiplicity, hit locations in XYZ, hits in XYZ weighted by energy, etc.
    h_nHits->Fill(tE.hitVect.size());
    for ( int iH = 0; iH < tE.hitVect.size(); ++iH ){

      //Gather hit information
      double hitX = tE.hitVect[iH].endX_mm;
      double hitY = tE.hitVect[iH].endY_mm;
      double hitZ = tE.hitVect[iH].endZ_mm;
      double startX = tE.hitVect[iH].startX_mm;
      double startY = tE.hitVect[iH].startY_mm;
      double startZ = tE.hitVect[iH].startZ_mm;
      double logEnergy_eV = TMath::Log10(tE.hitVect[iH].eDep_eV);


      
      //Plot hit information
      h_hitXY->Fill(hitX,hitY);
      h_hitYZ->Fill(hitY,hitZ);
      h_hitXZ->Fill(hitX,hitZ);
      h_eDep->Fill(logEnergy_eV);
      h_totalHitEnergyAtPrimaryXY->Fill(thePrim.X_mm,thePrim.Y_mm,tE.hitVect[iH].eDep_eV);
    }
  }

  //Post-processing division
  for( int iBX = 1; iBX <= h_totalHitEnergyAtPrimaryXY->GetNbinsX(); ++iBX ){
    for( int iBY = 1; iBY <= h_totalHitEnergyAtPrimaryXY->GetNbinsY(); ++iBY ){
      double num = h_totalHitEnergyAtPrimaryXY->GetBinContent(iBX,iBY);
      double denom = h_totalPrimaryEnergyAtPrimaryXY->GetBinContent(iBX,iBY);
      double pce = 0;
      if( denom != 0 ) pce = num/denom;
      h_pceVsXY->SetBinContent(iBX,iBY,pce);
    }
  }

  
  //Write stuff currently established
  fOut->Write();

}

void PrintPhononCollectionEfficiency()
{

  std::vector<int> al_vals = {100, 200, 300, 400, 500, 600};
  std::vector<int> nb_vals = {20, 40, 60, 80, 100};

  std::vector<std::string> hits_filenames;
  std::vector<std::string> primary_filenames;


  for (int al : al_vals) {
      for (int nb : nb_vals) {
          std::ostringstream ss;
          ss << "../G4Macros/data/Hits_Al" << al << "_Nb" << nb << ".txt";
          hits_filenames.push_back(ss.str());
      }
  }
  for (int al : al_vals) {
      for (int nb : nb_vals) {
          std::ostringstream ss;
          ss << "../G4Macros/data/Primary_Al" << al << "_Nb" << nb << ".txt";
          primary_filenames.push_back(ss.str());
      }
  }

  for (size_t i = 0; i < hits_filenames.size(); ++i) {

    const std::map<int,PrimaryInfo> primaryInfo = ParsePrimaryTextFileForPrimaries(primary_filenames[i]);
    const std::map<int,std::vector<Hit> > hitInfo = ParseHitTextFileForHits(hits_filenames[i]);

    double totalPrimaryEnergy_eV = 0.0;
    for (const auto& kv : primaryInfo) {
      totalPrimaryEnergy_eV += kv.second.energy_eV;
    }

    double totalHitEnergy_eV = 0.0;
    for (const auto& kv : hitInfo) {
      for (const Hit& h : kv.second) {
        totalHitEnergy_eV += h.eDep_eV;
      }
    }

    const double pce = (totalPrimaryEnergy_eV > 0.0) ? (totalHitEnergy_eV / totalPrimaryEnergy_eV) : 0.0;

    std::cout << "\n=== Phonon Collection Efficiency (PCE) ===\n"
              << "Primaries file: " << primary_filenames[i] << "\n"
              << "Hits file:      " << hits_filenames[i] << "\n"
              << "Total primary phonon energy: " << totalPrimaryEnergy_eV << " eV\n"
              << "Total hit phonon energy:     " << totalHitEnergy_eV << " eV\n"
              << "PCE = " << (100.0 * pce) << " %\n"
              << "========================================\n";
  }
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void PrintPhononCollectionEfficiencyAndPlot()
{
  std::vector<int> al_vals = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  std::vector<int> nb_vals = {20, 40, 60, 80, 100};

  const std::string baseDir = "../G4Macros/260324_run/lQPD_true";

  TFile* fOut = new TFile("PCE_Al_Nb.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file PCE_Al_Nb.root" << std::endl;
    return;
  }
  fOut->cd();

  // 2D PCE map
  TH2F* h_pce_al_nb = new TH2F("h_pce_al_nb",
                               "PCE vs Al and Nb thickness;Al thickness [nm];Nb thickness [nm];PCE [%]",
                               (int)al_vals.size(), 0, (int)al_vals.size(),
                               (int)nb_vals.size(), 0, (int)nb_vals.size());

  for (int i = 0; i < (int)al_vals.size(); ++i) {
    h_pce_al_nb->GetXaxis()->SetBinLabel(i + 1, std::to_string(al_vals[i]).c_str());
  }
  for (int j = 0; j < (int)nb_vals.size(); ++j) {
    h_pce_al_nb->GetYaxis()->SetBinLabel(j + 1, std::to_string(nb_vals[j]).c_str());
  }

  // Group histograms two ways:
  // 1) by Al (for fixed Al, varying Nb)
  // 2) by Nb (for fixed Nb, varying Al)
  std::vector<std::vector<TH1F*> > hists_by_al(al_vals.size());
  std::vector<std::vector<TH1F*> > hists_by_nb(nb_vals.size());

  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    for (int iN = 0; iN < (int)nb_vals.size(); ++iN) {

      const int al = al_vals[iA];
      const int nb = nb_vals[iN];

      std::ostringstream primPath, hitPath;
      primPath << baseDir << "/Primary_Al" << al << "_Nb" << nb << ".txt";
      hitPath  << baseDir << "/Hits_Al"    << al << "_Nb" << nb << ".txt";

      const std::string primariesFilename = primPath.str();
      const std::string hitsFilename      = hitPath.str();

      std::ifstream primFile(primariesFilename.c_str());
      std::ifstream hitFile(hitsFilename.c_str());

      if (!primFile.good() || !hitFile.good()) {
        std::cout << "Skipping Al=" << al << " Nb=" << nb
                  << " because one or both files are missing." << std::endl;
        continue;
      }

      primFile.close();
      hitFile.close();

      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primariesFilename);

      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);

      TString hName = TString::Format("h_eDep_Al%d_Nb%d", al, nb);
      TString hTitle = TString::Format("Hit EDeps for Al=%d nm, Nb=%d nm;log10(eDep[eV]);nEvents",
                                       al, nb);

      TH1F* h_eDep = new TH1F(hName, hTitle, 200, -6, 1);
      h_eDep->SetDirectory(fOut);
      h_eDep->SetLineWidth(2);
      h_eDep->SetStats(0);

      // Color by Nb for the fixed-Al overlays
      int color_nb = kBlack;
      if (iN == 0) color_nb = kBlack;
      else if (iN == 1) color_nb = kRed;
      else if (iN == 2) color_nb = kBlue;
      else if (iN == 3) color_nb = kGreen + 2;
      else if (iN == 4) color_nb = kMagenta + 1;

      h_eDep->SetLineColorAlpha(color_nb, 0.6);

      double totalPrimaryEnergy_eV = 0.0;
      for (const auto& kv : primaryInfo) {
        totalPrimaryEnergy_eV += kv.second.energy_eV;
      }

      double totalHitEnergy_eV = 0.0;
      for (const auto& kv : hitInfo) {
        for (const Hit& h : kv.second) {
          totalHitEnergy_eV += h.eDep_eV;
          if (h.eDep_eV > 0.0) {
            h_eDep->Fill(TMath::Log10(h.eDep_eV));
          }
        }
      }

      const double pce = (totalPrimaryEnergy_eV > 0.0)
                           ? (totalHitEnergy_eV / totalPrimaryEnergy_eV)
                           : 0.0;

      h_pce_al_nb->SetBinContent(iA + 1, iN + 1, 100.0 * pce);

      std::cout << "Al=" << al
                << " Nb=" << nb
                << "  Total primary energy=" << totalPrimaryEnergy_eV << " eV"
                << "  Total hit energy=" << totalHitEnergy_eV << " eV"
                << "  PCE=" << 100.0 * pce << " %" << std::endl;

      fOut->cd();
      h_eDep->Write();

      hists_by_al[iA].push_back(h_eDep);
      hists_by_nb[iN].push_back(h_eDep);
    }
  }

  // Write the 2D PCE map
  fOut->cd();
  h_pce_al_nb->Write();

  // Save the PCE map as a png
  TCanvas* c_pce = new TCanvas("c_pce_al_nb", "c_pce_al_nb", 900, 700);
  c_pce->SetRightMargin(0.18);
  h_pce_al_nb->SetStats(0);
  h_pce_al_nb->Draw("COLZ TEXT");
  c_pce->Write();
  c_pce->SaveAs("PCE_Al_Nb.png");

  // --------------------------------------------------------------------------
  // Fixed Al, varying Nb
  // --------------------------------------------------------------------------
  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    const int al_val = al_vals[iA];

    if (hists_by_al[iA].empty()) continue;

    TString cName  = TString::Format("c_eDep_overlay_Al%d", al_val);
    TString cTitle = TString::Format("Hit EDeps for Al=%d nm, varying Nb", al_val);
    TCanvas* c = new TCanvas(cName, cTitle, 900, 700);

    TLegend* leg = new TLegend(0.68, 0.68, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);

    double maxY = 0.0;
    for (TH1F* h : hists_by_al[iA]) {
      if (h->GetMaximum() > maxY) maxY = h->GetMaximum();
    }

    bool first = true;
    for (TH1F* h : hists_by_al[iA]) {
      h->SetTitle(TString::Format("Hit EDeps for Al=%d nm, varying Nb;log10(eDep[eV]);nEvents",
                                  al_val));
      h->SetMaximum(1.15 * maxY);

      TString histName = h->GetName();
      int parsed_al = -1;
      int nb_val = -1;
      sscanf(histName.Data(), "h_eDep_Al%d_Nb%d", &parsed_al, &nb_val);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }

      leg->AddEntry(h, TString::Format("Nb = %d nm", nb_val), "l");
    }

    leg->Draw();

    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("h_eDep_overlay_Al%d.png", al_val));

    delete leg;
    delete c;
  }

  // --------------------------------------------------------------------------
  // Fixed Nb, varying Al
  // --------------------------------------------------------------------------
  for (int iN = 0; iN < (int)nb_vals.size(); ++iN) {
    const int nb_val = nb_vals[iN];

    if (hists_by_nb[iN].empty()) continue;

    TString cName  = TString::Format("c_eDep_overlay_Nb%d", nb_val);
    TString cTitle = TString::Format("Hit EDeps for Nb=%d nm, varying Al", nb_val);
    TCanvas* c = new TCanvas(cName, cTitle, 900, 700);

    TLegend* leg = new TLegend(0.68, 0.58, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);

    double maxY = 0.0;
    for (TH1F* h : hists_by_nb[iN]) {
      if (h->GetMaximum() > maxY) maxY = h->GetMaximum();
    }

    bool first = true;
    int color_index = 0;

    for (TH1F* h : hists_by_nb[iN]) {
      h->SetTitle(TString::Format("Hit EDeps for Nb=%d nm, varying Al;log10(eDep[eV]);nEvents",
                                  nb_val));
      h->SetMaximum(1.15 * maxY);

      // Recolor by Al for these plots so each Al is visually distinct
      int color_al = kBlack;
      if (color_index == 0) color_al = kBlack;
      else if (color_index == 1) color_al = kRed;
      else if (color_index == 2) color_al = kBlue;
      else if (color_index == 3) color_al = kGreen + 2;
      else if (color_index == 4) color_al = kMagenta + 1;
      else if (color_index == 5) color_al = kOrange + 7;
      else if (color_index == 6) color_al = kCyan + 2;
      else if (color_index == 7) color_al = kViolet;
      else if (color_index == 8) color_al = kPink + 7;
      else if (color_index == 9) color_al = kAzure + 2;

      h->SetLineColorAlpha(color_al, 0.6);

      TString histName = h->GetName();
      int al_val = -1;
      int parsed_nb = -1;
      sscanf(histName.Data(), "h_eDep_Al%d_Nb%d", &al_val, &parsed_nb);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }

      leg->AddEntry(h, TString::Format("Al = %d nm", al_val), "l");
      ++color_index;
    }

    leg->Draw();

    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("h_eDep_overlay_Nb%d.png", nb_val));

    delete leg;
    delete c;
  }

  fOut->Write();
  fOut->Close();

  delete c_pce;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////  

//---------------------------------------------------------------------------------------
// Parsing function
//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,std::vector<Hit> > ParseHitTextFileForHits(std::string filename)
{
  //std::cout << "About to parse hits:    " << hitsFilename << std::endl;
  std::map<int,std::vector<Hit> > output;

  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cerr << "Error: could not open hit file: " << filename << std::endl;
    return output;
  }

  std::string wholeLine;
  int eventID = -1;
  int runID = -1;
  int counter = 0;

  while (std::getline(infile, wholeLine)) {
    // Skip empty / whitespace-only lines
    if (wholeLine.find_first_not_of(" \t\r\n") == std::string::npos) {
      continue;
    }

    // Skip header line(s)
    if (wholeLine.find("Run") != std::string::npos) {
      continue;
    }

    std::istringstream iss(wholeLine);
    Hit theHit;

    if (!(iss >> theHit.runID
              >> theHit.eventID
              >> theHit.trackID
              >> theHit.particleName
              >> theHit.startEnergy_eV
              >> theHit.startX_mm
              >> theHit.startY_mm
              >> theHit.startZ_mm
              >> theHit.startT_ns
              >> theHit.eDep_eV
              >> theHit.trackWeight
              >> theHit.endX_mm
              >> theHit.endY_mm
              >> theHit.endZ_mm
              >> theHit.endT_ns)) {
      std::cerr << "Warning: skipping malformed line in hit file: "
                << filename << std::endl
                << "  Line: " << wholeLine << std::endl;
      continue;
    }

    if (theHit.runID != runID || theHit.eventID != eventID) {
      output.emplace(theHit.eventID, std::vector<Hit>{});
      runID = theHit.runID;
      eventID = theHit.eventID;
      counter++;

      if (counter % 1000 == 0) {
        std::cout << "Done reading " << counter << " events for hits." << std::endl;
      }
    }

    output[eventID].push_back(theHit);
  }

  return output;
}


//---------------------------------------------------------------------------------------
// Parsing function
std::map<int,PrimaryInfo> ParsePrimaryTextFileForPrimaries(std::string filename)
{
  //std::cout << "About to parse primary: " << primariesFilename << std::endl;
  std::map<int,PrimaryInfo> output;

  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cerr << "Error: could not open primary file: " << filename << std::endl;
    return output;
  }

  std::string wholeLine;

  while (std::getline(infile, wholeLine)) {
    // Skip empty / whitespace-only lines
    if (wholeLine.find_first_not_of(" \t\r\n") == std::string::npos) {
      continue;
    }

    // Skip header line(s)
    if (wholeLine.find("Run") != std::string::npos) {
      continue;
    }

    std::istringstream iss(wholeLine);
    PrimaryInfo thePrim;

    // Your current primary text format appears to be:
    // runID eventID particleName energy X Y Z T
    //
    // So we do NOT read trackID here.
    if (!(iss >> thePrim.runID
              >> thePrim.eventID
              >> thePrim.particleName
              >> thePrim.energy_eV
              >> thePrim.X_mm
              >> thePrim.Y_mm
              >> thePrim.Z_mm
              >> thePrim.T_ns)) {
      std::cerr << "Warning: skipping malformed line in primary file: "
                << filename << std::endl
                << "  Line: " << wholeLine << std::endl;
      continue;
    }

    // This field exists in the struct but is not present in the file format.
    thePrim.trackID = -1;

    output[thePrim.eventID] = thePrim;
  }

  return output;
}



