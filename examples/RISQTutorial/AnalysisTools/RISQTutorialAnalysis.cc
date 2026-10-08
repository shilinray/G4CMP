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
#include <algorithm>
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
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"
#include "TRandom3.h"

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
  for( std::vector<Event>::size_type iE = 0; iE < eventList.size(); ++iE ){
    if( iE % 1000 == 0 ) std::cout << "Done with " << iE << " event histogram fills." << std::endl;
    
    //Get the event
    Event tE = eventList[iE];
    
    //Plot a number of hit-related things: 
    for ( std::vector<Hit>::size_type iH = 0; iH < tE.hitVect.size(); ++iH ){

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
  for( std::vector<Event>::size_type iE = 0; iE < eventList.size(); ++iE ){
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
    for ( std::vector<Hit>::size_type iH = 0; iH < tE.hitVect.size(); ++iH ){

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

  const std::string baseDir = "../G4Macros/260410_run/lQPD_Al_PF";
  // Create the ROOT output file where all histograms and canvases will be saved.
  TFile* fOut = new TFile("PCE_Al_Nb.root", "RECREATE");
  // Check that the file was created successfully.
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file PCE_Al_Nb.root" << std::endl;
    return;
  }
  // Make the output file the current ROOT directory so newly written objects go there.
  fOut->cd();

  // Create 2D PCE map
  TH2F* h_pce_al_nb = new TH2F("h_pce_al_nb",
                               "PCE vs Al and Nb thickness;Al thickness [nm];Nb thickness [nm];PCE [%]",
                               (int)al_vals.size(), 0, (int)al_vals.size(),
                                (int)nb_vals.size(), 0, (int)nb_vals.size());
  // labeling the x axis
  for (int i = 0; i < (int)al_vals.size(); ++i) {
    h_pce_al_nb->GetXaxis()->SetBinLabel(i + 1, std::to_string(al_vals[i]).c_str());
  }
  // labeling the y axis
  for (int j = 0; j < (int)nb_vals.size(); ++j) {
    h_pce_al_nb->GetYaxis()->SetBinLabel(j + 1, std::to_string(nb_vals[j]).c_str());
  }

  // Group histograms two ways:
  // hists_by_al[iA] = all histograms for one fixed Al value while Nb varies
  // hists_by_nb[iN] = all histograms for one fixed Nb value while Al varies
  std::vector<std::vector<TH1F*> > hists_by_al(al_vals.size());
  std::vector<std::vector<TH1F*> > hists_by_nb(nb_vals.size());

  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    for (int iN = 0; iN < (int)nb_vals.size(); ++iN) {

      const int al = al_vals[iA];
      const int nb = nb_vals[iN];
      // Construct the expected input filenames for this Al/Nb combination
      std::ostringstream primPath, hitPath;
      primPath << baseDir << "/Primary_Al" << al << "_Nb" << nb << ".txt";
      hitPath  << baseDir << "/Hits_Al"    << al << "_Nb" << nb << ".txt";

      const std::string primariesFilename = primPath.str();
      const std::string hitsFilename      = hitPath.str();

      std::ifstream primFile(primariesFilename.c_str());
      std::ifstream hitFile(hitsFilename.c_str());

      //if (!primFile.good() || !hitFile.good()) {
        //std::cout << "Skipping Al=" << al << " Nb=" << nb
                  //<< " because one or both files are missing." << std::endl;
        //continue;
      //}

      primFile.close();
      hitFile.close();
      // map is a sorted dictionary, key is event ID, value is vector with all the info
      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primariesFilename);

      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);

      TString hName = TString::Format("h_eDep_Al%d_Nb%d", al, nb);
      TString hTitle = TString::Format("Hit EDeps for Al=%d nm, Nb=%d nm;log10(eDep[eV]);nEvents",
                                       al, nb);
      // creates 1D histogram
      // 200 bins spanning 10^-6 eV to 10^1 eV in log10 space
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

      // Sum the total injected/primary energy over all events.
      double totalPrimaryEnergy_eV = 0.0;
      for (const auto& kv : primaryInfo) {
        totalPrimaryEnergy_eV += kv.second.energy_eV;
      }
      // Sum the total deposited hit energy over all events.
      // Also fill the histogram with log10(eDep_eV) for each positive-energy hit.
      double totalHitEnergy_eV = 0.0;
      for (const auto& kv : hitInfo) {
        //std::cout << "Event " << kv.first
          //<< " has " << kv.second.size()
          //<< " hits." << std::endl;
        for (const Hit& h : kv.second) {
          totalHitEnergy_eV += h.eDep_eV;
          if (h.eDep_eV > 0.0) {
            h_eDep->Fill(TMath::Log10(h.eDep_eV)); 
          }
          if (h.eDep_eV < 0.0) {
            std::cout << "There are negative energies." << std::endl;
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
      // Save pointers to the histogram in both grouping schemes
      // so we can make overlay plots later.
      hists_by_al[iA].push_back(h_eDep);
      hists_by_nb[iN].push_back(h_eDep);
    }
  }

  // Write the 2D PCE map to root file
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
  // For each Al thickness, overlay all corresponding histograms for different Nb values.
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


    TLine* line_Al1 = new TLine(TMath::Log10(0.34e-3), 0, TMath::Log10(0.34e-3), 1.15 * maxY);
    TLine* line_Al2 = new TLine(TMath::Log10(0.68e-3), 0, TMath::Log10(0.68e-3), 1.15 * maxY);
    TLine* line_Al3 = new TLine(TMath::Log10(1.02e-3), 0, TMath::Log10(1.02e-3), 1.15 * maxY);
    line_Al1->SetLineColor(kGreen+2);
    line_Al2->SetLineColor(kGreen+2);
    line_Al3->SetLineColor(kGreen+2);
    line_Al1->SetLineStyle(2);
    line_Al2->SetLineStyle(2);
    line_Al3->SetLineStyle(2);
    line_Al1->SetLineWidth(2);
    line_Al2->SetLineWidth(2);
    line_Al3->SetLineWidth(2);

    line_Al1->Draw();
    line_Al2->Draw();
    line_Al3->Draw();

    leg->AddEntry(line_Al1, "Al gap", "l");

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


    TLine* line_Al1 = new TLine(TMath::Log10(0.34e-3), 0, TMath::Log10(0.34e-3), 1.15 * maxY);
    TLine* line_Al2 = new TLine(TMath::Log10(0.68e-3), 0, TMath::Log10(0.68e-3), 1.15 * maxY);
    TLine* line_Al3 = new TLine(TMath::Log10(1.02e-3), 0, TMath::Log10(1.02e-3), 1.15 * maxY);
    line_Al1->SetLineColor(kGreen+2);
    line_Al2->SetLineColor(kGreen+2);
    line_Al3->SetLineColor(kGreen+2);
    line_Al1->SetLineStyle(2);
    line_Al2->SetLineStyle(2);
    line_Al3->SetLineStyle(2);
    line_Al1->SetLineWidth(2);
    line_Al2->SetLineWidth(2);
    line_Al3->SetLineWidth(2);

    line_Al1->Draw();
    line_Al2->Draw();
    line_Al3->Draw();


    leg->AddEntry(line_Al1, "Al gap", "l");

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
 
void PlotCumulativeHitEnergyVsEndTime_OverlayNb20()
{
  std::vector<int> al_vals = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  const int nb_val = 20;

  const std::string baseDir = "../G4Macros/260409_run/SQUAT_Al_Nb_PF";

  TFile* fOut = new TFile("CumulativeHitEnergyVsEndTime_Overlay_Nb20.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file." << std::endl;
    return;
  }

  std::vector<TH1F*> hists;
  double globalTMin = 0.0;
  double globalTMax = 0.0;
  bool firstTime = true;

  // --------------------------------------------------------------------------
  // First pass: find the global time range across all Al values
  // --------------------------------------------------------------------------
  for (int al : al_vals) {
    std::ostringstream hitPath;
    hitPath << baseDir << "/Hits_Al" << al << "_Nb" << nb_val << ".txt";

    const std::map<int, std::vector<Hit> > hitInfo = ParseHitTextFileForHits(hitPath.str());
    if (hitInfo.empty()) {
      std::cout << "Skipping empty or unreadable file: " << hitPath.str() << std::endl;
      continue;
    }

    for (const auto& kv : hitInfo) {
      for (const Hit& h : kv.second) {
        if (firstTime) {
          globalTMin = h.endT_ns;
          globalTMax = h.endT_ns;
          firstTime = false;
        } else {
          if (h.endT_ns < globalTMin) globalTMin = h.endT_ns;
          if (h.endT_ns > globalTMax) globalTMax = h.endT_ns;
        }
      }
    }
  }

  if (firstTime) {
    std::cerr << "Error: no valid hit files found for Nb = " << nb_val << std::endl;
    fOut->Close();
    delete fOut;
    return;
  }

  if (globalTMin == globalTMax) {
    globalTMin -= 0.5;
    globalTMax += 0.5;
  }

  const int nTimeBins = 200;
  double maxY = 0.0;

  // --------------------------------------------------------------------------
  // Second pass: build cumulative histograms
  // --------------------------------------------------------------------------
  for (size_t i = 0; i < al_vals.size(); ++i) {
    const int al = al_vals[i];

    std::ostringstream hitPath;
    hitPath << baseDir << "/Hits_Al" << al << "_Nb" << nb_val << ".txt";

    const std::map<int, std::vector<Hit> > hitInfo = ParseHitTextFileForHits(hitPath.str());
    if (hitInfo.empty()) continue;

    TString hName = TString::Format("h_cumulativeEnergyVsEndTime_Al%d_Nb%d", al, nb_val);
    TString hTitle = TString::Format("Cumulative Hit Energy vs End Time, Nb=%d nm;End Time [ns];Cumulative Deposited Energy [eV]",
                                     nb_val);

    TH1F* h_temp = new TH1F(hName, hTitle, nTimeBins, globalTMin, globalTMax);
    h_temp->SetDirectory(fOut);
    h_temp->SetStats(0);
    h_temp->SetLineWidth(2);

    // Give each Al curve a different color
    int color = kBlack;
    if      (i == 0) color = kBlack;
    else if (i == 1) color = kRed;
    else if (i == 2) color = kBlue;
    else if (i == 3) color = kGreen + 2;
    else if (i == 4) color = kMagenta + 1;
    else if (i == 5) color = kOrange + 7;
    else if (i == 6) color = kCyan + 2;
    else if (i == 7) color = kViolet;
    else if (i == 8) color = kPink + 7;
    else if (i == 9) color = kAzure + 2;

    h_temp->SetLineColor(color);

    // Fill deposited energy vs time
    for (const auto& kv : hitInfo) {
      for (const Hit& h : kv.second) {
        h_temp->Fill(h.endT_ns, h.eDep_eV);
      }
    }

    // Convert to cumulative
    double runningSum = 0.0;
    for (int iBin = 1; iBin <= h_temp->GetNbinsX(); ++iBin) {
      runningSum += h_temp->GetBinContent(iBin);
      h_temp->SetBinContent(iBin, runningSum);
    }

    if (h_temp->GetMaximum() > maxY) maxY = h_temp->GetMaximum();

    hists.push_back(h_temp);
    h_temp->Write();
  }

  if (hists.empty()) {
    std::cerr << "Error: no histograms were created." << std::endl;
    fOut->Close();
    delete fOut;
    return;
  }

  // --------------------------------------------------------------------------
  // Draw overlay
  // --------------------------------------------------------------------------
  TCanvas* c = new TCanvas("c_cumulativeOverlay_Nb20",
                           "Cumulative Hit Energy vs End Time, Nb=20",
                           1000, 700);

  TLegend* leg = new TLegend(0.62, 0.58, 0.88, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);

  bool firstDraw = true;
  for (TH1F* h : hists) {
    h->SetMaximum(1.08 * maxY);

    TString histName = h->GetName();
    int al_val = -1;
    int parsed_nb = -1;
    sscanf(histName.Data(), "h_cumulativeEnergyVsEndTime_Al%d_Nb%d", &al_val, &parsed_nb);

    if (firstDraw) {
      h->Draw("HIST");
      firstDraw = false;
    } else {
      h->Draw("HIST SAME");
    }

    leg->AddEntry(h, TString::Format("Al = %d nm", al_val), "l");
  }

  leg->Draw();

  fOut->cd();
  c->Write();
  c->SaveAs("CumulativeHitEnergyVsEndTime_Overlay_Nb20.png");

  fOut->Write();
  fOut->Close();

  delete leg;
  delete c;
  delete fOut;
}

void Ns_PrintPhononCollectionEfficiencyAndPlot()
{
  std::vector<int> al_vals = {100, 300, 500, 700, 900};
  std::vector<int> ns_vals = {1, 5, 10, 20, 40, 80, 100};

  const std::string baseDir = "../G4Macros/260419_run/SQUAT_Al_Nb_PF";

  // Al superconducting gap energy [eV]. Used to convert deposited energy to
  // total QPs: N_QP = 2 * E_dep / Delta_Al
  const double alGap_eV = 0.34e-3;

  // Create the ROOT output file where all histograms and canvases will be saved.
  TFile* fOut = new TFile("PCE_Al_ns.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file PCE_Al_ns.root" << std::endl;
    return;
  }
  fOut->cd();

  // Create 2D PCE map
  TH2F* h_pce_al_ns = new TH2F("h_pce_al_ns",
                               "PCE vs Al thickness and chip area;Al thickness [nm];Chip area [cm^{2}];PCE [%]",
                               (int)al_vals.size(), 0, (int)al_vals.size(),
                               (int)ns_vals.size(), 0, (int)ns_vals.size());

  // Label x axis with Al values
  for (int i = 0; i < (int)al_vals.size(); ++i) {
    h_pce_al_ns->GetXaxis()->SetBinLabel(i + 1, std::to_string(al_vals[i]).c_str());
  }

  // Label y axis with chip area values, where area = 1 cm^2 / ns
  for (int j = 0; j < (int)ns_vals.size(); ++j) {
    const double chipArea_cm2 = 1.0 / ns_vals[j];
    TString areaLabel = TString::Format("%.3f", chipArea_cm2);
    h_pce_al_ns->GetYaxis()->SetBinLabel(j + 1, areaLabel.Data());
  }

  // Create 2D total QP map: N_QP = 2 * totalHitEnergy / Delta_Al
  TH2F* h_nqp_al_ns = new TH2F("h_nqp_al_ns",
                               "Total QPs created vs Al thickness and chip area;Al thickness [nm];Chip area [cm^{2}];Total N_{QP}",
                               (int)al_vals.size(), 0, (int)al_vals.size(),
                               (int)ns_vals.size(), 0, (int)ns_vals.size());

  for (int i = 0; i < (int)al_vals.size(); ++i) {
    h_nqp_al_ns->GetXaxis()->SetBinLabel(i + 1, std::to_string(al_vals[i]).c_str());
  }
  for (int j = 0; j < (int)ns_vals.size(); ++j) {
    const double chipArea_cm2 = 1.0 / ns_vals[j];
    TString areaLabel = TString::Format("%.3f", chipArea_cm2);
    h_nqp_al_ns->GetYaxis()->SetBinLabel(j + 1, areaLabel.Data());
  }

  // Group histograms two ways:
  // hists_by_al[iA] = all histograms for one fixed Al value while ns varies
  // hists_by_ns[iS] = all histograms for one fixed ns value while Al varies
  std::vector<std::vector<TH1F*> > hists_by_al(al_vals.size());
  std::vector<std::vector<TH1F*> > hists_by_ns(ns_vals.size());

  // PCE [%] and binomial error tables, indexed [iA][iS], filled in the inner loop
  // and used after the loop to build the PCE vs Al line graph.
  // Binomial error: sigma_pce = 100 * sqrt(p*(1-p)/n_phonons).
  // n_phonons_pce is set dynamically from the number of events in Primary_Al100_ns1.txt.
  const std::string refPrimPath = baseDir + "/Primary_Al100_ns1.txt";
  const std::map<int, PrimaryInfo> refPrimInfo = ParsePrimaryTextFileForPrimaries(refPrimPath);
  if (refPrimInfo.empty()) {
    std::cerr << "Error: could not read primaries file for n_phonons_pce: " << refPrimPath << std::endl;
    fOut->Close();
    delete fOut;
    return;
  }
  const int n_phonons_pce = (int)refPrimInfo.size();
  std::cout << "n_phonons_pce = " << n_phonons_pce << " (from " << refPrimPath << ")" << std::endl;
  std::vector<std::vector<double> > pce_table(al_vals.size(),
                                              std::vector<double>(ns_vals.size(), 0.0));
  std::vector<std::vector<double> > pce_err_table(al_vals.size(),
                                                  std::vector<double>(ns_vals.size(), 0.0));

  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {

      const int al = al_vals[iA];
      const int ns = ns_vals[iS];

      // Construct the expected input filenames for this Al/ns combination
      std::ostringstream primPath, hitPath;
      primPath << baseDir << "/Primary_Al" << al << "_ns" << ns << ".txt";
      hitPath  << baseDir << "/Hits_Al"    << al << "_ns" << ns << ".txt";

      const std::string primariesFilename = primPath.str();
      const std::string hitsFilename      = hitPath.str();

      std::ifstream primFile(primariesFilename.c_str());
      std::ifstream hitFile(hitsFilename.c_str());

      // Uncomment if you want to skip missing files
      /*
      if (!primFile.good() || !hitFile.good()) {
        std::cout << "Skipping Al=" << al << " ns=" << ns
                  << " because one or both files are missing." << std::endl;
        primFile.close();
        hitFile.close();
        continue;
      }
      */

      primFile.close();
      hitFile.close();

      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primariesFilename);

      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);

      TString hName = TString::Format("h_eDep_Al%d_ns%d", al, ns);
      TString hTitle = TString::Format("Hit EDeps for Al=%d nm, ns=%d;eDep [eV];nEvents",
                                       al, ns);

      TH1F* h_eDep = new TH1F(hName, hTitle, 200, 0, 5e-3);
      h_eDep->SetDirectory(fOut);
      h_eDep->SetLineWidth(2);
      h_eDep->SetStats(0);

      // Color by ns for the fixed-Al overlays
      int color_ns = kBlack;
      if (iS == 0)  color_ns = kBlack;
      else if (iS == 1)  color_ns = kRed;
      else if (iS == 2)  color_ns = kBlue;
      else if (iS == 3)  color_ns = kGreen + 2;
      else if (iS == 4)  color_ns = kMagenta + 1;
      else if (iS == 5)  color_ns = kOrange + 7;
      else if (iS == 6)  color_ns = kCyan + 2;
      else if (iS == 7)  color_ns = kViolet;
      else if (iS == 8)  color_ns = kPink + 7;
      else if (iS == 9)  color_ns = kAzure + 2;
      else if (iS == 10) color_ns = kSpring + 5;

      h_eDep->SetLineColorAlpha(color_ns, 0.6);

      // Sum total primary energy
      double totalPrimaryEnergy_eV = 0.0;
      for (const auto& kv : primaryInfo) {
        totalPrimaryEnergy_eV += kv.second.energy_eV;
      }

      // Sum total hit energy and fill histogram
      double totalHitEnergy_eV = 0.0;
      for (const auto& kv : hitInfo) {
        for (const Hit& h : kv.second) {
          totalHitEnergy_eV += h.eDep_eV;

          if (h.eDep_eV > 0.0) {
            h_eDep->Fill(h.eDep_eV);
          }

          if (h.eDep_eV < 0.0) {
            std::cout << "There are negative energies." << std::endl;
          }
        }
      }

      const double pce = (totalPrimaryEnergy_eV > 0.0)
                           ? (totalHitEnergy_eV / totalPrimaryEnergy_eV)
                           : 0.0;

      // Total QPs = 2 * E_dep / Delta_Al (factor of 2: each broken Cooper pair yields 2 QPs)
      const double totalQP = 2.0 * totalHitEnergy_eV / alGap_eV;

      h_pce_al_ns->SetBinContent(iA + 1, iS + 1, 100.0 * pce);
      h_nqp_al_ns->SetBinContent(iA + 1, iS + 1, totalQP);

      // Store PCE and its binomial uncertainty for the line graph below
      pce_table[iA][iS] = 100.0 * pce;
      const double pce_err = (pce > 0.0 && pce < 1.0)
                               ? 100.0 * std::sqrt(pce * (1.0 - pce) / n_phonons_pce)
                               : 0.0;
      pce_err_table[iA][iS] = pce_err;

      std::cout << "Al=" << al
                << " ns=" << ns
                << "  Total primary energy=" << totalPrimaryEnergy_eV << " eV"
                << "  Total hit energy=" << totalHitEnergy_eV << " eV"
                << "  PCE=" << 100.0 * pce << " %"
                << "  Total QP=" << totalQP << std::endl;

      fOut->cd();
      h_eDep->Write();

      hists_by_al[iA].push_back(h_eDep);
      hists_by_ns[iS].push_back(h_eDep);
    }
  }

  // Write the 2D PCE map
  fOut->cd();
  h_pce_al_ns->Write();

  // Save the PCE map as a png
  TCanvas* c_pce = new TCanvas("c_pce_al_ns", "c_pce_al_ns", 900, 700);
  c_pce->SetRightMargin(0.18);
  h_pce_al_ns->SetStats(0);
  h_pce_al_ns->Draw("COLZ TEXT");
  c_pce->Write();
  c_pce->SaveAs("PCE_Al_ns.png");

  // Write the 2D total QP map
  fOut->cd();
  h_nqp_al_ns->Write();

  // Save the total QP map as a png
  TCanvas* c_nqp = new TCanvas("c_nqp_al_ns", "c_nqp_al_ns", 900, 700);
  c_nqp->SetRightMargin(0.18);
  h_nqp_al_ns->SetStats(0);
  h_nqp_al_ns->Draw("COLZ TEXT");
  c_nqp->Write();
  c_nqp->SaveAs("NQP_Al_ns.png");

  // --------------------------------------------------------------------------
  // PCE vs Al thickness: one line per chip area, with binomial error bars
  // --------------------------------------------------------------------------
  {
    const int nsColors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
                            kOrange+7, kCyan+2, kViolet, kPink+7, kAzure+2, kSpring+5};

    std::vector<double> xAl(al_vals.size());
    for (int iA = 0; iA < (int)al_vals.size(); ++iA)
      xAl[iA] = al_vals[iA];

    // TMultiGraph automatically scales axes to contain all constituent graphs
    TMultiGraph* mg_pce = new TMultiGraph("mg_pce_vs_Al",
        "PCE vs Al thickness;Al thickness [nm];PCE [%]");

    TCanvas* c_pceVsAl = new TCanvas("c_pce_vs_Al", "PCE vs Al thickness", 900, 700);
    TLegend* leg_pce = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_pce->SetBorderSize(1);
    leg_pce->SetFillStyle(0);

    for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
      const int ns = ns_vals[iS];
      const double chipArea_cm2 = 1.0 / ns;

      std::vector<double> yPCE(al_vals.size());
      std::vector<double> yErr(al_vals.size());
      std::vector<double> xErr(al_vals.size(), 0.0);
      for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
        yPCE[iA] = pce_table[iA][iS];
        yErr[iA] = pce_err_table[iA][iS];
      }

      TString gName = TString::Format("g_pce_vs_Al_ns%d", ns);
      TGraphErrors* g = new TGraphErrors((int)al_vals.size(), xAl.data(), yPCE.data(),
                                          xErr.data(), yErr.data());
      g->SetName(gName);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);

      const int color = nsColors[iS < 11 ? iS : 10];
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetFillColor(color);

      fOut->cd();
      g->Write();

      mg_pce->Add(g, "LP");
      leg_pce->AddEntry(g, TString::Format("%.3f cm^{2}", chipArea_cm2), "lp");
    }

    mg_pce->Draw("A");
    leg_pce->Draw();
    fOut->cd();
    mg_pce->Write();
    c_pceVsAl->Write();
    c_pceVsAl->SaveAs("PCE_vs_Al_byChipArea.png");

    delete leg_pce;
    delete c_pceVsAl;
  }

  // --------------------------------------------------------------------------
  // Normalized PCE vs Al thickness: each chip area's PCE is divided by its
  // value at 100 nm (the first Al thickness), so all curves start at 1.0.
  // Error bars use sigma_norm = sigma_pce[iA] / pce[iA=0]; the reference
  // point is defined to be exactly 1 (zero error bar at 100 nm).
  // --------------------------------------------------------------------------
  {
    const int nsColors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
                            kOrange+7, kCyan+2, kViolet, kPink+7, kAzure+2, kSpring+5};

    std::vector<double> xAl(al_vals.size());
    for (int iA = 0; iA < (int)al_vals.size(); ++iA)
      xAl[iA] = al_vals[iA];

    TMultiGraph* mg_norm = new TMultiGraph("mg_pce_vs_Al_norm",
        "Normalized PCE vs Al thickness (100 nm = 1);Al thickness [nm];Normalized PCE");

    TCanvas* c_pceNorm = new TCanvas("c_pce_vs_Al_norm", "Normalized PCE vs Al thickness", 900, 700);
    TLegend* leg_norm = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_norm->SetBorderSize(1);
    leg_norm->SetFillStyle(0);

    for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
      const int ns = ns_vals[iS];
      const double chipArea_cm2 = 1.0 / ns;

      const double ref = pce_table[0][iS];
      if (ref <= 0.0) continue;

      std::vector<double> yNorm(al_vals.size());
      std::vector<double> yErr(al_vals.size());
      std::vector<double> xErr(al_vals.size(), 0.0);
      for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
        yNorm[iA] = pce_table[iA][iS] / ref;
        // Reference point is defined as exactly 1: assign zero error there
        yErr[iA]  = (iA == 0) ? 0.0 : pce_err_table[iA][iS] / ref;
      }

      TString gName = TString::Format("g_pce_norm_vs_Al_ns%d", ns);
      TGraphErrors* g = new TGraphErrors((int)al_vals.size(), xAl.data(), yNorm.data(),
                                          xErr.data(), yErr.data());
      g->SetName(gName);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);

      const int color = nsColors[iS < 11 ? iS : 10];
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetFillColor(color);

      fOut->cd();
      g->Write();

      mg_norm->Add(g, "LP");
      leg_norm->AddEntry(g, TString::Format("%.3f cm^{2}", chipArea_cm2), "lp");
    }

    mg_norm->Draw("A");
    leg_norm->Draw();
    fOut->cd();
    mg_norm->Write();
    c_pceNorm->Write();
    c_pceNorm->SaveAs("PCE_vs_Al_byChipArea_norm.png");

    delete leg_norm;
    delete c_pceNorm;
  }

  // --------------------------------------------------------------------------
  // Fixed Al, varying ns
  // --------------------------------------------------------------------------
  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    const int al_val = al_vals[iA];

    if (hists_by_al[iA].empty()) continue;

    TString cName  = TString::Format("c_eDep_overlay_Al%d", al_val);
    TString cTitle = TString::Format("Hit EDeps for Al=%d nm, varying ns", al_val);
    TCanvas* c = new TCanvas(cName, cTitle, 900, 700);
    c->SetLogy();

    TLegend* leg = new TLegend(0.68, 0.55, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);

    double maxY = 0.0;
    for (TH1F* h : hists_by_al[iA]) {
      if (h->GetMaximum() > maxY) maxY = h->GetMaximum();
    }

    bool first = true;
    for (TH1F* h : hists_by_al[iA]) {
      h->SetTitle(TString::Format("Hit EDeps for Al=%d nm, varying ns;eDep [eV];nEvents",
                                  al_val));
      h->SetMaximum(1.15 * maxY);
      h->SetMinimum(0.5);

      TString histName = h->GetName();
      int parsed_al = -1;
      int ns_val = -1;
      sscanf(histName.Data(), "h_eDep_Al%d_ns%d", &parsed_al, &ns_val);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }

      leg->AddEntry(h, TString::Format("ns = %d", ns_val), "l");
    }

    TLine* line_Al1 = new TLine(0.34e-3, 0, 0.34e-3, 1.15 * maxY);
    TLine* line_Al2 = new TLine(0.68e-3, 0, 0.68e-3, 1.15 * maxY);
    TLine* line_Al3 = new TLine(1.02e-3, 0, 1.02e-3, 1.15 * maxY);
    line_Al1->SetLineColor(kGreen+2);
    line_Al2->SetLineColor(kGreen+2);
    line_Al3->SetLineColor(kGreen+2);
    line_Al1->SetLineStyle(2);
    line_Al2->SetLineStyle(2);
    line_Al3->SetLineStyle(2);
    line_Al1->SetLineWidth(2);
    line_Al2->SetLineWidth(2);
    line_Al3->SetLineWidth(2);

    line_Al1->Draw();
    line_Al2->Draw();
    line_Al3->Draw();

    leg->AddEntry(line_Al1, "Al gap", "l");

    leg->Draw();

    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("h_eDep_overlay_Al%d.png", al_val));

    delete leg;
    delete c;
  }

  // --------------------------------------------------------------------------
  // Fixed ns, varying Al
  // --------------------------------------------------------------------------
  for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
    const int ns_val = ns_vals[iS];

    if (hists_by_ns[iS].empty()) continue;

    TString cName  = TString::Format("c_eDep_overlay_ns%d", ns_val);
    TString cTitle = TString::Format("Hit EDeps for ns=%d, varying Al", ns_val);
    TCanvas* c = new TCanvas(cName, cTitle, 900, 700);
    c->SetLogy();

    TLegend* leg = new TLegend(0.68, 0.55, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);

    double maxY = 0.0;
    for (TH1F* h : hists_by_ns[iS]) {
      if (h->GetMaximum() > maxY) maxY = h->GetMaximum();
    }

    bool first = true;
    int color_index = 0;

    for (TH1F* h : hists_by_ns[iS]) {
      h->SetTitle(TString::Format("Hit EDeps for ns=%d, varying Al;eDep [eV];nEvents",
                                  ns_val));
      h->SetMaximum(1.15 * maxY);
      h->SetMinimum(0.5);

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
      int parsed_ns = -1;
      sscanf(histName.Data(), "h_eDep_Al%d_ns%d", &al_val, &parsed_ns);

      if (first) {
        h->Draw("HIST");
        first = false;
      } else {
        h->Draw("HIST SAME");
      }

      leg->AddEntry(h, TString::Format("Al = %d nm", al_val), "l");
      ++color_index;
    }

    TLine* line_Al1 = new TLine(0.34e-3, 0, 0.34e-3, 1.15 * maxY);
    TLine* line_Al2 = new TLine(0.68e-3, 0, 0.68e-3, 1.15 * maxY);
    TLine* line_Al3 = new TLine(1.02e-3, 0, 1.02e-3, 1.15 * maxY);
    line_Al1->SetLineColor(kGreen+2);
    line_Al2->SetLineColor(kGreen+2);
    line_Al3->SetLineColor(kGreen+2);
    line_Al1->SetLineStyle(2);
    line_Al2->SetLineStyle(2);
    line_Al3->SetLineStyle(2);
    line_Al1->SetLineWidth(2);
    line_Al2->SetLineWidth(2);
    line_Al3->SetLineWidth(2);

    line_Al1->Draw();
    line_Al2->Draw();
    line_Al3->Draw();

    leg->AddEntry(line_Al1, "Al gap", "l");

    leg->Draw();

    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("h_eDep_overlay_ns%d.png", ns_val));

    delete leg;
    delete c;
  }

  fOut->Write();
  fOut->Close();

  delete c_pce;
  delete c_nqp;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

// Ns_QuasiparticleAnalysis
//
// Estimates the peak number of concurrent quasiparticles (QPs) in the Al film
// as a function of Al thickness and chip area (encoded via the downsampling
// factor "ns", where chip area = 1/ns cm^2).
//
// The approach is a simplified "event-driven" QP population model:
//   1. Each phonon hit deposits energy E in the Al film. G4CMP reports this
//      as eDep_eV at the Final Time (endT_ns) -- i.e. the time the phonon
//      arrived at the Al surface and was absorbed via KaplanQP.
//   2. That deposited energy is converted to a number of QPs:
//        nQP = 2 * E / Delta_Al
//      where Delta_Al = 0.34 meV is the Al superconducting gap, and the
//      factor of 2 accounts for the pair of QPs broken per Cooper pair.
//   3. Each QP creation is modeled as an instantaneous step up at endT_ns,
//      followed by a step down of the same magnitude at endT_ns + qpLifetime_ns,
//      approximating an exponential decay by a rectangular pulse.
//   4. All step events are sorted in time, then scanned to build a running
//      total. The maximum of the running total is the peak concurrent QP count.
//   5. A TGraph of net QP vs time is saved for each (Al, ns) combination,
//      and the peak value is filled into a 2D summary histogram.
void Ns_QuasiparticleAnalysis()
{
  // Al film thickness values [nm] scanned in the simulation
  std::vector<int> al_vals = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

  // Downsampling factors: the simulation was run with 1/ns of the full chip
  // area, so ns=1 corresponds to the full 1 cm^2 chip, ns=2 to 0.5 cm^2, etc.
  std::vector<int> ns_vals = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  // Al superconducting gap energy [eV]. Each Cooper pair break produces 2 QPs,
  // each requiring at least Delta_Al of energy.
  const double alGap_eV = 0.34e-3;

  // QP lifetime [ns]: drawn per-hit from an exponential distribution with
  // mean 100 us (100000 ns). For an exponential, the mean equals the decay
  // constant tau, so rng.Exp(qpLifetimeMean_ns) gives the correct distribution.
  TRandom3 rng(0);
  const double qpLifetimeMean_ns = 100000.0;  // 100 us in ns (exponential mean = tau)

  // Directory containing the G4CMP hit output files for this scan
  const std::string baseDir = "../G4Macros/260419_run/lQPD_Al_PF";

  // Open the ROOT output file for all results
  TFile* fOut = new TFile("QP_Al_ns.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file QP_Al_ns.root" << std::endl;
    return;
  }
  fOut->cd();

  // 2D summary histogram: peak concurrent QPs for each (Al thickness, chip area) pair.
  // Each bin will be filled with the maximum QP population reached during the event.
  TH2F* h_qp_al_ns = new TH2F("h_qp_al_ns",
      "Peak concurrent quasiparticles vs Al thickness and chip area;Al thickness [nm];Chip area [cm^{2}];Peak N_{QP}",
      (int)al_vals.size(), 0, (int)al_vals.size(),
      (int)ns_vals.size(), 0, (int)ns_vals.size());

  // Label x-axis bins with the Al thickness values [nm]
  for (int i = 0; i < (int)al_vals.size(); ++i) {
    h_qp_al_ns->GetXaxis()->SetBinLabel(i + 1, std::to_string(al_vals[i]).c_str());
  }
  // Label y-axis bins with the corresponding chip area [cm^2] = 1/ns
  for (int j = 0; j < (int)ns_vals.size(); ++j) {
    const double chipArea_cm2 = 1.0 / ns_vals[j];
    TString areaLabel = TString::Format("%.3f", chipArea_cm2);
    h_qp_al_ns->GetYaxis()->SetBinLabel(j + 1, areaLabel.Data());
  }

  // Accumulates every qpLifetime_ns value actually used in the hit loop below.
  // Plotted as a histogram after the main loop.
  std::vector<double> usedLifetimes;

  // Storage for net-QP-vs-time graphs, indexed [iA][iS].
  // Filled in the main loop; used after the loop for overlay plots.
  std::vector<std::vector<TGraph*> > netQP_graphs(
      al_vals.size(), std::vector<TGraph*>(ns_vals.size(), nullptr));

  // Loop over all (Al thickness, ns) combinations
  for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
    for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {

      const int al = al_vals[iA];
      const int ns = ns_vals[iS];

      // Build the path to the hit output file for this (Al, ns) combination
      std::ostringstream hitPath;
      hitPath << baseDir << "/Hits_Al" << al << "_ns" << ns << ".txt";

      // Parse all hits from the file; keyed by event ID
      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitPath.str());

      // A QPEvent represents a discrete change in the QP population at a given time:
      //   deltaQP > 0: QPs are created (phonon absorbed into film)
      //   deltaQP < 0: QPs are removed (lifetime expiry)
      struct QPEvent {
        double time_ns;
        double deltaQP;
      };
      std::vector<QPEvent> qpEvents;

      // Running tally of all QPs created (before recombination), for diagnostics
      double totalCreatedQP = 0.0;

      // Loop over all events and hits to build the list of QP creation/removal events
      for (const auto& kv : hitInfo) {
        for (const Hit& h : kv.second) {
          if (h.eDep_eV > 0.0) {
            // Convert deposited energy to number of QPs: nQP = 2 * E / Delta_Al.
            // The factor of 2 is because each broken Cooper pair produces 2 QPs.
            // endT_ns is the time the phonon arrived at the film (Final Time from G4CMP).
            const int nQP = (int)std::round((h.eDep_eV / alGap_eV) * 2.0);
            totalCreatedQP += nQP;

            // Push a creation event at the phonon arrival time (endT_ns)
            qpEvents.push_back({h.endT_ns, (double)nQP});
            // Push one removal event per individual QP, each with its own exponentially-
            // drawn lifetime. This removes one QP at a time rather than all at once.
            for (int iQP = 0; iQP < nQP; ++iQP) {
              const double qpLifetime_ns = rng.Exp(qpLifetimeMean_ns);
              usedLifetimes.push_back(qpLifetime_ns);
              qpEvents.push_back({h.endT_ns + qpLifetime_ns, -1.0});
            }
          }
        }
      }

      // Sort all QP events chronologically. For events at exactly the same time,
      // process creation (+) before removal (-) to avoid spurious negative counts.
      std::sort(qpEvents.begin(), qpEvents.end(),
                [](const QPEvent& a, const QPEvent& b) {
                  if (a.time_ns == b.time_ns) return a.deltaQP > b.deltaQP;
                  return a.time_ns < b.time_ns;
                });

      // Scan through the sorted events to compute the running QP population vs time.
      // Record the peak population and build vectors for plotting.
      double runningQP = 0.0;
      double peakQP = 0.0;
      std::vector<double> times;
      std::vector<double> netQP;

      if (!qpEvents.empty()) {
        // Prepend a zero point just before the first event so the graph starts at 0
        times.push_back(qpEvents.front().time_ns);
        netQP.push_back(0.0);

        for (const auto& ev : qpEvents) {
          runningQP += ev.deltaQP;
          // Clamp to zero: QP count cannot go negative due to floating-point imprecision
          if (runningQP < 0.0) runningQP = 0.0;
          // Track the maximum concurrent QP population
          if (runningQP > peakQP) peakQP = runningQP;

          times.push_back(ev.time_ns);
          netQP.push_back(runningQP);
        }

        // Append a trailing zero if the final QP count is non-zero
        // (this shouldn't happen if all lifetimes are accounted for, but guards
        // against floating-point residuals in the graph)
        if (netQP.back() != 0.0) {
          times.push_back(qpEvents.back().time_ns);
          netQP.push_back(0.0);
        }
      }

      // Fill the 2D summary histogram and the per-ns storage table with the peak QP count
      const double finalQP = netQP.empty() ? 0.0 : netQP.back();
      h_qp_al_ns->SetBinContent(iA + 1, iS + 1, peakQP);

      std::cout << "Al=" << al
                << " ns=" << ns
                << " Total created QP=" << totalCreatedQP
                << " Peak QP=" << peakQP
                << " Final QP=" << finalQP
                << " Lifetime=" << qpLifetimeMean_ns << " ns (exponential mean)"
                << std::endl;

      // Store the TGraph for this (Al, ns) combination.
      // Overlay plots (fixed-ns and fixed-Al) are produced after the main loop.
      if (!times.empty()) {
        TString gName = TString::Format("g_netQP_Al%d_ns%d", al, ns);
        TGraph* g = new TGraph((int)times.size(), times.data(), netQP.data());
        g->SetName(gName);
        g->SetTitle(TString::Format(
            "Net QP vs Time, Al=%d nm, ns=%d;Time [ns];Net N_{QP}", al, ns));
        g->SetLineWidth(2);

        fOut->cd();
        g->Write();

        netQP_graphs[iA][iS] = g;
      }
    }
  }

  // --------------------------------------------------------------------------
  // Overlay plots: fixed ns (chip area), all Al thicknesses on one canvas
  // --------------------------------------------------------------------------
  {
    const int alColors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
                            kOrange+7, kCyan+2, kViolet, kPink+7, kAzure+2};

    for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
      const int ns = ns_vals[iS];
      const double chipArea_cm2 = 1.0 / ns;

      // Compute global axis range across all Al thicknesses for this ns
      double xMin = 1e99, xMax = -1e99, yMax = 0.0;
      for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
        TGraph* g = netQP_graphs[iA][iS];
        if (!g) continue;
        for (int ip = 0; ip < g->GetN(); ++ip) {
          double x, y;
          g->GetPoint(ip, x, y);
          if (x < xMin) xMin = x;
          if (x > xMax) xMax = x;
          if (y > yMax) yMax = y;
        }
      }
      if (xMin > xMax) continue;

      TString cName  = TString::Format("c_netQP_ns%d", ns);
      TString cTitle = TString::Format(
          "Net QP vs Time, chip area=%.3f cm^{2}", chipArea_cm2);
      TCanvas* c = new TCanvas(cName, cTitle, 900, 700);

      TLegend* leg = new TLegend(0.62, 0.55, 0.88, 0.88);
      leg->SetBorderSize(1);
      leg->SetFillStyle(0);

      bool firstDraw = true;
      for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
        TGraph* g = netQP_graphs[iA][iS];
        if (!g) continue;
        g->SetLineColor(alColors[iA < 10 ? iA : 9]);
        g->SetLineWidth(2);
        if (firstDraw) {
          g->Draw("AL");
          g->GetXaxis()->SetLimits(xMin, xMax);
          g->SetMinimum(0.0);
          g->SetMaximum(1.1 * yMax);
          firstDraw = false;
        } else {
          g->Draw("L SAME");
        }
        leg->AddEntry(g, TString::Format("Al = %d nm", al_vals[iA]), "l");
      }

      leg->Draw();
      fOut->cd();
      c->Write();
      c->SaveAs(TString::Format("netQP_ns%d_allAl.png", ns));
      delete leg;
      delete c;
    }
  }

  // --------------------------------------------------------------------------
  // Overlay plots: fixed Al thickness, all ns values (chip areas) on one canvas
  // --------------------------------------------------------------------------
  {
    const int nsColors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
                            kOrange+7, kCyan+2, kViolet, kPink+7, kAzure+2, kSpring+5};

    for (int iA = 0; iA < (int)al_vals.size(); ++iA) {
      const int al = al_vals[iA];

      // Compute global axis range across all ns values for this Al
      double xMin = 1e99, xMax = -1e99, yMax = 0.0;
      for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
        TGraph* g = netQP_graphs[iA][iS];
        if (!g) continue;
        for (int ip = 0; ip < g->GetN(); ++ip) {
          double x, y;
          g->GetPoint(ip, x, y);
          if (x < xMin) xMin = x;
          if (x > xMax) xMax = x;
          if (y > yMax) yMax = y;
        }
      }
      if (xMin > xMax) continue;

      TString cName  = TString::Format("c_netQP_Al%d", al);
      TString cTitle = TString::Format("Net QP vs Time, Al=%d nm", al);
      TCanvas* c = new TCanvas(cName, cTitle, 900, 700);

      TLegend* leg = new TLegend(0.62, 0.55, 0.88, 0.88);
      leg->SetBorderSize(1);
      leg->SetFillStyle(0);

      bool firstDraw = true;
      for (int iS = 0; iS < (int)ns_vals.size(); ++iS) {
        TGraph* g = netQP_graphs[iA][iS];
        if (!g) continue;
        g->SetLineColor(nsColors[iS < 11 ? iS : 10]);
        g->SetLineWidth(2);
        if (firstDraw) {
          g->Draw("AL");
          g->GetXaxis()->SetLimits(xMin, xMax);
          g->SetMinimum(0.0);
          g->SetMaximum(1.1 * yMax);
          firstDraw = false;
        } else {
          g->Draw("L SAME");
        }
        leg->AddEntry(g, TString::Format("%.3f cm^{2}", 1.0 / ns_vals[iS]), "l");
      }

      leg->Draw();
      fOut->cd();
      c->Write();
      c->SaveAs(TString::Format("netQP_Al%d_allNs.png", al));
      delete leg;
      delete c;
    }
  }

  // Write the 2D peak-QP summary histogram and save as a color-map PNG
  fOut->cd();
  h_qp_al_ns->Write();

  TCanvas* c_qp = new TCanvas("c_qp_al_ns", "c_qp_al_ns", 900, 700);
  c_qp->SetRightMargin(0.18);
  h_qp_al_ns->SetStats(0);
  h_qp_al_ns->Draw("COLZ TEXT");
  c_qp->Write();
  c_qp->SaveAs("QP_Al_ns.png");

  // Plot the distribution of QP lifetimes actually used in the simulation
  if (!usedLifetimes.empty()) {
    const double ltMin = *std::min_element(usedLifetimes.begin(), usedLifetimes.end());
    const double ltMax = *std::max_element(usedLifetimes.begin(), usedLifetimes.end());
    // Add a small margin so the outermost values are not on the bin edge
    const double margin = (ltMax - ltMin) * 0.05;
    TH1F* h_lt = new TH1F("h_qpLifetime",
        "QP lifetime distribution (exponential, #tau=100 #mus);Lifetime [ns];Counts",
        200, ltMin - margin, ltMax + margin);
    for (double v : usedLifetimes) h_lt->Fill(v);

    // Overlay an exponential reference curve: f(x) = N * binWidth * (1/tau) * exp(-x/tau)
    const double tau = qpLifetimeMean_ns;
    const double binWidth = h_lt->GetBinWidth(1);
    const double nEntries = (double)usedLifetimes.size();
    TF1* f_exp = new TF1("f_qpLifetimeExp", "[0]*exp(-x/[1])",
                         ltMin - margin, ltMax + margin);
    f_exp->SetParameter(0, nEntries * binWidth / tau);
    f_exp->SetParameter(1, tau);
    f_exp->SetLineColor(kRed);
    f_exp->SetLineWidth(2);

    TCanvas* c_lt = new TCanvas("c_qpLifetime", "QP lifetime distribution (exponential)", 800, 600);
    h_lt->SetLineColor(kBlue + 1);
    h_lt->SetLineWidth(2);
    h_lt->Draw("HIST");
    f_exp->Draw("SAME");

    TLegend* leg_lt = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg_lt->SetBorderSize(1);
    leg_lt->SetFillStyle(0);
    leg_lt->AddEntry(h_lt, "Sampled lifetimes", "l");
    leg_lt->AddEntry(f_exp, TString::Format("Exp(#tau=%.0f ns)", tau), "l");
    leg_lt->Draw();

    fOut->cd();
    h_lt->Write();
    f_exp->Write();
    c_lt->Write();
    c_lt->SaveAs("qpLifetime_distribution.png");
    delete leg_lt;
    delete c_lt;
  }

  fOut->Write();
  fOut->Close();

  delete c_qp;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

//---------------------------------------------------------------------------------------
// Helpers for Mems_PCEStudy: draw detector geometry overlay on a ROOT pad.
// Geometry is taken from RISQTutorialDetectorConstruction.cc and the
// single_squat STL files (ASCII format, coordinates in µm).
//---------------------------------------------------------------------------------------

// Scan an ASCII STL file and accumulate the bounding box in mm (coords * 1e-3).
static void UpdateSTLBounds(const std::string& path,
                             double& xmin, double& xmax,
                             double& ymin, double& ymax)
{
  std::ifstream f(path.c_str());
  if (!f.is_open()) return;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream ss(line);
    std::string tok; ss >> tok;
    if (tok != "vertex") continue;
    double x, y, z; ss >> x >> y >> z;
    x *= 1e-3; y *= 1e-3;   // µm → mm
    xmin = std::min(xmin, x); xmax = std::max(xmax, x);
    ymin = std::min(ymin, y); ymax = std::max(ymax, y);
  }
}

// Draw the boundary edges of the top face (|nz| > 0.9) of one ASCII STL part
// projected onto XY.  off_x / off_y are the placement offsets in mm.
static void DrawSTLTopFaceXY(const std::string& path,
                              double off_x, double off_y,
                              Color_t color, Width_t lw)
{
  std::ifstream f(path.c_str());
  if (!f.is_open()) {
    std::cerr << "DrawSTLTopFaceXY: cannot open " << path << std::endl;
    return;
  }

  // Quantize to 1 nm (1e-6 mm) for edge deduplication.
  typedef std::pair<long long, long long>  Pt2;
  typedef std::pair<Pt2, Pt2>             EdgeKey;

  std::map<EdgeKey, int>                     edgeCnt;
  std::map<EdgeKey, std::array<double, 4> >  edgeCoords;

  std::string line;
  double nx = 0, ny = 0, nz = 0;
  double vx[3], vy[3];
  int vi = -1;

  while (std::getline(f, line)) {
    std::istringstream ss(line);
    std::string tok; ss >> tok;
    if (tok == "facet") {
      std::string tmp; ss >> tmp;   // skip "normal"
      ss >> nx >> ny >> nz;
      vi = 0;
    } else if (tok == "vertex" && vi >= 0 && vi < 3) {
      double x, y, z; ss >> x >> y >> z;
      vx[vi] = x * 1e-3;
      vy[vi] = y * 1e-3;
      ++vi;
    } else if (tok == "endfacet") {
      if (vi == 3 && TMath::Abs(nz) > 0.9) {
        for (int j = 0; j < 3; ++j) {
          int k = (j + 1) % 3;
          Pt2 pa = std::make_pair((long long)(vx[j] * 1e6 + 0.5),
                                  (long long)(vy[j] * 1e6 + 0.5));
          Pt2 pb = std::make_pair((long long)(vx[k] * 1e6 + 0.5),
                                  (long long)(vy[k] * 1e6 + 0.5));
          if (pa > pb) std::swap(pa, pb);
          EdgeKey ek = std::make_pair(pa, pb);
          edgeCnt[ek]++;
          if (edgeCnt[ek] == 1) {
            std::array<double, 4> c = { {vx[j], vy[j], vx[k], vy[k]} };
            edgeCoords[ek] = c;
          }
        }
      }
      vi = -1;
    }
  }

  for (std::map<EdgeKey,int>::iterator it = edgeCnt.begin(); it != edgeCnt.end(); ++it) {
    if (it->second != 1) continue;   // interior edge
    const std::array<double, 4>& c = edgeCoords[it->first];
    TLine* l = new TLine(c[0] + off_x, c[1] + off_y, c[2] + off_x, c[3] + off_y);
    l->SetLineColor(color);
    l->SetLineWidth(lw);
    l->Draw();
  }
}

// Draw the SQUAT KID outline on the current pad.  Placement mirrors
// RISQTutorialDetectorConstruction: the combined bounding-box centre of all
// 5 parts is removed, then (targetX_mm, targetY_mm) is added.
static void DrawSquatOverlayXY(const std::string& stlDir,
                                double targetX_mm, double targetY_mm,
                                Color_t color, Width_t lw)
{
  double xmin =  1e18, xmax = -1e18;
  double ymin =  1e18, ymax = -1e18;
  for (int i = 1; i <= 5; ++i) {
    std::string path = stlDir + "/single_squat_BE" + std::to_string(i) + ".STL";
    UpdateSTLBounds(path, xmin, xmax, ymin, ymax);
  }
  if (xmin > xmax) return;
  const double offX = targetX_mm - 0.5 * (xmin + xmax);
  const double offY = targetY_mm - 0.5 * (ymin + ymax);
  for (int i = 1; i <= 5; ++i) {
    std::string path = stlDir + "/single_squat_BE" + std::to_string(i) + ".STL";
    DrawSTLTopFaceXY(path, offX, offY, color, lw);
  }
}

// Draw chip boundary, feedline, ground planes, and SQUAT KID on the current
// ROOT pad.  All dimensions computed from RISQTutorialDetectorConstruction.
// chipHalf_mm : Ge chip half-width (5.0 mm for numsensors = 1).
// stlDir      : directory containing single_squat_BE[1-5].STL.
static void DrawDetectorOverlayXY(double chipHalf_mm,
                                   const std::string& stlDir,
                                   Color_t color = kWhite, Width_t lw = 1)
{
  // Chip boundary
  TBox* chipBox = new TBox(-chipHalf_mm, -chipHalf_mm, chipHalf_mm, chipHalf_mm);
  chipBox->SetLineColor(color);
  chipBox->SetLineWidth(lw);
  chipBox->SetFillStyle(0);
  chipBox->Draw("l");

  // Feedline central conductor: full chip width, half-height = 4.5 µm
  const double flHY = 4.5e-3;
  TBox* flBox = new TBox(-chipHalf_mm, -flHY, chipHalf_mm, flHY);
  flBox->SetLineColor(color);
  flBox->SetLineWidth(lw);
  flBox->SetFillStyle(0);
  flBox->Draw("l");

  // Upper ground plane: centre +56.5 µm, half-height 50 µm
  TBox* ugpBox = new TBox(-chipHalf_mm, 56.5e-3 - 50.0e-3,
                           chipHalf_mm, 56.5e-3 + 50.0e-3);
  ugpBox->SetLineColor(color);
  ugpBox->SetLineWidth(lw);
  ugpBox->SetFillStyle(0);
  ugpBox->Draw("l");

  // Lower ground plane: centre -8.5 µm, half-height 2 µm
  TBox* lgpBox = new TBox(-chipHalf_mm, -8.5e-3 - 2.0e-3,
                           chipHalf_mm, -8.5e-3 + 2.0e-3);
  lgpBox->SetLineColor(color);
  lgpBox->SetLineWidth(lw);
  lgpBox->SetFillStyle(0);
  lgpBox->Draw("l");

  // SQUAT KID placed at (0, -200 µm) per RISQTutorialDetectorConstruction
  DrawSquatOverlayXY(stlDir, 0.0, -0.2, color, lw);
}

//---------------------------------------------------------------------------------------

// Mems_PCEStudy
//
// Reads G4CMP output from the MEMS grid scan: phonons launched from each
// point of a grid_size x grid_size grid across the 1 cm x 1 cm chip.
// For each of the four configurations (Al-only vs Al+Nb) x (SW vs PF loss),
// computes the phonon collection efficiency (PCE = total hit energy / total
// primary energy) and produces a 2D histogram of PCE with X position on the
// x-axis and Y position on the y-axis.
//
// Expected file naming (from mems.sh / run_mems.sh):
//   ../G4Macros/260519_run/<config>/Hits_Al<Al>_Nb<Nb>_ix<ix>_iy<iy>.txt
//   ../G4Macros/260519_run/<config>/Primary_Al<Al>_Nb<Nb>_ix<ix>_iy<iy>.txt
void Mems_PCEStudy(int grid_size = 30, bool drawOverlay = false, std::string stlDir = "../../single_squat")
{
  const int al = 30;
  const int nb = 70;
  const std::string baseRunDir = "../../../../261005_run_mems";

  // The four run configurations produced by run_mems.sh
  struct Config {
    std::string dirName;
    std::string label;
  };
  const std::vector<Config> configs = {
    //{"MEMS_Al_SW",    "Al only, SW loss"},
    //{"MEMS_Al_PF",    "Al only, PF loss"},
    {"MEMS_Al_Nb_SW", "Al+Nb, SW loss"},
    {"MEMS_Al_Nb_PF", "Al+Nb, PF loss"},
  };

  TFile* fOut = new TFile("Mems_PCE.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create output file Mems_PCE.root" << std::endl;
    return;
  }
  fOut->cd();

  // Line plots combining both loss configurations: PCE along the X axis at
  // Y=0, and PCE along the Y axis at X=0.
  const int configColors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};

  TMultiGraph* mg_pceVsX_atYZero = new TMultiGraph(
      "mg_pceVsX_atYZero", "PCE vs X Position (Y=0);X [mm];PCE");
  TCanvas* c_pceVsX_atYZero = new TCanvas("c_pceVsX_atYZero", "PCE vs X (Y=0)", 900, 700);
  TLegend* leg_pceVsX_atYZero = new TLegend(0.68, 0.75, 0.88, 0.88);
  leg_pceVsX_atYZero->SetBorderSize(1);
  leg_pceVsX_atYZero->SetFillStyle(0);

  TMultiGraph* mg_pceVsY_atXZero = new TMultiGraph(
      "mg_pceVsY_atXZero", "PCE vs Y Position (X=0);Y [mm];PCE");
  TCanvas* c_pceVsY_atXZero = new TCanvas("c_pceVsY_atXZero", "PCE vs Y (X=0)", 900, 700);
  TLegend* leg_pceVsY_atXZero = new TLegend(0.68, 0.75, 0.88, 0.88);
  leg_pceVsY_atXZero->SetBorderSize(1);
  leg_pceVsY_atXZero->SetFillStyle(0);

  TMultiGraph* mg_pceVsDiag_xEqY = new TMultiGraph(
      "mg_pceVsDiag_xEqY", "PCE along x=y Diagonal;Position along x=y diagonal [mm];PCE");
  TCanvas* c_pceVsDiag_xEqY = new TCanvas("c_pceVsDiag_xEqY", "PCE vs diagonal (x=y)", 900, 700);
  TLegend* leg_pceVsDiag_xEqY = new TLegend(0.68, 0.75, 0.88, 0.88);
  leg_pceVsDiag_xEqY->SetBorderSize(1);
  leg_pceVsDiag_xEqY->SetFillStyle(0);

  for (size_t iCfg = 0; iCfg < configs.size(); ++iCfg) {
    const Config& cfg = configs[iCfg];
    const std::string baseDir = baseRunDir + "/" + cfg.dirName;
    const TString tag = TString(cfg.dirName.c_str());

    // 2D PCE map: x-axis = X position [mm], y-axis = Y position [mm], z = PCE
    TH2F* h_pce = new TH2F(
        TString::Format("h_pce_%s", tag.Data()),
        TString::Format("PCE Map (%s);X [mm];Y [mm];PCE", cfg.label.c_str()),
        grid_size, -5.0, 5.0,
        grid_size, -5.0, 5.0);
    h_pce->SetDirectory(fOut);
    h_pce->SetStats(0);

    for (int ix = 0; ix < grid_size; ++ix) {
      for (int iy = 0; iy < grid_size; ++iy) {

        std::ostringstream primPath, hitPath;
        primPath << baseDir << "/Primary_Al" << al << "_Nb" << nb
                 << "_ix" << ix << "_iy" << iy << ".txt";
        hitPath  << baseDir << "/Hits_Al"    << al << "_Nb" << nb
                 << "_ix" << ix << "_iy" << iy << ".txt";

        const std::map<int, PrimaryInfo> primaryInfo =
            ParsePrimaryTextFileForPrimaries(primPath.str());
        const std::map<int, std::vector<Hit> > hitInfo =
            ParseHitTextFileForHits(hitPath.str());

        // All primaries in this file share the same launch position;
        // use the accumulated position from the first entry for the grid coordinates.
        double primX_mm = 0.0;
        double primY_mm = 0.0;
        double totalPrimaryEnergy_eV = 0.0;
        for (const auto& kv : primaryInfo) {
          primX_mm = kv.second.X_mm;
          primY_mm = kv.second.Y_mm;
          totalPrimaryEnergy_eV += kv.second.energy_eV;
        }

        double totalHitEnergy_eV = 0.0;
        for (const auto& kv : hitInfo) {
          for (const Hit& h : kv.second) {
            totalHitEnergy_eV += h.eDep_eV;
          }
        }

        const double pce = (totalPrimaryEnergy_eV > 0.0)
                             ? (totalHitEnergy_eV / totalPrimaryEnergy_eV)
                             : 0.0;

        const int binX = h_pce->GetXaxis()->FindBin(primX_mm);
        const int binY = h_pce->GetYaxis()->FindBin(primY_mm);
        h_pce->SetBinContent(binX, binY, pce);

        std::cout << cfg.dirName
                  << "  ix=" << ix << " iy=" << iy
                  << "  pos=(" << primX_mm << ", " << primY_mm << ") mm"
                  << "  PCE=" << 100.0 * pce << " %"
                  << std::endl;
      }
    }

    fOut->cd();
    h_pce->Write();

    TCanvas* c = new TCanvas(
        TString::Format("c_pce_%s", tag.Data()),
        TString::Format("PCE Map: %s", cfg.label.c_str()),
        900, 700);
    c->SetRightMargin(0.18);
    h_pce->Draw("COLZ");
    if (drawOverlay) DrawDetectorOverlayXY(5.0, stlDir);
    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("Mems_PCE_%s.png", tag.Data()));
    delete c;

    const int color = configColors[iCfg < 4 ? iCfg : 3];

    // Slice through the PCE map at the Y bin nearest Y=0
    {
      const int binY0 = h_pce->GetYaxis()->FindBin(0.0);
      std::vector<double> xVals, pceVals;
      for (int ix = 1; ix <= h_pce->GetNbinsX(); ++ix) {
        xVals.push_back(h_pce->GetXaxis()->GetBinCenter(ix));
        pceVals.push_back(h_pce->GetBinContent(ix, binY0));
      }
      TGraph* g = new TGraph((int)xVals.size(), xVals.data(), pceVals.data());
      g->SetName(TString::Format("g_pceVsX_atYZero_%s", tag.Data()));
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      fOut->cd();
      g->Write();
      mg_pceVsX_atYZero->Add(g, "LP");
      leg_pceVsX_atYZero->AddEntry(g, cfg.label.c_str(), "lp");
    }

    // Slice through the PCE map at the X bin nearest X=0
    {
      const int binX0 = h_pce->GetXaxis()->FindBin(0.0);
      std::vector<double> yVals, pceVals;
      for (int iy = 1; iy <= h_pce->GetNbinsY(); ++iy) {
        yVals.push_back(h_pce->GetYaxis()->GetBinCenter(iy));
        pceVals.push_back(h_pce->GetBinContent(binX0, iy));
      }
      TGraph* g = new TGraph((int)yVals.size(), yVals.data(), pceVals.data());
      g->SetName(TString::Format("g_pceVsY_atXZero_%s", tag.Data()));
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      fOut->cd();
      g->Write();
      mg_pceVsY_atXZero->Add(g, "LP");
      leg_pceVsY_atXZero->AddEntry(g, cfg.label.c_str(), "lp");
    }

    // Slice through the PCE map along the diagonal x=y
    {
      const int nBins = std::min(h_pce->GetNbinsX(), h_pce->GetNbinsY());
      std::vector<double> diagVals, pceVals;
      for (int i = 1; i <= nBins; ++i) {
        const double x = h_pce->GetXaxis()->GetBinCenter(i);
        // Signed distance from the origin along the x=y diagonal
        diagVals.push_back(x * std::sqrt(2.0));
        pceVals.push_back(h_pce->GetBinContent(i, i));
      }
      TGraph* g = new TGraph((int)diagVals.size(), diagVals.data(), pceVals.data());
      g->SetName(TString::Format("g_pceVsDiag_xEqY_%s", tag.Data()));
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      fOut->cd();
      g->Write();
      mg_pceVsDiag_xEqY->Add(g, "LP");
      leg_pceVsDiag_xEqY->AddEntry(g, cfg.label.c_str(), "lp");
    }
  }

  c_pceVsX_atYZero->cd();
  c_pceVsX_atYZero->SetGrid();
  mg_pceVsX_atYZero->Draw("A");
  leg_pceVsX_atYZero->Draw();
  fOut->cd();
  mg_pceVsX_atYZero->Write();
  c_pceVsX_atYZero->Write();
  c_pceVsX_atYZero->SaveAs("Mems_PCE_vs_X_atYZero.png");
  delete leg_pceVsX_atYZero;
  delete c_pceVsX_atYZero;

  c_pceVsY_atXZero->cd();
  c_pceVsY_atXZero->SetGrid();
  mg_pceVsY_atXZero->Draw("A");
  leg_pceVsY_atXZero->Draw();
  fOut->cd();
  mg_pceVsY_atXZero->Write();
  c_pceVsY_atXZero->Write();
  c_pceVsY_atXZero->SaveAs("Mems_PCE_vs_Y_atXZero.png");
  delete leg_pceVsY_atXZero;
  delete c_pceVsY_atXZero;

  c_pceVsDiag_xEqY->cd();
  c_pceVsDiag_xEqY->SetGrid();
  mg_pceVsDiag_xEqY->Draw("A");
  leg_pceVsDiag_xEqY->Draw();
  fOut->cd();
  mg_pceVsDiag_xEqY->Write();
  c_pceVsDiag_xEqY->Write();
  c_pceVsDiag_xEqY->SaveAs("Mems_PCE_vs_Diag_xEqY.png");
  delete leg_pceVsDiag_xEqY;
  delete c_pceVsDiag_xEqY;

  fOut->Write();
  fOut->Close();
  delete fOut;
}

//---------------------------------------------------------------------------------------
// Plot PCE for the sensor-count scan submitted by G4Macros/run_numsensors.sh.
//---------------------------------------------------------------------------------------
void NumSensors_PCEStudy()
{
  const int al = 100;
  const std::vector<int> numSensors = {1, 2, 5, 10, 20, 30, 40, 50,
                                       60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250, 275, 300, 330, 360, 390, 420};
  const std::string baseRunDir = "../../../../261006_run_numsensors";

  struct Config {
    std::string directory;
    std::string label;
    Color_t color;
  };
  const std::vector<Config> configs = {
    {"SQUAT_Al_Nb_SW", "Side-wall loss", kBlue + 1},
    {"SQUAT_Al_Nb_PF", "Polished-wall loss", kRed + 1},
  };

  // Resolution on energy absorbed by a single sensor [eV]: 100 meV to 1 eV in 200 meV steps
  std::vector<double> resAbsorbed_eV;
  for (double r = 0.1; r <= 1.0 + 1e-9; r += 0.2) resAbsorbed_eV.push_back(r);

  TFile* fOut = new TFile("PCE_vs_NumSensors.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create PCE_vs_NumSensors.root" << std::endl;
    return;
  }

  TCanvas* canvas = new TCanvas("c_pce_vs_numSensors",
                                "PCE vs chip area", 900, 700);
  TMultiGraph* multigraph = new TMultiGraph(
      "mg_pce_vs_numSensors",
      "Phonon Collection Efficiency vs Chip Area;Chip area [cm^{2}];PCE [%]");
  TLegend* legend = new TLegend(0.72, 0.80, 0.88, 0.88);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.025);

  // Same plot, but with number of sensors on the x-axis instead of chip area
  TCanvas* canvas_ns = new TCanvas("c_pce_vs_numSensors_bySensorCount",
                                   "PCE vs number of sensors", 900, 700);
  TMultiGraph* multigraph_ns = new TMultiGraph(
      "mg_pce_vs_numSensors_bySensorCount",
      "Phonon Collection Efficiency vs Number of Sensors;Number of sensors;PCE [%]");
  TLegend* legend_ns = new TLegend(0.72, 0.80, 0.88, 0.88);
  legend_ns->SetBorderSize(1);
  legend_ns->SetFillStyle(0);
  legend_ns->SetTextSize(0.025);

  for (const Config& config : configs) {
    std::vector<double> sensorValues;
    // Chip area [cm^2] = 1/numsensors, per RISQTutorialDetectorConstruction (geHalf = 0.5cm/sqrt(numsensors))
    std::vector<double> areaValues_cm2;
    std::vector<double> pceValues;

    for (int ns : numSensors) {
      const std::string tag = "Al" + std::to_string(al) + "_ns" + std::to_string(ns);
      const std::string runDir = baseRunDir + "/" + config.directory;
      const std::string primaryFilename = runDir + "/Primary_" + tag + ".txt";
      const std::string hitsFilename = runDir + "/Hits_" + tag + ".txt";

      std::ifstream primaryFile(primaryFilename.c_str());
      std::ifstream hitsFile(hitsFilename.c_str());
      if (!primaryFile.good() || !hitsFile.good()) {
        std::cerr << "Skipping ns=" << ns << " for " << config.label
                  << ": missing primary or hit file." << std::endl;
        continue;
      }

      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primaryFilename);
      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);

      double totalPrimaryEnergy_eV = 0.0;
      for (const auto& entry : primaryInfo) {
        totalPrimaryEnergy_eV += entry.second.energy_eV;
      }

      double totalHitEnergy_eV = 0.0;
      for (const auto& entry : hitInfo) {
        for (const Hit& hit : entry.second) {
          totalHitEnergy_eV += hit.eDep_eV;
        }
      }

      if (totalPrimaryEnergy_eV <= 0.0) {
        std::cerr << "Skipping ns=" << ns << " for " << config.label
                  << ": primary energy is zero." << std::endl;
        continue;
      }

      const double pcePercent = 100.0 * totalHitEnergy_eV / totalPrimaryEnergy_eV;
      sensorValues.push_back(ns);
      areaValues_cm2.push_back(1.0 / ns);
      pceValues.push_back(pcePercent);
      std::cout << config.label << "  ns=" << ns
                << "  PCE=" << pcePercent << " %" << std::endl;
    }

    if (sensorValues.empty()) continue;

    // Table of number-of-sensors vs PCE for this loss configuration
    {
      std::cout << "\n=== PCE vs Number of Sensors (" << config.label << ") ===\n"
                << "NumSensors\tPCE[%]\n";
      std::ofstream tableOut(TString::Format("PCE_table_%s.txt", config.directory.c_str()).Data());
      tableOut << "NumSensors\tPCE[%]\n";
      for (int iN = 0; iN < (int)sensorValues.size(); ++iN) {
        std::cout << (int)sensorValues[iN] << "\t\t" << pceValues[iN] << "\n";
        tableOut << (int)sensorValues[iN] << "\t" << pceValues[iN] << "\n";
      }
      std::cout << "========================================\n";
    }

    // One line per resolution-on-energy-absorbed value: x = number of sensors N,
    // y = resolution on deposited energy, where
    //   res_deposited(N) = res_absorbed_1sensor * sqrt(N) / PCE
    {
      const int resColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kCyan + 2, kViolet, kPink + 7, kAzure + 2};

      TMultiGraph* mg_resDeposited = new TMultiGraph(
          TString::Format("mg_resDeposited_%s", config.directory.c_str()),
          TString::Format("Resolution on Deposited Energy (%s);Chip area [cm^{2}];Resolution on deposited energy [eV]",
                          config.label.c_str()));

      TCanvas* c_res = new TCanvas(TString::Format("c_resDeposited_%s", config.directory.c_str()),
          TString::Format("Resolution on Deposited Energy: %s", config.label.c_str()), 900, 700);
      TLegend* leg_res = new TLegend(0.76, 0.62, 0.88, 0.88);
      leg_res->SetBorderSize(1);
      leg_res->SetFillStyle(0);
      leg_res->SetTextSize(0.022);

      for (int iR = 0; iR < (int)resAbsorbed_eV.size(); ++iR) {
        const double resAbsorbed_eV_val = resAbsorbed_eV[iR];

        std::vector<double> yResDeposited(sensorValues.size());
        for (int iN = 0; iN < (int)sensorValues.size(); ++iN) {
          const double N = sensorValues[iN];
          const double pceFraction = pceValues[iN] / 100.0;
          yResDeposited[iN] = (pceFraction > 0.0)
                                ? (resAbsorbed_eV_val * std::sqrt(N) / pceFraction)
                                : 0.0;
        }

        TString gName = TString::Format("g_resDeposited_%s_res%dmeV",
                                        config.directory.c_str(), (int)(resAbsorbed_eV_val * 1000.0));
        TGraph* g = new TGraph((int)areaValues_cm2.size(), areaValues_cm2.data(), yResDeposited.data());
        g->SetName(gName);
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);

        const int color = resColors[iR < 10 ? iR : 9];
        g->SetLineColor(color);
        g->SetMarkerColor(color);

        fOut->cd();
        g->Write();

        mg_resDeposited->Add(g, "LP");
        leg_res->AddEntry(g, TString::Format("%.0f meV", resAbsorbed_eV_val * 1000.0), "lp");

        // Highlight the lowest (best) resolution point on this line
        int minIdx = 0;
        for (int iN = 1; iN < (int)yResDeposited.size(); ++iN) {
          if (yResDeposited[iN] < yResDeposited[minIdx]) minIdx = iN;
        }
        TGraph* gMin = new TGraph(1, &areaValues_cm2[minIdx], &yResDeposited[minIdx]);
        gMin->SetName(TString::Format("%s_min", gName.Data()));
        gMin->SetMarkerStyle(29);
        gMin->SetMarkerSize(2.2);
        gMin->SetMarkerColor(color);
        fOut->cd();
        gMin->Write();
        mg_resDeposited->Add(gMin, "P");
      }

      c_res->SetGrid();
      c_res->SetLogx();
      c_res->SetLogy();
      mg_resDeposited->Draw("A");
      leg_res->Draw();
      fOut->cd();
      mg_resDeposited->Write();
      c_res->Write();
      c_res->SaveAs(TString::Format("ResDeposited_%s.png", config.directory.c_str()));

      delete leg_res;
      delete c_res;
    }

    // Same resolution plot, but with number of sensors on the x-axis instead of chip area
    {
      const int resColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kCyan + 2, kViolet, kPink + 7, kAzure + 2};

      TMultiGraph* mg_resDeposited_ns = new TMultiGraph(
          TString::Format("mg_resDeposited_bySensorCount_%s", config.directory.c_str()),
          TString::Format("Resolution on Deposited Energy (%s);Number of sensors;Resolution on deposited energy [eV]",
                          config.label.c_str()));

      TCanvas* c_res_ns = new TCanvas(TString::Format("c_resDeposited_bySensorCount_%s", config.directory.c_str()),
          TString::Format("Resolution on Deposited Energy: %s", config.label.c_str()), 900, 700);
      TLegend* leg_res_ns = new TLegend(0.76, 0.62, 0.88, 0.88);
      leg_res_ns->SetBorderSize(1);
      leg_res_ns->SetFillStyle(0);
      leg_res_ns->SetTextSize(0.022);

      for (int iR = 0; iR < (int)resAbsorbed_eV.size(); ++iR) {
        const double resAbsorbed_eV_val = resAbsorbed_eV[iR];

        std::vector<double> yResDeposited(sensorValues.size());
        for (int iN = 0; iN < (int)sensorValues.size(); ++iN) {
          const double N = sensorValues[iN];
          const double pceFraction = pceValues[iN] / 100.0;
          yResDeposited[iN] = (pceFraction > 0.0)
                                ? (resAbsorbed_eV_val * std::sqrt(N) / pceFraction)
                                : 0.0;
        }

        TString gName = TString::Format("g_resDeposited_bySensorCount_%s_res%dmeV",
                                        config.directory.c_str(), (int)(resAbsorbed_eV_val * 1000.0));
        TGraph* g = new TGraph((int)sensorValues.size(), sensorValues.data(), yResDeposited.data());
        g->SetName(gName);
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);

        const int color = resColors[iR < 10 ? iR : 9];
        g->SetLineColor(color);
        g->SetMarkerColor(color);

        fOut->cd();
        g->Write();

        mg_resDeposited_ns->Add(g, "LP");
        leg_res_ns->AddEntry(g, TString::Format("%.0f meV", resAbsorbed_eV_val * 1000.0), "lp");

        // Highlight the lowest (best) resolution point on this line
        int minIdx = 0;
        for (int iN = 1; iN < (int)yResDeposited.size(); ++iN) {
          if (yResDeposited[iN] < yResDeposited[minIdx]) minIdx = iN;
        }
        TGraph* gMin = new TGraph(1, &sensorValues[minIdx], &yResDeposited[minIdx]);
        gMin->SetName(TString::Format("%s_min", gName.Data()));
        gMin->SetMarkerStyle(29);
        gMin->SetMarkerSize(2.2);
        gMin->SetMarkerColor(color);
        fOut->cd();
        gMin->Write();
        mg_resDeposited_ns->Add(gMin, "P");
      }

      c_res_ns->SetGrid();
      c_res_ns->SetLogx();
      c_res_ns->SetLogy();
      mg_resDeposited_ns->Draw("A");
      leg_res_ns->Draw();
      fOut->cd();
      mg_resDeposited_ns->Write();
      c_res_ns->Write();
      c_res_ns->SaveAs(TString::Format("ResDeposited_bySensorCount_%s.png", config.directory.c_str()));

      delete leg_res_ns;
      delete c_res_ns;
    }

    TString graphName = TString::Format("g_pce_vs_numSensors_%s", config.directory.c_str());
    TGraph* graph = new TGraph((int)areaValues_cm2.size(), areaValues_cm2.data(), pceValues.data());
    graph->SetName(graphName);
    graph->SetLineColor(config.color);
    graph->SetMarkerColor(config.color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.9);
    fOut->cd();
    graph->Write();
    multigraph->Add(graph, "LP");
    legend->AddEntry(graph, config.label.c_str(), "lp");

    TString graphName_ns = TString::Format("g_pce_vs_numSensors_bySensorCount_%s", config.directory.c_str());
    TGraph* graph_ns = new TGraph((int)sensorValues.size(), sensorValues.data(), pceValues.data());
    graph_ns->SetName(graphName_ns);
    graph_ns->SetLineColor(config.color);
    graph_ns->SetMarkerColor(config.color);
    graph_ns->SetLineWidth(2);
    graph_ns->SetMarkerStyle(20);
    graph_ns->SetMarkerSize(0.9);
    fOut->cd();
    graph_ns->Write();
    multigraph_ns->Add(graph_ns, "LP");
    legend_ns->AddEntry(graph_ns, config.label.c_str(), "lp");
  }

  if (multigraph->GetListOfGraphs()->GetSize() == 0) {
    std::cerr << "Error: no valid num-sensor scan files were found." << std::endl;
    fOut->Close();
    delete legend;
    delete canvas;
    delete legend_ns;
    delete canvas_ns;
    delete fOut;
    return;
  }

  canvas->cd();
  canvas->SetGrid();
  multigraph->Draw("A");
  legend->Draw();
  fOut->cd();
  multigraph->Write();
  canvas->Write();
  canvas->SaveAs("PCE_vs_NumSensors.png");

  canvas_ns->cd();
  canvas_ns->SetGrid();
  multigraph_ns->Draw("A");
  legend_ns->Draw();
  fOut->cd();
  multigraph_ns->Write();
  canvas_ns->Write();
  canvas_ns->SaveAs("PCE_vs_NumSensors_bySensorCount.png");

  fOut->Write();
  fOut->Close();

  delete legend;
  delete canvas;
  delete legend_ns;
  delete canvas_ns;
  delete fOut;
}

//---------------------------------------------------------------------------------------
// Analyze the initial-phonon-energy scan submitted by G4Macros/run_energy.sh
// (run_energy.sh -> energy.sh -> pceStudy.mac). Produces, for the side-wall (SW) and
// polished-face (PF) loss configurations:
//   1. PCE vs initial phonon energy, SW and PF as separate lines on one plot.
//   2. Histogram of the phonon energy absorbed per event, one canvas per initial
//      energy with SW and PF overlaid.
//   3. Estimated net quasiparticle population vs time, one canvas per loss
//      configuration with one line per initial energy.
//   4. Net QP vs time at 3 meV for QP lifetimes of 200 us to 2 ms (200 us steps),
//      one canvas per loss configuration with one line per lifetime.
//   5. Rise time, fall time, peak height and average of the net QP curves vs initial
//      energy and vs QP lifetime (SW and PF overlaid), 8 canvases in total, each with a
//      matching .txt data file. Rise time is the time from t = 0 to the peak; fall time
//      is the time from the peak until the population returns to zero (last QP removed);
//      the average is the time average of Nqp from t = 0 to that last removal. All are
//      taken from the exact event list, not the plot grid.
// The Nqp(t) lines are rainbow-colored from red (lowest) to violet (highest).
// The QP model is the same as Ns_QuasiparticleAnalysis: nQP = 2*E/Delta_Al per hit,
// created at endT_ns and removed one QP at a time with exponential lifetimes.
//---------------------------------------------------------------------------------------
void Energy_PCEAndQPStudy(std::string baseRunDir = "../../../../261005_run_energy")
{
  const int al = 40;
  const int ns = 1;
  // Same strings as run_energy.sh, since they are part of the output filenames [eV]
  const std::vector<std::string> energyStrings = {"0.0002", "0.0003", "0.0004", "0.0005", 
    "0.0006", "0.0008", "0.001", "0.0015", "0.002", "0.0025", "0.003", "0.0035", "0.004", 
    "0.0045", "0.005", "0.0075", "0.01"};

  struct Config {
    std::string directory;
    std::string label;
    Color_t color;
  };
  const std::vector<Config> configs = {
    {"SQUAT_Al_Nb_SW", "Side-wall loss", kBlue + 1},
    {"SQUAT_Al_Nb_PF", "Polished-wall loss", kRed + 1},
  };

  const double alGap_eV = 0.34e-3;
  const double qpLifetimeMean_ns = 300000.0;  // 300 us
  const int nTimeSamples = 2000;  // Net-QP curves are sampled on a uniform time grid

  // Rainbow color for index i of n, running from red (i = 0) to violet (i = n - 1)
  auto rainbowColor = [](int i, int n) {
    const double hue = (n > 1) ? 270.0 * i / (n - 1) : 0.0;
    Float_t r, g, b;
    TColor::HSV2RGB(hue, 1.0, 1.0, r, g, b);
    return TColor::GetColor(r, g, b);
  };

  TRandom3 rng(0);

  // QP lifetime scan, made for one initial energy only: 200 us to 2 ms in 200 us steps
  const std::string lifetimeScanEnergy = "0.003";
  std::vector<double> scanLifetimes_us;
  for (int tau_us = 200; tau_us <= 2000; tau_us += 200) scanLifetimes_us.push_back(tau_us);

  struct QPEvent {
    double time_ns;
    double deltaQP;
  };
  struct QPCurve {
    std::vector<double> times, netQP;
    double peakQP = 0.0;
    double riseTime_ns = -1.0;  // t = 0 to the peak; negative if undefined
    double fallTime_ns = -1.0;  // peak to zero population; negative if undefined
    double avgQP = 0.0;  // time average over [0, last QP removed]
  };

  // Net QP population vs time for the given hits and mean lifetime: nQP = 2*E/Delta_Al
  // per hit is created at endT_ns, and each QP is removed after an exponential lifetime.
  // Creation events sort before removals on ties; the running total is sampled on a
  // uniform time grid. Peak height, rise time and fall time come from the exact event list.
  auto computeNetQP = [&](const std::map<int, std::vector<Hit> >& hits, double tau_ns) {
    std::vector<QPEvent> qpEvents;
    for (const auto& kv : hits) {
      for (const Hit& hit : kv.second) {
        if (hit.eDep_eV <= 0.0) continue;
        const int nQP = (int)std::round((hit.eDep_eV / alGap_eV) * 2.0);
        qpEvents.push_back({hit.endT_ns, (double)nQP});
        for (int iQP = 0; iQP < nQP; ++iQP) {
          qpEvents.push_back({hit.endT_ns + rng.Exp(tau_ns), -1.0});
        }
      }
    }

    QPCurve curve;
    if (qpEvents.empty()) return curve;
    std::sort(qpEvents.begin(), qpEvents.end(),
              [](const QPEvent& a, const QPEvent& b) {
                if (a.time_ns == b.time_ns) return a.deltaQP > b.deltaQP;
                return a.time_ns < b.time_ns;
              });

    double total = 0.0, integral = 0.0;
    size_t iPeak = 0;
    for (size_t i = 0; i < qpEvents.size(); ++i) {
      if (i > 0) integral += total * (qpEvents[i].time_ns - qpEvents[i - 1].time_ns);
      total += qpEvents[i].deltaQP;
      if (total > curve.peakQP) {
        curve.peakQP = total;
        iPeak = i;
      }
    }
    if (curve.peakQP > 0.0) {
      curve.riseTime_ns = qpEvents[iPeak].time_ns;
      curve.fallTime_ns = qpEvents.back().time_ns - qpEvents[iPeak].time_ns;
      curve.avgQP = integral / qpEvents.back().time_ns;
    }

    const double dt = qpEvents.back().time_ns / nTimeSamples;
    double runningQP = 0.0;
    size_t iEv = 0;
    for (int iT = 0; iT <= nTimeSamples; ++iT) {
      const double t = iT * dt;
      while (iEv < qpEvents.size() && qpEvents[iEv].time_ns <= t) {
        runningQP += qpEvents[iEv].deltaQP;
        if (runningQP < 0.0) runningQP = 0.0;
        ++iEv;
      }
      curve.times.push_back(t);
      curve.netQP.push_back(runningQP);
    }
    return curve;
  };

  TFile* fOut = new TFile("PCE_vs_Energy.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create PCE_vs_Energy.root" << std::endl;
    return;
  }

  TMultiGraph* mg_pce = new TMultiGraph("mg_pce_vs_energy",
      "Phonon Collection Efficiency vs Initial Phonon Energy;Initial phonon energy [meV];PCE [%]");
  TLegend* leg_pce = new TLegend(0.72, 0.80, 0.88, 0.88);
  leg_pce->SetBorderSize(1);
  leg_pce->SetFillStyle(0);
  leg_pce->SetTextSize(0.025);

  // Absorbed-energy histograms, indexed [iE][iC], for the overlay canvases below
  std::vector<std::vector<TH1F*> > h_absorbed(
      energyStrings.size(), std::vector<TH1F*>(configs.size(), nullptr));

  // PCE [%], indexed [iE][iC]; negative means the point is missing
  std::vector<std::vector<double> > pceTable(
      energyStrings.size(), std::vector<double>(configs.size(), -1.0));

  // Rise time [us], fall time [us], peak height and average of the net QP curves vs a
  // scan variable
  struct QPMetrics {
    std::vector<double> x, rise_us, fall_us, peak, avg;
  };
  std::vector<QPMetrics> energyMetrics(configs.size()), lifetimeMetrics(configs.size());

  for (int iC = 0; iC < (int)configs.size(); ++iC) {
    const Config& config = configs[iC];

    std::vector<double> energies_meV, pceValues, pceErrors;

    TCanvas* c_qp = new TCanvas(TString::Format("c_netQP_%s", config.directory.c_str()),
        TString::Format("Net QP vs Time: %s", config.label.c_str()), 900, 700);
    c_qp->SetLogy();
    TLegend* leg_qp = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_qp->SetBorderSize(1);
    leg_qp->SetFillStyle(0);
    bool firstQPDraw = true;
    std::vector<TGraph*> qpGraphs;
    std::vector<TGraph*> lifetimeGraphs;
    TLegend* leg_lt = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_lt->SetBorderSize(1);
    leg_lt->SetFillStyle(0);

    for (int iE = 0; iE < (int)energyStrings.size(); ++iE) {
      const double energy_eV = std::atof(energyStrings[iE].c_str());
      const std::string tag = "Al" + std::to_string(al) + "_ns" + std::to_string(ns)
                              + "_E" + energyStrings[iE] + "eV";
      const std::string runDir = baseRunDir + "/" + config.directory;
      const std::string primaryFilename = runDir + "/Primary_" + tag + ".txt";
      const std::string hitsFilename = runDir + "/Hits_" + tag + ".txt";

      std::ifstream primaryFile(primaryFilename.c_str());
      std::ifstream hitsFile(hitsFilename.c_str());
      if (!primaryFile.good() || !hitsFile.good()) {
        std::cerr << "Skipping E=" << energyStrings[iE] << " eV for " << config.label
                  << ": missing primary or hit file." << std::endl;
        continue;
      }

      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primaryFilename);
      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);

      double totalPrimaryEnergy_eV = 0.0;
      for (const auto& entry : primaryInfo) totalPrimaryEnergy_eV += entry.second.energy_eV;
      if (totalPrimaryEnergy_eV <= 0.0) {
        std::cerr << "Skipping E=" << energyStrings[iE] << " eV for " << config.label
                  << ": primary energy is zero." << std::endl;
        continue;
      }

      // Histogram of the total energy absorbed per event
      TString hName = TString::Format("h_absorbed_%s_E%seV", config.directory.c_str(),
                                      energyStrings[iE].c_str());
      TH1F* h = new TH1F(hName,
          TString::Format("Absorbed phonon energy, E_{init} = %.2f meV;Absorbed energy [eV];nEvents",
                          1e3 * energy_eV),
          100, 0, 1.05 * energy_eV);
      h->SetDirectory(fOut);
      h->SetStats(0);
      h->SetLineWidth(2);
      h->SetLineColorAlpha(config.color, 0.6);

      double totalHitEnergy_eV = 0.0;

      for (const auto& kv : hitInfo) {
        double eventAbsorbed_eV = 0.0;
        for (const Hit& hit : kv.second) {
          if (hit.eDep_eV <= 0.0) continue;
          eventAbsorbed_eV += hit.eDep_eV;
        }
        totalHitEnergy_eV += eventAbsorbed_eV;
        if (eventAbsorbed_eV > 0.0) h->Fill(eventAbsorbed_eV);
      }

      const double pce = totalHitEnergy_eV / totalPrimaryEnergy_eV;
      const double nPrimaries = (double)primaryInfo.size();
      const double pceErr = (pce > 0.0 && pce < 1.0)
                              ? 100.0 * std::sqrt(pce * (1.0 - pce) / nPrimaries) : 0.0;
      energies_meV.push_back(1e3 * energy_eV);
      pceValues.push_back(100.0 * pce);
      pceErrors.push_back(pceErr);
      pceTable[iE][iC] = 100.0 * pce;
      std::cout << config.label << "  E=" << energy_eV << " eV  PCE=" << 100.0 * pce
                << " %" << std::endl;

      fOut->cd();
      h->Write();
      h_absorbed[iE][iC] = h;

      // Same hits, many QP lifetimes, overlaid on one plot
      if (energyStrings[iE] == lifetimeScanEnergy) {
        for (int iL = 0; iL < (int)scanLifetimes_us.size(); ++iL) {
          const double tau_us = scanLifetimes_us[iL];
          const QPCurve curve = computeNetQP(hitInfo, 1e3 * tau_us);
          if (curve.times.empty()) continue;
          std::cout << "    tau=" << tau_us << " us  Peak QP=" << curve.peakQP
                    << "  Average QP=" << curve.avgQP
                    << "  rise=" << curve.riseTime_ns / 1e3 << " us  fall="
                    << curve.fallTime_ns / 1e3 << " us" << std::endl;
          if (curve.peakQP > 0.0 && curve.fallTime_ns >= 0.0) {
            lifetimeMetrics[iC].x.push_back(tau_us);
            lifetimeMetrics[iC].rise_us.push_back(curve.riseTime_ns / 1e3);
            lifetimeMetrics[iC].fall_us.push_back(curve.fallTime_ns / 1e3);
            lifetimeMetrics[iC].peak.push_back(curve.peakQP);
            lifetimeMetrics[iC].avg.push_back(curve.avgQP);
          }

          TGraph* gl = new TGraph((int)curve.times.size(), curve.times.data(),
                                  curve.netQP.data());
          gl->SetName(TString::Format("g_netQP_lifetimeScan_%s_tau%dus",
                                      config.directory.c_str(), (int)tau_us));
          gl->SetTitle(TString::Format(
              "Net QP vs Time, E = %.2f meV, %s;Time [ns];Net N_{QP}",
              1e3 * energy_eV, config.label.c_str()));
          gl->SetLineWidth(2);
          gl->SetLineColor(rainbowColor(iL, (int)scanLifetimes_us.size()));
          fOut->cd();
          gl->Write();
          lifetimeGraphs.push_back(gl);
          leg_lt->AddEntry(gl, TString::Format("#tau = %.0f #mus", tau_us), "l");
        }
      }

      const QPCurve curve = computeNetQP(hitInfo, qpLifetimeMean_ns);
      if (curve.times.empty()) continue;
      const std::vector<double>& times = curve.times;
      const std::vector<double>& netQP = curve.netQP;
      std::cout << "    Peak QP=" << curve.peakQP << "  Average QP=" << curve.avgQP
                << "  rise=" << curve.riseTime_ns / 1e3 << " us  fall="
                << curve.fallTime_ns / 1e3 << " us" << std::endl;
      if (curve.peakQP > 0.0 && curve.fallTime_ns >= 0.0) {
        energyMetrics[iC].x.push_back(1e3 * energy_eV);
        energyMetrics[iC].rise_us.push_back(curve.riseTime_ns / 1e3);
        energyMetrics[iC].fall_us.push_back(curve.fallTime_ns / 1e3);
        energyMetrics[iC].peak.push_back(curve.peakQP);
        energyMetrics[iC].avg.push_back(curve.avgQP);
      }

      TGraph* g = new TGraph((int)times.size(), times.data(), netQP.data());
      g->SetName(TString::Format("g_netQP_%s_E%seV", config.directory.c_str(),
                                 energyStrings[iE].c_str()));
      g->SetTitle(TString::Format("Net QP vs Time, %s;Time [ns];Net N_{QP}",
                                  config.label.c_str()));
      g->SetLineWidth(2);
      g->SetLineColor(rainbowColor(iE, (int)energyStrings.size()));
      fOut->cd();
      g->Write();
      qpGraphs.push_back(g);
      leg_qp->AddEntry(g, TString::Format("E = %.2f meV", 1e3 * energy_eV), "l");
    }

    // Net QP vs time, one line per initial energy
    if (!qpGraphs.empty()) {
      double xMax = 0.0, yMax = 0.0;
      for (TGraph* g : qpGraphs) {
        for (int ip = 0; ip < g->GetN(); ++ip) {
          double x, y;
          g->GetPoint(ip, x, y);
          if (x > xMax) xMax = x;
          if (y > yMax) yMax = y;
        }
      }
      for (TGraph* g : qpGraphs) {
        if (firstQPDraw) {
          g->Draw("AL");
          g->GetXaxis()->SetLimits(0.0, xMax);
          g->SetMinimum(1.0);
          g->SetMaximum(3.0 * yMax);
          firstQPDraw = false;
        } else {
          g->Draw("L SAME");
        }
      }
      leg_qp->Draw();
      fOut->cd();
      c_qp->Write();
      c_qp->SaveAs(TString::Format("netQP_vs_Time_%s.png", config.directory.c_str()));
    }
    delete leg_qp;
    delete c_qp;

    if (!lifetimeGraphs.empty()) {
      TCanvas* c_lt = new TCanvas(
          TString::Format("c_netQP_lifetimeScan_%s", config.directory.c_str()),
          TString::Format("Net QP vs Time, lifetime scan: %s", config.label.c_str()), 900, 700);
      c_lt->SetLogy();
      double xMax = 0.0, yMax = 0.0;
      for (TGraph* g : lifetimeGraphs) {
        for (int ip = 0; ip < g->GetN(); ++ip) {
          double x, y;
          g->GetPoint(ip, x, y);
          if (x > xMax) xMax = x;
          if (y > yMax) yMax = y;
        }
      }
      for (int iL = 0; iL < (int)lifetimeGraphs.size(); ++iL) {
        TGraph* g = lifetimeGraphs[iL];
        if (iL == 0) {
          g->Draw("AL");
          g->GetXaxis()->SetLimits(0.0, xMax);
          g->SetMinimum(1.0);
          g->SetMaximum(3.0 * yMax);
        } else {
          g->Draw("L SAME");
        }
      }
      leg_lt->Draw();
      fOut->cd();
      c_lt->Write();
      c_lt->SaveAs(TString::Format("netQP_lifetimeScan_E%seV_%s.png",
                                   lifetimeScanEnergy.c_str(), config.directory.c_str()));
      delete c_lt;
    }
    delete leg_lt;

    if (energies_meV.empty()) continue;

    std::vector<double> xErr(energies_meV.size(), 0.0);
    TGraphErrors* g_pce = new TGraphErrors((int)energies_meV.size(), energies_meV.data(),
                                           pceValues.data(), xErr.data(), pceErrors.data());
    g_pce->SetName(TString::Format("g_pce_vs_energy_%s", config.directory.c_str()));
    g_pce->SetLineColor(config.color);
    g_pce->SetMarkerColor(config.color);
    g_pce->SetFillColor(config.color);
    g_pce->SetLineWidth(2);
    g_pce->SetMarkerStyle(20);
    g_pce->SetMarkerSize(0.9);
    fOut->cd();
    g_pce->Write();
    mg_pce->Add(g_pce, "LP");
    leg_pce->AddEntry(g_pce, config.label.c_str(), "lp");
  }

  // Rise time, fall time, peak height and average vs initial energy and vs QP lifetime,
  // SW and PF overlaid on each canvas
  auto drawMetric = [&](const std::vector<QPMetrics>& metrics,
                        std::vector<double> QPMetrics::*field, const std::string& name,
                        const std::string& title, const std::string& xTitle,
                        const std::string& yTitle, bool logX) {
    TMultiGraph* mg = new TMultiGraph(("mg_" + name).c_str(),
        (title + ";" + xTitle + ";" + yTitle).c_str());
    TLegend* leg = new TLegend(0.15, 0.78, 0.40, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    for (int iC = 0; iC < (int)configs.size(); ++iC) {
      const QPMetrics& m = metrics[iC];
      if (m.x.empty()) continue;
      TGraph* g = new TGraph((int)m.x.size(), m.x.data(), (m.*field).data());
      g->SetName(("g_" + name + "_" + configs[iC].directory).c_str());
      g->SetLineColor(configs[iC].color);
      g->SetMarkerColor(configs[iC].color);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.9);
      mg->Add(g, "LP");
      leg->AddEntry(g, configs[iC].label.c_str(), "lp");
    }

    // Matching data table: one row per scan value, one column per configuration
    std::map<double, std::vector<double> > rows;
    for (int iC = 0; iC < (int)configs.size(); ++iC) {
      const QPMetrics& m = metrics[iC];
      for (size_t i = 0; i < m.x.size(); ++i) {
        std::vector<double>& row = rows.emplace(
            m.x[i], std::vector<double>(configs.size(), -1.0)).first->second;
        row[iC] = (m.*field)[i];
      }
    }
    auto plainText = [](std::string str) {
      for (const auto& sub : std::vector<std::pair<std::string, std::string> >{
               {"#mu", "u"}, {"N_{QP}", "N_QP"}}) {
        for (size_t pos = str.find(sub.first); pos != std::string::npos;
             pos = str.find(sub.first, pos + sub.second.size())) {
          str.replace(pos, sub.first.size(), sub.second);
        }
      }
      return str;
    };
    std::ofstream dataFile((name + ".txt").c_str());
    if (!dataFile.good()) {
      std::cerr << "Error: could not create " << name << ".txt" << std::endl;
    } else {
      dataFile << "# " << plainText(title) << "\n# "
               << TString::Format("%-16s", plainText(xTitle).c_str());
      for (const Config& config : configs) {
        dataFile << TString::Format("%-34s",
            (config.directory + " " + plainText(yTitle)).c_str());
      }
      dataFile << "\n";
      for (const auto& row : rows) {
        dataFile << "  " << TString::Format("%-16.6g", row.first);
        for (double value : row.second) {
          if (value < 0.0) dataFile << TString::Format("%-34s", "NA");
          else dataFile << TString::Format("%-34.6g", value);
        }
        dataFile << "\n";
      }
    }

    if (mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 0) {
      TCanvas* c = new TCanvas(("c_" + name).c_str(), title.c_str(), 900, 700);
      c->SetGrid();
      if (logX) c->SetLogx();
      mg->Draw("A");
      leg->Draw();
      fOut->cd();
      mg->Write();
      c->Write();
      c->SaveAs((name + ".png").c_str());
      delete c;
    }
    delete leg;
  };

  drawMetric(energyMetrics, &QPMetrics::rise_us, "QPRiseTime_vs_Energy",
             "QP rise time (0 to peak) vs initial phonon energy",
             "Initial phonon energy [meV]", "Rise time [#mus]", true);
  drawMetric(energyMetrics, &QPMetrics::fall_us, "QPFallTime_vs_Energy",
             "QP fall time (peak to zero) vs initial phonon energy",
             "Initial phonon energy [meV]", "Fall time [#mus]", true);
  drawMetric(energyMetrics, &QPMetrics::peak, "QPPeak_vs_Energy",
             "Peak N_{QP} vs initial phonon energy",
             "Initial phonon energy [meV]", "Peak N_{QP}", true);
  drawMetric(lifetimeMetrics, &QPMetrics::rise_us, "QPRiseTime_vs_Lifetime",
             "QP rise time (0 to peak) vs QP lifetime",
             "QP lifetime [#mus]", "Rise time [#mus]", false);
  drawMetric(lifetimeMetrics, &QPMetrics::fall_us, "QPFallTime_vs_Lifetime",
             "QP fall time (peak to zero) vs QP lifetime",
             "QP lifetime [#mus]", "Fall time [#mus]", false);
  drawMetric(lifetimeMetrics, &QPMetrics::peak, "QPPeak_vs_Lifetime",
             "Peak N_{QP} vs QP lifetime",
             "QP lifetime [#mus]", "Peak N_{QP}", false);
  drawMetric(energyMetrics, &QPMetrics::avg, "QPAvg_vs_Energy",
             "Average N_{QP} (0 to last QP removed) vs initial phonon energy",
             "Initial phonon energy [meV]", "Average N_{QP}", true);
  drawMetric(lifetimeMetrics, &QPMetrics::avg, "QPAvg_vs_Lifetime",
             "Average N_{QP} (0 to last QP removed) vs QP lifetime",
             "QP lifetime [#mus]", "Average N_{QP}", false);

  // PCE vs initial energy as a text table, one row per energy
  {
    std::ofstream tableFile("PCE_vs_Energy.txt");
    if (!tableFile.good()) {
      std::cerr << "Error: could not create PCE_vs_Energy.txt" << std::endl;
    } else {
      tableFile << "# PCE vs initial phonon energy\n";
      tableFile << "# " << TString::Format("%-14s", "E_init[meV]");
      for (const Config& config : configs) {
        tableFile << TString::Format("%-24s", (config.directory + "_PCE[%]").c_str());
      }
      tableFile << "\n";
      for (int iE = 0; iE < (int)energyStrings.size(); ++iE) {
        tableFile << "  " << TString::Format("%-14.4g", 1e3 * std::atof(energyStrings[iE].c_str()));
        for (int iC = 0; iC < (int)configs.size(); ++iC) {
          if (pceTable[iE][iC] < 0.0) {
            tableFile << TString::Format("%-24s", "NA");
          } else {
            tableFile << TString::Format("%-24.4f", pceTable[iE][iC]);
          }
        }
        tableFile << "\n";
      }
    }
  }

  // PCE vs initial energy, SW and PF on the same plot
  if (mg_pce->GetListOfGraphs() && mg_pce->GetListOfGraphs()->GetSize() > 0) {
    TCanvas* c_pce = new TCanvas("c_pce_vs_energy", "PCE vs initial phonon energy", 900, 700);
    c_pce->SetGrid();
    c_pce->SetLogx();
    mg_pce->Draw("A");
    leg_pce->Draw();
    fOut->cd();
    mg_pce->Write();
    c_pce->Write();
    c_pce->SaveAs("PCE_vs_Energy.png");
    delete c_pce;
  } else {
    std::cerr << "Error: no valid energy scan files were found in " << baseRunDir << std::endl;
  }

  // Absorbed-energy histograms: one canvas per initial energy, SW and PF overlaid,
  // with the Al gap multiples marked as in the other eDep plots
  for (int iE = 0; iE < (int)energyStrings.size(); ++iE) {
    double maxY = 0.0;
    for (TH1F* h : h_absorbed[iE]) {
      if (h && h->GetMaximum() > maxY) maxY = h->GetMaximum();
    }
    if (maxY <= 0.0) continue;

    const double energy_eV = std::atof(energyStrings[iE].c_str());
    TCanvas* c = new TCanvas(TString::Format("c_absorbed_E%seV", energyStrings[iE].c_str()),
        TString::Format("Absorbed energy, E=%.2f meV", 1e3 * energy_eV), 900, 700);
    c->SetLogy();
    TLegend* leg = new TLegend(0.68, 0.70, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);

    bool first = true;
    for (int iC = 0; iC < (int)configs.size(); ++iC) {
      TH1F* h = h_absorbed[iE][iC];
      if (!h) continue;
      h->SetMaximum(1.15 * maxY);
      h->SetMinimum(0.5);
      h->Draw(first ? "HIST" : "HIST SAME");
      first = false;
      leg->AddEntry(h, configs[iC].label.c_str(), "l");
    }

    std::vector<TLine*> gapLines;
    for (int k = 1; k * alGap_eV <= 1.05 * energy_eV; ++k) {
      TLine* line = new TLine(k * alGap_eV, 0, k * alGap_eV, 1.15 * maxY);
      line->SetLineColor(kGreen + 2);
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw();
      gapLines.push_back(line);
    }
    if (!gapLines.empty()) leg->AddEntry(gapLines[0], "Al gap", "l");
    leg->Draw();

    fOut->cd();
    c->Write();
    c->SaveAs(TString::Format("h_absorbed_E%seV.png", energyStrings[iE].c_str()));
    delete leg;
    delete c;
  }

  fOut->Write();
  fOut->Close();
  delete fOut;
}


//---------------------------------------------------------------------------------------
// Same energy-scan files and QP model as Energy_PCEAndQPStudy, but the phonons are injected
// at a constant average rate (Poisson process) over duration_s instead of all at t = 0.
// Each of the nPhonons primaries is given a start time drawn uniformly in [0, duration_s]
// (a Poisson process conditioned on nPhonons events), and that offset is added to the
// endT_ns of all of its hits. Phonons do not interact, so this is exact for a
// superposition of independent events. If a run has a different number of primaries than
// nPhonons, events are re-used (sampled with replacement). Produces, for the side-wall (SW)
// and polished-face (PF) configurations:
//   1. Net QP vs time, one canvas per configuration, one line per initial energy.
//   2. Net QP vs time at 3 meV for QP lifetimes of 200 us to 2 ms, one canvas per
//      configuration, one line per lifetime.
//   3. Time-averaged net QP over [0, duration_s] vs initial energy and vs QP lifetime,
//      SW and PF overlaid, each with a matching .txt file.
// Nqp(t) lines are rainbow-colored from red (lowest) to violet (highest). Output files
// carry the "ConstRate" tag so they do not overwrite the t = 0 study.
//---------------------------------------------------------------------------------------
void Energy_ConstantRateQPStudy(std::string baseRunDir = "../../../../261005_run_energy",
                                int nPhonons = 100000, double duration_s = 1.0)
{
  const int al = 40;
  const int ns = 1;
  // Same strings as run_energy.sh, since they are part of the output filenames [eV]
  const std::vector<std::string> energyStrings = {"0.0002", "0.0003", "0.0004", "0.0005",
    "0.0006", "0.0008", "0.001", "0.0015", "0.002", "0.0025", "0.003", "0.0035", "0.004",
    "0.0045", "0.005", "0.0075", "0.01"};

  struct Config {
    std::string directory;
    std::string label;
    Color_t color;
  };
  const std::vector<Config> configs = {
    {"SQUAT_Al_Nb_SW", "Side-wall loss", kBlue + 1},
    {"SQUAT_Al_Nb_PF", "Polished-wall loss", kRed + 1},
  };

  const double alGap_eV = 0.34e-3;
  const double qpLifetimeMean_ns = 300000.0;  // 300 us
  const double duration_ns = duration_s * 1e9;
  const int nTimeSamples = 5000;  // Net-QP curves are sampled on a uniform time grid

  TRandom3 rng(0);

  // QP lifetime scan, made for one initial energy only: 200 us to 2 ms in 200 us steps
  const std::string lifetimeScanEnergy = "0.003";
  std::vector<double> scanLifetimes_us;
  for (int tau_us = 200; tau_us <= 2000; tau_us += 200) scanLifetimes_us.push_back(tau_us);

  // Rainbow color for index i of n, running from red (i = 0) to violet (i = n - 1)
  auto rainbowColor = [](int i, int n) {
    const double hue = (n > 1) ? 270.0 * i / (n - 1) : 0.0;
    Float_t r, g, b;
    TColor::HSV2RGB(hue, 1.0, 1.0, r, g, b);
    return TColor::GetColor(r, g, b);
  };

  struct Injection {
    int eventID;
    double t0_ns;
  };
  struct QPEvent {
    double time_ns;
    double deltaQP;
  };
  struct QPCurve {
    std::vector<double> times_s, netQP;
    double avgQP = 0.0;  // time average of the net QP over [0, duration_s]
  };
  // Time-averaged net QP of the curves vs a scan variable
  struct QPMetrics {
    std::vector<double> x, avg;
  };

  // Start time for each of the nPhonons phonons, constant average rate over duration_ns
  auto makeInjections = [&](const std::map<int, PrimaryInfo>& primaryInfo) {
    std::vector<int> ids;
    for (const auto& kv : primaryInfo) ids.push_back(kv.first);
    std::vector<Injection> injections;
    injections.reserve(nPhonons);
    for (int i = 0; i < nPhonons; ++i) {
      const int id = ((int)ids.size() == nPhonons) ? ids[i] : ids[rng.Integer(ids.size())];
      injections.push_back({id, rng.Uniform(duration_ns)});
    }
    return injections;
  };

  // Net QP population vs time: nQP = 2*E/Delta_Al per hit is created at endT_ns plus the
  // phonon's injection time, and each QP is removed after an exponential lifetime.
  // Creation events sort before removals on ties. The average over the injection window
  // [0, duration_ns] comes from the exact event list; the plotted curve is sampled on a
  // uniform time grid.
  auto computeNetQP = [&](const std::map<int, std::vector<Hit> >& hits,
                          const std::vector<Injection>& injections, double tau_ns) {
    std::vector<QPEvent> qpEvents;
    for (const Injection& inj : injections) {
      const auto it = hits.find(inj.eventID);
      if (it == hits.end()) continue;
      for (const Hit& hit : it->second) {
        if (hit.eDep_eV <= 0.0) continue;
        const int nQP = (int)std::round((hit.eDep_eV / alGap_eV) * 2.0);
        const double tCreate_ns = hit.endT_ns + inj.t0_ns;
        qpEvents.push_back({tCreate_ns, (double)nQP});
        for (int iQP = 0; iQP < nQP; ++iQP) {
          qpEvents.push_back({tCreate_ns + rng.Exp(tau_ns), -1.0});
        }
      }
    }

    QPCurve curve;
    if (qpEvents.empty()) return curve;
    std::sort(qpEvents.begin(), qpEvents.end(),
              [](const QPEvent& a, const QPEvent& b) {
                if (a.time_ns == b.time_ns) return a.deltaQP > b.deltaQP;
                return a.time_ns < b.time_ns;
              });

    // Integrate the step-function net QP over [0, duration_ns]
    double total = 0.0, integral = 0.0, tPrev = 0.0;
    for (const QPEvent& ev : qpEvents) {
      const double tClipped = std::min(ev.time_ns, duration_ns);
      if (tClipped > tPrev) {
        integral += total * (tClipped - tPrev);
        tPrev = tClipped;
      }
      total += ev.deltaQP;
    }
    if (tPrev < duration_ns) integral += total * (duration_ns - tPrev);
    curve.avgQP = integral / duration_ns;

    const double dt = qpEvents.back().time_ns / nTimeSamples;
    double runningQP = 0.0;
    size_t iEv = 0;
    for (int iT = 0; iT <= nTimeSamples; ++iT) {
      const double t = iT * dt;
      while (iEv < qpEvents.size() && qpEvents[iEv].time_ns <= t) {
        runningQP += qpEvents[iEv].deltaQP;
        if (runningQP < 0.0) runningQP = 0.0;
        ++iEv;
      }
      curve.times_s.push_back(t * 1e-9);
      curve.netQP.push_back(runningQP);
    }
    return curve;
  };

  TFile* fOut = new TFile("QP_ConstantRate.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "Error: could not create QP_ConstantRate.root" << std::endl;
    return;
  }

  std::vector<QPMetrics> energyMetrics(configs.size()), lifetimeMetrics(configs.size());

  // Draws one Nqp(t) canvas from a set of graphs and legend, log y from 1 up
  auto drawNetQPCanvas = [&](TCanvas* canvas, const std::vector<TGraph*>& graphs,
                             TLegend* legend, const std::string& pngName) {
    double xMax = 0.0, yMax = 0.0;
    for (TGraph* g : graphs) {
      for (int ip = 0; ip < g->GetN(); ++ip) {
        double x, y;
        g->GetPoint(ip, x, y);
        if (x > xMax) xMax = x;
        if (y > yMax) yMax = y;
      }
    }
    canvas->cd();
    for (size_t i = 0; i < graphs.size(); ++i) {
      if (i == 0) {
        graphs[i]->Draw("AL");
        graphs[i]->GetXaxis()->SetLimits(0.0, xMax);
        graphs[i]->SetMinimum(1.0);
        graphs[i]->SetMaximum(3.0 * yMax);
      } else {
        graphs[i]->Draw("L SAME");
      }
    }
    legend->Draw();
    fOut->cd();
    canvas->Write();
    canvas->SaveAs(pngName.c_str());
  };

  for (int iC = 0; iC < (int)configs.size(); ++iC) {
    const Config& config = configs[iC];

    std::vector<TGraph*> qpGraphs, lifetimeGraphs;
    TLegend* leg_qp = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_qp->SetBorderSize(1);
    leg_qp->SetFillStyle(0);
    TLegend* leg_lt = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg_lt->SetBorderSize(1);
    leg_lt->SetFillStyle(0);

    for (int iE = 0; iE < (int)energyStrings.size(); ++iE) {
      const double energy_eV = std::atof(energyStrings[iE].c_str());
      const std::string tag = "Al" + std::to_string(al) + "_ns" + std::to_string(ns)
                              + "_E" + energyStrings[iE] + "eV";
      const std::string runDir = baseRunDir + "/" + config.directory;
      const std::string primaryFilename = runDir + "/Primary_" + tag + ".txt";
      const std::string hitsFilename = runDir + "/Hits_" + tag + ".txt";

      std::ifstream primaryFile(primaryFilename.c_str());
      std::ifstream hitsFile(hitsFilename.c_str());
      if (!primaryFile.good() || !hitsFile.good()) {
        std::cerr << "Skipping E=" << energyStrings[iE] << " eV for " << config.label
                  << ": missing primary or hit file." << std::endl;
        continue;
      }

      const std::map<int, PrimaryInfo> primaryInfo =
          ParsePrimaryTextFileForPrimaries(primaryFilename);
      const std::map<int, std::vector<Hit> > hitInfo =
          ParseHitTextFileForHits(hitsFilename);
      if (primaryInfo.empty()) {
        std::cerr << "Skipping E=" << energyStrings[iE] << " eV for " << config.label
                  << ": no primaries." << std::endl;
        continue;
      }

      // The same injection times are used for every lifetime, so only tau changes
      const std::vector<Injection> injections = makeInjections(primaryInfo);

      if (energyStrings[iE] == lifetimeScanEnergy) {
        for (int iL = 0; iL < (int)scanLifetimes_us.size(); ++iL) {
          const double tau_us = scanLifetimes_us[iL];
          const QPCurve curve = computeNetQP(hitInfo, injections, 1e3 * tau_us);
          if (curve.times_s.empty()) continue;
          std::cout << config.label << "  tau=" << tau_us << " us  Average QP="
                    << curve.avgQP << std::endl;
          lifetimeMetrics[iC].x.push_back(tau_us);
          lifetimeMetrics[iC].avg.push_back(curve.avgQP);

          TGraph* gl = new TGraph((int)curve.times_s.size(), curve.times_s.data(),
                                  curve.netQP.data());
          gl->SetName(TString::Format("g_netQP_ConstRate_lifetimeScan_%s_tau%dus",
                                      config.directory.c_str(), (int)tau_us));
          gl->SetTitle(TString::Format(
              "Net QP vs Time, %d phonons over %.3g s, E = %.2f meV, %s;Time [s];Net N_{QP}",
              nPhonons, duration_s, 1e3 * energy_eV, config.label.c_str()));
          gl->SetLineWidth(2);
          gl->SetLineColor(rainbowColor(iL, (int)scanLifetimes_us.size()));
          fOut->cd();
          gl->Write();
          lifetimeGraphs.push_back(gl);
          leg_lt->AddEntry(gl, TString::Format("#tau = %.0f #mus", tau_us), "l");
        }
      }

      const QPCurve curve = computeNetQP(hitInfo, injections, qpLifetimeMean_ns);
      if (curve.times_s.empty()) continue;
      std::cout << config.label << "  E=" << energy_eV << " eV  Average QP="
                << curve.avgQP << std::endl;
      energyMetrics[iC].x.push_back(1e3 * energy_eV);
      energyMetrics[iC].avg.push_back(curve.avgQP);

      TGraph* g = new TGraph((int)curve.times_s.size(), curve.times_s.data(),
                             curve.netQP.data());
      g->SetName(TString::Format("g_netQP_ConstRate_%s_E%seV", config.directory.c_str(),
                                 energyStrings[iE].c_str()));
      g->SetTitle(TString::Format("Net QP vs Time, %d phonons over %.3g s, %s;Time [s];Net N_{QP}",
                                  nPhonons, duration_s, config.label.c_str()));
      g->SetLineWidth(2);
      g->SetLineColor(rainbowColor(iE, (int)energyStrings.size()));
      fOut->cd();
      g->Write();
      qpGraphs.push_back(g);
      leg_qp->AddEntry(g, TString::Format("E = %.2f meV", 1e3 * energy_eV), "l");
    }

    if (!qpGraphs.empty()) {
      TCanvas* c_qp = new TCanvas(
          TString::Format("c_netQP_ConstRate_%s", config.directory.c_str()),
          TString::Format("Net QP vs Time, constant rate: %s", config.label.c_str()), 900, 700);
      c_qp->SetLogy();
      drawNetQPCanvas(c_qp, qpGraphs, leg_qp,
                      "netQP_ConstantRate_vs_Time_" + config.directory + ".png");
      delete c_qp;
    }
    delete leg_qp;

    if (!lifetimeGraphs.empty()) {
      TCanvas* c_lt = new TCanvas(
          TString::Format("c_netQP_ConstRate_lifetimeScan_%s", config.directory.c_str()),
          TString::Format("Net QP vs Time, constant rate, lifetime scan: %s",
                          config.label.c_str()), 900, 700);
      c_lt->SetLogy();
      drawNetQPCanvas(c_lt, lifetimeGraphs, leg_lt,
                      "netQP_ConstantRate_lifetimeScan_E" + lifetimeScanEnergy + "eV_"
                      + config.directory + ".png");
      delete c_lt;
    }
    delete leg_lt;
  }

  // Average net QP vs initial energy and vs QP lifetime, SW and PF overlaid on each
  // canvas, with a matching text file
  auto drawMetric = [&](const std::vector<QPMetrics>& metrics,
                        std::vector<double> QPMetrics::*field, const std::string& name,
                        const std::string& title, const std::string& xTitle,
                        const std::string& yTitle, const std::string& xTitlePlain,
                        const std::string& yTitlePlain, bool logX) {
    TMultiGraph* mg = new TMultiGraph(("mg_" + name).c_str(),
        (title + ";" + xTitle + ";" + yTitle).c_str());
    TLegend* leg = new TLegend(0.15, 0.78, 0.40, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    for (int iC = 0; iC < (int)configs.size(); ++iC) {
      const QPMetrics& m = metrics[iC];
      if (m.x.empty()) continue;
      TGraph* g = new TGraph((int)m.x.size(), m.x.data(), (m.*field).data());
      g->SetName(("g_" + name + "_" + configs[iC].directory).c_str());
      g->SetLineColor(configs[iC].color);
      g->SetMarkerColor(configs[iC].color);
      g->SetLineWidth(2);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.9);
      mg->Add(g, "LP");
      leg->AddEntry(g, configs[iC].label.c_str(), "lp");
    }

    // Matching data table: one row per scan value, one column per configuration
    std::map<double, std::vector<double> > rows;
    for (int iC = 0; iC < (int)configs.size(); ++iC) {
      const QPMetrics& m = metrics[iC];
      for (size_t i = 0; i < m.x.size(); ++i) {
        std::vector<double>& row = rows.emplace(
            m.x[i], std::vector<double>(configs.size(), -1.0)).first->second;
        row[iC] = (m.*field)[i];
      }
    }
    std::ofstream dataFile((name + ".txt").c_str());
    if (!dataFile.good()) {
      std::cerr << "Error: could not create " << name << ".txt" << std::endl;
    } else {
      dataFile << "# " << yTitlePlain << " vs " << xTitlePlain << ", " << nPhonons
               << " phonons at a constant rate over " << duration_s << " s\n# "
               << TString::Format("%-16s", xTitlePlain.c_str());
      for (const Config& config : configs) {
        dataFile << TString::Format("%-34s", (config.directory + " " + yTitlePlain).c_str());
      }
      dataFile << "\n";
      for (const auto& row : rows) {
        dataFile << "  " << TString::Format("%-16.6g", row.first);
        for (double value : row.second) {
          if (value < 0.0) dataFile << TString::Format("%-34s", "NA");
          else dataFile << TString::Format("%-34.6g", value);
        }
        dataFile << "\n";
      }
    }

    if (mg->GetListOfGraphs() && mg->GetListOfGraphs()->GetSize() > 0) {
      TCanvas* c = new TCanvas(("c_" + name).c_str(), title.c_str(), 900, 700);
      c->SetGrid();
      if (logX) c->SetLogx();
      mg->Draw("A");
      leg->Draw();
      fOut->cd();
      mg->Write();
      c->Write();
      c->SaveAs((name + ".png").c_str());
      delete c;
    }
    delete leg;
  };

  const std::string energyX = "Initial phonon energy [meV]";
  const std::string energyXPlain = "E_init[meV]";
  const std::string lifetimeX = "QP lifetime [#mus]";
  const std::string lifetimeXPlain = "tau[us]";
  drawMetric(energyMetrics, &QPMetrics::avg, "QPConstRateAvg_vs_Energy",
             "Average N_{QP}, constant-rate injection, vs initial phonon energy",
             energyX, "Average N_{QP}", energyXPlain, "average_N_QP", true);
  drawMetric(lifetimeMetrics, &QPMetrics::avg, "QPConstRateAvg_vs_Lifetime",
             "Average N_{QP}, constant-rate injection, vs QP lifetime",
             lifetimeX, "Average N_{QP}", lifetimeXPlain, "average_N_QP", false);

  fOut->Write();
  fOut->Close();
  delete fOut;
}

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



