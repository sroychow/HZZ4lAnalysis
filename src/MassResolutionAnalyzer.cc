#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <typeinfo>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "HZZ4lUtil.h"
#include "AnaUtil.h"
#include "MassResolutionAnalyzer.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;
using std::setw;
using std::setiosflags;
using std::resetiosflags;
using std::ios;
using std::stoi;
using namespace vhtm;
// -----------
// Constructor
// -----------
MassResolutionAnalyzer::MassResolutionAnalyzer()
  : PhysicsObjSelector()
{
}
// ----------
// Destructor
// ----------
MassResolutionAnalyzer::~MassResolutionAnalyzer() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MassResolutionAnalyzer::beginJob() 
{ 
  //open the output tree
  tfout_ = TFile::Open(treeFileName_.c_str(),"RECREATE");
  tfout_->cd();
  massRtree_ = new TTree("mRestree","mass resolution study tree");
  massRtree_->Branch("massZ",&mr.m);
  massRtree_->Branch("massZErr", &mr.merr);
  massRtree_->Branch("eta1", &mr.eta1);
  massRtree_->Branch("eta2",&mr.eta2);
  massRtree_->Branch("pT1",&mr.pT1);
  massRtree_->Branch("pT2",&mr.pT2);
  massRtree_->Branch("weight",&mr.w);

  PhysicsObjSelector::beginJob();
  // Open the output ROOT file
  histf()->cd();
  PhysicsObjSelector::bookHistograms();
  bookHistograms();
  return true;
  mr.m = -1.;
  mr.merr = -1.;
  mr.eta1 = -1.;
  mr.eta2 = -1.;
  mr.pT1 = -1.;
  mr.pT2 = -1.;
  mr.w = 1.;
}
// ---------------
// Book Common histograms
// ---------------
void MassResolutionAnalyzer::bookHistograms() 
{
  histf()->cd();
  histf()->mkdir("Event");
  histf()->cd("Event");
  new TH1D("nGoodvtx","Number of Vertices",50,0,50);
  //Rho Values
  new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
  histf()->mkdir("MassResolution");
  histf()->cd("MassResolution");
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);
  new TH1F("nMu","number of Muons",20,-0.5,19.5);
  new TH1F("nEle","number of Electrons",20,-0.5,19.5);
  new TH1F("muPt", " Muon pt",200,-0.5,199.5);
  new TH1F("elePt"," Electron pt",200,-0.5,199.5);
}
//-----------------------------------------------------------------------------------
// -------------------
// The main event loop
// -------------------
void MassResolutionAnalyzer::clearLists() {
 PhysicsObjSelector::clear();
 tightMuVec_.clear();
 tightEleVec_.clear();
 vtxList.clear();
}
void MassResolutionAnalyzer::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = max(10000, nEvents()/1000);
  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  //long int nRecoMuons=0,nGoodRecoMuons=0;
  std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
  for (int ev = 0; ev < nEvents(); ++ev) {

    clearEvent();
    clearLists();
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    //For data or for MC without pileup
    puevWt_ = 1;
    /*
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
    */
    // Show status of the run
    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << chain()->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

    op.verbose = (logOption() >> 1 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
   
    histf()->cd();
    histf()->cd("MassResolution");
    AnaUtil::fillHist1D("evtCutFlow", 0, puevWt_);
    if( vtxList.empty() )             continue;
    AnaUtil::fillHist1D("evtCutFlow", 1, puevWt_);
    int ngoodVtx=vtxList.size();
    //to check triggered and non-triggered cases
    if (useTrigger() && !isTriggered(true, false))      continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, puevWt_);

    // main analysis object selection
    findObjects(puevWt_);

    // access selected objects (Tight, isolated leptons)
    //const auto& elePhotonPairVec = getTightIsoElePhotonPairList();
    tightMuVec_ = getTightIsoMuList();
    histf()->cd("MassResolution");
    AnaUtil::fillHist1D("nMu", tightMuVec_.size(), puevWt_);
    if(tightMuVec_.size() < 2)                      continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, puevWt_);
    getZmmMassInfo();
  }
  // Analysis is over
  endJob();
}
//from tight iso muons for Z
void MassResolutionAnalyzer::getZmmMassInfo() {
  for (unsigned int i = 0; i < tightMuVec_.size(); ++i) {
    const auto& ip = tightMuVec_[i];
    TLorentzVector lep1P4, lep1P4_p; 
    lep1P4.SetPtEtaPhiE(ip.smearedKpt,ip.eta,ip.phi,ip.energy);
    lep1P4_p.SetPtEtaPhiE(ip.smearedKpt + ip.smearedKpterr,ip.eta,ip.phi,ip.energy);
    for (unsigned int j = i+1; j < tightMuVec_.size(); ++j) {
      const auto& jp = tightMuVec_[j];
      //std::cout << "Ip charge=" << ip.charge << "\tJp charge=" << jp.charge << std::endl;
      // opposite charge
      if ((ip.charge + jp.charge) != 0) continue;
      TLorentzVector lep2P4, lep2P4_p;
      lep2P4.SetPtEtaPhiE(jp.smearedKpt,jp.eta,jp.phi,jp.energy);
      lep2P4_p.SetPtEtaPhiE(jp.smearedKpt + jp.smearedKpterr,jp.eta,jp.phi,jp.energy);
      double moriginal = (lep1P4 + lep2P4).M();
      double dm1 = (lep1P4_p + lep2P4).M() - moriginal;
      double dm2 = (lep1P4 + lep2P4_p).M() - moriginal;
      double delm = TMath::Sqrt(TMath::Power(dm1,2) + TMath::Power(dm2,2));
      mr.m = moriginal;
      mr.merr = delm;
      mr.eta1 = ip.eta;
      mr.eta2 = jp.phi;
      mr.pT1 = ip.smearedKpt;
      mr.pT2 = jp.smearedKpt;
      massRtree_->Fill();
    }
  } 
}

bool MassResolutionAnalyzer::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    string key = tokens[0];
    string value = tokens[1];
    if(key == "masstreeFile")
      treeFileName_ = value;
    tokens.clear();
  }
  // Close the file
  fin.close();
  printJob();
  return true;
}
void MassResolutionAnalyzer::endJob() {
  closeFiles();
  histf()->cd();
  histf()->Write();
  histf()->Close();
  std::cout << "Mtree" << massRtree_ << std::endl;
  std::cout << "tfout" << tfout_ << std::endl;
  //massRtree_->AutoSave();
  tfout_->Write();
  tfout_->Close();
  delete histf();
  //delete massRtree_; 
}
