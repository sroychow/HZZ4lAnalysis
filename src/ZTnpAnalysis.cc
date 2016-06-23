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
#include <utility> 
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
#include "TH3.h"
#include "TProfile.h"

#include "AnaUtil.h"
#include "HZZ4lUtil.h"
#include "ZTnpAnalysis.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
ZTnpAnalysis::ZTnpAnalysis()
  : PhysicsObjSelector(),
    treeFname_("zTnPisoTree.root")
{
}
// ----------
// Destructor
// ----------
ZTnpAnalysis::~ZTnpAnalysis() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ZTnpAnalysis::beginJob() 
{
  outTreeFile_ = TFile::Open(treeFname_.c_str(),"recreate");
  outTreeFile_->cd();
  outTree_ = new TTree("tree","muTnP tree for iso");
  outTree_->Branch("TnP_eta",&zcand.TnP_eta);
  outTree_->Branch("TnP_phi",&zcand.TnP_phi);
  outTree_->Branch("TnP_mass",&zcand.TnP_mass);
  outTree_->Branch("TnP_hasFSR",&zcand.TnP_hasFSR);
  outTree_->Branch("TnP_mll",&zcand.TnP_mll);
  outTree_->Branch("TnP_l1_pdgId",&zcand.TnP_l1_pdgId);
  outTree_->Branch("TnP_l1_pt",&zcand.TnP_l1_pt);
  outTree_->Branch("TnP_l1_eta",&zcand.TnP_l1_eta);
  outTree_->Branch("TnP_l1_phi",&zcand.TnP_l1_phi);
  outTree_->Branch("TnP_l1_mass",&zcand.TnP_l1_mass);
  outTree_->Branch("TnP_l1_charge",&zcand.TnP_l1_charge);
  outTree_->Branch("TnP_l1_tightId",&zcand.TnP_l1_tightId);
  outTree_->Branch("TnP_l1_looseId",&zcand.TnP_l1_looseId);
  outTree_->Branch("TnP_l1_dxy",&zcand.TnP_l1_dxy);
  outTree_->Branch("TnP_l1_dz",&zcand.TnP_l1_dz);
  outTree_->Branch("TnP_l1_edxy",&zcand.TnP_l1_edxy);
  outTree_->Branch("TnP_l1_edz",&zcand.TnP_l1_edz);
  outTree_->Branch("TnP_l1_ip3d",&zcand.TnP_l1_ip3d);
  outTree_->Branch("TnP_l1_sip3d",&zcand.TnP_l1_sip3d);
  outTree_->Branch("TnP_l1_ptErr",&zcand.TnP_l1_ptErr);
  outTree_->Branch("TnP_l1_lostHits",&zcand.TnP_l1_lostHits);
  outTree_->Branch("TnP_l1_trackerLayers",&zcand.TnP_l1_trackerLayers);
  outTree_->Branch("TnP_l1_pixelLayers",&zcand.TnP_l1_pixelLayers);
  outTree_->Branch("TnP_l1_etaSc",&zcand.TnP_l1_etaSc);
  outTree_->Branch("TnP_l1_isGap",&zcand.TnP_l1_isGap);
  outTree_->Branch("TnP_l1_r9",&zcand.TnP_l1_r9);
  outTree_->Branch("TnP_l1_convVeto",&zcand.TnP_l1_convVeto);
  outTree_->Branch("TnP_l1_mvaIdSpring15",&zcand.TnP_l1_mvaIdSpring15);
  outTree_->Branch("TnP_l1_relIsoAfterFSR",&zcand.TnP_l1_relIsoAfterFSR);
  outTree_->Branch("TnP_l1_chargedHadIso03",&zcand.TnP_l1_chargedHadIso03);
  outTree_->Branch("TnP_l1_hasOwnFSR",&zcand.TnP_l1_hasOwnFSR);
  outTree_->Branch("TnP_l1_hlt1L",&zcand.TnP_l1_hlt1L);
  outTree_->Branch("TnP_l1_p4WithFSR_pt",&zcand.TnP_l1_p4WithFSR_pt);
  outTree_->Branch("TnP_l1_p4WithFSR_eta",&zcand.TnP_l1_p4WithFSR_eta);
  outTree_->Branch("TnP_l1_p4WithFSR_phi",&zcand.TnP_l1_p4WithFSR_phi);
  outTree_->Branch("TnP_l1_p4WithFSR_mass",&zcand.TnP_l1_p4WithFSR_mass);
  outTree_->Branch("TnP_l2_pdgId",&zcand.TnP_l2_pdgId);
  outTree_->Branch("TnP_l2_pt",&zcand.TnP_l2_pt);
  outTree_->Branch("TnP_l2_eta",&zcand.TnP_l2_eta);
  outTree_->Branch("TnP_l2_phi",&zcand.TnP_l2_phi);
  outTree_->Branch("TnP_l2_mass",&zcand.TnP_l2_mass);
  outTree_->Branch("TnP_l2_charge",&zcand.TnP_l2_charge);
  outTree_->Branch("TnP_l2_tightId",&zcand.TnP_l2_tightId);
  outTree_->Branch("TnP_l2_looseId",&zcand.TnP_l2_looseId);
  outTree_->Branch("TnP_l2_dxy",&zcand.TnP_l2_dxy);
  outTree_->Branch("TnP_l2_dz",&zcand.TnP_l2_dz);
  outTree_->Branch("TnP_l2_edxy",&zcand.TnP_l2_edxy);
  outTree_->Branch("TnP_l2_edz",&zcand.TnP_l2_edz);
  outTree_->Branch("TnP_l2_ip3d",&zcand.TnP_l2_ip3d);
  outTree_->Branch("TnP_l2_sip3d",&zcand.TnP_l2_sip3d);
  outTree_->Branch("TnP_l2_ptErr",&zcand.TnP_l2_ptErr);
  outTree_->Branch("TnP_l2_lostHits",&zcand.TnP_l2_lostHits);
  outTree_->Branch("TnP_l2_trackerLayers",&zcand.TnP_l2_trackerLayers);
  outTree_->Branch("TnP_l2_pixelLayers",&zcand.TnP_l2_pixelLayers);
  outTree_->Branch("TnP_l2_etaSc",&zcand.TnP_l2_etaSc);
  outTree_->Branch("TnP_l2_isGap",&zcand.TnP_l2_isGap);
  outTree_->Branch("TnP_l2_r9",&zcand.TnP_l2_r9);
  outTree_->Branch("TnP_l2_convVeto",&zcand.TnP_l2_convVeto);
  outTree_->Branch("TnP_l2_mvaIdSpring15",&zcand.TnP_l2_mvaIdSpring15);
  outTree_->Branch("TnP_l2_relIsoAfterFSR",&zcand.TnP_l2_relIsoAfterFSR);
  outTree_->Branch("TnP_l2_chargedHadIso03",&zcand.TnP_l2_chargedHadIso03);
  outTree_->Branch("TnP_l2_hasOwnFSR",&zcand.TnP_l2_hasOwnFSR);
  outTree_->Branch("TnP_l2_hlt1L",&zcand.TnP_l2_hlt1L);
  outTree_->Branch("TnP_l2_p4WithFSR_pt",&zcand.TnP_l2_p4WithFSR_pt);
  outTree_->Branch("TnP_l2_p4WithFSR_eta",&zcand.TnP_l2_p4WithFSR_eta);
  outTree_->Branch("TnP_l2_p4WithFSR_phi",&zcand.TnP_l2_p4WithFSR_phi);
  outTree_->Branch("TnP_l2_p4WithFSR_mass",&zcand.TnP_l2_p4WithFSR_mass);
  
  
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();
  
  histf()->cd();
  histf()->mkdir("ZTnpAnalysis");
  bookHistograms();
  
  return true;
}
// ---------------
// Book histograms
// ---------------
void ZTnpAnalysis::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("ZTnpAnalysis");
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);
  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodmuon", "Number of Good muons(with selection cuts) per event", 20, 0, 20);
  new TH1D("nGoodelectron", "Number of Good electrons(with selection cuts) per event", 20, 0, 20);
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void ZTnpAnalysis::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
  tightMuon_.clear();
  fsrPVec_.clear();
  tObj_.clear();
}
// -------------------
// The main event loop
// -------------------
void ZTnpAnalysis::eventLoop()
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = std::max(10000, nEvents()/1000);
  
  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int fevt = (firstEvent() > -1) ? firstEvent() : 0;
  int levt = (lastEvent() > -1) ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent(); // reset tree variables
    clearLists();
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
    const Event& evt = eventColl()->at(0);
    
    histf()->cd();
    
    // For data or for MC without pileup
    puevWt_ = 1;
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
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
    if (ev%nPrint == 0 || firstEvent() > -1)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << run
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> "
	   << " Events proc. " << setw(8) << ((firstEvent() > -1) ? ev - firstEvent() : ev)
	   << endl;
    
    histf()->cd();
    histf()->cd("ZTnpAnalysis");
    AnaUtil::fillHist1D("evtCutFlow", 0, puevWt_);
    
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    double ngoodVtx = vtxList_.size();
    
    // crucial
    histf()->cd();
    histf()->cd("ZTnpAnalysis");

    AnaUtil::fillHist1D("evtCutFlow", 1, puevWt_);
    // is event triggered?
    if (useTrigger() && !isTriggered(true, false)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, puevWt_);
    
    // at least 1 good PV
    if (ngoodVtx < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, puevWt_);
    
    // main analysis object selection
    findObjects(puevWt_);
    histf()->cd("ZTnpAnalysis");
    tightMuon_ = getTightMuPhotonPairList();
    if(tightMuon_.size() < 2)         continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);   
    findTriggerObjectInfo(tObj_);
    //dumpTriggerObjectInfo(tObj_);
    fsrPVec_ = getFSRPhotonVec();
    getZtnPpair();
  }
  // Analysis is over
  endJob();
}
void ZTnpAnalysis::getZtnPpair() {
  for (unsigned int i = 0; i < tightMuon_.size(); ++i) {
    const auto& tagmu = tightMuon_[i].first;
    TLorentzVector tagP4 = HZZ4lUtil::getP4(tagmu);
    TLorentzVector tagP4wfsr = tagP4;
    bool taghasfsr = false;
    if(!tightMuon_[i].second.empty()) { 
      tagP4wfsr += HZZ4lUtil::getP4(tightMuon_[i].second.at(0));
      taghasfsr = true;
    }
    if(tagP4.Pt() <= 20.)    continue;  
    double tagiso = HZZ4lUtil::computeMuonReliso(tagmu, fsrPVec_, 0.01, 0.3);   
    if (tagiso >= 0.35) continue; // it is already scaled by pT
    int matchedidx;
    bool trigMatched = matchTriggerObject(tObj_, tagP4, "HLT_IsoMu20_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoTkMu20_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoMu22_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoTkMu22_v", 0, 30,matchedidx) < 0.02;
    if (!trigMatched) continue; // it is already scaled by pT

    for (unsigned int j = i+1; j < tightMuon_.size(); ++j) {
      const auto& probemu = tightMuon_[j].first;
      // opposite charge
      if ((tagmu.charge + probemu.charge) != 0) continue; 
      //std::cout << "Tag & Probe pair found!!" << std::endl;
      TLorentzVector probeP4 = HZZ4lUtil::getP4(probemu); 
      TLorentzVector probeP4wfsr = probeP4;
      bool probehasfsr = false;
      if(!tightMuon_[j].second.empty()) {
        probeP4wfsr += HZZ4lUtil::getP4(tightMuon_[j].second.at(0));
        probehasfsr = true;
      }
      double probeiso = HZZ4lUtil::computeMuonReliso(probemu, fsrPVec_, 0.01, 0.3); 
      //compute TnP pair P4
      TLorentzVector tnpP4 = tagP4 + probeP4;
      TLorentzVector tnpP4wfsr = tagP4wfsr + probeP4wfsr;
      //Fill TnP pair property
      zcand.flavour = HZZ4lUtil::ZType::mumu;
      zcand.TnP_eta = tnpP4.Eta();   
      zcand.TnP_phi = tnpP4.Phi();   
      zcand.TnP_mass = tnpP4wfsr.M();   //with fsr
      zcand.TnP_hasFSR = taghasfsr || probehasfsr;   
      zcand.TnP_mll = tnpP4.M();   //without fsr
      //tag property
      zcand.TnP_l1_pdgId = (tagmu.charge > 0) ? 13 : -13;   
      zcand.TnP_l1_pt = tagP4.Pt();   
      zcand.TnP_l1_eta = tagP4.Eta();   
      zcand.TnP_l1_phi = tagP4.Phi();   
      zcand.TnP_l1_mass = tagP4.M();   
      zcand.TnP_l1_charge = tagmu.charge;   
      zcand.TnP_l1_tightId = 1;   
      zcand.TnP_l1_looseId = 1;   
      zcand.TnP_l1_dxy = tagmu.dxyPV;   
      zcand.TnP_l1_dz  = tagmu.dzPV;   
      zcand.TnP_l1_edxy = -1;   //not saved in our tree 
      zcand.TnP_l1_edz = -1;//not saved in out tree
      zcand.TnP_l1_ip3d = tagmu.dB3D;   //check
      zcand. TnP_l1_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      zcand.TnP_l1_ptErr = -1;//NS   
      zcand.TnP_l1_lostHits = 0;//NS   
      zcand.TnP_l1_trackerLayers = tagmu.trkHits;   
      zcand.TnP_l1_pixelLayers = tagmu.pixHits;   
      //zcand.TnP_l1_etaSc = -1;   
      //zcand.TnP_l1_isGap = -1;   
      //zcand.TnP_l1_r9 = -1;   
      //zcand.TnP_l1_convVeto = -1;   
      //zcand.TnP_l1_mvaIdSpring15 = -1;;   
      zcand.TnP_l1_relIsoAfterFSR = tagiso;   
      zcand.TnP_l1_chargedHadIso03 = tagmu.pfChargedHadIsoR03;   
      zcand.TnP_l1_hasOwnFSR  = 1; //doubt
      zcand.TnP_l1_hlt1L = 1; //doubt  
      zcand.TnP_l1_p4WithFSR_pt = tagP4wfsr.Pt();   
      zcand.TnP_l1_p4WithFSR_eta = tagP4wfsr.Eta();   
      zcand.TnP_l1_p4WithFSR_phi = tagP4wfsr.Phi();   
      zcand.TnP_l1_p4WithFSR_mass = tagP4wfsr.M();   
      //probe property
      zcand.TnP_l2_pdgId = (probemu.charge > 0) ? 13 : -13;   
      zcand.TnP_l2_pt = probeP4.Pt();   
      zcand.TnP_l2_eta = probeP4.Eta();   
      zcand.TnP_l2_phi = probeP4.Phi();   
      zcand.TnP_l2_mass = probeP4.M();   
      zcand.TnP_l2_charge = probemu.charge;   
      zcand.TnP_l2_tightId = 1;   
      zcand.TnP_l2_looseId = 1;   
      zcand.TnP_l2_dxy = probemu.dxyPV;   
      zcand.TnP_l2_dz  = probemu.dzPV;   
      zcand.TnP_l2_edxy = -1;   //correct? 
      zcand.TnP_l2_edz = -1;
      zcand.TnP_l2_ip3d = tagmu.dB3D;   //check
      zcand. TnP_l2_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      zcand.TnP_l2_ptErr = -1;   
      zcand.TnP_l2_lostHits = 0;   
      zcand.TnP_l2_trackerLayers = probemu.trkHits;   
      zcand.TnP_l2_pixelLayers = probemu.pixHits;   
      //zcand.TnP_l2_etaSc = -1;   
      //zcand.TnP_l2_isGap = -1;   
      //zcand.TnP_l2_r9 = -1;   
      //zcand.TnP_l2_convVeto = -1;   
      //zcand.TnP_l2_mvaIdSpring15 = -1;;   
      zcand.TnP_l2_relIsoAfterFSR = probeiso;   
      zcand.TnP_l2_chargedHadIso03 = probemu.pfChargedHadIsoR03;   
      zcand.TnP_l2_hasOwnFSR  = 1; //doubt
      zcand.TnP_l2_hlt1L = matchTriggerObject(tObj_, probeP4, "HLT_IsoMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoMu22_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu22_v", 0, 30,matchedidx) < 0.02; 
      zcand.TnP_l2_p4WithFSR_pt = probeP4wfsr.Pt();   
      zcand.TnP_l2_p4WithFSR_eta = probeP4wfsr.Eta();   
      zcand.TnP_l2_p4WithFSR_phi = probeP4wfsr.Phi();   
      zcand.TnP_l2_p4WithFSR_mass = probeP4wfsr.M();  
      outTree_->Fill();
    }
  }
}
void ZTnpAnalysis::endJob() {
  closeFiles();
  histf()->cd();
  histf()->Write();
  histf()->Close();
  outTree_->AutoSave();
  //outTreeFile_->Write();
  outTreeFile_->Close();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool ZTnpAnalysis::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;
  
  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), std::ios::in);    
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
    if (key == "outputTreeFilename")
      treeFname_ = value;
    tokens.clear();
  }
  // Close the file
  fin.close();
  
  printJob();
  
  return true;
}
void ZTnpAnalysis::printJob(ostream& os) const
{
  AnaBase::printJob(os);
}
