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
  nTnP = 0;
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
  outTree_->Branch("run",&run,"run/I");
  outTree_->Branch("lumi",&lumi,"lumi/I");
  outTree_->Branch("event",&event,"event/I");
  outTree_->Branch("isData",&isData,"isData/I");
  outTree_->Branch("nTnP",&nTnP,"nTnP/I");
  
  outTree_->Branch("TnP_pt",TnP_pt,"TnP_pt[nTnP]/F");
  outTree_->Branch("TnP_eta",TnP_eta,"TnP_eta[nTnP]/F");
  outTree_->Branch("TnP_phi",TnP_phi,"TnP_phi[nTnP]/F");
  outTree_->Branch("TnP_mass",TnP_mass,"TnP_mass[nTnP]/F");
  outTree_->Branch("TnP_hasFSR",TnP_hasFSR,"TnP_hasFSR[nTnP]/I");
  outTree_->Branch("TnP_mll",TnP_mll,"TnP_mll[nTnP]/F");
  outTree_->Branch("TnP_l1_pdgId",TnP_l1_pdgId,"TnP_l1_pdgId[nTnP]/I");
  outTree_->Branch("TnP_l1_pt",TnP_l1_pt,"TnP_l1_pt[nTnP]/F");
  outTree_->Branch("TnP_l1_eta",TnP_l1_eta,"TnP_l1_eta[nTnP]/F");
  outTree_->Branch("TnP_l1_phi",TnP_l1_phi,"TnP_l1_phi[nTnP]/F");
  outTree_->Branch("TnP_l1_mass",TnP_l1_mass,"TnP_l1_mass[nTnP]/F");
  outTree_->Branch("TnP_l1_charge",TnP_l1_charge,"TnP_l1_charge[nTnP]/I");
  outTree_->Branch("TnP_l1_tightId",TnP_l1_tightId,"TnP_l1_tightId[nTnP]/I");
  outTree_->Branch("TnP_l1_looseId",TnP_l1_looseId,"TnP_l1_looseId[nTnP]/I");
  outTree_->Branch("TnP_l1_dxy",TnP_l1_dxy,"TnP_l1_dxy[nTnP]/F");
  outTree_->Branch("TnP_l1_dz",TnP_l1_dz,"TnP_l1_dz[nTnP]/F");
  outTree_->Branch("TnP_l1_edxy",TnP_l1_edxy,"TnP_l1_edxy[nTnP]/F");
  outTree_->Branch("TnP_l1_edz",TnP_l1_edz,"TnP_l1_edz[nTnP]/F");
  outTree_->Branch("TnP_l1_ip3d",TnP_l1_ip3d,"TnP_l1_ip3d[nTnP]/F");
  outTree_->Branch("TnP_l1_sip3d",TnP_l1_sip3d,"TnP_l1_sip3d[nTnP]/F");
  outTree_->Branch("TnP_l1_ptErr",TnP_l1_ptErr,"TnP_l1_ptErr[nTnP]/F");
  outTree_->Branch("TnP_l1_lostHits",TnP_l1_lostHits,"TnP_l1_lostHits[nTnP]/I");
  outTree_->Branch("TnP_l1_trackerLayers",TnP_l1_trackerLayers,"TnP_l1_trackerLayer[nTnP]/I");
  outTree_->Branch("TnP_l1_pixelLayers",TnP_l1_pixelLayers,"TnP_l1_pixelLayers[nTnP]/I");
  outTree_->Branch("TnP_l1_etaSc",TnP_l1_etaSc,"TnP_l1_etaSc[nTnP]/F");
  outTree_->Branch("TnP_l1_isGap",TnP_l1_isGap,"TnP_l1_isGap[nTnP]/F");
  outTree_->Branch("TnP_l1_r9",TnP_l1_r9,"TnP_l1_r9[nTnP]/F");
  outTree_->Branch("TnP_l1_convVeto",TnP_l1_convVeto,"TnP_l1_convVeto[nTnP]/F");
  outTree_->Branch("TnP_l1_mvaIdSpring15",TnP_l1_mvaIdSpring15,"TnP_l1_mvaIdSpring15[nTnP]/F");
  outTree_->Branch("TnP_l1_relIsoAfterFSR",TnP_l1_relIsoAfterFSR,"TnP_l1_relIsoAfterFSR[nTnP]/F");
  outTree_->Branch("TnP_l1_chargedHadIso03",TnP_l1_chargedHadIso03,"TnP_l1_chargedHadIso03[nTnP]/F");
  outTree_->Branch("TnP_l1_hasOwnFSR",TnP_l1_hasOwnFSR,"TnP_l1_hasOwnFSR[nTnP]/I");
  if(isMC()) {
    outTree_->Branch("TnP_l1_mcMatchId",TnP_l1_mcMatchId,"TnP_l1_mcMatchId[nTnP]/I");
    outTree_->Branch("TnP_l1_mcMatchAny",TnP_l1_mcMatchAny,"TnP_l1_mcMatchAny[nTnP]/I");
    outTree_->Branch("TnP_l1_mcPt",TnP_l1_mcPt,"TnP_l1_mcPt[nTnP]/F");
    outTree_->Branch("TnP_l1_mcPt1",TnP_l1_mcPt1,"TnP_l1_mcPt1[nTnP]/F");
  }

  outTree_->Branch("TnP_l1_hlt1L",TnP_l1_hlt1L,"TnP_l1_hlt1L[nTnP]/I");
  outTree_->Branch("TnP_l1_p4WithFSR_pt",TnP_l1_p4WithFSR_pt,"TnP_l1_p4WithFSR_pt[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_eta",TnP_l1_p4WithFSR_eta,"TnP_l1_p4WithFSR_eta[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_phi",TnP_l1_p4WithFSR_phi,"TnP_l1_p4WithFSR_phi[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_mass",TnP_l1_p4WithFSR_mass,"TnP_l1_p4WithFSR_mass[nTnP]/F");

  outTree_->Branch("TnP_l2_pdgId",TnP_l2_pdgId,"TnP_l2_pdgId[nTnP]/I");
  outTree_->Branch("TnP_l2_pt",TnP_l2_pt,"TnP_l2_pt[nTnP]/F");
  outTree_->Branch("TnP_l2_eta",TnP_l2_eta,"TnP_l2_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_phi",TnP_l2_phi,"TnP_l2_phi[nTnP]/F");
  outTree_->Branch("TnP_l2_mass",TnP_l2_mass,"TnP_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_charge",TnP_l2_charge,"TnP_l2_charge[nTnP]/I");
  outTree_->Branch("TnP_l2_tightId",TnP_l2_tightId,"TnP_l2_tightId[nTnP]/I");
  outTree_->Branch("TnP_l2_looseId",TnP_l2_looseId,"TnP_l2_looseId[nTnP]/I");
  outTree_->Branch("TnP_l2_dxy",TnP_l2_dxy,"TnP_l2_dxy[nTnP]/F");
  outTree_->Branch("TnP_l2_dz",TnP_l2_dz,"TnP_l2_dz[nTnP]/F");
  outTree_->Branch("TnP_l2_edxy",TnP_l2_edxy,"TnP_l2_edxy[nTnP]/F");
  outTree_->Branch("TnP_l2_edz",TnP_l2_edz,"TnP_l2_edz[nTnP]/F");
  outTree_->Branch("TnP_l2_ip3d",TnP_l2_ip3d,"TnP_l2_ip3d[nTnP]/F");
  outTree_->Branch("TnP_l2_sip3d",TnP_l2_sip3d,"TnP_l2_sip3d[nTnP]/F");
  outTree_->Branch("TnP_l2_ptErr",TnP_l2_ptErr,"TnP_l2_ptErr[nTnP]/F");
  outTree_->Branch("TnP_l2_lostHits",TnP_l2_lostHits,"TnP_l2_lostHits[nTnP]/I");
  outTree_->Branch("TnP_l2_trackerLayers",TnP_l2_trackerLayers,"TnP_l2_trackerLayers[nTnP]/I");
  outTree_->Branch("TnP_l2_pixelLayers",TnP_l2_pixelLayers,"TnP_l2_pixelLayers[nTnP]/I");
  outTree_->Branch("TnP_l2_etaSc",TnP_l2_etaSc,"TnP_l2_etaSc[nTnP]/F");
  outTree_->Branch("TnP_l2_isGap",TnP_l2_isGap,"TnP_l2_isGap[nTnP]/I");
  outTree_->Branch("TnP_l2_r9",TnP_l2_r9,"TnP_l2_r9[nTnP]/F");
  outTree_->Branch("TnP_l2_convVeto",TnP_l2_convVeto,"TnP_l2_convVeto[nTnP]/I");
  outTree_->Branch("TnP_l2_mvaIdSpring15",TnP_l2_mvaIdSpring15,"TnP_l2_mvaIdSpring15[nTnP]/F");
  outTree_->Branch("TnP_l2_relIsoAfterFSR",TnP_l2_relIsoAfterFSR,"TnP_l2_relIsoAfterFSR[nTnP]/F");
  outTree_->Branch("TnP_l2_chargedHadIso03",TnP_l2_chargedHadIso03,"TnP_l2_chargedHadIso03[nTnP]/F");
  outTree_->Branch("TnP_l2_hasOwnFSR",TnP_l2_hasOwnFSR,"TnP_l2_hasOwnFSR[nTnP]/I");
  if(isMC()) {
    outTree_->Branch("TnP_l2_mcMatchId",TnP_l2_mcMatchId,"TnP_l2_mcMatchId[nTnP]/I");
    outTree_->Branch("TnP_l2_mcMatchAny",TnP_l2_mcMatchAny,"TnP_l2_mcMatchAny[nTnP]/I");
    outTree_->Branch("TnP_l2_mcPt",TnP_l2_mcPt,"TnP_l2_mcPt[nTnP]/F");
    outTree_->Branch("TnP_l2_mcPt1",TnP_l2_mcPt1,"TnP_l2_mcPt1[nTnP]/F");
  }
  outTree_->Branch("TnP_l2_hlt1L",TnP_l2_hlt1L,"TnP_l2_hlt1L[nTnP]/I");
  outTree_->Branch("TnP_l2_p4WithFSR_pt",TnP_l2_p4WithFSR_pt,"TnP_l2_p4WithFSR_pt[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_eta",TnP_l2_p4WithFSR_eta,"TnP_l2_p4WithFSR_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_phi",TnP_l2_p4WithFSR_phi,"TnP_l2_p4WithFSR_phi[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_mass",TnP_l2_p4WithFSR_mass,"TnP_l2_p4WithFSR_mass[nTnP]/F");
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
  genObj_.clear();
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
    ///std::cout << "P1\n";
    histf()->cd("ZTnpAnalysis");
    tightMuon_ = getTightMuPhotonPairList();
    if(tightMuon_.size() < 2)         continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);   
    findTriggerObjectInfo(tObj_);
    //std::cout << "P2\n";
    //dumpTriggerObjectInfo(tObj_);
    fsrPVec_ = getFSRPhotonVec();
    nTnP = 0;
    run = eventColl()->at(0).run;
    event = eventColl()->at(0).event;
    lumis = eventColl()->at(0).lumis;
    if(isMC())   {
      genObj_ = *genParticleColl();
      isData = 0;
    }
    else isData = 1;
    getZtnPpair();
    //std::cout << "P3\n";
    outTree_->Fill();
    //std::cout << "P4\n"; 
  }
  // Analysis is over
  endJob();
}
void ZTnpAnalysis::getZtnPpair() {
  int n = 0;
  std::vector<vhtm::ZtnP>  zvec;
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
      vhtm::ZtnP ztemp;
      ztemp.TnP_pt = tnpP4.Pt();
      //Fill TnP pair property
      ztemp.TnP_pt = tnpP4.Pt();
      ztemp.TnP_eta = tnpP4.Eta();   
      ztemp.TnP_phi = tnpP4.Phi();   
      ztemp.TnP_mass = tnpP4wfsr.M();   //with fsr
      ztemp.TnP_hasFSR = taghasfsr || probehasfsr;   
      ztemp.TnP_mll = tnpP4.M();   //without fsr
      //tag property
      ztemp.TnP_l1_pdgId = (tagmu.charge > 0) ? -13 : 13;   
      ztemp.TnP_l1_pt = tagP4.Pt();   
      ztemp.TnP_l1_eta = tagP4.Eta();   
      ztemp.TnP_l1_phi = tagP4.Phi();   
      ztemp.TnP_l1_mass = tagP4.M();   
      ztemp.TnP_l1_charge = tagmu.charge;   
      ztemp.TnP_l1_tightId = 1;   
      ztemp.TnP_l1_looseId = 1;   
      ztemp.TnP_l1_dxy = tagmu.dxyPV;   
      ztemp.TnP_l1_dz  = tagmu.dzPV;   
      ztemp.TnP_l1_edxy = -1;   //not saved in our tree 
      ztemp.TnP_l1_edz = -1;//not saved in out tree
      ztemp.TnP_l1_ip3d = tagmu.dB3D;   //check
      ztemp.TnP_l1_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      ztemp.TnP_l1_ptErr = -1;//NS   
      ztemp.TnP_l1_lostHits = 0;//NS   
      ztemp.TnP_l1_trackerLayers = tagmu.trkHits;   
      ztemp.TnP_l1_pixelLayers = tagmu.pixHits;   
      ztemp.TnP_l1_etaSc = -1;   
      ztemp.TnP_l1_isGap = -1;   
      ztemp.TnP_l1_r9 = -1;   
      ztemp.TnP_l1_convVeto = -1;   
      ztemp.TnP_l1_mvaIdSpring15 = -1;;   
      ztemp.TnP_l1_relIsoAfterFSR = tagiso;   
      ztemp.TnP_l1_chargedHadIso03 = tagmu.pfChargedHadIsoR03;   
      ztemp.TnP_l1_hasOwnFSR  = taghasfsr; //doubt
      if(isMC()) {
        TLorentzVector genP4;
        int mid = 0;
        int genId = GenLevelMatching(tagP4, genObj_, genP4, mid);
        if(genId == -1) {
          //std::cout << "NOGEN" << std::endl;
          ztemp.TnP_l1_mcMatchId = 0;
          ztemp.TnP_l1_mcMatchAny = 0;
          ztemp.TnP_l1_mcPt = 0.;
          ztemp.TnP_l1_mcPt1 = 0.;
        } else {
          //if(abs(mid) != 23)   std::cout << "MIDMatched=" << mid << std::endl;
          ztemp.TnP_l1_mcMatchId = mid;
          ztemp.TnP_l1_mcMatchAny = 1;
          ztemp.TnP_l1_mcPt = genP4.Pt();
          ztemp.TnP_l1_mcPt1 = genP4.Pt();
        }
      }
      ztemp.TnP_l1_hlt1L = 1; //doubt  t 
      ztemp.TnP_l1_p4WithFSR_pt = tagP4wfsr.Pt();   
      ztemp.TnP_l1_p4WithFSR_eta = tagP4wfsr.Eta();   
      ztemp.TnP_l1_p4WithFSR_phi = tagP4wfsr.Phi();   
      ztemp.TnP_l1_p4WithFSR_mass = tagP4wfsr.M();   
      //probe property
      ztemp.TnP_l2_pdgId = (probemu.charge > 0) ? -13 : 13;   
      ztemp.TnP_l2_pt = probeP4.Pt();   
      ztemp.TnP_l2_eta = probeP4.Eta();   
      ztemp.TnP_l2_phi = probeP4.Phi();   
      ztemp.TnP_l2_mass = probeP4.M();   
      ztemp.TnP_l2_charge = probemu.charge;   
      ztemp.TnP_l2_tightId = 1;   
      ztemp.TnP_l2_looseId = 1;   
      ztemp.TnP_l2_dxy = probemu.dxyPV;   
      ztemp.TnP_l2_dz  = probemu.dzPV;   
      ztemp.TnP_l2_edxy = -1;   //correct? 
      ztemp.TnP_l2_edz = -1;
      ztemp.TnP_l2_ip3d = tagmu.dB3D;   //check
      ztemp.TnP_l2_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      ztemp.TnP_l2_ptErr = -1;   
      ztemp.TnP_l2_lostHits = 0;   
      ztemp.TnP_l2_trackerLayers = probemu.trkHits;   
      ztemp.TnP_l2_pixelLayers = probemu.pixHits;   
      ztemp.TnP_l2_etaSc = -1;   
      ztemp.TnP_l2_isGap = -1;   
      ztemp.TnP_l2_r9 = -1;   
      ztemp.TnP_l2_convVeto = -1;   
      ztemp.TnP_l2_mvaIdSpring15 = -1;;   
      ztemp.TnP_l2_relIsoAfterFSR = probeiso;
      //if(probeiso > 0.35) std::cout << "FAILINGPROBE" << std::endl;      
      ztemp.TnP_l2_chargedHadIso03 = probemu.pfChargedHadIsoR03;   
      ztemp.TnP_l2_hasOwnFSR  = probehasfsr; //doubt
      if(isMC()) {
        TLorentzVector genP4;
        int mid = 0;
        int genId = GenLevelMatching(probeP4, genObj_, genP4, mid);
        if(genId == -1) {
          //std::cout << "NOGEN" << std::endl;
          ztemp.TnP_l2_mcMatchId = 0;
          ztemp.TnP_l2_mcMatchAny = 0;
          ztemp.TnP_l2_mcPt = 0.;
          ztemp.TnP_l2_mcPt1 = 0.;
        } else {
          //if(abs(mid) != 23)   std::cout << "MIDMatched=" << mid << std::endl;
          ztemp.TnP_l2_mcMatchId = mid;
          ztemp.TnP_l2_mcMatchAny = 1;
          ztemp.TnP_l2_mcPt = genP4.Pt();
          ztemp.TnP_l2_mcPt1 = genP4.Pt();
        }
      }
      ztemp.TnP_l2_hlt1L = matchTriggerObject(tObj_, probeP4, "HLT_IsoMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoMu22_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu22_v", 0, 30,matchedidx) < 0.02; 
      ztemp.TnP_l2_p4WithFSR_pt = probeP4wfsr.Pt();   
      ztemp.TnP_l2_p4WithFSR_eta = probeP4wfsr.Eta();   
      ztemp.TnP_l2_p4WithFSR_phi = probeP4wfsr.Phi();   
      ztemp.TnP_l2_p4WithFSR_mass = probeP4wfsr.M(); 
      n++;
      zvec.push_back(ztemp);
    }
  }
  nTnP = n;
  for(int i = 0; i<n; i++) {
    if(i>=10)   continue;
    TnP_pt[i] = zvec[i].TnP_pt;
    
    TnP_eta[i] = zvec[i].TnP_eta;   
    TnP_phi[i] = zvec[i].TnP_phi;  
    TnP_mass[i] = zvec[i].TnP_mass;   
    TnP_hasFSR[i] = zvec[i].TnP_hasFSR;   
    TnP_mll[i] = zvec[i].TnP_mll;   
    TnP_l1_pdgId[i] = zvec[i].TnP_l1_pdgId;   
    TnP_l1_pt[i] = zvec[i].TnP_l1_pt;   
    TnP_l1_eta[i] = zvec[i].TnP_l1_eta;   
    TnP_l1_phi[i] = zvec[i].TnP_l1_phi;   
    TnP_l1_mass[i] = zvec[i].TnP_l1_mass;   
    TnP_l1_charge[i] = zvec[i].TnP_l1_charge;   
    TnP_l1_tightId[i] = zvec[i].TnP_l1_tightId;   
    TnP_l1_looseId[i] = zvec[i].TnP_l1_looseId;   
    TnP_l1_dxy[i] = zvec[i].TnP_l1_dxy;   
    TnP_l1_dz[i] = zvec[i].TnP_l1_dz;   
    TnP_l1_edxy[i] = zvec[i].TnP_l1_edxy;   
    TnP_l1_edz[i] = zvec[i].TnP_l1_edz;   
    TnP_l1_ip3d[i] = zvec[i].TnP_l1_ip3d;   
    TnP_l1_sip3d[i] = zvec[i].TnP_l1_sip3d;   
    TnP_l1_ptErr[i] = zvec[i].TnP_l1_ptErr;   
    TnP_l1_lostHits[i] = zvec[i].TnP_l1_lostHits;   
    TnP_l1_trackerLayers[i] = zvec[i].TnP_l1_trackerLayers;   
    TnP_l1_pixelLayers[i] = zvec[i].TnP_l1_pixelLayers;   
    TnP_l1_etaSc[i] = zvec[i].TnP_l1_etaSc;   
    TnP_l1_isGap[i] = zvec[i].TnP_l1_isGap;   
    TnP_l1_r9[i] = zvec[i].TnP_l1_r9;   
    TnP_l1_convVeto[i] = zvec[i].TnP_l1_convVeto;   
    TnP_l1_mvaIdSpring15[i] = zvec[i].TnP_l1_mvaIdSpring15;   
    TnP_l1_relIsoAfterFSR[i] = zvec[i].TnP_l1_relIsoAfterFSR;   
    TnP_l1_chargedHadIso03[i] = zvec[i].TnP_l1_chargedHadIso03;   
    TnP_l1_hasOwnFSR[i] = zvec[i].TnP_l1_hasOwnFSR;  
    if(isMC()) {
      TnP_l1_mcMatchId[i] = zvec[i].TnP_l1_mcMatchId;
      TnP_l1_mcMatchAny[i] = zvec[i].TnP_l1_mcMatchAny;
      TnP_l1_mcPt[i] = zvec[i].TnP_l1_mcPt;
      TnP_l1_mcPt1[i] = zvec[i].TnP_l1_mcPt1;
    } 
    TnP_l1_hlt1L[i] = zvec[i].TnP_l1_hlt1L;   
    TnP_l1_p4WithFSR_pt[i] = zvec[i].TnP_l1_p4WithFSR_pt;   
    TnP_l1_p4WithFSR_eta[i] = zvec[i].TnP_l1_p4WithFSR_eta;   
    TnP_l1_p4WithFSR_phi[i] = zvec[i].TnP_l1_p4WithFSR_phi;   
    TnP_l1_p4WithFSR_mass[i] = zvec[i].TnP_l1_p4WithFSR_mass;   
    TnP_l2_pdgId[i] = zvec[i].TnP_l2_pdgId;   
    TnP_l2_pt[i] = zvec[i].TnP_l2_pt;   
    TnP_l2_eta[i] = zvec[i].TnP_l2_eta;   
    TnP_l2_phi[i] = zvec[i].TnP_l2_phi;   
    TnP_l2_mass[i] = zvec[i].TnP_l2_mass;   
    TnP_l2_charge[i] = zvec[i].TnP_l2_charge;   
    TnP_l2_tightId[i] = zvec[i].TnP_l2_tightId;   
    TnP_l2_looseId[i] = zvec[i].TnP_l2_looseId;   
    TnP_l2_dxy[i] = zvec[i].TnP_l2_dxy;   
    TnP_l2_dz[i] = zvec[i].TnP_l2_dz;   
    TnP_l2_edxy[i] = zvec[i].TnP_l2_edxy;   
    TnP_l2_edz[i] = zvec[i].TnP_l2_edz;   
    TnP_l2_ip3d[i] = zvec[i].TnP_l2_ip3d;   
    TnP_l2_sip3d[i] = zvec[i].TnP_l2_sip3d;   
    TnP_l2_ptErr[i] = zvec[i].TnP_l2_ptErr;   
    TnP_l2_lostHits[i] = zvec[i].TnP_l2_lostHits;   
    TnP_l2_trackerLayers[i] = zvec[i].TnP_l2_trackerLayers;   
    TnP_l2_pixelLayers[i] = zvec[i].TnP_l2_pixelLayers;   
    TnP_l2_etaSc[i] = zvec[i].TnP_l2_etaSc;   
    TnP_l2_isGap[i] = zvec[i].TnP_l2_isGap;   
    TnP_l2_r9[i] = zvec[i].TnP_l2_r9;   
    TnP_l2_convVeto[i] = zvec[i].TnP_l2_convVeto;   
    TnP_l2_mvaIdSpring15[i] = zvec[i].TnP_l2_mvaIdSpring15;   
    TnP_l2_relIsoAfterFSR[i] = zvec[i].TnP_l2_relIsoAfterFSR;   
    TnP_l2_chargedHadIso03[i] = zvec[i].TnP_l2_chargedHadIso03;   
    TnP_l2_hasOwnFSR[i] = zvec[i].TnP_l2_hasOwnFSR;   
    if(isMC()) {
      TnP_l2_mcMatchId[i] = zvec[i].TnP_l2_mcMatchId;
      TnP_l2_mcMatchAny[i] = zvec[i].TnP_l2_mcMatchAny;
      TnP_l2_mcPt[i] = zvec[i].TnP_l2_mcPt;
      TnP_l2_mcPt1[i] = zvec[i].TnP_l2_mcPt1; 
    }
    TnP_l2_hlt1L[i] = zvec[i].TnP_l2_hlt1L;   
    TnP_l2_p4WithFSR_pt[i] = zvec[i].TnP_l2_p4WithFSR_pt;   
    TnP_l2_p4WithFSR_eta[i] = zvec[i].TnP_l2_p4WithFSR_eta;   
    TnP_l2_p4WithFSR_phi[i] = zvec[i].TnP_l2_p4WithFSR_phi;   
    TnP_l2_p4WithFSR_mass[i] = zvec[i].TnP_l2_p4WithFSR_mass;
  }
  //std::cout << "NTnP=" << n << std::endl;
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
