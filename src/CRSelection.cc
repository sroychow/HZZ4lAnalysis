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
#include "PhysicsObjects.h"
#include "CRSelection.h"

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
using std::stoi;
using std::ostream;

using namespace vhtm;

// -----------
// Constructor
// -----------
CRSelection::CRSelection()
  : PhysicsObjSelector(),
    dumpGenInfo_(false),
    useEventList_(false),
    selectEvType_(false),
    evtype_(-1)
{
}
// ----------
// Destructor
// ----------
CRSelection::~CRSelection() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool CRSelection::beginJob() 
{ 
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();

  syncDumpf_.open(syncDumpFile_);

  histf()->cd();
  histf()->mkdir("CRSelection");
  histf()->mkdir("GeneratorLevel");
  bookHistograms();
  
  return true;
}
// ---------------
// Book histograms
// ---------------
void CRSelection::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("CRSelection");

  // Object PLots
  new TH1D("nGoodmuon", "Number of Good muons(with selection cuts) per event", 20, 0, 20);
  new TH1D("nGoodelectron", "Number of Good electrons(with selection cuts) per event", 20, 0, 20);

  // Z and Hplots
  new TH1D("nZcand", "Number of selected Zcandidates per event", 20, -0.5, 19.5);
  new TH1D("nSSllcand", "Number of selected same sign ll pair candidate per event", 20, -0.5, 19.5);
  new TH1D("nOSllcand", "Number of selected opp sign ll pair candidate per event before iso cut", 20, -0.5, 19.5);
  new TH1D("nOSllcand_iso", "Number of selected opp sign ll pair candidate per event after iso cut", 20, -0.5, 19.5);
    
  new TH1D("nZZcand", "Number of selected ZZ candidates per event", 20, 0, 20);
  new TH1F("massnZcand", "Mass of selected Zcandidates", 200, 0., 200.);
  new TH1F("massZ1", "Mass of selected Z1", 200, 0., 200.);
  new TH1F("massZ2", "Mass of selected Z2", 200, 0., 200.);
  new TH1F("mass4l", "4lepton mass", 200, 0., 200.);
  
  new TH1D("evtCutFlow", "Event CutFlow", 8, -0.5, 7.5);
  new TH1D("crSelCutFlowSS", "CR Selection CutFlow (SS)", 9, -0.5, 8.5);
  new TH1D("crSelCutFlowOS", "CR Selection CutFlow (OS)", 10, -0.5, 9.5);
  
  new TH1D("dRlepZal1Zal2", "dRlepZal1Zal2", 100, 0, 5);
  new TH1D("dRlepZbl1Zbl2", "dRlepZbl1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal1Zbl1", "dRlepZal1Zbl1", 100, 0, 5);
  new TH1D("dRlepZal1Zbl2", "dRlepZal1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal2Zbl1", "dRlepZal2Zbl1", 100, 0, 5);
  new TH1D("dRlepZal2Zbl2", "dRlepZal2Zbl2", 100, 0, 5);
  
  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void CRSelection::clearLists() {
  PhysicsObjSelector::clear();

  vtxList_.clear();
  ZCandList_.clear();
  SSllpairCandList_.clear();
  OSllpairCandList_.clear();
  ZZPairVec_.clear();
  ZssllPairVec_.clear();
  ZosllPairVec_.clear();
  evtype_ = -1;
}
// -------------------
// The main event loop
// -------------------
void CRSelection::eventLoop()
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
  for (int ev = 0; ev < nEvents(); ++ev) {
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
    if (ev%nPrint == 0)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << run
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> "
	   << " Events proc. " << setw(8) << ev
	   << endl;
    
    if (useEventList_ && eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;  // not found
    }
    
    histf()->cd();
    histf()->cd("CRSelection");
    AnaUtil::fillHist1D("evtCutFlow", 0, puevWt_);

    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    double ngoodVtx = vtxList_.size();

    // crucial
    histf()->cd();
    histf()->cd("CRSelection");
    
    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt_); 
    
    // is event triggered?
    if (useTrigger() && !isTriggered(true, false)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1, puevWt_);
    
    // at least 1 good PV
    if (ngoodVtx < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, puevWt_);
    
    // main analysis object selection
    findObjects(puevWt_);

    // dump the event and skip the rest
    if (0) {
      showEventNumber();
      dumpEvent();
    }
    
    // access selected objects 
    const auto& tightElePhotonPairVec = getTightIsoElePhotonPairList();
    const auto& tightMuPhotonPairVec  = getTightIsoMuPhotonPairList();
    
    const auto& looseElePhotonPairVec = getLooseElePhotonPairList();
    const auto& looseMuPhotonPairVec  = getLooseMuPhotonPairList();
    
    histf()->cd();
    histf()->cd("CRSelection");
    AnaUtil::fillHist1D("nGoodmuon", looseMuPhotonPairVec.size(), puevWt_);
    AnaUtil::fillHist1D("nGoodelectron", looseElePhotonPairVec.size(), puevWt_);

    // We need at least 4 loose leptons
    if ((looseElePhotonPairVec.size() + looseMuPhotonPairVec.size()) < 4) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, puevWt_);
    
    // Must have at least 2 OSSF leptons
    if (tightElePhotonPairVec.size() < 2 && tightMuPhotonPairVec.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);
    
    // Find the real Z candidate(s)
    if (tightElePhotonPairVec.size() >= 2) ZSelector<vhtm::Electron>(tightElePhotonPairVec, ZCandList_);
    if (tightMuPhotonPairVec.size() >= 2)  ZSelector<vhtm::Muon>(tightMuPhotonPairVec, ZCandList_);
    AnaUtil::fillHist1D("nZcand", ZCandList_.size(), puevWt_);

    // add lepton isolation to the Z candidates found above
    addLeptonIsolation(ZCandList_, tightElePhotonPairVec, tightMuPhotonPairVec);
      
    // At least one real Z candidate found, now find the other SS/OS lepton pair
    if (ZCandList_.size() > 0) {
      if (0) {
	for (auto& v: ZCandList_)
	  HZZ4lUtil::printZCandidate(v, " >> Potential Z Candidate");
      }
      AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
      
      int crtype = static_cast<int>(AnaUtil::cutValue(evselCutMap(), "CRtype"));
      if (crtype == 0) { // SS
	for (auto& v: ZCandList_) {
	  if (looseMuPhotonPairVec.size() >= 2)
	    leptonPairSelector<vhtm::Muon>(looseMuPhotonPairVec, v, false, SSllpairCandList_, ZssllPairVec_);
	  if (looseElePhotonPairVec.size() >= 2)
	    leptonPairSelector<vhtm::Electron>(looseElePhotonPairVec, v, false, SSllpairCandList_, ZssllPairVec_);
	}
	// add lepton isolation to the ll candidates found above
	addLeptonIsolation(SSllpairCandList_, looseElePhotonPairVec, looseMuPhotonPairVec);

	if (ZssllPairVec_.size() > 0) {
	  AnaUtil::fillHist1D("nSSllcand", ZssllPairVec_.size(), puevWt_);
	  finalZllSelector(ZssllPairVec_, false, run, lumis, event);
	}
      } 
      else if (crtype == 1) { // OS
	for (auto& v: ZCandList_) {
	  if (looseMuPhotonPairVec.size() >= 2)
	    leptonPairSelector<vhtm::Muon>(looseMuPhotonPairVec, v, true, OSllpairCandList_, ZosllPairVec_);
	  if (looseElePhotonPairVec.size() >= 2)
	    leptonPairSelector<vhtm::Electron>(looseElePhotonPairVec, v, true, OSllpairCandList_, ZosllPairVec_);
	}
	// add lepton isolation to the ll candidates found above
	addLeptonIsolation(OSllpairCandList_, looseElePhotonPairVec, looseMuPhotonPairVec);

	if (ZosllPairVec_.size() > 0) {
	  AnaUtil::fillHist1D("nOSllcand", ZosllPairVec_.size(), puevWt_);
	  finalZllSelector(ZosllPairVec_, true, run, lumis, event);
	}
      }
      else 
	cout << "Wrong parameter! CRtype = " << crtype << endl;
    }
  }
  // Analysis is over
  endJob();
}
void CRSelection::finalZllSelector(std::vector<std::pair<ZCandidate, ZCandidate> >& objPairList, bool studyOSPair, int run, int lumi, int event) {
  AnaUtil::fillHist1D("evtCutFlow", 6, puevWt_);
  unsigned int evtHasCR = 0;
  ZCandidate evZ, evll;
  cout << setprecision(3);
  for (unsigned int i = 0; i < objPairList.size(); ++i) {
    auto& p = objPairList[i];
    bool decision = ((studyOSPair) ? CRSelectorZOSll(p.first, p.second) : CRSelectorZSSll(p.first, p.second));
    if (decision) {
      evtHasCR++;
      if (evtHasCR == 1 || p.first.massDiff < evZ.massDiff) {
	evZ = p.first;
	evll = p.second;
      }
      else if (p.first.massDiff == evZ.massDiff) {
	if ( (p.second.l1P4.Pt() + p.second.l2P4.Pt())/2. > (evll.l1P4.Pt() + evll.l2P4.Pt())/2. )
	  evll = p.second;
      }
    }
  }
  if (evtHasCR) {
    AnaUtil::fillHist1D("evtCutFlow", 7, puevWt_);
    double mass4l = (evZ.l1P4 + evZ.l2P4 + evZ.fsrPhoP4 + evll.l1P4 + evll.l2P4 + evll.fsrPhoP4).M();

    // dump event
    cout << "---- New Event" << endl;
    showEventNumber(false);
    cout << " --- #Z Candidates: " << ZCandList_.size() 
	 << ", #Zll Candidates: " << objPairList.size() 
	 << ", #Zll Candidates(final): " << evtHasCR 
	 << ", flavour: "  << (HZZ4lUtil::sameFlavourZPair(evZ, evll) ? "same" : "different") 
	 << endl;
    cout << "--- mass4l: " << mass4l << " GeV" 
	 << endl;
    cout << "--- fGridRhoFastjetAll: " << getEventGridRho() << endl;
    HZZ4lUtil::printZCandidate(evZ, " >> Z Candidate");
    HZZ4lUtil::printZCandidate(evll, " >> ll Candidate");
    dumpEvent();
    
    const auto& jetVec = getLeptonCleanedLooseJetList();
    double jet1Pt = -1.0, jet2Pt = -1.0;
    if (jetVec.size()) jet1Pt = jetVec.at(0).pt;
    if (jetVec.size() > 1) jet2Pt = jetVec.at(1).pt;
    HZZ4lUtil::syncDumper(run, lumi, event, evZ, evll, jetVec.size(), jet1Pt, jet2Pt, syncDumpf_);
  }
}
bool CRSelection::CRSelectorZSSll(ZCandidate& Z, ZCandidate& ssll, bool verbose) {
  histf()->cd();
  histf()->cd("CRSelection");
  AnaUtil::fillHist1D("crSelCutFlowSS", 0, puevWt_);

  // -- from twiki --
  // dR(eta,phi)>0.02 between each of the four leptons (to remove ghosts)
  // and candidates where the sam eleptons contribute more than once 
  double dra1a2 = Z.l1P4.DeltaR(Z.l2P4);
  double drb1b2 = ssll.l1P4.DeltaR(ssll.l2P4);
  double dra1b1 = Z.l1P4.DeltaR(ssll.l1P4);
  double dra1b2 = Z.l1P4.DeltaR(ssll.l2P4);
  double dra2b1 = Z.l2P4.DeltaR(ssll.l1P4);
  double dra2b2 = Z.l2P4.DeltaR(ssll.l2P4);
  AnaUtil::fillHist1D("dRlepZal1Zal2", dra1a2, puevWt_);
  AnaUtil::fillHist1D("dRlepZbl1Zbl2", drb1b2, puevWt_);
  AnaUtil::fillHist1D("dRlepZal1Zbl1", dra1b1, puevWt_);
  AnaUtil::fillHist1D("dRlepZal1Zbl2", dra1b2, puevWt_);
  AnaUtil::fillHist1D("dRlepZal2Zbl1", dra2b1, puevWt_);
  AnaUtil::fillHist1D("dRlepZal2Zbl2", dra2b2, puevWt_);
  bool dRlep = dra1a2 > 0.02 &&
    drb1b2 > 0.02 &&
    dra1b1 > 0.02 &&
    dra1b2 > 0.02 &&
    dra2b1 > 0.02 &&
    dra2b2 > 0.02;
  if (!dRlep) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 1, puevWt_);

  // Lepton pT
  vector<TLorentzVector> lepP4List;
  lepP4List.push_back(Z.l1P4);
  lepP4List.push_back(Z.l2P4);
  lepP4List.push_back(ssll.l1P4);
  lepP4List.push_back(ssll.l2P4);
  std::sort(lepP4List.begin(), lepP4List.end(), PtComparatorTL<TLorentzVector>());
  
  // -- from twiki --
  // any two leptons of the four have pt > 20/10
  if (lepP4List[0].Pt() <= 20) return false;;
  AnaUtil::fillHist1D("crSelCutFlowSS", 2, puevWt_);
  
  if (lepP4List[1].Pt() <= 10) return false;;
  AnaUtil::fillHist1D("crSelCutFlowSS", 3, puevWt_);

  if (Z.mass <= 40. || Z.mass >= 120.) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 4, puevWt_);
  
  if (ssll.mass <= 12. || ssll.mass >= 120.) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 5, puevWt_);

  // QCD suppression: m(ll) > 4 GeV on all OS lepton pairs (only three in this case given that the Z2 is SS)
  TLorentzVector ZaP4, ZbP4;
  if (Z.l1Charge + ssll.l1Charge == 0) {
    ZaP4 = Z.l1P4 + ssll.l1P4;
    ZbP4 = Z.l2P4 + ssll.l2P4;
  }
  else {
    ZaP4 = Z.l1P4 + ssll.l2P4;
    ZbP4 = Z.l2P4 + ssll.l1P4;
  }
  if (ZaP4.M() <= 4. || ZbP4.M() <= 4.) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 6, puevWt_);
  
  // Smart Cut
  bool smartcutFlag = false;
  if (HZZ4lUtil::sameFlavourZPair(Z, ssll)) {
    TLorentzVector ZSSlepP4, ZOSlepP4;
    TLorentzVector ZSSlepFsrP4, ZOSlepFsrP4;
    if (Z.l1Charge == ssll.l1Charge) {
      ZSSlepP4 = Z.l1P4;
      ZSSlepFsrP4 = Z.l1FsrP4;

      ZOSlepP4 = Z.l2P4;
      ZOSlepFsrP4 = Z.l2FsrP4;
    } 
    else {
      ZSSlepP4 = Z.l2P4;
      ZSSlepFsrP4 = Z.l2FsrP4;

      ZOSlepP4 = Z.l1P4;
      ZOSlepFsrP4 = Z.l1FsrP4;
    }
    
    TLorentzVector alt1Z1P4 = ZOSlepP4 + ssll.l1P4 + ZOSlepFsrP4 + ssll.l1FsrP4;
    TLorentzVector alt1Z2P4 = ZSSlepP4 + ssll.l2P4 + ZSSlepFsrP4 + ssll.l2FsrP4;
    
    TLorentzVector alt2Z1P4 = ZOSlepP4 + ssll.l2P4 + ZOSlepFsrP4 + ssll.l2FsrP4;
    TLorentzVector alt2Z2P4 = ZSSlepP4 + ssll.l1P4 + ZSSlepFsrP4 + ssll.l1FsrP4;
    
    if (0)
      cout << "Smart Cut combination 1" << endl
	   << "--- SmartCut Leptons: Z1: (" << Z.l1Index << "," << Z.l2Index
	   << "), Z2: (" << ssll.l1Index << "," << ssll.l2Index << ")" << endl
	   << "--- SmartCut: (" << Z.flavour << ", " << ssll.flavour
	   << ", " << Z.mass << ", " << alt1Z1P4.M()
	   << ", " << ssll.mass << ", " << alt1Z2P4.M() << ")"
	   << endl
	   << "Smart Cut combination 2" << endl
	   << "--- SmartCut Leptons: Z1: (" << Z.l1Index << "," << Z.l2Index
	   << "), Z2: (" << ssll.l1Index << "," << ssll.l2Index << ")" << endl
	   << "--- SmartCut: (" << Z.flavour << ", " << ssll.flavour
	   << ", " << Z.mass << ", " << alt2Z1P4.M()
	   << ", " << ssll.mass << ", " << alt2Z2P4.M() << ")"
	   << endl;
    if ( (std::fabs(alt1Z1P4.M() - HZZ4lUtil::MZnominal) < Z.massDiff && alt1Z2P4.M() < 12.) ||
	 (std::fabs(alt2Z1P4.M() - HZZ4lUtil::MZnominal) < Z.massDiff && alt2Z2P4.M() < 12.) ) {
      smartcutFlag = true;
    }
  }
  if (smartcutFlag) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 7, puevWt_);
  
  double mass4l = (Z.l1P4 + Z.l2P4 + Z.fsrPhoP4 + ssll.l1P4 + ssll.l2P4 + ssll.fsrPhoP4).M();
  if (mass4l <= 70.) return false;
  AnaUtil::fillHist1D("crSelCutFlowSS", 8, puevWt_);

  return true;
}
bool CRSelection::CRSelectorZOSll(ZCandidate& Z, ZCandidate& osll, bool verbose) {
  AnaUtil::fillHist1D("crSelCutFlowOS", 0, puevWt_);
  
  // -- from twiki --
  // dR(eta,phi)>0.02 between each of the four leptons (to remove ghosts)
  // and candidates where the sam eleptons contribute more than once 
  double dra1a2 = Z.l1P4.DeltaR(Z.l2P4);
  double drb1b2 = osll.l1P4.DeltaR(osll.l2P4);
  double dra1b1 = Z.l1P4.DeltaR(osll.l1P4);
  double dra1b2 = Z.l1P4.DeltaR(osll.l2P4);
  double dra2b1 = Z.l2P4.DeltaR(osll.l1P4);
  double dra2b2 = Z.l2P4.DeltaR(osll.l2P4);
  AnaUtil::fillHist1D("dRlepZal1Zal2", dra1a2, puevWt_);
  AnaUtil::fillHist1D("dRlepZbl1Zbl2", drb1b2, puevWt_);
  AnaUtil::fillHist1D("dRlepZal1Zbl1", dra1b1, puevWt_);
  AnaUtil::fillHist1D("dRlepZal1Zbl2", dra1b2, puevWt_);
  AnaUtil::fillHist1D("dRlepZal2Zbl1", dra2b1, puevWt_);
  AnaUtil::fillHist1D("dRlepZal2Zbl2", dra2b2, puevWt_);
  bool dRlep = dra1a2 > 0.02 &&
    drb1b2 > 0.02 &&
    dra1b1 > 0.02 &&
    dra1b2 > 0.02 &&
    dra2b1 > 0.02 &&
    dra2b2 > 0.02;
  if (!dRlep) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 1, puevWt_);

  // Lepton pT
  vector<TLorentzVector> lepP4List;
  lepP4List.push_back(Z.l1P4);
  lepP4List.push_back(Z.l2P4);
  lepP4List.push_back(osll.l1P4);
  lepP4List.push_back(osll.l2P4);
  std::sort(lepP4List.begin(), lepP4List.end(), PtComparatorTL<TLorentzVector>());
  
  // -- from twiki --
  // any two leptons of the four have pt > 20/10
  if (lepP4List[0].Pt() <= 20) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 2, puevWt_);

  if (lepP4List[1].Pt() <= 10) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 3, puevWt_);

  if (Z.mass <= 40. || Z.mass >= 120.) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 4, puevWt_);
  
  if (osll.mass <= 12. || osll.mass >= 120.) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 5, puevWt_);

  // Tight + isolation veto
  bool oslepIdIso1 = false,
    oslepIdIso2 = false;
  if (osll.flavour == HZZ4lUtil::ZType::mumu) {
    const auto& muPhotonPairVec = getLooseMuPhotonPairList();
    const vhtm::Muon& mu1 = muPhotonPairVec.at(osll.l1Index).first;
    const vhtm::Muon& mu2 = muPhotonPairVec.at(osll.l2Index).first;

    // made by requiring exactly one ("3P1F") or two ("2P2F") of the additional
    // leptons to pass loose cuts (loose ID + SIP) but fail tight ID and isolation.
    if (!mu1.isPFMuon || osll.l1Isolation >= 0.4) oslepIdIso1 = true;
    if (!mu2.isPFMuon || osll.l2Isolation >= 0.4) oslepIdIso2 = true;
  }
  else if (osll.flavour == HZZ4lUtil::ZType::ee) {
    const auto& elePhotonPairVec = getLooseElePhotonPairList();
    const vhtm::Electron& ele1 = elePhotonPairVec.at(osll.l1Index).first;
    const vhtm::Electron& ele2 = elePhotonPairVec.at(osll.l2Index).first;

    if (!HZZ4lUtil::electronBDT(ele1) || osll.l1Isolation >= 0.5) oslepIdIso1 = true;
    if (!HZZ4lUtil::electronBDT(ele2) || osll.l2Isolation >= 0.5) oslepIdIso2 = true;
  }
  
  // 0=3P1F && 1=2P2F
  bool oslepIdIso = false;
  int os_crtype = static_cast<int>(AnaUtil::cutValue(evselCutMap(), "osCRtype"));
  if (os_crtype == 0) { // 3P1F
    if ((!oslepIdIso1 && oslepIdIso2) || (oslepIdIso1 && !oslepIdIso2)) oslepIdIso = true;
  }
  else if (os_crtype == 1) { // 2P2F
    if (oslepIdIso1 && oslepIdIso2) oslepIdIso = true;
  }
  else 
    cout << "Wrong parameter! osCRType: " << os_crtype << endl;
  
  if (!oslepIdIso) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 6, puevWt_);

  // QCD suppression: m(ll)>4 GeV on all OS lepton pairs (only three in this case given that the Z2 is SS)
  TLorentzVector ZaP4, ZbP4;
  if (Z.l1Charge + osll.l1Charge == 0) {
    ZaP4 = Z.l1P4 + osll.l1P4;
    ZbP4 = Z.l2P4 + osll.l2P4;
  }
  else {
    ZaP4 = Z.l1P4 + osll.l2P4;
    ZbP4 = Z.l2P4 + osll.l1P4;
  }
  if (Z.mass <= 4 || osll.mass <= 4 || ZaP4.M() <= 4. || ZbP4.M() <= 4.) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 7, puevWt_);

  // Smart Cut
  bool smartcutFlag = false;
  if (HZZ4lUtil::sameFlavourZPair(Z, osll)) {
    TLorentzVector altZ1P4, altZ2P4;
    if (Z.l1Charge + osll.l1Charge == 0) {
      altZ1P4 = Z.l1P4 + osll.l1P4 + Z.l1FsrP4 + osll.l1FsrP4;
      altZ2P4 = Z.l2P4 + osll.l2P4 + Z.l2FsrP4 + osll.l2FsrP4;
    }
    else {
      altZ1P4 = Z.l1P4 + osll.l2P4 + Z.l1FsrP4 + osll.l2FsrP4;
      altZ2P4 = Z.l2P4 + osll.l1P4 + Z.l2FsrP4 + osll.l1FsrP4;
    }
    if (std::fabs(altZ2P4.M() - HZZ4lUtil::MZnominal) < std::fabs(altZ1P4.M() - HZZ4lUtil::MZnominal)) {
      TLorentzVector lv = altZ1P4;
      altZ1P4 = altZ2P4;
      altZ2P4 = lv;
    }
    if (0) 
      cout << "--- SmartCut Leptons: Z1: (" << Z.l1Index << "," << Z.l2Index
	   << "), Z2: (" << osll.l1Index << "," << osll.l2Index << ")" << endl
	   << "--- SmartCut: (" << Z.flavour << ", " << osll.flavour
	   << ", " << Z.mass << ", " << altZ1P4.M() 
	   << ", " << osll.mass << ", " << altZ2P4.M() << ")"
	   << endl;
    // Finally apply the Smart cut
    if (std::fabs(altZ1P4.M() - HZZ4lUtil::MZnominal) < Z.massDiff && altZ2P4.M() < 12) {
      if (0) cout << " -- SmartCut: skip Z+ll Candidate" << endl;
      smartcutFlag = true;
    }
  }
  if (smartcutFlag) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 8, puevWt_);  
    
  double mass4l = (Z.l1P4 + Z.l2P4 + Z.fsrPhoP4 + osll.l1P4 + osll.l2P4 + osll.fsrPhoP4).M();
  if (mass4l <= 70) return false;
  AnaUtil::fillHist1D("crSelCutFlowOS", 9, puevWt_);
  
  return true;
}
void CRSelection::endJob() {
  syncDumpf_.close();
  closeFiles();
  
  histf()->cd();
  histf()->cd("CRSelection");
  string evlabels[] = {
    "Events processed",
    "Events passing trigger",
    "Events with > 0 good vertex",
    "4 leptons",
    "# tight leptons >= 2",
    "# of Z Candidates > 0",
    "Candidate Phase 1",
    "Candidate Phase 2"
  };
  HZZ4lUtil::showEfficiency("evtCutFlow", evlabels, "Event Selection");  

  string crlabels [] = {
    "All",
    "dRlep",
    "lep pT > 20 GeV",
    "lep pT > 10 GeV",
    "Zmasscut",
    "llmasscut",
    "tight+iso",
    "qcdsupressflag",
    "smartcutFlag",
    "mass4lCut"
  };
  HZZ4lUtil::showEfficiency("crSelCutFlow", crlabels, "CR Event Selection");  

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
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
bool CRSelection::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  std::ifstream fin(jobFile.c_str(), std::ios::in);    
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
    if (key == "useEventList")
      useEventList_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "syncDumpFile")
      syncDumpFile_ = value;

    tokens.clear();
  }
  // Close the file
  fin.close();

  selectEvType_ = (static_cast<int>(AnaUtil::cutValue(evselCutMap(), "selectEvType")) > 0) ? true : false;
  printJob();

  return true;
}
void CRSelection::printJob(ostream& os) const
{
  AnaBase::printJob(os);
}
