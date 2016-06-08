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
#include "DataMcComparator.h"

#include "MEMCalculators.h"
using namespace MEMNames;

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
DataMcComparator::DataMcComparator()
  : PhysicsObjSelector(),
    checkGen_(true),
    dumpGenInfo_(false),
    useEventList_(false),
    selectEvType_(false),
    evtype_(-1),
    doKDcalc_(false),
    dumpFilename_("syncDumpFile.txt")
{
}
// ----------
// Destructor
// ----------
DataMcComparator::~DataMcComparator() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool DataMcComparator::beginJob() 
{ 
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();

  histf()->cd();
  histf()->mkdir("DataMcComparator");
  bookHistograms();
  
  return true;
}
// ---------------
// Book histograms
// ---------------
void DataMcComparator::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("DataMcComparator");

  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodmuon", "Number of Good muons(with selection cuts) per event", 20, 0, 20);
  new TH1D("nGoodelectron", "Number of Good electrons(with selection cuts) per event", 20, 0, 20);

  // Z and Hplots
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);
  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);
  new TH1D("nZcand", "Number of selected Zcandidates per event", 20, 0, 20);
  new TH1F("massZcand", "Mass of selected Zcandidates with tag + loose probe", 200, 0., 200.);
  new TH1F("massZcand_SIPpass", "Mass of selected Zcandidates with tag + probe(SIP)", 200, 0., 200.);
  new TH1F("massZcand_SIPISOpass", "Mass of selected Zcandidates with tag + probe(SIP+ISO)", 200, 0., 200.);
  new TH1F("relIso_tag_ele", "Relative Isolation of tag muon", 100, 0., 1.);
  new TH1F("relIso_probe_ele", "Relative Isolation of probe muon", 100, 0., 1.);
  new TH1F("relIso_Zpx_ele", "Relative Isolation of loose muons not in Z", 100, 0., 1.);
  new TH1F("relIso_tag_mu", "Relative Isolation of tag muon", 100, 0., 1.);
  new TH1F("relIso_probe_mu", "Relative Isolation of probe muon", 100, 0., 1.);
  new TH1F("relIso_Zpx_mu", "Relative Isolation of loose muons not in Z", 100, 0., 1.);
  new TH1F("sip_tag_ele", "sip of tag muon", 100, 0., 50.);
  new TH1F("sip_probe_ele", "sip of probe muon", 100, 0., 50.);
  new TH1F("sip_Zpx_ele", "SIP of loose muons not in Z", 100, 0., 50.);
  new TH1F("sip_tag_mu", "SIP of tag muon", 100, 0., 50.);
  new TH1F("sip_probe_mu", "SIP of probe muon", 100, 0., 50.);
  new TH1F("sip_Zpx_mu", "SIP of loose muons not in Z", 100, 0., 50.);
  
  
  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void DataMcComparator::clearLists() {
  PhysicsObjSelector::clear();

  vtxList_.clear();
  genZList_.clear();
  ZCandList_.clear();
  ZZPairVec_.clear();
  evtype_ = -1;
}
// -------------------
// The main event loop
// -------------------
void DataMcComparator::eventLoop()
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
    
    if (useEventList_ && eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
      //if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue;
    }
    
    histf()->cd();
    histf()->cd("DataMcComparator");
    AnaUtil::fillHist1D("evtCutFlow", 0, puevWt_);
    
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);

//    cout << "point 1" << endl;
    findVtxInfo(vtxList_, op, fLog());
    double ngoodVtx = vtxList_.size();      // crucial
  
//    cout << "point 2" << endl;

    histf()->cd();
    histf()->cd("DataMcComparator");
    
    AnaUtil::fillHist1D("evtCutFlow", 1, puevWt_);

    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt_); 
    
    // is event triggered?
    if (useTrigger() && !isTriggered(true, false)) continue;

    AnaUtil::fillHist1D("evtCutFlow", 2, puevWt_);
//    cout << "point 3" << endl;
    
    // at least 1 good PV
    if (ngoodVtx < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, puevWt_);
//   cout << "point 4" << endl;

    // main analysis object selection
    findObjects(puevWt_);

//    cout << "point 5" << endl;

    // access selected objects 
    const auto& tightEleVector = getTightEleList();
    const auto& tightMuVector  = getTightMuList();

//   cout << "point 6" << endl;

    if( tightEleVector.empty() && tightMuVector.empty() )     continue;

    histf()->cd();
    histf()->cd("DataMcComparator");

    AnaUtil::fillHist1D("nGoodmuon", tightMuVector.size(), puevWt_);
    AnaUtil::fillHist1D("nGoodelectron", tightEleVector.size(), puevWt_);
    AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);
    
    //mmmm = 0, eeee = 1, eemm = 2
    
    if( static_cast<int>(AnaUtil::cutValue(evselCutMap(), "evType")) == EventType::eeee ) {
      if ( tightEleVector.size() ) {
        histf()->cd();
        histf()->cd("DataMcComparator");
        AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
        int tagidx = getTagElectronidx(tightEleVector);
        if( tagidx < 0 )    continue;
        AnaUtil::fillHist1D("evtCutFlow", 6, puevWt_);
        const auto& tagEle = tightEleVector[tagidx];
        const auto& loosepreSIPEleVector = getPreSIPLooseEleList();
	ZSelector<vhtm::Electron>(tagEle, loosepreSIPEleVector, ZCandList_ );
        if( ZCandList_.empty() )   continue;
        AnaUtil::fillHist1D("evtCutFlow", 7, puevWt_);
        int probeidx = ZCandList_[0].l2Index;
        const auto& probeEle = loosepreSIPEleVector[probeidx];
        MassComparatorElectron(tagEle,probeEle);
      }
    } 
    else if( static_cast<int>(AnaUtil::cutValue(evselCutMap(), "evType")) == EventType::mmmm ) {
      if ( tightMuVector.size() ) {
        histf()->cd();
        histf()->cd("DataMcComparator");
        AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
        int tagidx = getTagMuonidx(tightMuVector);
        if( tagidx < 0 )    continue;
        AnaUtil::fillHist1D("evtCutFlow", 6, puevWt_);
        const auto& tagMu = tightMuVector[tagidx];
        const auto& loosepreSIPMuVector = getPreSIPLooseMuList();
	ZSelector<vhtm::Muon>(tagMu,loosepreSIPMuVector, ZCandList_);
        if( ZCandList_.empty() )   continue;
        AnaUtil::fillHist1D("evtCutFlow", 7, puevWt_);
        int probeidx = ZCandList_[0].l2Index;
        const auto& probeMu = loosepreSIPMuVector[probeidx];
        MassComparatorMuon(tagMu,probeMu);
      }
    }
    
    // At least one Z candidate found 
  }	
  // Analysis is over
  endJob();
}
// create unique lepton pair combination giving a Z statisfying Z2 mass cuts
// and push them into a vector
template <typename T>
void DataMcComparator::ZSelector( const T& tagLep, const std::vector<T>& leptonVec, 
                                  std::vector<ZCandidate>& candList) {
  for(unsigned int i = 0; i < leptonVec.size(); ++i) {
    const auto& probeLep = leptonVec[i];
    if( tagLep.charge + probeLep.charge != 0 ) continue; 
    const TLorentzVector& tagLepP4 = HZZ4lUtil::getP4(tagLep);
    const TLorentzVector& probeLepP4 = HZZ4lUtil::getP4(probeLep);
    if( tagLepP4.DeltaR(probeLepP4) < 0.02 )  continue;
    double Zmass  = (tagLepP4 + probeLepP4).M();
    if( Zmass < 60. || Zmass > 120. )    continue;
    ZCandidate ztmp;
    ztmp.l1P4 = tagLepP4;
    ztmp.l2P4 = probeLepP4;
    ztmp.l1Charge = tagLep.charge;
    ztmp.l2Charge = probeLep.charge;
    ztmp.l2Index = i;
    ztmp.mass = Zmass;
    ztmp.massDiff = std::fabs(Zmass - HZZ4lUtil::MZnominal);
    candList.push_back(ztmp);
    break;
  }
}
//T&Pfuntions
void DataMcComparator::MassComparatorElectron(const vhtm::Electron& tagLep, const vhtm::Electron& probeLep) {
  histf()->cd();
  histf()->cd("DataMcComparator");
  double Zmass = ZCandList_[0].mass;

  AnaUtil::fillHist1D("massZcand", Zmass, puevWt_);
  
  double tagLepsip = std::fabs(tagLep.dB3D/tagLep.edB3D); 
  double probeLepsip = std::fabs(probeLep.dB3D/probeLep.edB3D); 
  AnaUtil::fillHist1D("sip_tag_ele", tagLepsip, puevWt_);
  AnaUtil::fillHist1D("sip_probe_ele", probeLepsip, puevWt_);

  for( const auto& lmu : getPreSIPLooseMuList() ) {
    TLorentzVector lp4 = HZZ4lUtil::getP4(lmu);
    if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
      continue;
    AnaUtil::fillHist1D("sip_Zpx_mu", std::fabs(lmu.dB3D/lmu.edB3D), puevWt_);
  }
  for( const auto& lele : getPreSIPLooseEleList() ) {
    TLorentzVector lp4 = HZZ4lUtil::getP4(lele);
    if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
      continue;
    AnaUtil::fillHist1D("sip_Zpx_ele", std::fabs(lele.dB3D/lele.edB3D), puevWt_);
  }

  if( probeLepsip < AnaUtil::cutValue(electronCutMap(), "SIP3D"))  {
    AnaUtil::fillHist1D("massZcand_SIPpass", Zmass, puevWt_);
 
    double taglepIso = HZZ4lUtil::pfiso( tagLep, getEventGridRho(), 0.);
    double probelepIso = HZZ4lUtil::pfiso( probeLep, getEventGridRho(), 0.);
    AnaUtil::fillHist1D("relIso_tag_ele", taglepIso/tagLep.pt, puevWt_);
    AnaUtil::fillHist1D("relIso_probe_ele", probelepIso/probeLep.pt, puevWt_);
    for( const auto& lmu : getLooseMuList() ) {
      TLorentzVector lp4 = HZZ4lUtil::getP4(lmu);
      if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
        continue;
      double liso = HZZ4lUtil::pfiso( lmu, 0.);
      AnaUtil::fillHist1D("relIso_Zpx_mu", liso/lmu.pt, puevWt_);
    }
    for( const auto& lele : getLooseEleList() ) {
      TLorentzVector lp4 = HZZ4lUtil::getP4(lele);
      if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
        continue;
      double liso = HZZ4lUtil::pfiso( lele, 0.);
      AnaUtil::fillHist1D("relIso_Zpx_ele", liso/lele.pt, puevWt_);
    }
    if( probelepIso/probeLep.pt < 0.5 )  
      AnaUtil::fillHist1D("massZcand_SIPISOpass", Zmass, puevWt_);
  }
}
void DataMcComparator::MassComparatorMuon(const vhtm::Muon& tagLep, const vhtm::Muon& probeLep ) {
  histf()->cd();
  histf()->cd("DataMcComparator");
  double Zmass = ZCandList_[0].mass;
  double taglepIso = HZZ4lUtil::pfiso( probeLep, 0.);
  double probelepIso = HZZ4lUtil::pfiso( probeLep, 0.);

  AnaUtil::fillHist1D("relIso_tag_mu", taglepIso/tagLep.pt, puevWt_);
  AnaUtil::fillHist1D("relIso_probe_mu", probelepIso/probeLep.pt, puevWt_);
  AnaUtil::fillHist1D("massZcand", Zmass, puevWt_);

  double tagLepsip = std::fabs(tagLep.dB3D/tagLep.edB3D); 
  double probeLepsip = std::fabs(probeLep.dB3D/probeLep.edB3D); 
  AnaUtil::fillHist1D("sip_tag_mu", tagLepsip, puevWt_);
  AnaUtil::fillHist1D("sip_probe_mu", probeLepsip, puevWt_);

  for( const auto& lmu : getPreSIPLooseMuList() ) {
    TLorentzVector lp4 = HZZ4lUtil::getP4(lmu);
    if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
      continue;
    AnaUtil::fillHist1D("sip_Zpx_mu", std::fabs(lmu.dB3D/lmu.edB3D), puevWt_);
  }
  for( const auto& lele : getPreSIPLooseEleList() ) {
    TLorentzVector lp4 = HZZ4lUtil::getP4(lele);
    if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
      continue;
    AnaUtil::fillHist1D("sip_Zpx_ele", std::fabs(lele.dB3D/lele.edB3D), puevWt_);
  }

  if( std::fabs(probeLep.dB3D/probeLep.edB3D) < AnaUtil::cutValue(electronCutMap(), "SIP3D")) {
    AnaUtil::fillHist1D("massZcand_SIPpass", Zmass, puevWt_);
    for( const auto& lmu : getLooseMuList() ) {
      TLorentzVector lp4 = HZZ4lUtil::getP4(lmu);
      if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
        continue;
      double liso = HZZ4lUtil::pfiso( lmu, 0.);
      AnaUtil::fillHist1D("relIso_Zpx_mu", liso/lmu.pt, puevWt_);
    }
    for( const auto& lele : getLooseEleList() ) {
      TLorentzVector lp4 = HZZ4lUtil::getP4(lele);
      if( lp4 == HZZ4lUtil::getP4(tagLep) || lp4 == HZZ4lUtil::getP4(probeLep) )
        continue;
      double liso = HZZ4lUtil::pfiso( lele, 0.);
      AnaUtil::fillHist1D("relIso_Zpx_ele", liso/lele.pt, puevWt_);
    }

    if( probelepIso/probeLep.pt < 0.4 )       
      AnaUtil::fillHist1D("massZcand_SIPISOpass", Zmass, puevWt_);
   }
}

int DataMcComparator::getTagElectronidx(const std::vector<vhtm::Electron>& elVec) {
 int idx = -1;
 for( unsigned int i = 0; i<elVec.size(); i++ ) {
   if( HZZ4lUtil::pfiso( elVec[i] , getEventGridRho(), 0.) < 0.5 ) {
     idx =  i;
     break;
   }
 }
 return idx;
}
int DataMcComparator::getTagMuonidx(const std::vector<vhtm::Muon>& muVec) {
 int idx = -1;
 for( unsigned int i = 0; i<muVec.size(); i++ ) {
   if( HZZ4lUtil::pfiso( muVec[i] ,  0.) < 0.4 ) {
     idx =  i;
     break;
   }
 }
 return idx;
}

void DataMcComparator::endJob() {
  syncDumpf_.close();
  closeFiles();

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
bool DataMcComparator::readJob(const string& jobFile, int& nFiles)
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
    if (key == "useEventList")
      useEventList_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "checkGenInfo")
      checkGen_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "doKD")
      doKDcalc_ = (stoi(value.c_str()) > 0) ? true : false;
    else if (key == "syncDumpFile")
      dumpFilename_ = value.c_str();

    tokens.clear();
  }
  // Close the file
  fin.close();

  syncDumpf_.open(dumpFilename_, ios::out);
  if (!syncDumpf_) {
    cerr << "Output File: " << dumpFilename_ << " could not be opened!" << endl;
    return false;
  }

  selectEvType_ = (static_cast<int>(AnaUtil::cutValue(evselCutMap(), "selectEvType")) > 0) ? true : false;
  printJob();

  return true;
}
void DataMcComparator::printJob(ostream& os) const
{
  AnaBase::printJob(os);
}
