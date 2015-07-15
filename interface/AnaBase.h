#ifndef __AnaBase__hh
#define __AnaBase__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TVector3.h"

#include "PhysicsObjects.h"

typedef unsigned int uint;
typedef unsigned long ulong;

// REQUIRED, most probably to solve circular dependence problem!!!
class AnaBase;
class TChain;
class TFile;

template <class T>
class PtComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.pt > b.pt;
  }
};

template <class T>
class PtComparatorTL {
public:
  bool operator()(const T &a, const T &b) const {
    return a.Pt() > b.Pt();
  }
};

template <class T>
class MassComparator {
public:
  bool operator()(const T &a, const T &b) const {
    TLorentzVector l1,l2;
    l1.SetPtEtaPhiE(a.pt, a.eta, a.phi, a.energy);
    l2.SetPtEtaPhiE(b.pt, b.eta, b.phi, b.energy);
    return l1.M() > l2.M();
  }
};

typedef struct  
{
  bool verbose;
  bool usesbit;
  bool printselected;
} Options;

class AnaBase {
    
public:

  AnaBase();
  virtual ~AnaBase();
    
  virtual void eventLoop() = 0;  // the main analysis
  virtual bool beginJob();
  virtual void endJob() = 0;

  virtual void selectEvent() = 0;

  int setInputFile(const std::string& fname);
  bool branchFound(const std::string& b);
  int getEntries() const;
  void setData(int val);  
  int getRunNumber() const;
  virtual bool openFiles();
  virtual void closeFiles(); 
  virtual bool readJob(const std::string& jobFile, int& nFiles);
  virtual void printJob(std::ostream& os=std::cout) const;
  void showEventNumber(bool addNewline=true, std::ostream& os=std::cout) const {
    const vhtm::Event& evt = eventList_->at(0);
    std::cout << ">>> Run " << evt.run
	      << " Lumis " << evt.lumis 
	      << " Event " << evt.event;
    if (addNewline) std::cout << std::endl;
  }
  bool isTriggered(bool check_prescale=true, bool verbose=false) const;
  virtual void dumpTriggerPaths(std::ostream& os=std::cout, bool check_prescale=true) const;
  void dumpTriggerObjectInfo(const std::vector<vhtm::TriggerObject>& list, std::ostream& os=std::cout) const;
  double wtPileUp(int& nPU) const;
  bool readPileUpHist();
  bool matchTriggerPath(const std::vector<std::string>& v, const std::string& path) const;
  double matchTriggerObject(const std::vector<vhtm::TriggerObject>& trigObjList, 
                            const TLorentzVector& obj, 
                            const std::string& trigpath, 
                            int trig_skip, 
                            double maxPtDiff, 
                            int& trig_indx) const;
  
  void clearEvent();
  void enableBranches();
   int getEntry(int lflag) const;
  void setAddresses(); 

  void dumpGenInfo(std::ostream& os=std::cout) const;
  int getMotherId(const vhtm::GenParticle& gp, int& mmid) const;

  void findVtxInfo(std::vector<vhtm::Vertex>& list, Options& op, std::ostream& os=std::cout);
  void findElectronInfo(std::vector<vhtm::Electron>& list, double vz, Options& op, std::ostream& os=std::cout);
  void findMuonInfo(std::vector<vhtm::Muon>& list, double vz, Options& op, std::ostream& os=std::cout);
  void findTauInfo(std::vector<vhtm::Tau>& list, double vz, Options& op, std::ostream& os=std::cout);
  void findJetInfo(std::vector<vhtm::Jet>& list, Options& op, std::ostream& os=std::cout);
  void findTriggerObjectInfo(std::vector<vhtm::TriggerObject>& list);
  void dumpEvent(const char* optstr, std::ostream& os=std::cout, bool ps=false);
  TVector3 findLeptonVtx(int index, bool& isGoodVtx);
  void findGenInfo(int leptonId, std::vector<vhtm::GenParticle>& genLepList, std::vector<vhtm::GenParticle>& genTauList);
  void findLLTGenInfo( std::vector<vhtm::GenParticle>& genMuList, std::vector<vhtm::GenParticle>& genElList, 
		       std::vector<vhtm::GenParticle>& genTauList);
  int vetoMuon(double vetoPtCut, double dzTauCut);
  int vetoElectron(double vetoPtCut, double dzTauCut);
  int vetoElectronSpl(double vetoPtCut);

  bool electronMVA(const vhtm::Electron& electron);
  bool eleId(const vhtm::Electron& ele, int bx);
  bool DYtoeeVeto(TLorentzVector TOS, TLorentzVector TSS, double tauOSemfrac, vhtm::Electron ele, int elIndx);
  int GenLevelMatching(const TLorentzVector& DetObj, const std::vector<vhtm::GenParticle>& genList);
  void isGenMatchedDy(int leptonId, std::vector<vhtm::GenParticle>&, bool&);

  const std::vector<vhtm::Event>* eventColl() const {return eventList_;}
  const std::vector<vhtm::Vertex>* vertexColl() const {return vertexList_;}
  const std::vector<vhtm::GenEvent>* genEventColl() const {return genEventList_;}
  const std::vector<vhtm::Tau>* tauColl() const {return tauList_;}
  const std::vector<vhtm::Electron>* electronColl() const {return electronList_;}
  const std::vector<vhtm::Muon>* muonColl() const {return muonList_;}
  const std::vector<vhtm::Photon>* photonColl() const {return photonList_;}
  const std::vector<vhtm::PackedPFCandidate>* packedPFCandidateColl() const {return packedPFCandidateList_;}
  const std::vector<vhtm::Jet>* jetColl() const {return jetList_;}
  const std::vector<vhtm::MET>* metColl() const {return metList_;}
  const std::vector<vhtm::GenParticle>* genParticleColl() const {return genParticleList_;}
  const std::vector<vhtm::GenJet>* genJetColl() const {return genJetList_;}
  const std::vector<vhtm::GenMET>* genMetColl() const {return genMetList_;}
  const std::vector<vhtm::TriggerObject>* triggerObjColl() {return triggerObjList_;}

  const std::vector<int>* l1physbits() const {return l1physbits_;}
  const std::vector<int>* l1techbits() const {return l1techbits_;}
  const std::vector<std::string>* hltpaths() const {return hltpaths_;}
  const std::vector<int>* hltresults() const {return hltresults_;}
  const std::vector<int>* hltprescales() const {return hltprescales_;}

  int nvertex() const {return vertexList_->size();}
  int nelectron() const {return electronList_->size();}
  int nmuon() const {return muonList_->size();}
  int nphoton() const {return photonList_->size();}
  int npackedPFCandidate() const { return packedPFCandidateList_->size(); }
  int ntau() const {return tauList_->size();}
  int njet() const {return jetList_->size();}
  int nmet() const {return metList_->size();}
  int ngenparticle() const {return genParticleList_->size();}
  int ntriggerobj() const {return triggerObjList_->size();}
  int ngenjet() const {return genJetList_->size();}
  int ngenmet() const {return genMetList_->size();}

  TChain* chain() const {return chain_;}
  TFile* histf() const {return histf_;}

  int nEvents() const {return nEvents_;}
  int firstEvent() const {return firstEvt_;}
  int lastEvent() const {return lastEvt_;}
  ofstream& fLog() {return fLog_;}
  ofstream& evLog() {return evLog_;}

  bool isMC() const {return isMC_;}
  bool isSignal() const {return isSignal_;}
  int logOption() const {return logOption_;}
  bool useTrigger() const {return useTrigger_;}
  bool usePUWt() const {return usePUWt_;}
  const std::vector<std::string>& trigPathList() const {return trigPathList_;}
  bool useTrueNInt() const {return useTrueNInt_;}

  const std::map<std::string, double>& vtxCutMap() const {return vtxCutMap_;}
  const std::map<std::string, double>& muonCutMap() const {return muonCutMap_;}
  const std::map<std::string, double>& photonCutMap() const {return photonCutMap_;}
  const std::map<std::string, double>& packedPFCandidateCutMap() const {return packedPFCandidateCutMap_;}
  const std::map<std::string, double>& electronCutMap() const {return electronCutMap_;}
  const std::map<std::string, double>& tauCutMap() const {return tauCutMap_;}
  const std::map<std::string, double>& bjetCutMap() const {return bjetCutMap_;}
  const std::map<std::string, double>& jetCutMap() const {return jetCutMap_;}
  const std::map<std::string, double>& evselCutMap() const {return evselCutMap_;}
  const std::map<std::string, int>& eventIdMap() const {return eventIdMap_;}
  int bunchCrossing() const {return bunchCrossing_;}
  
public:
  double puevWt_;
  bool emt, eet, mmt;

private:
  TChain* chain_;      // chain contains a list of root files containing the same tree
  TFile* histf_;       // The output file with histograms

  // The tree branches
  std::vector<vhtm::Event>* eventList_;
  std::vector<vhtm::Vertex>* vertexList_;
  std::vector<vhtm::GenEvent>* genEventList_;
  std::vector<vhtm::Tau>* tauList_;
  std::vector<vhtm::Electron>* electronList_;
  std::vector<vhtm::Muon>* muonList_;
  std::vector<vhtm::Photon>* photonList_;
  std::vector<vhtm::PackedPFCandidate>* packedPFCandidateList_;
  std::vector<vhtm::Jet>* jetList_;
  std::vector<vhtm::MET>* metList_;
  std::vector<vhtm::GenParticle>* genParticleList_;
  std::vector<vhtm::GenJet>* genJetList_;
  std::vector<vhtm::GenMET>* genMetList_;
  std::vector<vhtm::TriggerObject>* triggerObjList_;

  std::vector<int>* l1physbits_;
  std::vector<int>* l1techbits_;
  std::vector<std::string>* hltpaths_;
  std::vector<int>* hltresults_;
  std::vector<int>* hltprescales_;

  std::vector<std::string> brList_;
  std::vector<double> puWtList_;

  int nEvents_;

  ofstream fLog_;   
  ofstream evLog_;   

  bool isMC_;
  bool isSignal_;
  bool readTrk_;
  bool readTrigObject_;
  std::vector<std::string> fileList_;
  int logOption_;
  bool useTrigger_;
  bool usePUWt_;
  std::vector<std::string> trigPathList_;

  std::string histFile_;
  std::string puHistFile_;
  bool useTrueNInt_;
  std::string logFile_;
  std::string evFile_;
  int maxEvt_;
  int bunchCrossing_;
  int firstEvt_;
  int lastEvt_;

  std::map<std::string, double> vtxCutMap_;
  std::map<std::string, double> muonCutMap_;
  std::map<std::string, double> photonCutMap_;
  std::map<std::string, double> packedPFCandidateCutMap_;
  std::map<std::string, double> electronCutMap_;
  std::map<std::string, double> tauCutMap_;
  std::map<std::string, double> bjetCutMap_;
  std::map<std::string, double> jetCutMap_;
  std::map<std::string, double> evselCutMap_;

  std::map<std::string, int> eventIdMap_;
};
#endif
