#ifndef __ZTnpAnalysis__hh
#define __ZTnpAnalysis__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "PhysicsObjSelector.h"
#include "ZCandidate.h"

using std::ofstream;
using std::ifstream;
using std::ostream;
class ZTnpAnalysis: public PhysicsObjSelector {
    
public:
  enum EventType {
    mmmm = 0, eeee, eemm
  };
  enum MassRange {
    low = 0, high, fullblind, unblind 
  };
  ZTnpAnalysis();
  virtual ~ZTnpAnalysis();
  
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;

  void selectEvent() {};
  virtual void bookHistograms();
  void getZtnPpair();  
public:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >   tightMuon_;
  std::vector<vhtm::PackedPFCandidate>  fsrPVec_;
  std::vector<vhtm::TriggerObject> tObj_;
  std::string treeFname_;
  TTree* outTree_;
  TFile* outTreeFile_;
  vhtm::ZtnP  zcand;
};
#endif
