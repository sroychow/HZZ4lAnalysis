#ifndef __MassResolutionAnalyzer__hh 
#define __MassResolutionAnalyzer__hh
#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

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
#include "TTree.h"
#include "TFile.h"

using std::string;
struct mresVar {
  double m;
  double merr;
  double eta1;
  double eta2;
  double pT1;
  double pT2;
  double w;
};
class MassResolutionAnalyzer : public PhysicsObjSelector {
    
public:

  MassResolutionAnalyzer();
  virtual ~MassResolutionAnalyzer();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  bool readJob(const string& jobFile, int& nFiles);  
  void clearLists();
  void selectEvent(){};
  virtual void bookHistograms();
  void getZmmMassInfo();
  void getZeeMassInfo();

  std::vector<vhtm::Vertex> vtxList;  
  bool _dumpEvent;
  //bool useTightLepton_;
  std::vector<ZCandidate> Z;
  std::vector<vhtm::Muon> tightMuVec_;
  std::vector<vhtm::Electron> tightEleVec_;
  TTree* massRtree_;
  mresVar mr;
  TFile* tfout_;
  std::string treeFileName_;
};
#endif
