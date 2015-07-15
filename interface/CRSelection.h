#ifndef __CRSelection__hh
#define __CRSelection__hh

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
#include "AnaBase.h"
#include "PhysicsObjSelector.h"
#include "ZCandidate.h"

class CRSelection: public PhysicsObjSelector {
    
public:
  enum EventType {
    mmmm = 0, eeee, eemm
  };
  CRSelection();
  virtual ~CRSelection();
  
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;

  void selectEvent() {};
  virtual void bookHistograms();
  
  template <class T>
    void ZSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec);
  
  template <class T> 
    void leptonPairSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >&, ZCandidate& Z, 
			    bool studyOSPair, std::vector<ZCandidate>& candList, std::vector<std::pair<ZCandidate, ZCandidate> >& objPairList);
 
  void ZZselector();
  void finalZllSelector(std::vector<std::pair<ZCandidate, ZCandidate> >& objPairList, bool studyOSPair, int run, int lumi, int event);

  bool CRSelectorZSSll(ZCandidate& Z, ZCandidate& ssll, bool verbose=false);
  bool CRSelectorZOSll(ZCandidate& Z, ZCandidate& osll, bool verbose=false);

public:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<std::pair<ZCandidate, ZCandidate> > ZZPairVec_;
  std::vector<std::pair<ZCandidate, ZCandidate> > ZssllPairVec_;
  std::vector<std::pair<ZCandidate, ZCandidate> > ZosllPairVec_;

  std::vector<ZCandidate> ZCandList_;
  std::vector<ZCandidate> SSllpairCandList_;
  std::vector<ZCandidate> OSllpairCandList_;
  
  bool dumpGenInfo_;
  bool useEventList_;
  bool selectEvType_;
  int evtype_;
  ofstream syncDumpf_;
};
#endif

