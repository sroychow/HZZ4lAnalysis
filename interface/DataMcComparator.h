#ifndef __DataMcComparator__hh
#define __DataMcComparator__hh

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
class DataMcComparator: public PhysicsObjSelector {
    
public:
  enum EventType {
    mmmm = 0, eeee, eemm
  };
  DataMcComparator();
  virtual ~DataMcComparator();
  
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;

  void selectEvent() {};
  virtual void bookHistograms();
    
  template <typename T>
    void ZSelector(const T& tagLep, const std::vector<T>& leptonVec, std::vector<ZCandidate>& candList);
    void MassComparatorElectron(const vhtm::Electron& tagLep, const vhtm::Electron& probeLep);
    void MassComparatorMuon(const vhtm::Muon& tagLep, const vhtm::Muon& probeLep );
    int getTagElectronidx(const std::vector<vhtm::Electron>& elVec);
    int getTagMuonidx(const std::vector<vhtm::Muon>& muVec);
public:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<std::pair<ZCandidate, ZCandidate> > ZZPairVec_;

  std::vector<ZCandidate> ZCandList_;
  std::vector<vhtm::GenParticle> genZList_;
  bool checkGen_;
  bool dumpGenInfo_;
  bool useEventList_;
  bool selectEvType_;
  int evtype_;
  ofstream syncDumpf_;
  bool doKDcalc_;
  std::string dumpFilename_;
};
#endif
