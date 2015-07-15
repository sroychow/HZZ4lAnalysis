#ifndef __EventSelection__hh
#define __EventSelection__hh

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

class EventSelection: public PhysicsObjSelector {
    
public:
  enum EventType {
    mmmm = 0, eeee, eemm
  };
  EventSelection();
  virtual ~EventSelection();
  
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;

  void selectEvent() {};
  virtual void bookHistograms();
    
  template <typename T>
    void ZSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec, std::vector<ZCandidate>& candList);

  void ZZselector();
  void finalZZSelector(int run, int lumi, int event);
  int findExtraLeptons(const ZCandidate& Z1, const ZCandidate& Z2);
  int findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
			const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, bool verbose=true);
  // Functions for gen Level
  int getGenDauPgd(const vhtm::GenParticle& gp);
  double getHmassfromZdau(const vhtm::GenParticle& Z1, const vhtm::GenParticle& Z2);
  bool genOk();
  static bool hasJetPair(const std::vector<vhtm::Jet>& jetList);
  static void ZZkbg(const ZCandidate& Z1, const ZCandidate& Z2, const TLorentzVector& jet1P4, const TLorentzVector& jet2P4,
		    int nJets, std::map<std::string, double>& kd);

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
