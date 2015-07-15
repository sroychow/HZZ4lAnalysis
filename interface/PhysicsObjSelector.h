#ifndef __PhysicsObjSelector__hh
#define __PhysicsObjSelector__hh

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
#include "AnaBase.h"
#include "ZCandidate.h"

class PhysicsObjSelector: public AnaBase {
 public:
  PhysicsObjSelector();
  virtual ~PhysicsObjSelector() {}
  virtual void bookHistograms();
  virtual bool beginJob();
  
  const std::vector<vhtm::Muon>& getLooseMuList() const {return looseMuVec_;}
  const std::vector<vhtm::Muon>& getTightMuList() const {return tightMuVec_;}
  const std::vector<vhtm::Electron>& getLooseEleList() const {return looseEleVec_;}
  const std::vector<vhtm::Electron>& getTightEleList() const {return tightEleVec_;}
  const std::vector<vhtm::Jet>& getLeptonCleanedLooseJetList() const {return looseJetVec_;}
  const std::vector<vhtm::Jet>& getLeptonCleanedTightJetList() const {return tightJetVec_;}
  const std::vector<vhtm::PackedPFCandidate>& getFSRPhotonList() const {return photonVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getLooseElePhotonPairList() const {return looseElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getTightElePhotonPairList() const {return tightElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getLooseMuPhotonPairList() const {return looseMuPhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getTightMuPhotonPairList() const {return tightMuPhotonPairVec_;}
  void setTightElePhotonPairList(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& list) {tightElePhotonPairVec_ = list;}
  int getNLooseJets() const {return looseJetVec_.size();}
  int getNTightJets() const {return tightJetVec_.size();}
  int getNbJets() const {return nbJets_;} 
  
  void findObjects(double wt=1);
  void muonSelector(double wt=1);
  void electronSelector(double wt=1);
  void photonSelector();
  void jetSelector(double wt=1);
  bool jetLeptonCleaning(const vhtm::Jet& jet, double dR) const;
  void setEventGridRho(double evRho) {fGridRhoFastjetAll_ = evRho;}
  double getEventGridRho() const {return fGridRhoFastjetAll_;}
  void clear();

  bool passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  double findClosestLepton(const vhtm::PackedPFCandidate& photon, int& muindx, int& elindx) const;
  void leptonCrossCleaning();
  bool crossCleaned(const vhtm::Electron& electron) const;
  bool ZmumuIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.01, double isoCone=0.4, bool verbose=false) const;
  bool ZeeIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.08, double isoCone=0.4, bool verbose=false) const;
  void ZZisoMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec);
  void dumpEvent(bool dumpGen=false, bool showEvent=false, ostream& os=std::cout) const;

 private:
  bool dumpEvent_;
  std::vector<vhtm::Muon> looseMuVec_, tightMuVec_;
  std::vector<vhtm::Electron> looseEleVec_, tightEleVec_;
  std::vector<vhtm::PackedPFCandidate> photonVec_;
  std::vector<vhtm::Jet> looseJetVec_, tightJetVec_;
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > looseElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > tightElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > looseMuPhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > tightMuPhotonPairVec_; 
  double fGridRhoFastjetAll_;
  bool searchedEle_, 
    searchedMu_, 
    searchedPhoton_;
  int nbJets_;
};
#endif
