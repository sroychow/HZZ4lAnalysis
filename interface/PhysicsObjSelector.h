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
#include "HZZ4lUtil.h"

class PhysicsObjSelector: public AnaBase {
 public:
  PhysicsObjSelector();
  virtual ~PhysicsObjSelector() {}
  virtual void bookHistograms();
  virtual bool beginJob();

  const std::vector<vhtm::PackedPFCandidate>& getFSRPhotonVec() const {return fsrPhotonVec_;}
  const std::vector<vhtm::Muon>& getPreSIPLooseMuList() const {return preSIPlooseMuVec_;}
  const std::vector<vhtm::Electron>& getPreSIPLooseEleList() const {return preSIPlooseEleVec_;}
  
  const std::vector<vhtm::Muon>& getLooseMuList() const {return looseMuVec_;}
  const std::vector<vhtm::Muon>& getTightMuList() const {return tightMuVec_;}
  const std::vector<vhtm::Muon>& getTightIsoMuList() const {return tightIsoMuVec_;}
  
  const std::vector<vhtm::Electron>& getLooseEleList() const {return looseEleVec_;}
  const std::vector<vhtm::Electron>& getTightEleList() const {return tightEleVec_;}
  const std::vector<vhtm::Electron>& getTightIsoEleList() const {return tightIsoEleVec_;}
  
  const std::vector<vhtm::Jet>& getLeptonCleanedLooseJetList() const {return looseJetVec_;}
  const std::vector<vhtm::Jet>& getLeptonCleanedTightJetList() const {return tightJetVec_;}
  
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getLooseElePhotonPairList() const {return looseElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getTightElePhotonPairList() const {return tightElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getTightIsoElePhotonPairList() const {return tightIsoElePhotonPairVec_;}
  
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getLooseMuPhotonPairList() const {return looseMuPhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getTightMuPhotonPairList() const {return tightMuPhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getTightIsoMuPhotonPairList() const {return tightIsoMuPhotonPairVec_;}
  
  void setTightElePhotonPairList(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate>>>& list) {tightElePhotonPairVec_ = list;}
  int getNLooseJets() const {return looseJetVec_.size();}
  int getNTightJets() const {return tightJetVec_.size();}
  int getNbJets() const {return nbJets_;} 
  
  void findObjects(double wt=1);
  void muonSelector(double wt=1);
  void electronSelector(double wt=1);
  void photonSelector();
  void isoLeptonSelector();
  
  void jetSelector(double wt=1);
  bool jetLeptonCleaning(const vhtm::Jet& jet, double dR) const;
  void setEventGridRho(double evRho) {fGridRhoFastjetAll_ = evRho;}
  double getEventGridRho() const {return fGridRhoFastjetAll_;}
  void clear();
  
  bool passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  bool passedSuperClusterVetobyReference(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  double findClosestLepton(const vhtm::PackedPFCandidate& photon, int& muindx, int& elindx) const;
  void leptonCrossCleaning();
  bool crossCleaned(const vhtm::Electron& electron) const;
  void ZZMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec);
  void dumpEvent(bool dumpGen=false, bool showEvent=false, std::ostream& os=std::cout) const;
  
  void addLeptonIsolation(std::vector<ZCandidate>& candList, 
			  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& elePhotonPairVec,  
			  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& muPhotonPairVec);
  // create unique lepton pair combination giving a Z and add the candidate into a vector
  template <typename T>
    void ZSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec, std::vector<ZCandidate>& candList) {
    for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
      const auto& ip = lepPhotonPairVec[i];
      TLorentzVector lep1P4 = HZZ4lUtil::getP4(ip.first);
      
      TLorentzVector lep1FsrP4;
      bool lep1HasFsr = HZZ4lUtil::fsrPhotonP4(ip.second, lep1FsrP4); 
      
      for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
	const auto& jp = lepPhotonPairVec[j];
	
	// opposite charge
	if ((ip.first.charge + jp.first.charge) != 0) continue; 
	
	TLorentzVector lep2P4 = HZZ4lUtil::getP4(jp.first);
	TLorentzVector lep2FsrP4;
	bool lep2HasFsr = HZZ4lUtil::fsrPhotonP4(jp.second, lep2FsrP4); 
	
	ZCandidate ztmp;
	if (typeid(jp.first) == typeid(vhtm::Muon)) {
	  ztmp.flavour = HZZ4lUtil::ZType::mumu;
	}
	else if (typeid(jp.first) == typeid(vhtm::Electron)) {
	  ztmp.flavour = HZZ4lUtil::ZType::ee;
	}
	else 
	  ztmp.flavour = -1;
	
	ztmp.l1Index = i;
	ztmp.l1P4 = lep1P4;
	ztmp.l1Charge = ip.first.charge;
	ztmp.l1FsrP4 = lep1FsrP4;
	
	ztmp.l2Index = j;
	ztmp.l2P4 = lep2P4;
	ztmp.l2Charge = jp.first.charge;
	ztmp.l2FsrP4 = lep2FsrP4;
	
	int whichLep = 0;
	if (lep1HasFsr && lep2HasFsr) {
	  whichLep = 3;
	}
	else if (lep1HasFsr) {
	  whichLep = 1;
	}
	else
	  whichLep = 2;
	ztmp.fsrWithLep = whichLep; 
	ztmp.fsrPhoP4 = lep1FsrP4 + lep2FsrP4;
	
	double Zmass = (lep1P4 + lep2P4 + lep1FsrP4 + lep2FsrP4).M();
	ztmp.mass = Zmass;
	ztmp.massDiff = std::fabs(Zmass - HZZ4lUtil::MZnominal);
	
	candList.push_back(ztmp);
      }
    }
  }
  // Z+ll (SS) CR, (method AA), made by relaxing selection requirement on the two additional leptons; 
  // they are required to be same-sign but to pass only loose ID (see above) and SIP
  // (i.e. no cut on iso, no tight muon/tight electron requirements) 
  template <class T> 
    void leptonPairSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec,
			    ZCandidate& Z, bool studyOSPair, std::vector<ZCandidate>& candList, 
			    std::vector<std::pair<ZCandidate, ZCandidate> >& objPairList) {
    for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
      const auto& ip = lepPhotonPairVec[i];
      TLorentzVector lep1P4 = HZZ4lUtil::getP4(ip.first);
      if (lep1P4 == Z.l1P4 || lep1P4 == Z.l2P4) continue;  // Keep aside the lepton that forms the real Z
      
      TLorentzVector lep1FsrP4;
      bool lep1HasFsr = HZZ4lUtil::fsrPhotonP4(ip.second, lep1FsrP4); 
      
      for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
	const auto& jp = lepPhotonPairVec[j];
	
	if (studyOSPair) {
	  if (ip.first.charge + jp.first.charge != 0) continue; // opposite charge
	}
	else {
	  if (ip.first.charge != jp.first.charge) continue;     // same charge
	}
	TLorentzVector lep2P4 = HZZ4lUtil::getP4(jp.first);
	if (lep2P4 == Z.l1P4 || lep2P4 == Z.l2P4) continue;     // Keep aside the lepton that forms the real Z
	
	TLorentzVector lep2FsrP4;
	bool lep2HasFsr = HZZ4lUtil::fsrPhotonP4(jp.second, lep2FsrP4); 
	
	ZCandidate ztmp;
	if (typeid(jp.first) == typeid(vhtm::Muon))
	  ztmp.flavour = HZZ4lUtil::ZType::mumu;
	else if (typeid(jp.first) == typeid(vhtm::Electron))
	  ztmp.flavour = HZZ4lUtil::ZType::ee;
	else 
	  ztmp.flavour = -1;
	
	ztmp.l1Index = i;
	ztmp.l1P4 = lep1P4;
	ztmp.l1Charge = ip.first.charge;
	ztmp.l1FsrP4 = lep1FsrP4;
	
	ztmp.l2Index = j;
	ztmp.l2P4 = lep2P4;
	ztmp.l2Charge = jp.first.charge;
	ztmp.l2FsrP4 = lep2FsrP4;
	
	int whichLep = 0;
	if (lep1HasFsr && lep2HasFsr) {
	  whichLep = 3;
	}
	else if (lep1HasFsr) {
	  whichLep = 1;
	}
	else
	  whichLep = 2;
	ztmp.fsrWithLep = whichLep; 
	ztmp.fsrPhoP4 = lep1FsrP4 + lep2FsrP4;
	double Zmass = (lep1P4 + lep2P4 + lep1FsrP4 + lep2FsrP4).M();
	ztmp.mass = Zmass;
	ztmp.massDiff = std::fabs(Zmass - HZZ4lUtil::MZnominal);
	
	candList.push_back(ztmp);
	objPairList.push_back({Z, ztmp});
      }
    }
  }
 private:
  bool dumpEvent_;
  std::vector<vhtm::Muon> preSIPlooseMuVec_, looseMuVec_, tightMuVec_,tightIsoMuVec_;
  std::vector<vhtm::Electron> preSIPlooseEleVec_, looseEleVec_, tightEleVec_,tightIsoEleVec_;
  std::vector<vhtm::PackedPFCandidate> fsrPhotonVec_;
  std::vector<vhtm::Jet> looseJetVec_, tightJetVec_;
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > looseElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > tightElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > tightIsoElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > looseMuPhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > tightMuPhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > tightIsoMuPhotonPairVec_; 
  double fGridRhoFastjetAll_;
  bool searchedEle_, 
    searchedMu_, 
    searchedPhoton_;
  int nbJets_;
};
#endif
