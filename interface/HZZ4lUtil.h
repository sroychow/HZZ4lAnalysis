#ifndef __HZZ4lUtil__hh
#define __HZZ4lUtil__hh

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"
#include "ZCandidate.h"

namespace HZZ4lUtil {
  const double MZnominal = 91.1876;

  enum ZType {
    mumu = 0, ee, wrong
  };
  
  // must be defined inside the header to be effective
  template <class T> 
  TLorentzVector getP4(const T& obj) {
    TLorentzVector lv;
    lv.SetPtEtaPhiE(obj.pt, obj.eta, obj.phi, obj.energy);
    return lv;
  }
  double getEleRhoEffectiveArea(double eta);
  double computeMuonIso(const vhtm::Muon& mu, const ZCandidate& Za, const TLorentzVector& lP4,
			const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.01, double isoCone=0.4, bool verbose=false);
  double computeElectronIso(const vhtm::Electron& ele, const ZCandidate& Za, const TLorentzVector& lP4,
			    const std::vector<TLorentzVector>& fsrVec, double eventRho, 
			    double vetoCone=0.08, double isoCone=0.4, bool verbose=false);
  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum=0.0);
  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum=0.0);
  bool jetpuMVAid(const vhtm::Jet& jet);
  void printP4(const TLorentzVector& lv, const std::string& tag);
  bool sameFlavourZPair(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand);
  bool electronBDT(const vhtm::Electron& electron);
  void printZCandidate(const ZCandidate& Za, const std::string& tag, std::ostream& os=std::cout);
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate& Z2, int nJets,
		  double jet1Pt, double jet2Pt, std::ostream& os=std::cout);
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate& Z2, int nJets,
		  double jet1Pt, double jet2Pt, int category, std::ostream& os=std::cout); 
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, int nJets,
		  double jet1Pt, double jet2Pt,  const std::map<std::string, double>& kd, int category, std::ostream& os=std::cout);
  void selectFSRPhoforLepton(const TLorentzVector& lep1P4, const TLorentzVector& lep2P4,
			     const std::vector<vhtm::PackedPFCandidate>& lep1prePhoVec,
			     std::vector<vhtm::PackedPFCandidate>& lep1PhoVec);
  int selectBestFSRforLepton(const TLorentzVector& lepP4, const std::vector<vhtm::PackedPFCandidate>& lepPhoVec);
  void addFSRtoAltZ(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, 
		    int l1indx, int l2indx, TLorentzVector& altZP4, const std::string& tag, bool verbose=false);
  void showEfficiency(const std::string& hname, const std::string slist[], const std::string& tag);
  bool isLooseJet(const vhtm::Jet& jet);
  bool isTightJet(const vhtm::Jet& jet);
}
#endif
