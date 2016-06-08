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
  double getEleRhoEffectiveArea(double etax);
  double getEleRhoEffectiveArea03(double etax);
  double computeMuonReliso(const vhtm::Muon& mu, const std::vector<vhtm::PackedPFCandidate>& fsrP4Vec, 
			   double vetoCone=0.01, double isoCone=0.4, bool verbose=false);
  double computeElectronReliso(const vhtm::Electron& ele, const std::vector<vhtm::PackedPFCandidate>& fsrP4Vec, double eventRho, 
			       double vetoCone=0.08, double isoCone=0.4, bool verbose=false);
  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum=0.0);
  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum=0.0);
  double pfiso(const vhtm::PackedPFCandidate& cand);
  
  double pfiso03(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum=0.0);
  double pfiso03(const vhtm::Muon& mu, double fsrPhotonEtSum=0.0);

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
  void showEfficiency(const std::string& hname, const std::string slist[], const std::string& tag);
  bool isLooseJet(const vhtm::Jet& jet);
  bool isTightJet(const vhtm::Jet& jet);
  bool fsrPhotonP4(const std::vector<vhtm::PackedPFCandidate>& fsrList, TLorentzVector& p4);
}
#endif
