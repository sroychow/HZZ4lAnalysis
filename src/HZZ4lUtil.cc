#include "configana.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 
#include <typeinfo>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AnaUtil.h"
#include "HZZ4lUtil.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::ostream;

namespace HZZ4lUtil {
  double getEleRhoEffectiveArea(double etax) {
    double area;
    double eta = std::fabs(etax);
    if (eta >= 0.0 && eta < 0.8) area = 0.1830;
    else if (eta >= 0.8 && eta < 1.3) area = 0.1734;
    else if (eta >= 1.3 && eta < 2.0) area = 0.1077;
    else if (eta >= 2.0 && eta < 2.2) area = 0.1565;
    else if (eta >= 2.2) area = 0.2680;
    return area;
  }

  double getEleRhoEffectiveArea03(double etax) {
    double area;
     double eta = std::fabs(etax);
    if (eta >= 0.0000 && 1.0000) area = 0.1752;
    else if (eta >= 1.0000 && 1.4790) area = 0.1862;
    else if (eta >= 1.4790 && 2.0000) area = 0.1411;
    else if (eta >= 2.0000 && 2.2000) area = 0.1534;
    else if (eta >= 2.2000 && 2.3000) area = 0.1903;
    else if (eta >= 2.3000 && 2.4000) area = 0.2243;
    else if (eta >= 2.4000 && 5.0000) area = 0.2687;
    return area;
  }
  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum) {
    return (ele.chargedHadronIso + std::max(0., ele.neutralHadronIso + ele.photonIso - fsrPhotonEtSum
  					    - getEleRhoEffectiveArea(std::fabs(ele.eta)) * eventRho));
  }


  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum) {
    return (mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - fsrPhotonEtSum - 0.5 * mu.sumPUPt));
  }
  double pfiso(const vhtm::PackedPFCandidate& cand) {
    return (cand.isolationMap.at("c30").at(0) 
	    + cand.isolationMap.at("c30").at(2) 
	    + cand.isolationMap.at("c30").at(3) 
	    + cand.isolationMap.at("c30").at(4));
  }
  //***********************************not EG pog now uses scETA*********************************
  double pfiso03(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum) {
    return ( ele.sumChargedHadronPt + 
             std::max(0., ele.sumNeutralHadronEt + ele.sumPhotonEt - fsrPhotonEtSum
                          - getEleRhoEffectiveArea03(std::fabs(ele.scEta)) * eventRho));
  }

  double pfiso03(const vhtm::Muon& mu, double fsrPhotonEtSum) {
    return mu.pfChargedHadIsoR03 + std::max(0.0, mu.pfNeutralHadIsoR03 + mu.pfPhotonIso03 
                                              - fsrPhotonEtSum - 0.5*mu.sumPUPt03);
  }
	
  // Isolation in new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone of all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons 
  double computeElectronReliso(const vhtm::Electron& ele, const std::vector<vhtm::PackedPFCandidate>& fsrVec, double eventRho,
			       double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = getP4(ele);
    for (const auto& v: fsrVec) {
      TLorentzVector fsrP4 = getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if ((std::fabs(ele.scEta) < 1.479 || dR > vetoCone) && dR < isoCone)
  	phoEtSum += fsrP4.Et();
    }
    double iso = 9999.;
    if( isoCone == 0.3 ) iso = pfiso03(ele, eventRho, phoEtSum);
    else iso = pfiso(ele, eventRho, phoEtSum);

    if (verbose) {
      cout << "electron isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt  effArea eventRho" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
  	   << setw(9) << ele.chargedHadronIso
  	   << setw(9) << ele.neutralHadronIso 
  	   << setw(9) << ele.photonIso
  	   << setw(9) << phoEtSum
  	   << setw(9) << getEleRhoEffectiveArea(std::fabs(lP4.Eta()))
           << setw(9) << eventRho
  	   << endl;
    }
    return iso/lP4.Pt();
  }
  // Isolation with new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // dR > 0.01 for muons
  double computeMuonReliso(const vhtm::Muon& mu, const std::vector<vhtm::PackedPFCandidate>& fsrP4Vec, 
			   double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = getP4(mu);
    for (const auto& v: fsrP4Vec) {
      TLorentzVector fsrP4 = getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if (dR > vetoCone && dR < isoCone)
  	phoEtSum += fsrP4.Et();
    }
    double iso = 9999.;
    if( isoCone == 0.3 ) iso = pfiso03(mu, phoEtSum);
    else iso = pfiso(mu, phoEtSum);

    if (verbose) {
      cout << "muon isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt     PUPt" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
  	   << setw(9) << mu.sumChargedHadronPt
  	   << setw(9) << mu.sumNeutralHadronEt
  	   << setw(9) << mu.sumPhotonEt
  	   << setw(9) << phoEtSum
  	   << setw(9) << mu.sumPUPt
  	   << endl;
    }
    
    return iso/lP4.Pt();
  }  
  bool jetpuMVAid(const vhtm::Jet& jet) {
    float jpumva = jet.jpumva;
    double pt = jet.pt;
    double eta = std::fabs(jet.eta);
    bool passPU = true;
    if (pt > 20) {
      if (eta > 3.) {
  	if (jpumva <= -0.45) passPU = false;
      }
      else if (eta > 2.75) {
  	if (jpumva <= -0.55) passPU = false;
      }
      else if (eta > 2.5) {
  	if (jpumva <= -0.6) passPU = false;
      }
      else if (jpumva <= -0.63) passPU = false;
    }
    else {
      if (eta > 3.) {
  	if (jpumva <= -0.95) passPU = false;
      }
      else if (eta > 2.75) {
  	if (jpumva <= -0.94) passPU = false;
      }
      else if (eta > 2.5) {
  	if (jpumva <= -0.96) passPU = false;
      }
      else if (jpumva <= -0.95) passPU = false;
    }
    return passPU;
  }
  bool isLooseJet(const vhtm::Jet& jet) {
    bool centralCut = (std::fabs(jet.eta) <= 2.4) 
      ? (jet.chargedHadronEnergyFraction > 0 && 
  	 jet.chargedMultiplicity > 0 && 
  	 jet.chargedEmEnergyFraction < 0.99)
      : true;
    
    return (jet.neutralHadronEnergyFraction < 0.99 && 
  	    jet.neutralEmEnergyFraction < 0.99 &&
  	    (jet.chargedMultiplicity + jet.neutralMultiplicity) > 1 &&
  	    jet.muonEnergyFraction < 0.8 &&
  	    centralCut);
  }
  bool isTightJet(const vhtm::Jet& jet) {
    bool centralCut = (std::fabs(jet.eta) <= 2.4)
      ? (jet.chargedHadronEnergyFraction > 0 && 
  	 jet.chargedMultiplicity > 0 && 
  	 jet.chargedEmEnergyFraction < 0.9)
      : true;
    
    return (jet.neutralHadronEnergyFraction < 0.9 && 
  	    jet.neutralEmEnergyFraction < 0.9 && 
  	    (jet.chargedMultiplicity + jet.neutralMultiplicity) > 1 &&
  	    jet.muonEnergyFraction < 0.8 && 
  	    centralCut);
  }
  bool electronBDT(const vhtm::Electron& electron) {
    double scEta = std::fabs(electron.scEta);
    double elePt = electron.pt;
    double BDT = electron.BDT; //preComp;
    bool isBDT = (elePt <= 10 && ((scEta < 0.8                   && BDT > -0.265)   ||
                                  (scEta >= 0.8 && scEta < 1.479 && BDT > -0.556)   ||
                                  (scEta >= 1.479                && BDT > -0.551))) ||
                 (elePt > 10  && ((scEta < 0.8                   && BDT > -0.072)   ||
                                  (scEta >= 0.8 && scEta < 1.479 && BDT > -0.286)   || 
                                  (scEta >= 1.479                && BDT > -0.267)));
    return isBDT;
  }
  void printZCandidate(const ZCandidate& Za, const std::string& tag, std::ostream& os) {
    os << std::setprecision(3);
    os << tag << endl;
    os << " -- Leptons: " << endl;
    os << "  indx      pT     eta     phi  energy  charge  reliso   fsrPt"
       << endl; 
    os << setw(6) << 1
       << setw(8) << Za.l1P4.Pt()
       << setw(8) << Za.l1P4.Eta()
       << setw(8) << Za.l1P4.Phi()
       << setw(8) << Za.l1P4.Energy()
       << setw(8) << Za.l1Charge
       << setw(8) << Za.l1Isolation
       << setw(8) << Za.l1FsrP4.Pt()
       << endl;
    os << setw(6) << 2
       << setw(8) << Za.l2P4.Pt()
       << setw(8) << Za.l2P4.Eta()
       << setw(8) << Za.l2P4.Phi()
       << setw(8) << Za.l2P4.Energy()
       << setw(8) << Za.l2Charge
       << setw(8) << Za.l2Isolation
       << setw(8) << Za.l2FsrP4.Pt()
       << endl;
    
    if (Za.fsrPhoP4.Et() > 2) {
      os << " -- FSR Photons (combined for both leptons): " << endl;
      os << "       eT     eta     phi   energy" 
  	 << endl; 
      os << setw(9) << Za.fsrPhoP4.Et()
  	 << setw(8) << Za.fsrPhoP4.Eta()
  	 << setw(8) << Za.fsrPhoP4.Phi()
  	 << setw(9) << Za.fsrPhoP4.Energy()
  	 << endl;
    }
    os << " -- Z Properties:" << endl;
    os << "  flavour    mass  massDiff" 
       << endl; 
    os << setw(9) << ((!Za.flavour) ? "muon" : "electron") 
       << setw(8) << Za.mass
       << setw(10) << Za.massDiff
       << endl << endl;
  }
  bool sameFlavourZPair(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
    if ( (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) ||
  	 (Z1Cand.flavour == ZType::ee   && Z2Cand.flavour == ZType::mumu) ) return false;
    return true;  
  }
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate& Z2, int nJets,
  			   double jet1Pt, double jet2Pt, std::ostream& os) {
    //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1:.2fX}:{mZ2:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
    //{njets30:dX}: {jet1pt:.2fX}:{jet2pt:.2fX}:{category}
    os << std::fixed << setprecision(2);
    os << run << ":"
       << lumi << ":"
       << event << ":"
       << (Z1.l1P4 + Z1.l2P4 + Z1.fsrPhoP4 + Z2.l1P4 + Z2.l2P4 + Z2.fsrPhoP4).M() << ":"
       << Z1.mass << ":"
       << Z2.mass << ":"
       << 0 << ":"
       << nJets << ":"
       << jet1Pt << ":"
       << jet2Pt
       << endl;
  }
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate& Z2, int nJets,
  			   double jet1Pt, double jet2Pt, const int category, std::ostream& os) {
    //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1:.2fX}:{mZ2:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
    //{njets30:dX}: {jet1pt:.2fX}:{jet2pt:.2fX}:{category}
    os << std::fixed << setprecision(2);
    os << run << ":"
       << lumi << ":"
       << event << ":"
       << (Z1.l1P4 + Z1.l2P4 + Z1.fsrPhoP4 + Z2.l1P4 + Z2.l2P4 + Z2.fsrPhoP4).M() << ":"
       << Z1.mass << ":"
       << Z2.mass << ":"
       << 0 << ":"
       << nJets << ":"
       << jet1Pt << ":"
       << jet2Pt << ":"
       << category
       << endl;
  }
  void syncDumper(int run, int lumi, int event, const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, int nJets,
  			   double jet1Pt, double jet2Pt, const std::map<std::string, double>& kd, int category, std::ostream& os) {
    //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1Cand:.2fX}:{mZ2Cand:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
    //{njets30:dX}:  {jet1pt:.2fX}:{jet2pt:.2fX}:{category} 
    
    os << std::fixed << setprecision(2);
    os << run << ":"
       << lumi << ":"
       << event << ":"
       << (Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 +  Z2Cand.fsrPhoP4).M() << ":"
       << Z1Cand.mass << ":"
       << Z2Cand.mass << ":";
    os << std::fixed << setprecision(3);
    os << kd.find("D_bkg_kin")->second << ":"
       << kd.find("D_bkg")->second  << ":"
       << kd.find("Dgg10_VAMCFM")->second << ":"
       << kd.find("Djet_VAJHU")->second << ":"
       << kd.find("D_g4")->second << ":";
    os << std::fixed << setprecision(2);
    os << nJets <<  ":"
       << jet1Pt << ":"
       << jet2Pt << ":"
       << category 
       << endl;
  }
  void showEfficiency(const string& hname, const string slist[], const string& tag) {
    cout << ">>> " << tag << " Efficiency" << endl;
    TH1 *h = AnaUtil::getHist1D(hname);
    if (h) {
      cout << setw(64) << "CutFlow"
  	   << setw(10) << "Events"
  	   << setw(10) << "AbsEff"
  	   << setw(10) << "RelEff"
  	   << endl;
      cout << setprecision(3);
      int nbins = h->GetNbinsX();
      for (int i = 1; i <= nbins; ++i)
  	cout << setw(64) << slist[i-1]
  	     << setprecision(1) 
  	     << setw(10) << int(h->GetBinContent(i))
  	     << setprecision(5) 
  	     << setw(10) << ((h->GetBinContent(1)>0) ? h->GetBinContent(i)/h->GetBinContent(1) : 0.0)
  	     << setw(10) << ( i == 1 ? 1.0 : (h->GetBinContent(i-1)>0) ? h->GetBinContent(i)/h->GetBinContent(i-1) : 0.0)
  	     << endl;
    }
  }
  void printP4(const TLorentzVector& lv, const string& tag) {
    cout << setprecision(3);
    cout << tag << " = (" 
  	 << setw(7) << lv.Pt()  << "," 
  	 << setw(7) << lv.Eta() << "," 
  	 << setw(7) << lv.Phi() << "," 
  	 << setw(7) << lv.Energy() << ")" 
  	 << endl;
  }
  bool fsrPhotonP4(const vector<vhtm::PackedPFCandidate>& fsrList, TLorentzVector& p4) {
    bool hasFsr = false;
    p4.SetPtEtaPhiE(0., 0., 0., 0);  
    if (!fsrList.empty()) {
      const vhtm::PackedPFCandidate& cand = fsrList[0]; // only the first one is important
      p4.SetPtEtaPhiE(cand.pt, cand.eta, cand.phi, cand.energy);
      hasFsr = true;
    }
    return hasFsr;
  }
}
