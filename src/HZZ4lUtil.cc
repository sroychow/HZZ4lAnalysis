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

namespace HZZ4lUtil {
  double getEleRhoEffectiveArea(double eta) {
    double area;
    if (eta >= 0.0 && eta < 0.8) area = 0.1830;
    else if (eta >= 0.8 && eta < 1.3) area = 0.1734;
    else if (eta >= 1.3 && eta < 2.0) area = 0.1077;
    else if (eta >= 2.0 && eta < 2.2) area = 0.1565;
    else if (eta >= 2.2) area = 0.2680;
    return area;
  }
  double pfiso(const vhtm::Electron& ele, double eventRho, double fsrPhotonEtSum) {
    return (ele.chargedHadronIso + std::max(0., ele.neutralHadronIso + ele.photonIso  - fsrPhotonEtSum
  					    - getEleRhoEffectiveArea(std::fabs(ele.eta)) * eventRho));
  }
  double pfiso(const vhtm::Muon& mu, double fsrPhotonEtSum) {
    return (mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - fsrPhotonEtSum - 0.5 * mu.sumPUPt));
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
    double BDT = electron.BDT;
    bool isBDT = (elePt <= 10 && ((scEta < 0.8 && BDT > -0.202)                     ||
  				  ((scEta >= 0.8 && scEta < 1.479) && BDT > -0.444) ||
  				  (scEta >= 1.479 && BDT > 0.264)))                 ||
                 (elePt >  10 && ((scEta < 0.8 && BDT > -0.110)                     ||
  				  ((scEta >= 0.8 && scEta < 1.479) && BDT > -0.284) ||
  				  (scEta >= 1.479 && BDT > -0.212)));
    return isBDT;
  }
  void printZCandidate(const ZCandidate& Za, const std::string& tag, std::ostream& os) {
    os << std::setprecision(3);
    os << tag << endl;
    os << " -- Leptons: " << endl;
    os << "  indx      pT     eta     phi  energy  charge  reliso"
       << endl; 
    os << setw(6) << 1
       << setw(8) << Za.l1P4.Pt()
       << setw(8) << Za.l1P4.Eta()
       << setw(8) << Za.l1P4.Phi()
       << setw(8) << Za.l1P4.Energy()
       << setw(8) << Za.l1Charge
       << setw(8) << Za.l1Isolation
       << endl;
    os << setw(6) << 2
       << setw(8) << Za.l2P4.Pt()
       << setw(8) << Za.l2P4.Eta()
       << setw(8) << Za.l2P4.Phi()
       << setw(8) << Za.l2P4.Energy()
       << setw(8) << Za.l2Charge
       << setw(8) << Za.l2Isolation
       << endl;
    
    if (Za.fsrPhoP4.Et() > 2) {
      os << " -- FSR Photon: " << endl;
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
  // -- from twiki --
  // for each photon attached to a lepton of each Z, test if the ll(g) candidate would satisfy the two criteria:
  // 4 < m(llg))) < 100 GeV
  // |m(llg)))-91.1876| < |m(ll)-91.1876| 
  void selectFSRPhoforLepton(const TLorentzVector& lep1P4, const TLorentzVector& lep2P4,
  					const vector<vhtm::PackedPFCandidate>& lepprePhoVec,
  					vector<vhtm::PackedPFCandidate>& lepPhoVec) 
  {
    for (const auto& v: lepprePhoVec) {
      const TLorentzVector& pfPhoP4 = getP4(v);
      double mass = (lep1P4 + lep2P4 + pfPhoP4).M();
      if ((mass > 4. && mass < 100.)
  	  && std::fabs(mass - MZnominal) < std::fabs((lep1P4 + lep2P4).M() - MZnominal))
  	lepPhoVec.push_back(v);
    }
    if (lepPhoVec.size() > 1) 
      std::sort(lepPhoVec.begin(), lepPhoVec.end(), PtComparator<vhtm::PackedPFCandidate>());
  }
  // -- from twiki --
  // if within the same Z multiple photons satisfy these criteria, pick the best one according to this logic:
  //    - if there's at least one photon with pT > 4 GeV, pick the one with highest pT
  //    - if all photons have pT < 4 GeV, pick the one that has the smallest dR to its closest lepton. 
  int selectBestFSRforLepton(const TLorentzVector& lepP4, const vector<vhtm::PackedPFCandidate>& lepPhoVec) {
    int index = -1;
    if (lepPhoVec[0].pt > 4.) index = 0;
    else {
      double dRmin = 999.;
      for (unsigned int i = 0; i < lepPhoVec.size(); ++i) {
  	const TLorentzVector& phoP4 = getP4(lepPhoVec[i]);
  	double dR = phoP4.DeltaR(lepP4);
  	if (dR < dRmin) {
  	  index = i;
  	  dRmin = dR;
  	}
      }
    }
    return index;
  }
  void addFSRtoAltZ(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, 
  		    int l1indx, int l2indx, TLorentzVector& altZP4, const std::string& tag, bool verbose) {
    if (Z1Cand.fsrWithLep == l1indx) {
      if (std::fabs((altZP4+Z1Cand.fsrPhoP4).M() - MZnominal) < std::fabs(altZP4.M() - MZnominal)) {
  	if (verbose) cout << "--- SmartCut: adding Photon, l1indx: " << l1indx 
  			  << ", Photon pT: " << Z1Cand.fsrPhoP4.Pt() << tag 
  			  << endl;
  	altZP4 += Z1Cand.fsrPhoP4;
      }
    }
    if (Z2Cand.fsrWithLep == l2indx) {
      if (std::fabs((altZP4+Z2Cand.fsrPhoP4).M() - MZnominal) < std::fabs(altZP4.M() - MZnominal)) {
  	if (verbose) cout << "--- SmartCut: adding Photon, l2indx: " << l2indx 
  			  << ", Photon pT: " << Z2Cand.fsrPhoP4.Pt() << tag 
  			  << endl;
  	altZP4 += Z2Cand.fsrPhoP4;
      }
    }
  }
  // -- from twiki --
  // For each FSR photon that was selected, exclude that photon from the photon isolation sum of all four leptons of both Zs
  // (Note that in Run I we were instead doing this only for the Z the photon is attached to) if it was in the isolation cone
  // and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons;
  // note: these requirements should probably be rechecked for consistency with isolation algorithms)
  double computeElectronIso(const vhtm::Electron& ele, const ZCandidate& Za, const TLorentzVector& lP4,
  				       const std::vector<TLorentzVector>& fsrVec, double eventRho,
  				       double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    for (const auto& v: fsrVec) {
      if ((std::fabs(ele.scEta) < 1.479 || lP4.DeltaR(v) > vetoCone) && lP4.DeltaR(v) < isoCone)
  	phoEtSum += v.Et();
    }
    double iso = pfiso(ele, eventRho, phoEtSum);
    if (verbose)
      cout << "lepton isolation: " << iso << ", " << lP4.Pt() << ", "
  	   << ele.chargedHadronIso << ", "
  	   << ele.neutralHadronIso << ", "
  	   << ele.photonIso << ", "
  	   << phoEtSum << ", " 
  	   << getEleRhoEffectiveArea(std::fabs(lP4.Eta())) << ", " 
  	   << eventRho
  	   << endl;
    
    return iso/lP4.Pt();
  }
  // -- from twiki --
  // For each FSR photon that was selected, exclude that photon from the photon isolation sum of all four leptons of both Zs
  // (Note that in Run I we were instead doing this only for the Z the photon is attached to) if it was in the isolation cone
  // and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons;
  // note: these requirements should probably be rechecked for consistency with isolation algorithms)
  double computeMuonIso(const vhtm::Muon& mu, const ZCandidate& Za, const TLorentzVector& lP4,
  			const std::vector<TLorentzVector>& fsrVec, 
  			double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    for (const auto& v: fsrVec) {
      if (lP4.DeltaR(v) > vetoCone && lP4.DeltaR(v) < isoCone)
  	phoEtSum += v.Et();
    }
    double iso = pfiso(mu, phoEtSum);
    if (verbose) 
      cout << "l1 isolation: " << iso << ", " << lP4.Pt() << ", "
  	   << mu.sumChargedHadronPt << ", "
  	   << mu.sumNeutralHadronEt << ", "
  	   << mu.sumPhotonEt << ", "
  	   << phoEtSum << ", " 
  	   << mu.sumPUPt
  	   << endl;
    
    return iso/lP4.Pt();
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
  	     << setw(10) << ( i == 1?1.0:(h->GetBinContent(i-1)>0) ? h->GetBinContent(i)/h->GetBinContent(i-1) : 0.0)
  	     << endl;
    }
  }
}
