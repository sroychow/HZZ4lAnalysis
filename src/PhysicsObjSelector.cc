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

#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector2.h"

#include "AnaUtil.h"
#include "HZZ4lUtil.h"
#include "PhysicsObjSelector.h"

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

PhysicsObjSelector::PhysicsObjSelector()
  : AnaBase(),
    searchedEle_(false),
    searchedMu_(false),
    searchedPhoton_(false),
    nbJets_(0)
{
}
bool PhysicsObjSelector::beginJob() {
  AnaBase::beginJob();
  histf()->cd();
  histf()->mkdir("ObjectSelection");
  return true;
}
void PhysicsObjSelector::bookHistograms() {
  histf()->cd();
  histf()->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 10, -0.5, 10.5);
  new TH1D("eleCutFlow", "Electron Cut Flow", 9, -0.5, 8.5);
  new TH1D("jetCutFlow", "Jet Cut Flow", 7, -0.5, 6.5);
  new TH1D("nJets", "Number of jets cleaned w.r.t tight leptons passing iso per event", 20, 0, 20);
  new TH1D("btagfactor", "bDiscriminator", 50, 0., 2.);
  new TH1D("nbJets", "Number of b-jets per event", 20, 0, 20);

  new TH1F("muPt", "Muon p_{T}", 200, 0., 200.);
  new TH1F("elePt", "Electron p_{T}", 200, 0., 200.);
  new TH1F("phoPt", "Photon p_{T}", 200, 0., 200.);
  new TH1F("jetPt", "Jet p_{T}", 200, 0., 200.);
}
// muon selection
void PhysicsObjSelector::muonSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  for (const auto& muon: *muonColl()) {
    AnaUtil::fillHist1D("muCutFlow", 0, wt);

    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt"))                             continue;
    AnaUtil::fillHist1D("muCutFlow", 1, wt);

    if (std::fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                continue;
    AnaUtil::fillHist1D("muCutFlow", 2, wt);

    if (std::fabs(muon.dxyPV) >= AnaUtil::cutValue(muonCutMap(), "dxyPV") )           continue;
    AnaUtil::fillHist1D("muCutFlow", 3, wt);

    if (std::fabs(muon.dzPV) >= AnaUtil::cutValue(muonCutMap(), "dzPV") )             continue;
    AnaUtil::fillHist1D("muCutFlow", 4, wt);

    bool muType = muon.isGlobalMuon || (muon.isTrackerMuon && muon.matches > 0);
    if (!muType)                                                                      continue;
    AnaUtil::fillHist1D("muCutFlow", 5, wt);

    if (muon.muonBestTrackType == 2)                                                  continue;
    AnaUtil::fillHist1D("muCutFlow", 6, wt);

    if (!muon.isghostCleaned)                                                         continue;
    AnaUtil::fillHist1D("muCutFlow", 7, wt);

    preSIPlooseMuVec_.push_back(muon);

    if (std::fabs(muon.dB3D/muon.edB3D) >= AnaUtil::cutValue(muonCutMap(), "SIP3D"))  continue;
    AnaUtil::fillHist1D("muCutFlow", 8, wt);

    // loose muon 
    looseMuVec_.push_back(muon);

    // tight muon 
    bool isTight = muon.isPFMuon || (muon.passTrackerhighPtid && muon.pt > 200.);
    if (isTight) {
      AnaUtil::fillHist1D("muCutFlow", 9, wt);
      tightMuVec_.push_back(muon);
    }
  }
  if (preSIPlooseMuVec_.size() > 1)
    std::sort(preSIPlooseMuVec_.begin(), preSIPlooseMuVec_.end(), PtComparator<vhtm::Muon>());
  if (looseMuVec_.size() > 1)
    std::sort(looseMuVec_.begin(), looseMuVec_.end(), PtComparator<vhtm::Muon>());

  // attach an empty vector<pair<vhtm::Muon, vector<vhtm::PackedPFCandidate>>> for each loose muon
  for (auto const& muon: looseMuVec_) {
    std::vector<vhtm::PackedPFCandidate> phov;
    looseMuPhotonPairVec_.push_back({muon, phov});
  }

  if (tightMuVec_.size() > 1)           
    std::sort(tightMuVec_.begin(), tightMuVec_.end(), PtComparator<vhtm::Muon>());

  searchedMu_ = true;
}
// electron selecton
void PhysicsObjSelector::electronSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  if (!searchedMu_) muonSelector();

  for (const auto& electron: *electronColl()) {
    AnaUtil::fillHist1D("eleCutFlow", 0, wt);
    
    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt"))                                continue;
    AnaUtil::fillHist1D("eleCutFlow", 1, wt);

    if (std::fabs(electron.eta) >= AnaUtil::cutValue(electronCutMap(), "eta"))                   continue;
    AnaUtil::fillHist1D("eleCutFlow", 2, wt);

    if (std::fabs(electron.dxyPV) >= AnaUtil::cutValue(electronCutMap(), "dxyPV"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 3, wt);

    if (std::fabs(electron.dzPV) >= AnaUtil::cutValue(electronCutMap(), "dzPV"))                 continue;
    AnaUtil::fillHist1D("eleCutFlow", 4, wt);

    // we don't use the missingHits cut anymore since this variable is now part of the BDT inputs. 
    //if (electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits"))             continue;
    AnaUtil::fillHist1D("eleCutFlow", 5, wt);

    // Cross cleaning
    if (!crossCleaned(electron))                                                                 continue;
    AnaUtil::fillHist1D("eleCutFlow", 6, wt);

    preSIPlooseEleVec_.push_back(electron);
    
    if (std::fabs(electron.dB3D/electron.edB3D) >= AnaUtil::cutValue(electronCutMap(), "SIP3D")) continue;
    AnaUtil::fillHist1D("eleCutFlow", 7, wt);

    // loose electrons
    looseEleVec_.push_back(electron);
    
    // tight electrons
    if (HZZ4lUtil::electronBDT(electron)) {
      AnaUtil::fillHist1D("eleCutFlow", 8, wt);
      tightEleVec_.push_back(electron);
    }
  }
  if (preSIPlooseEleVec_.size() > 1)
    std::sort(preSIPlooseEleVec_.begin(), preSIPlooseEleVec_.end(), PtComparator<vhtm::Electron>());
  if (looseEleVec_.size() > 1)
    std::sort(looseEleVec_.begin(), looseEleVec_.end(), PtComparator<vhtm::Electron>());

  // attach an empty vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate>>> for each loose electron
  for (auto const& electron: looseEleVec_) {
    std::vector<vhtm::PackedPFCandidate> phov;
    looseElePhotonPairVec_.push_back({electron, phov});
  }

  if (tightEleVec_.size() > 1) 
    std::sort(tightEleVec_.begin(), tightEleVec_.end(), PtComparator<vhtm::Electron>());

  searchedEle_ = true;
}
// photon from PFCandidates
void PhysicsObjSelector::photonSelector() {
  if (!searchedMu_)  muonSelector();
  if (!searchedEle_) electronSelector();

  // Pre-selection of PF Photons
  if (looseEleVec_.size() || looseMuVec_.size()) {
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4
      if (pfcand.pdgId != 22 || pfcand.pt <= 2. || std::fabs(pfcand.eta) >= 2.4) continue;

      // Photon PF relative isolation less than 1.8.
      // The PF isolation is computed using a cone of 0.3, a threshold of 0.2 GeV on charged 
      // hadrons with a veto cone of 0.0001, and 0.5 GeV on neutral hadrons and photons with 
      // a veto cone of 0.01, including also the contribution from PU vertices (same radius 
      // and threshold as per charged isolation) 
      double iso = HZZ4lUtil::pfiso(pfcand);
      if (iso/pfcand.pt >= 1.8)                                                  continue;
      //new in 2016
      if (!passedSuperClusterVetobyReference(pfcand, false))                                continue;
      
      // Photons are associated to the closest lepton in the event among all those passing loose ID + SIP cut
      // It is constructed in such a way that we must check elindx first 
      // Discard photons that do not satisfy the cuts dR(pho,l)/ET_pho^2 < 0.012, and dR(pho,l) < 0.5. 
      int muindx = -1, elindx = -1;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (muindx < 0 && elindx < 0) continue;

      TLorentzVector pfcandP4 = HZZ4lUtil::getP4(pfcand); 
      double dRovEt2 = dRmin/pfcandP4.Et2();
      if (dRovEt2 >= 0.012 || dRmin >= 0.5)                                       continue; 

      // If more than one photon is associated to the same lepton, the lowest dR(pho,l)/ET_pho^2 is selected. 
      if (elindx > -1) {  // check electron index first
        if (looseElePhotonPairVec_.at(elindx).second.empty())
	  looseElePhotonPairVec_.at(elindx).second.push_back(pfcand);
        else {
          TLorentzVector prephoP4 = HZZ4lUtil::getP4(looseElePhotonPairVec_.at(elindx).second.at(0));
          TLorentzVector eleP4 = HZZ4lUtil::getP4(looseElePhotonPairVec_.at(elindx).first);
          if (eleP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
            looseElePhotonPairVec_.at(elindx).second.clear();
            looseElePhotonPairVec_.at(elindx).second.push_back(pfcand);
          }
        }
      }
      else if (muindx > -1) {
        if (looseMuPhotonPairVec_.at(muindx).second.empty())
	  looseMuPhotonPairVec_.at(muindx).second.push_back(pfcand);
        else {
          TLorentzVector prephoP4 = HZZ4lUtil::getP4(looseMuPhotonPairVec_.at(muindx).second.at(0));
          TLorentzVector muP4 = HZZ4lUtil::getP4(looseMuPhotonPairVec_.at(muindx).first);
          if (muP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
            looseMuPhotonPairVec_.at(muindx).second.clear();
            looseMuPhotonPairVec_.at(muindx).second.push_back(pfcand);
          }
        }
      }
    }
  }

  // now fill the vectors with tight leptons
  for (const auto& elem: looseElePhotonPairVec_) {
    if (!elem.second.empty())
      fsrPhotonVec_.push_back(elem.second.at(0));

    if (HZZ4lUtil::electronBDT(elem.first)) 
      tightElePhotonPairVec_.push_back(elem);
  }
  for (const auto& elem: looseMuPhotonPairVec_) {
    if (!elem.second.empty()) 
      fsrPhotonVec_.push_back(elem.second.at(0));

    if (elem.first.isPFMuon) 
      tightMuPhotonPairVec_.push_back(elem);
  }
}
void PhysicsObjSelector::leptonCrossCleaning() {
  vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate> > > list;
  for (const auto& elem: tightElePhotonPairVec_) 
    if (crossCleaned(elem.first)) list.push_back(elem);

  tightElePhotonPairVec_.clear();
  tightElePhotonPairVec_ = list;
}
// from twiki
// Lepton cross cleaning: Remove electrons which are 
// within dR(eta,phi)<0.05 of a muon passing tight ID && SIP<4 
bool PhysicsObjSelector::crossCleaned(const vhtm::Electron& electron) const {
  bool flag = true;
  for (const auto& mu: tightMuVec_) {
    if (HZZ4lUtil::getP4(electron).DeltaR(HZZ4lUtil::getP4(mu)) < 0.05) {
      flag = false;
      break;
    }
  }
  return flag;
}
void PhysicsObjSelector::jetSelector(double wt) {
  for (const auto& jet: *jetColl()) {
    AnaUtil::fillHist1D("jetCutFlow", 0, wt);
    
    if (jet.pt <= AnaUtil::cutValue(jetCutMap(), "pt"))                   continue;
    AnaUtil::fillHist1D("jetCutFlow", 1, wt);

    if (std::fabs(jet.eta) >= AnaUtil::cutValue(jetCutMap(), "eta"))      continue;
    AnaUtil::fillHist1D("jetCutFlow", 2, wt);

    // from twiki
    // PU jet ID is currently buggy; a fix will arrive soon. In the meanwhile, it should not be applied.
    //if (!HZZ4lUtil::jetpuMVAid(jet))                                      continue; 
    AnaUtil::fillHist1D("jetCutFlow", 3, wt);

    if (!jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep"))) continue;
    AnaUtil::fillHist1D("jetCutFlow", 4, wt);

    if (!HZZ4lUtil::isLooseJet(jet))                                      continue;   
    AnaUtil::fillHist1D("jetCutFlow", 5, wt);
    looseJetVec_.push_back(jet);

    AnaUtil::fillHist1D("btagfactor", jet.pfCombinedInclusiveSecondaryVertexV2BJetTags, wt);
    if (jet.pfCombinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor"))
      nbJets_++;
    AnaUtil::fillHist1D("jetPt", jet.pt, wt);
    if (HZZ4lUtil::isTightJet(jet)) {
      tightJetVec_.push_back(jet);
      AnaUtil::fillHist1D("jetCutFlow", 6, wt);
    }
  }
    
  int nJets = looseJetVec_.size();
  AnaUtil::fillHist1D("nJets", nJets, wt);
  AnaUtil::fillHist1D("nbJets", nbJets_, wt);

  if (nJets > 1)
    std::sort(looseJetVec_.begin(), looseJetVec_.end(), PtComparator<vhtm::Jet>());
  if (tightJetVec_.size() > 1)
    std::sort(tightJetVec_.begin(), tightJetVec_.end(), PtComparator<vhtm::Jet>());
}
bool PhysicsObjSelector::jetLeptonCleaning(const vhtm::Jet& jet, double dR) const {
  TLorentzVector jetP4 = HZZ4lUtil::getP4(jet);
  // Muons
  for (const auto& mu: tightMuVec_) {
    double reliso = HZZ4lUtil::computeMuonReliso(mu, fsrPhotonVec_, 0.01,0.3);
    //76X sync, cone 0.3, riso<0.35
    if (reliso >= 0.35) continue;
    if (jetP4.DeltaR(HZZ4lUtil::getP4(mu)) <= dR) return false;
  }
  // Electrons
  //for 76X sync 03 cone is used
  for (const auto& ele: tightEleVec_) {
    //for 04 cone
    //double reliso = HZZ4lUtil::computeElectronReliso(ele, fsrPhotonVec_, getEventGridRho());
    //for cone 03
    double reliso = HZZ4lUtil::computeElectronReliso(ele, fsrPhotonVec_, getEventGridRho(),0.08, 0.3);
    //for cone size 0.4
    //if (reliso >= 0.5) continue;
    //for cone size 03, cut < 0.35
    if (reliso >= 0.35) continue;
    if (jetP4.DeltaR(HZZ4lUtil::getP4(ele)) <= dR) return false;
  }
  // clean against FSR also
  for (const auto& v: fsrPhotonVec_) {
    // how about photon isolation?
    if (jetP4.DeltaR(HZZ4lUtil::getP4(v)) <= dR) return false;
  }
  return true;
}
// clear lists
void PhysicsObjSelector::clear() {
  preSIPlooseMuVec_.clear();
  looseMuVec_.clear();
  tightMuVec_.clear();
  tightIsoMuVec_.clear();

  preSIPlooseEleVec_.clear();
  looseEleVec_.clear();
  tightEleVec_.clear();
  tightIsoEleVec_.clear();

  looseJetVec_.clear();
  tightJetVec_.clear();
  nbJets_ = 0;

  looseElePhotonPairVec_.clear();
  tightElePhotonPairVec_.clear();
  tightIsoElePhotonPairVec_.clear();

  looseMuPhotonPairVec_.clear();
  tightMuPhotonPairVec_.clear();
  tightIsoMuPhotonPairVec_.clear();

  fsrPhotonVec_.clear();

  searchedMu_ = false; 
  searchedEle_ = false; 
  searchedPhoton_ = false;
}

bool PhysicsObjSelector::passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose) const {
  // Supercluster veto: remove all PF photons that match with any electron passing loose ID and SIP cuts; 
  // matching is according to (|deta| < 2, |dphi| < 0.05) OR (dR < 0.15), with respect to the electron's supercluster. 
  bool passedVeto = true;
  if (verbose && looseEleVec_.size())
    cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << endl;
  for (const auto& ele: looseEleVec_) {
    double deta = ele.scEta - pfcand.eta;
    double dphi = TVector2::Phi_mpi_pi(ele.scPhi - pfcand.phi);
    double dR = std::sqrt(deta * deta + dphi * dphi);
    if (verbose) {
      cout << setprecision(3);
      cout << setw(8) << pfcand.pt
	   << setw(8) << pfcand.eta
	   << setw(8) << pfcand.phi
	   << setw(8) << ele.pt
	   << setw(8) << ele.scEta
	   << setw(8) << ele.scPhi
	   << setw(8) << deta
	   << setw(8) << dphi
	   << setw(8) << dR
	   << endl;
    }
    if ( (std::fabs(deta) < 0.05 && std::fabs(dphi) < 2.0) || dR < 0.15 ) {
      passedVeto = false;
      break;
    }
  }
  return passedVeto;
}
bool PhysicsObjSelector::passedSuperClusterVetobyReference(const vhtm::PackedPFCandidate& pfcand, bool verbose) const {
  // Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, as returned by the method 
  bool passedVeto = true;
  TLorentzVector pfcandP4 = HZZ4lUtil::getP4(pfcand);
  if (verbose && looseEleVec_.size())
    cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << endl;
  for (const auto& ele: looseEleVec_) {
    for(const auto& prefP4 : ele.associatedPackedPFCandidatesP4 ) {
      if(pfcandP4 == prefP4) {
        passedVeto = false;
        break;
      }
    }
    if(!passedVeto)   break;
    if (verbose) {
      cout << setprecision(3);
      cout << setw(8) << pfcand.pt
	   << setw(8) << pfcand.eta
	   << setw(8) << pfcand.phi
	   << setw(8) << ele.pt
	   << setw(8) << ele.scEta
	   << setw(8) << ele.scPhi
	   << endl;
    }
  }
  return passedVeto;
}
double PhysicsObjSelector::findClosestLepton(const vhtm::PackedPFCandidate& pfPho, int& muindx, int& elindx) const {
  TLorentzVector phoP4 = HZZ4lUtil::getP4(pfPho);
  double dRmin = 999;
  muindx = -1;
  // First consider loose+SIP muons
  for (unsigned int i = 0; i < looseMuVec_.size(); ++i) {
    const vhtm::Muon& mu = looseMuVec_[i];
    TLorentzVector muP4 = HZZ4lUtil::getP4(mu);
    double dR = muP4.DeltaR(phoP4);
    //double dRovEt2 = dR/phoP4.Et2();
    //if (dRovEt2 >= 0.012 || dR >= 0.5) continue; 
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  // Then consider loose+SIP electrons
  elindx = -1;
  for (unsigned int i = 0; i < looseEleVec_.size(); ++i) {
    const vhtm::Electron& ele = looseEleVec_[i];
    TLorentzVector eleP4 = HZZ4lUtil::getP4(ele);
    double dR = eleP4.DeltaR(phoP4);
    //double dRovEt2 = dR/phoP4.Et2();
    //if (dRovEt2 >= 0.012 || dR >= 0.5) continue; 
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  return dRmin;
}
void PhysicsObjSelector::findObjects(double wt) {
  // order of execution is crucial!
  const vhtm::Event& evt = eventColl()->at(0);
  fGridRhoFastjetAll_ = evt.fGridRhoFastjetAll;

  // muonSelector must precede electronSelector
  muonSelector(wt);
  electronSelector(wt);
  
  // After electrons and muons are found, find photons
  photonSelector();

  // Now find isolated leptons after removing contribution from the associated photons
  isoLeptonSelector();

  // Jets
  jetSelector(wt);
}
void PhysicsObjSelector::isoLeptonSelector() {
  for (const auto& elem: tightMuPhotonPairVec_) {
    double iso = HZZ4lUtil::computeMuonReliso(elem.first, fsrPhotonVec_, 0.01, 0.3);
    if (iso >= 0.35) continue; // it is already scaled by pT
    tightIsoMuVec_.push_back(elem.first);
    tightIsoMuPhotonPairVec_.push_back(elem);
  }
  for (const auto& elem: tightElePhotonPairVec_) {
    //double iso = HZZ4lUtil::computeElectronReliso(elem.first, fsrPhotonVec_, getEventGridRho());
    double iso = HZZ4lUtil::computeElectronReliso(elem.first, fsrPhotonVec_, getEventGridRho(),0.08,0.3);
    //if (iso >= 0.5) continue; // it is already scaled by pT
    if (iso >= 0.31) continue; // it is already scaled by pT
    tightIsoEleVec_.push_back(elem.first);
    tightIsoElePhotonPairVec_.push_back(elem);
  }
}
void PhysicsObjSelector::ZZMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec) {
  // -- from twiki --
  // Require both Z candidate masses (computed including the FSR photons, if present) to be 12 < m(ll(g)) < 120 GeV
  // leptons are already isolated in the new FSR scheme
  // Z pair pass the mass cut, we have a ZZ Candidate
  bool ZZmasscond = (Za.mass > 12 && Za.mass < 120) && (Zb.mass > 12 && Zb.mass < 120);
  
  // -- from twiki --
  // define the Z1 as the one with mass closest to the nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2.
  if (ZZmasscond) {
    if (Za.massDiff < Zb.massDiff) {
      if (Za.mass > 40.)
  	ZZVec.push_back({Za, Zb});
    } 
    else {
      if (Zb.mass > 40.)
	ZZVec.push_back({Zb, Za});
    }
  }
}
void PhysicsObjSelector::addLeptonIsolation(std::vector<ZCandidate>& candList, 
					    const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& elePhotonPairVec, 
					    const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& muPhotonPairVec) 
{
  for (auto& v: candList) {
    if (v.flavour == HZZ4lUtil::ZType::mumu) {
      const vhtm::Muon& mu1 = muPhotonPairVec.at(v.l1Index).first;
      v.l1Isolation = HZZ4lUtil::computeMuonReliso(mu1, fsrPhotonVec_);

      const vhtm::Muon& mu2 = muPhotonPairVec.at(v.l2Index).first;
      v.l2Isolation = HZZ4lUtil::computeMuonReliso(mu2, fsrPhotonVec_);
    }
    else if (v.flavour == HZZ4lUtil::ZType::ee) {
      const vhtm::Electron& ele1 = elePhotonPairVec.at(v.l1Index).first;
      v.l1Isolation = HZZ4lUtil::computeElectronReliso(ele1, fsrPhotonVec_, getEventGridRho());

      const vhtm::Electron& ele2 = elePhotonPairVec.at(v.l2Index).first;
      v.l2Isolation = HZZ4lUtil::computeElectronReliso(ele2, fsrPhotonVec_, getEventGridRho());
    }
  }
}
void PhysicsObjSelector::dumpEvent(bool dumpGen, bool showEvent, ostream& os) const {
  os << std::setprecision(3);

  // Event
  if (showEvent) showEventNumber(os);

  // Muons
  if (muonColl()->size()) {
    os << " -- # Muons: " << muonColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  global tracker      PF  nMatch  Type       SIP ghostCleaned  reliso"
       << endl;
    int indx = 0;
    for (auto const& muon: *muonColl()) {
      os << setw(6) << ++indx 
	 << setw(8) << muon.pt 
	 << setw(8) << muon.eta
	 << setw(8) << muon.phi
         << setw(8) << muon.charge
	 << setw(8) << muon.dxyPV
	 << setw(8) << muon.dzPV
	 << setw(8) << muon.isGlobalMuon
	 << setw(8) << muon.isTrackerMuon
	 << setw(8) << muon.isPFMuon
	 << setw(8) << muon.matches
	 << setw(6) << muon.muonBestTrackType
	 << setw(10) << muon.dB3D/muon.edB3D
	 << setw(13) << (muon.isghostCleaned ? "yes" : "no")
         << setw(8) << HZZ4lUtil::computeMuonReliso(muon, fsrPhotonVec_, 0.01, 0.4, false)
	 << endl;
    }
  }
  // Electrons
  if (electronColl()->size()) {
    os << " -- # Electrons: " << electronColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  misHit       SIP   isBDT  reliso crossCleaned      BDT BDTpreComp"
       << endl; 
    int indx = 0;
    for (auto const& electron: *electronColl()) {
      os << setw(6)  << ++indx
	 << setw(8)  << electron.pt 
	 << setw(8)  << electron.eta
	 << setw(8)  << electron.phi
         << setw(8)  << electron.charge
	 << setw(8)  << electron.dxyPV
	 << setw(8)  << electron.dzPV
	 << setw(8)  << electron.missingHits
	 << setw(10) << electron.dB3D/electron.edB3D
	 << setw(8)  << HZZ4lUtil::electronBDT(electron)
         << setw(8)  << HZZ4lUtil::computeElectronReliso(electron, fsrPhotonVec_, getEventGridRho(), 0.08, 0.4, false)
         << setw(13) << crossCleaned(electron)
         << setw(9) << electron.BDT
         << setw(11) << electron.BDTpreComp
	 << endl;
    }
  }
  // Photons
  int npho = 0;
  for (const auto& pfcand: *packedPFCandidateColl()) {
    // Preselection: pT > 2 GeV, |eta| < 2.4 
    if (pfcand.pdgId != 22 || pfcand.pt <= 2. || std::fabs(pfcand.eta) >= 2.4) continue;
    ++npho;
  }
  if (npho) {
    os << " -- Photons: " << endl;
    os << "  indx      pT     eta     phi  reliso  passedSCVeto closestLepton   dRmin        Et2 dRovEt2(%)"
       << endl; 
    int indx = 0;
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt <= 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      int muindx = -1, elindx = -1;
      ostringstream ltype;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (elindx > -1) ltype << "electron(" << elindx << ")";
      else if (muindx > -1) ltype << "muon(" << muindx << ")";
      else ltype << "none";
      double iso = HZZ4lUtil::pfiso(pfcand);
      if (0) os << "Photon iso: "
	        << setw(9) << pfcand.isolationMap.at("c30").at(0)
		<< setw(9) << pfcand.isolationMap.at("c30").at(2)
		<< setw(9) << pfcand.isolationMap.at("c30").at(3)
		<< setw(9) << pfcand.isolationMap.at("c30").at(4)
		<< endl;
      os << setw(6) << ++indx
	 << setw(8) << pfcand.pt 
	 << setw(8) << pfcand.eta
	 << setw(8) << pfcand.phi
	 << setw(8) << iso/pfcand.pt
	 << setw(14) << (passedSuperClusterVetobyReference(pfcand) ? "yes" : "no")
	 << setw(14) << ltype.str()
	 << setw(8) << dRmin
         << setprecision(2)
	 << setw(11) << HZZ4lUtil::getP4(pfcand).Et2()
         << setprecision(3)
	 << setw(11) << dRmin*100/HZZ4lUtil::getP4(pfcand).Et2()
	 << endl;
    }
  }
  // Jets
  if (jetColl()->size()) {
    os << " -- # Jets: " << jetColl()->size() << endl;
    os << "  indx       pT      eta      phi NumConst   CHM      CHF     CEMF      NHF     NEMF      MUF     puId   bDiscrim lepCleaned  looseId  tightId"
       << endl;
    int indx = 0;
    for (auto const& jet: *jetColl()) {
      os << setw(6) << ++indx
	 << setw(9) << jet.pt
	 << setw(9) << jet.eta
	 << setw(9) << jet.phi
         << setw(9) << (jet.chargedMultiplicity + jet.neutralMultiplicity)
         << setw(6) << jet.chargedMultiplicity
         << setw(9) << jet.chargedHadronEnergyFraction
         << setw(9) << jet.chargedEmEnergyFraction
         << setw(9) << jet.neutralHadronEnergyFraction
         << setw(9) << jet.neutralEmEnergyFraction
         << setw(9) << jet.muonEnergyFraction
	 << setw(9) << HZZ4lUtil::jetpuMVAid(jet)
         << setw(11) << jet.pfCombinedInclusiveSecondaryVertexV2BJetTags 
	 << setw(11) << (jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep")) ? "yes" : "no")
         << setw(9) << HZZ4lUtil::isLooseJet(jet)
         << setw(9) << HZZ4lUtil::isTightJet(jet)
	 << endl;
    }
  }
  // Selected vertices
  if (vertexColl()->size()) {
    os << " -- # Vertices: " << vertexColl()->size() << endl;
    os << "  indx      ndf     rho    chi2     dxy       z  isfake" 
       << endl; 
    int indx = 0;
    for (auto const& vtx: *vertexColl()) {
      double dxy = std::sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
      os << setw(6) << ++indx
	 << setw(9) << vtx.ndf
	 << setw(8) << vtx.rho
	 << setw(8) << vtx.chi2
	 << setw(8) << dxy
	 << setw(8) << vtx.z
	 << setw(8) << (vtx.isfake ? "yes" : "no")
	 << endl;
    }
  }
  if (isMC() && dumpGen) dumpGenInfo(os); 
}
