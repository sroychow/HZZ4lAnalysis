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

    if (std::fabs(muon.dB3D/muon.edB3D) >= AnaUtil::cutValue(muonCutMap(), "SIP3D"))  continue;
    AnaUtil::fillHist1D("muCutFlow", 7, wt);

    if (!muon.isghostCleaned)                                                         continue;
    AnaUtil::fillHist1D("muCutFlow", 8, wt);

    // loose muon 
    looseMuVec_.push_back(muon);

    // attach an empty PackedPFCandidate vector to each loose Muon
    std::vector<vhtm::PackedPFCandidate> phov;
    looseMuPhotonPairVec_.push_back({muon, phov});

    // tight muon 
    if (muon.isPFMuon) {
      AnaUtil::fillHist1D("muCutFlow", 9, wt);
      tightMuVec_.push_back(muon);
    }
  }
  if (looseMuVec_.size() > 1)
    std::sort(looseMuVec_.begin(), looseMuVec_.end(), PtComparator<vhtm::Muon>());
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

    if (electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 5, wt);

    if (std::fabs(electron.dB3D/electron.edB3D) >= AnaUtil::cutValue(electronCutMap(), "SIP3D")) continue;
    AnaUtil::fillHist1D("eleCutFlow", 6, wt);
    
    // Cross cleaning
    if (!crossCleaned(electron))  continue;
    AnaUtil::fillHist1D("eleCutFlow", 7, wt);

    // loose electrons
    looseEleVec_.push_back(electron);
    
    // attach an empty PackedPFCandidate vector to each loose electron
    std::vector<vhtm::PackedPFCandidate> phov;
    looseElePhotonPairVec_.push_back({electron, phov});

    // tight electrons
    if (HZZ4lUtil::electronBDT(electron)) {
      AnaUtil::fillHist1D("eleCutFlow", 8, wt);
      tightEleVec_.push_back(electron);
    }
  }
  if (looseEleVec_.size() > 1)
    std::sort(looseEleVec_.begin(), looseEleVec_.end(), PtComparator<vhtm::Electron>());
  if (tightEleVec_.size() > 1) 
    std::sort(tightEleVec_.begin(), tightEleVec_.end(), PtComparator<vhtm::Electron>());

  searchedEle_ = true;
}
// photon from PFCandidates
void PhysicsObjSelector::photonSelector() {
  if (!searchedMu_)  muonSelector();
  if (!searchedEle_) electronSelector();

  // -- from twiki --
  // FSR Photon selection: for each photon, consider the closest lepton:
  // Note: Multiple photons can be attached to a lepton
  // Muons

  // Pre-selection of PF Photons
  if (looseEleVec_.size() || looseMuVec_.size()) {
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt <= 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      
      // Supercluster veto
      if (!passedSuperClusterVeto(pfcand, false)) continue;
      
      int muindx = -1, elindx = -1;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      bool leptonMatched = false;
      if (dRmin < 0.07) {
	if (pfcand.pt > 2.) leptonMatched = true;
      }
      else if (dRmin < 0.5) {
	if (pfcand.pt > 4.) {
	  double iso = pfcand.isolationMap.at("c30").at(0)
	    + pfcand.isolationMap.at("c30").at(2) 
	    + pfcand.isolationMap.at("c30").at(3) 
	    + pfcand.isolationMap.at("c30").at(4);
	  if (iso/pfcand.pt < 1.) leptonMatched = true;
	}
      } 
      if (!leptonMatched) continue;
      photonVec_.push_back(pfcand);
      
      // Now update vector of pair
      if (elindx > -1)
	looseElePhotonPairVec_.at(elindx).second.push_back(pfcand);
      else if (muindx > -1)
	looseMuPhotonPairVec_.at(muindx).second.push_back(pfcand);
    }
  }
  // now fill the vectors with tight leptons
  for (const auto& el: looseElePhotonPairVec_)
    if (HZZ4lUtil::electronBDT(el.first)) tightElePhotonPairVec_.push_back(el);

  for (const auto& el: looseMuPhotonPairVec_)
    if (el.first.isPFMuon) tightMuPhotonPairVec_.push_back(el);
}
void PhysicsObjSelector::leptonCrossCleaning() {
  vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate> > > list;
  for (const auto& it: tightElePhotonPairVec_) 
    if (crossCleaned(it.first)) list.push_back(it);

  tightElePhotonPairVec_.clear();
  tightElePhotonPairVec_ = list;
}
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

    if (!HZZ4lUtil::jetpuMVAid(jet))                                       continue; 
    AnaUtil::fillHist1D("jetCutFlow", 3, wt);

    if (!jetLeptonCleaning(jet, AnaUtil::cutValue(jetCutMap(), "dRlep"))) continue;
    AnaUtil::fillHist1D("jetCutFlow", 4, wt);

    if (!HZZ4lUtil::isLooseJet(jet))                                       continue;   
    AnaUtil::fillHist1D("jetCutFlow", 5, wt);
    looseJetVec_.push_back(jet);

    AnaUtil::fillHist1D("btagfactor", jet.combinedInclusiveSecondaryVertexV2BJetTags, wt);
    if (jet.combinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor"))
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
  for (const auto& mu: tightMuVec_) {
    if (HZZ4lUtil::pfiso(mu)/mu.pt >= 0.4) continue;
    if (jetP4.DeltaR(HZZ4lUtil::getP4(mu)) <= dR) return false;
  }
  for (const auto& ele: tightEleVec_) {
    if (HZZ4lUtil::pfiso(ele, getEventGridRho())/ele.pt >= 0.5) continue;
    if (jetP4.DeltaR(HZZ4lUtil::getP4(ele)) <= dR) return false;
  }
  return true;
}
// clear lists
void PhysicsObjSelector::clear() {
  looseMuVec_.clear();
  tightMuVec_.clear();
  looseEleVec_.clear();
  tightEleVec_.clear();
  photonVec_.clear();
  looseJetVec_.clear();
  tightJetVec_.clear();
  nbJets_ = 0;

  looseElePhotonPairVec_.clear();
  tightElePhotonPairVec_.clear();
  looseMuPhotonPairVec_.clear();
  tightMuPhotonPairVec_.clear();

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
    float dR = std::sqrt(deta * deta + dphi * dphi);
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
double PhysicsObjSelector::findClosestLepton(const vhtm::PackedPFCandidate& pfPho, int& muindx, int& elindx) const {
  TLorentzVector phoP4 = HZZ4lUtil::getP4(pfPho);
  double dRmin = 999;
  muindx = -1;
  // First consider loose muons
  for (unsigned int i = 0; i < looseMuVec_.size(); ++i) {
    const vhtm::Muon& mu = looseMuVec_[i];
    TLorentzVector muP4 = HZZ4lUtil::getP4(mu);
    double dR = muP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  // Then consider loose electron
  elindx = -1;
  for (unsigned int i = 0; i < looseEleVec_.size(); ++i) {
    const vhtm::Electron& ele = looseEleVec_[i];
    TLorentzVector eleP4 = HZZ4lUtil::getP4(ele);
    double dR = eleP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  return dRmin;
}
void PhysicsObjSelector::findObjects(double wt) {
  const vhtm::Event& evt = eventColl()->at(0);
  fGridRhoFastjetAll_ = evt.fGridRhoFastjetAll;

  // muonSelector must precede electronSelector
  muonSelector(wt);
  electronSelector(wt);
  
  // after electrons and muons are found, find photons
  photonSelector();

  // Jets
  jetSelector(wt);
}
bool PhysicsObjSelector::ZeeIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone, double isoCone, bool verbose) const {
  const auto& elePhotonPairVec = getTightElePhotonPairList();
  
  const vhtm::Electron& ele1 = elePhotonPairVec.at(Za.l1Index).first;
  double l1iso = HZZ4lUtil::computeElectronIso(ele1, Za, Za.l1P4, fsrVec, getEventGridRho(), vetoCone, isoCone, verbose);
  Za.l1Isolation = l1iso;

  const vhtm::Electron& ele2 = elePhotonPairVec.at(Za.l2Index).first;
  double l2iso = HZZ4lUtil::computeElectronIso(ele2, Za, Za.l2P4, fsrVec, getEventGridRho(), vetoCone, isoCone, verbose);
  Za.l2Isolation = l2iso;
  
  if (Za.l1Isolation < 0.5 && Za.l2Isolation < 0.5) return true;
  return false;
}
bool PhysicsObjSelector::ZmumuIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone, double isoCone, bool verbose) const {
  const auto& muPhotonPairVec = getTightMuPhotonPairList();

  const vhtm::Muon& mu1 = muPhotonPairVec.at(Za.l1Index).first;
  double l1iso = HZZ4lUtil::computeMuonIso(mu1, Za, Za.l1P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l1Isolation = l1iso;
  
  const vhtm::Muon& mu2 = muPhotonPairVec.at(Za.l2Index).first;
  double l2iso = HZZ4lUtil::computeMuonIso(mu2, Za, Za.l2P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l2Isolation = l2iso;
  
  if (Za.l1Isolation < 0.4 && Za.l2Isolation < 0.4) return true;
  return false;
}
// Function to compute the isolation for the ZZ leptons and check the Z masses
// If they pass, form a candidate ZZ pair
void PhysicsObjSelector::ZZisoMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec) {
  std::vector<TLorentzVector> fsrVec;
  fsrVec.push_back(Za.fsrPhoP4);
  fsrVec.push_back(Zb.fsrPhoP4);
  // -- from twiki --
  // Require all four leptons to pass the isolation requirement
  
  // First Z Candidate
  bool ZalepIso = false; 
  if (Za.flavour == HZZ4lUtil::ZType::mumu)
    ZalepIso = ZmumuIso(Za, fsrVec);
  else if (Za.flavour == HZZ4lUtil::ZType::ee)
    ZalepIso = ZeeIso(Za, fsrVec);
  
  // Second Z candidate 
  bool ZblepIso = false;
  if (Zb.flavour == HZZ4lUtil::ZType::mumu)
    ZblepIso = ZmumuIso(Zb, fsrVec);
  else if (Zb.flavour == HZZ4lUtil::ZType::ee)
    ZblepIso = ZeeIso(Zb, fsrVec);
  
  // -- from twiki --
  // Require both Z candidate masses (computed including the FSR photons, if present) to be 12 < m(ll(g)) < 120 GeV
  bool ZZmasscond = (Za.mass > 12 && Za.mass < 120) && (Zb.mass > 12 && Zb.mass < 120);
  
  // if all the leptons pass isolation and Z pair pass the mass cut, we have a ZZ Candidate
  if (ZalepIso && ZblepIso && ZZmasscond) {
    // -- from twiki --
    // define the Z1 as the one with mass closest to the nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2.
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
         << setw(8) << HZZ4lUtil::pfiso(muon)/muon.pt
	 << endl;
    }
  }
  // Electrons
  if (electronColl()->size()) {
    os << " -- # Electrons: " << electronColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  misHit       SIP   isBDT  reliso crossCleaned"
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
         << setw(8)  << HZZ4lUtil::pfiso(electron, getEventGridRho())/electron.pt
         << setw(13)  << crossCleaned(electron)
	 << endl;
    }
  }
  // Photons
  int npho = 0;
  for (const auto& pfcand: *packedPFCandidateColl()) {
    // Preselection: pT > 2 GeV, |eta| < 2.4 
    if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
    ++npho;
  }
  if (npho) {
    os << " -- Photons: " << endl;
    os << "  indx      pT     eta     phi  reliso  passedSCVeto closestLepton   dRmin"
       << endl; 
    int indx = 0;
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      int muindx = -1, elindx = -1;
      string ltype = "none";
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (elindx > -1) ltype = "electron";
      else if (muindx > -1) ltype = "muon";
      double iso = pfcand.isolationMap.at("c30").at(0) 
	+ pfcand.isolationMap.at("c30").at(2) 
	+ pfcand.isolationMap.at("c30").at(3) 
	+ pfcand.isolationMap.at("c30").at(4);
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
	 << setw(14) << (passedSuperClusterVeto(pfcand) ? "yes" : "no")
	 << setw(14) << ltype
	 << setw(8) << dRmin
	 << endl;
    }
  }
  // Jets
  if (jetColl()->size()) {
    os << " -- # Jets: " << jetColl()->size() << endl;
    os << "  indx       pT      eta      phi NumConst   CHM      CHF     CEMF      NHF     NEMF      MUF     puId bDiscrim lepCleaned  looseId  tightId"
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
	 << setw(9) << jet.combinedInclusiveSecondaryVertexV2BJetTags
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
