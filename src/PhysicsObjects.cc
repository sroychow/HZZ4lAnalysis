#include "PhysicsObjects.h"

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

ClassImp(vhtm::Event)
ClassImp(vhtm::GenEvent)
ClassImp(vhtm::Electron)
ClassImp(vhtm::GenParticle)
ClassImp(vhtm::GenJet)
ClassImp(vhtm::GenMET)
ClassImp(vhtm::MET)
ClassImp(vhtm::Tau)
ClassImp(vhtm::Muon)
ClassImp(vhtm::Jet)
ClassImp(vhtm::Vertex)
ClassImp(vhtm::TriggerObject)
ClassImp(vhtm::Candidate)
ClassImp(vhtm::Photon)
ClassImp(vhtm::PackedPFCandidate)
ClassImp(vhtm::ZtnP)
ClassImp(vhtm::SelectedEvent)

vhtm::Candidate::Candidate():
  pt(-999), eta(-999), phi(-999) {} vhtm::Candidate::Candidate(float _pt, float _eta, float _phi):
  pt(_pt), eta(_eta), phi(_phi) {} vhtm::PackedPFCandidate::PackedPFCandidate():
  pt(-999.),
  eta(-999.), 
  phi(-999.),
  energy(-999.),
  pdgId(0),
  charge(-999),
  vx(-999.),
  vy(-999.),
  vz(-999.),
  fromPV(-999),
  dxy(-999.),
  dz(-999.),
  dxyError(-999.),
  dzError(-999.)
{}

vhtm::Event::Event():
  run(0),
  event(0),
  lumis(0),
  bunch(0),
  orbit(0),
  time(-1),
  isdata(false),
  isPhysDeclared(false),
  isBPTX0(false),
  isBSCMinBias(false),
  isBSCBeamHalo(false),
  isPrimaryVertex(false),
  //isBeamScraping(false),
  rho(-1),
  rhoNeutral(-1),
  nvtx(0),
  fGridRhoAll(-1),
  fGridRhoFastjetAll(-1),
  fGridRhoFastjetAllCalo(-1),
  fGridRhoFastjetCentralCalo(-1),
  fGridRhoFastjetCentralChargedPileUp(-1),
  fGridRhoFastjetCentralNeutral(-1)
{
  nPU.clear();
  bunchCrossing.clear();
  trueNInt.clear();
}

vhtm::GenEvent::GenEvent():
  processID(0),
  ptHat(-999)
{
  pdfWeights.clear();
}

vhtm::Electron::Electron():
  eta(-999),
  phi(-999),
  pt(-999),
  ecalDriven(false),
  isGap(false),
  hasGsfTrack(false),
  trackPt(-999),
  energy(-999),
  caloEnergy(-999),
  charge(-999),
  pixHits(-1),
  trkHits(-1),
  nValidHits(-1),
  trkD0(-999),
  trkDz(-999),
  hoe(-999),
  eop(-999),
  sigmaEtaEta(-999),
  sigmaIEtaIEta(-999),
  deltaPhiTrkSC(-999),
  deltaEtaTrkSC(-999),
  classif(-1),
  vx(-999),
  vy(-999),
  vz(-999),
  scEn(-999),
  scEta(-999),
  scPhi(-999),
  scET(-999),
  scRawEnergy(-999),
  BDT(-999),
  BDTpreComp(-999),
  dxyPV(-999),
  dzPV(-999),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  relIso(-999),
  pfRelIso(-999),
  chargedHadronIso(-999),
  neutralHadronIso(-999),
  photonIso(-999),
  sumChargedHadronPt(-999),
  sumNeutralHadronEt(-999.),
  sumPhotonEt(-999.),
  sumPUPt(-999),
  missingHits(-1),
  dB(-999),
  edB(-999),
  dB3D(-999),
  edB3D(-999),
  nBrems(-1),
  fbrem(-999),
  hasMatchedConv(false),
  mvaPreselection(false),
  isTriggerElectron(false),
  fidFlag(0),
  selbit(0),
  passMediumId(false),
  passTightId(false),
  mvaCategory(false),
  ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values(-999.)
{
  idmap.clear();
}

vhtm::GenParticle::GenParticle():
  eta(-999),
  phi(-999),
  p(-999),
  px(-999),
  py(-999),
  pz(-999),
  pt(-999),
  energy(-999),
  pdgId(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  status(-999),
  charge(-999),
  numDaught(-1),
  numMother(-1),
  motherIndex(-1)
{
  motherIndices.clear();
  daughtIndices.clear();
}

vhtm::GenJet::GenJet():
  eta(-999),
  phi(-999),
  p(-999),
  pt(-999),
  energy(-999),
  emf(-999),
  hadf(-999)
{}

vhtm::MET::MET():
  met(-999),
  metphi(-999),
  sumet(-999),
  metuncorr(-999),
  metphiuncorr(-999),
  sumetuncorr(-999)
{}

vhtm::GenMET::GenMET():
  met(-999),
  metphi(-999),
  sumet(-999)
{}

vhtm::Tau::Tau():
  eta(-999),
  phi(-999),
  pt(-999),
  energy(-999),
  charge(-999),
  mass(-999),
  dxyPV(-999),
  dzPV(-999),
  vtxIndex(-1),
  vtxDxy(-999),
  vtxDz(-999),
  leadChargedParticlePt(-999),
  leadNeutralParticlePt(-999),
  leadParticlePt(-999),
  ptSumChargedHadronsIsoCone(-999),
  ptSumNeutralHadronsIsoCone(-999),
  ptSumPhotonsIsoCone(-999),
  decayModeFinding(-1),
  decayModeFindingNewDMs(-1),
  decayModeFindingOldDMs(-1),
  againstMuonLoose(-1),
  againstMuonMedium(-1),
  againstMuonTight(-1),
  againstMuonLoose3(-1),
  againstMuonTight3(-1),
  againstElectronLoose(-1),
  againstElectronMedium(-1),
  againstElectronTight(-1),
  //againstElectronMVA(-1),
  againstElectronLooseMVA5(-1),
  againstElectronMediumMVA5(-1),
  againstElectronTightMVA5(-1),
  byLooseCombinedIsolationDeltaBetaCorr3Hits(-1),
  byMediumCombinedIsolationDeltaBetaCorr3Hits(-1),
  byTightCombinedIsolationDeltaBetaCorr3Hits(-1),
  byCombinedIsolationDeltaBetaCorrRaw3Hits(-1),
  chargedIsoPtSum(-1),
  neutralIsoPtSum(-1),
  puCorrPtSum(-1),
  jetPt(-999),
  jetEta(-999),
  jetPhi(-999),
  emFraction(-999),
  vx(-999), vy(-999), vz(-999),
  zvertex(-999), 
  dxySig(-999),
  selbit(0)
{
  sigChHadList.clear();
  sigNeHadList.clear();
  sigGammaList.clear();
  isoChHadList.clear();
  isoNeHadList.clear();
  isoGammaList.clear();
}
vhtm::Muon::Muon():
  isGlobalMuon(false),
  isTrackerMuon(false),
  isPFMuon(false),
  isghostCleaned(false),
  passTrackerhighPtid(false),
  eta(-999),
  phi(-999),
  pt(-999),
  p(-999),
  energy(-999),
  charge(-999),
  trkD0(-999),
  trkDz(-999),
  muonBestTrackType(-999),
  globalChi2(9999.),
  tkNChi2(9999.),
  trkIso(-999),
  ecalIso(-999),
  hcalIso(-999),
  hoIso(-999),
  pfChargedIsoR03(-999),
  pfChargedHadIsoR03(999.),
  pfNeutralHadIsoR03(999.),
  pfPhotonIso03(999.),
  sumPUPt03(-999),
  pfRelIso03(-1),
  sumChargedParticlePt(999.),
  sumChargedHadronPt(999.),
  sumNeutralHadronEt(999.),
  sumPhotonEt(999.),
  sumPUPt(999.),   
  pfRelIso04(-1),
  passID(false),
  dxyPV(-999),
  dzPV(-999),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  pixHits(-1),
  trkHits(-1),
  muoHits(-1),
  matches(-1),
  trackerLayersWithMeasurement(-1),
  vx(-999),
  vy(-999),
  vz(-999),
  dB(-999),
  edB(-999),
  dB3D(-999),
  edB3D(-999),
  isGlobalMuonPromptTight(false),
  isAllArbitrated(false),
  nChambers(-1),
  nMatches(-1),
  nMatchedStations(-1),
  stationMask(0),
  stationGapMaskDistance(0),
  stationGapMaskPull(0),
  muonID(false),
  nSegments(-1),
  selbit(0)
{}

vhtm::Jet::Jet():
  eta(-999),
  phi(-999),
  pt(-999),
  pt_raw(-999),
  energy(-999),
  energy_raw(-999),
  jecUnc(-999),
  resJEC(-999),
  partonFlavour(-1),
  chargedEmEnergyFraction(-999),
  chargedHadronEnergyFraction(-999),
  chargedMuEnergyFraction(-999),
  electronEnergyFraction(-999),
  muonEnergyFraction(-999),
  neutralEmEnergyFraction(-999),
  neutralHadronEnergyFraction(-999),
  photonEnergyFraction(-999),
  chargedHadronMultiplicity(-1),
  chargedMultiplicity(-1),
  electronMultiplicity(-1),
  muonMultiplicity(-1),
  neutralHadronMultiplicity(-1),
  neutralMultiplicity(-1),
  photonMultiplicity(-1),
  nConstituents(-1),
  //simpleSecondaryVertexHighEffBTag(-999),
  //simpleSecondaryVertexHighPurBTag(-999),
  combinedSecondaryVertexBTag(-999),
  //combinedSecondaryVertexMVABTag(-999),
  combinedInclusiveSecondaryVertexBTag(-999),
  pfCombinedInclusiveSecondaryVertexV2BJetTags(-999.),
  //combinedMVABTag(-999),
  jpumva(9999.),
  passLooseID(-1),
  passTightID(-1),
  selbit(0) 
{
  discrimap.clear();
}

vhtm::Vertex::Vertex():
  x(-999),
  y(-999),
  z(-999),
  xErr(-999),
  yErr(-999),
  zErr(-999),
  rho(-999),
  chi2(999.),
  ndf(-1),
  //ntracks(-1),
  //ntracksw05(-1),
  isfake(true),
  isvalid(false),
  //sumPt(-999),
  selbit(0) {}

vhtm::TriggerObject::TriggerObject():
  energy(-999),
  pt(-999),
  eta(-999),
  phi(-999)
{
  pathList.clear();
}
vhtm::Photon::Photon():
  et(-999),
  eta(-999),
  clusterEta(-999),
  phi(-999),
  clusterPhi(-999),
  energy(-999),
  theta(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  scEnergy(-999),
  scEta(-999),
  scPhi(-999),
  scSize(-999),
  scEtaWidth(-999),
  scPhiWidth(-999),
  scEt(-999),
  scRawEnergy(-999),
  scx(-999),
  scy(-999),
  scz(-999),
  isoEcalRecHit03(-999),
  isoHcalRecHit03(-999),
  isoSolidTrkCone03(-999),
  isoHollowTrkCone03(-999),
  nTrkSolidCone03(-1),
  nTrkHollowCone03(-1),
  isoEcalRecHit04(-999),
  isoHcalRecHit04(-999),
  isoSolidTrkCone04(-999),
  isoHollowTrkCone04(-999),
  nTrkSolidCone04(-1),
  nTrkHollowCone04(-1),
  isEB(false),
  isEE(false),
  isEBGap(false),
  isEEGap(false),
  isEBEEGap(false),
  fidFlag(0),
  passElectronVeto(false),
  hasPixelSeed(false),
  ecalIso(-9999),
  hcalIso(-999),
  trackIso(-999),
  chargedHadIso(-999),
  neutralHadIso(-999),
  photonIso(-999),
  r9(-999),
  hoe(-999),
  sigmaEtaEta(-999),
  sigmaIEtaIEta(-999),
  e1x5(-999),
  e2x5(-999), 
  e3x3(-999),
  e5x5(-999), 
  r1x5(-999), 
  r2x5(-999), 
  maxEnergyXtal(-999),
  hasConversionTracks(false),
  nTracks(-1),
  isConverted(false),
  pairInvMass(-999),
  pairCotThetaSeparation(-999),
  pairPx(-999),
  pairPy(-999),
  pairPz(-999),
  conv_vx(-999),
  conv_vy(-999),
  conv_vz(-999),
  eovp(-999),
  zpv(-999),
  distOfMinApproach(-999),
  dPhiTracksAtVtx(-999),
  dPhiTracksAtEcal(-999),
  dEtaTracksAtEcal(-999),
  selbit(0) {}


vhtm::SelectedEvent::SelectedEvent() : 
      run(0),
      lumi(0),
      event(0),
      mass4l(0.),
      mZ1(0.),
      mZ2(0.),
      nJets(0),
      jet1pt(0.),
      jet2pt(0.),
      category(0),
      m4lRefit(0.),
      m4lRefitError(0.),
      weight(0.),
      flavour(4)	
{
      kd["Dgg10_VAMCFM"] = -1.; 
      kd["D_bkg_kin"] = -1.; 
      kd["D_bkg"] = -1.; 
      kd["D_g4"] = -1.; 
      kd["Djet_VAJHU"] = -1.;
}

void vhtm::SelectedEvent::reset() 
{
      run = 0;
      lumi = 0;
      event = 0;
      mass4l = 0.;
      mZ1 = 0.;
      mZ2 = 0.;
      nJets = 0.;
      jet1pt = 0.;
      jet2pt = 0.;
      category = 0;
      m4lRefit = 0.;
      m4lRefitError = 0.;
      weight = 0.;
      kd["Dgg10_VAMCFM"] = -1.; 
      kd["D_bkg_kin"] = -1.; 
      kd["D_bkg"] = -1.; 
      kd["D_g4"] = -1.; 
      kd["Djet_VAJHU"] = -1.;
      flavour = 4;
}

vhtm::ZtnP::ZtnP() 
{

}


