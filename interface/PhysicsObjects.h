#ifndef __PhysicsObjects_h
#define __PhysicsObjects_h

#include <vector>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include "TObject.h"

namespace vhtm {
  class Candidate: public TObject {
  public:
    Candidate();
    Candidate(float pt, float eta, float phi);
    virtual ~Candidate() {}

    float pt;
    float eta;
    float phi;
   
    ClassDef(Candidate,1)
  };

  class PackedPFCandidate: public TObject {
  public:
    PackedPFCandidate();
    //Candidate(float pt, float eta, float phi);
    virtual ~PackedPFCandidate() {}

    float pt;
    float eta;
    float phi;
    float energy;
    
    int pdgId;
    int charge;
    
    double vx;
    double vy;
    double vz;
   
    int fromPV;
    //w.r.t PV
    float dxy;
    float dz;
    float dxyError;
    float dzError;

    std::map< std::string,std::vector<double> > isolationMap;

    ClassDef(PackedPFCandidate,1)
  };
 
  class Event: public TObject {
  public:
    Event();
    virtual ~Event() {}
  
    unsigned int run;
    unsigned int event;
    unsigned int lumis;
    unsigned int bunch;
    unsigned int orbit;
    double time;
    bool isdata;
  
    bool isPhysDeclared;
    bool isBPTX0;
    bool isBSCMinBias;
    bool isBSCBeamHalo;
    bool isPrimaryVertex;
    //bool isBeamScraping;
    bool passHBHENoiseFilter;
  
    double rho;
    double rhoNeutral;
    int nvtx; 

    double fGridRhoAll;
    double fGridRhoFastjetAll;
    double fGridRhoFastjetAllCalo;
    double fGridRhoFastjetCentralCalo;
    double fGridRhoFastjetCentralChargedPileUp;
    double fGridRhoFastjetCentralNeutral;
    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;

    ClassDef(Event,1)
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    virtual ~GenEvent() {}
  
    unsigned int processID;
    double ptHat;
    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent,1)
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}
    double eta;
    double phi;
    double pt;
    bool ecalDriven;
    bool isGap;
    bool hasGsfTrack;
    double trackPt;
    double energy;
    double caloEnergy;
    int charge;
    int pixHits;
    int trkHits;
    int nValidHits;
    double trkD0;
    double trkDz;
    // ID variables
    double hoe;
    double eop;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double deltaPhiTrkSC;
    double deltaEtaTrkSC;
    int classif;
  
    // Vertex
    double vx;
    double vy;
    double vz;

    // SC associated with electron
    double scEn;
    double scEta;
    double scPhi;
    double scET;
    double scRawEnergy;
    float BDT;  
    float BDTpreComp;  
    // Vertex association variables
    double dxyPV;
    double dzPV;
    double vtxDist3D;
    int vtxIndex;
    double vtxDistZ;
    float relIso;
    float pfRelIso;
    // PFlow isolation variable
    float chargedHadronIso;
    float neutralHadronIso;
    float photonIso;

    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;

    int missingHits;

    double dB;
    double edB;    

    double dB3D;
    double edB3D;    

    int nBrems;
    float fbrem;
  
    float hasMatchedConv;
    bool mvaPreselection;
    bool isTriggerElectron;
    int fidFlag;
    std::map<std::string, float> idmap;
    int selbit;
    bool passMediumId;
    bool passTightId;
    int mvaCategory;
    float ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values;
    std::map< std::string,std::vector<double> > isolationMap;
    std::vector<TLorentzVector> associatedPackedPFCandidatesP4;

    ClassDef(Electron, 1)
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
  
    double eta;
    double phi;
    double p;
    double px;
    double py;
    double pz;
    double pt;
    double energy;
    int pdgId;
    double vx;
    double vy;
    double vz;
    int status;
    double charge;
    int numDaught;
    int numMother;
    int motherIndex;
    std::vector<int> motherIndices;
    std::vector<int> daughtIndices;
  
    ClassDef(GenParticle,1)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    double eta;
    double phi;
    double p;
    double pt;
    double energy;
    double emf;
    double hadf;
  
    ClassDef(GenJet,1)
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    double met;
    double metphi;
    double sumet;
    double metuncorr;
    double metphiuncorr;
    double sumetuncorr;
  
    ClassDef(MET,1)
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    double eta;
    double phi;
    double pt;
    double energy;
    int charge;
    double mass;
  
    double dxyPV;
    double dzPV;
    int vtxIndex;
    double vtxDxy;
    double vtxDz;

    // Leading particle pT
    double leadChargedParticlePt;
    double leadNeutralParticlePt;
    double leadParticlePt;

    std::vector<vhtm::Candidate> sigChHadList;
    std::vector<vhtm::Candidate> sigNeHadList;
    std::vector<vhtm::Candidate> sigGammaList;
    std::vector<vhtm::Candidate> isoChHadList;
    std::vector<vhtm::Candidate> isoNeHadList;
    std::vector<vhtm::Candidate> isoGammaList;

    float ptSumChargedHadronsIsoCone;
    float ptSumNeutralHadronsIsoCone;
    float ptSumPhotonsIsoCone;

     // tau id. discriminators
     float decayModeFinding;
     float decayModeFindingNewDMs;
     float decayModeFindingOldDMs;

     // discriminators against electrons/muons
     float againstMuonLoose;
     float againstMuonMedium;
     float againstMuonTight;

     float againstMuonLoose3;
     float againstMuonTight3;

     float againstElectronLoose;
     float againstElectronMedium;
     float againstElectronTight;
     //float againstElectronMVA;
  
     float againstElectronLooseMVA5;
     float againstElectronMediumMVA5;
     float againstElectronTightMVA5;

     float byLooseCombinedIsolationDeltaBetaCorr3Hits;
     float byMediumCombinedIsolationDeltaBetaCorr3Hits;
     float byTightCombinedIsolationDeltaBetaCorr3Hits;
     float byCombinedIsolationDeltaBetaCorrRaw3Hits;
     float chargedIsoPtSum;
     float neutralIsoPtSum;
     float puCorrPtSum;

     // kinematic variables for PFJet associated to PFTau
    double jetPt;
    double jetEta;
    double jetPhi;
     float emFraction;
    double vx;
    double vy;
    double vz;
  
    double zvertex;
    double dxySig;
       int selbit;
  
    ClassDef(Tau, 1)
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
    bool isGlobalMuon;
    bool isTrackerMuon;
    bool isPFMuon;
    bool isghostCleaned;
    bool passTrackerhighPtid;

    double eta;
    double phi;
    double pt;

    double p;
    double energy;
    int charge;
    double trkD0;

    double trkDz;
    int muonBestTrackType;  
   
    double globalChi2;
    double tkNChi2;

    float trkIso;
    float ecalIso;
    float hcalIso;
    float hoIso;
    float pfChargedIsoR03;
    float pfChargedHadIsoR03;
    float pfNeutralHadIsoR03;
    float pfPhotonIso03;    
    float sumPUPt03;
    float pfRelIso03;

    float sumChargedParticlePt;
    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;
    float pfRelIso04;

    int passID;
    double dxyPV;
    double dzPV;
    double vtxDist3D;
    int vtxIndex;
    double vtxDistZ;
    int pixHits;
    int trkHits;
    int muoHits;
    int matches;
    int trackerLayersWithMeasurement;

    double vx;
    double vy;
    double vz;

    double dB;//2D
    double edB;    
  
    double dB3D;
    double edB3D;    

    // UW Recommendation
    bool isGlobalMuonPromptTight;
    bool isAllArbitrated;
    int nChambers;
    int nMatches;
    int nMatchedStations;

    unsigned int stationMask;
    unsigned int stationGapMaskDistance;
    unsigned int stationGapMaskPull;

    bool muonID;
   
    int nSegments;

    int selbit;
    std::map< std::string,std::vector<double> > isolationMap;
    ClassDef(Muon, 1)
//    ClassDef(Muon, 2)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
    int partonFlavour;

    float chargedEmEnergyFraction;
    float chargedHadronEnergyFraction;
    float chargedMuEnergyFraction;
    float electronEnergyFraction;
    float muonEnergyFraction;
    float neutralEmEnergyFraction;
    float neutralHadronEnergyFraction;
    float photonEnergyFraction;
    int chargedHadronMultiplicity;
    int chargedMultiplicity;
    int electronMultiplicity;
    int muonMultiplicity;
    int neutralHadronMultiplicity;
    int neutralMultiplicity;
    int photonMultiplicity;
    int nConstituents;

    //float simpleSecondaryVertexHighEffBTag;
    //float simpleSecondaryVertexHighPurBTag;

    float combinedSecondaryVertexBTag;
    //double combinedSecondaryVertexMVABTag;
    float combinedInclusiveSecondaryVertexBTag;
    float combinedInclusiveSecondaryVertexV2BJetTags;
    float pfCombinedInclusiveSecondaryVertexV2BJetTags;
    //double combinedMVABTag;
    std::map<std::string, float> discrimap;
    
    float jpumva;

    int passLooseID;
    int passTightID;
  
    int selbit;

    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    virtual ~Vertex() {}
  
    double x;
    double y;
    double z;
    double xErr;
    double yErr;
    double zErr;
    double rho;
    double chi2;
    double ndf;
    //int ntracks;
    //int ntracksw05;
    bool isfake;
    bool isvalid;
    //double sumPt; // vector sum
  
    int selbit;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    virtual ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    virtual ~TriggerObject() {}
  
    double energy;
    double pt;
    double eta;
    double phi;
  
    std::map<std::string, unsigned int> pathList;
  
    ClassDef(TriggerObject, 1)
  };
  class Photon : public TObject {
  public:
    Photon();
    virtual ~Photon() {}
    
    double et;
    double eta;
    double clusterEta;
    double phi;
    double clusterPhi;
    double energy;
    double theta; 
    double vx;
    double vy;
    double vz;
    
    double scEnergy;
    double scEta;
    double scPhi;
    double scSize;
    double scEtaWidth;
    double scPhiWidth;
    double scEt;
    double scRawEnergy;
    double scx;
    double scy;
    double scz; 
    double isoEcalRecHit03;
    double isoHcalRecHit03;
    double isoSolidTrkCone03;
    double isoHollowTrkCone03;
    int nTrkSolidCone03;
    int nTrkHollowCone03;
    
    double isoEcalRecHit04;
    double isoHcalRecHit04;
    double isoSolidTrkCone04;
    double isoHollowTrkCone04;
    int nTrkSolidCone04;
    int nTrkHollowCone04;
    
    bool isEB;
    bool isEE; 
    bool isEBGap;
    bool isEEGap;
    bool isEBEEGap;
    int fidFlag;
    
    bool passElectronVeto;
    bool hasPixelSeed;
    double ecalIso;
    double hcalIso;
    double trackIso;
    double chargedHadIso;
    double neutralHadIso;
    double photonIso;
    
    double r9;
    double hoe;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double e1x5;
    double e2x5; 
    double e3x3;
    double e5x5;
    double r1x5;
    double r2x5;
    double maxEnergyXtal;
    
    bool hasConversionTracks;
    int nTracks;
    bool isConverted;
    double pairInvMass;
    double pairCotThetaSeparation;
    double pairPx;
    double pairPy;
    double pairPz;
    double conv_vx;
    double conv_vy;
    double conv_vz;
    double eovp;
    double zpv;
    double distOfMinApproach;
    double dPhiTracksAtVtx;
    double dPhiTracksAtEcal;
    double dEtaTracksAtEcal;  
    
    int selbit;
  
    std::map< std::string,std::vector<double> > isolationMap;

    ClassDef(Photon, 1)
  };

class SelectedEvent : public TObject {
  public:
    SelectedEvent();
    virtual ~SelectedEvent(){}
    unsigned long int run;
    unsigned long int lumi;
    unsigned long int event;
    double mass4l;
    double mZ1;
    double mZ2;
    std::map<std::string,double>  kd;//D_bkg^kin,D_bkg,D_gg,D_HJJ^VBF,D_0
    int nJets;
    double jet1pt;
    double jet2pt;
    int category;
    double m4lRefit;
    double m4lRefitError;
    double weight;
    int flavour;//0=4mu,1=4e,2=2e2mu,3=2mu2e,4=wrong
    void reset();
    ClassDef(SelectedEvent,1)
};
class ZtnP : public TObject {
   public :
  ZtnP();
  virtual ~ZtnP(){}

  int  flavour;//0=mumu,1=ee,2=wrong   
  float         TnP_pt;   //[nTnP]
  float         TnP_eta;   //[nTnP]
  float         TnP_phi;   //[nTnP]
  float         TnP_mass;   //[nTnP]
  int           TnP_hasFSR;   //[nTnP]
  float         TnP_mll;   //[nTnP]
  int           TnP_l1_pdgId;   //[nTnP]
  float         TnP_l1_pt;   //[nTnP]
  float         TnP_l1_eta;   //[nTnP]
  float         TnP_l1_phi;   //[nTnP]
  float         TnP_l1_mass;   //[nTnP]
  int           TnP_l1_charge;   //[nTnP]
  int           TnP_l1_tightId;   //[nTnP]
  int           TnP_l1_looseId;   //[nTnP]
  float         TnP_l1_dxy;   //[nTnP]
  float         TnP_l1_dz;   //[nTnP]
  float         TnP_l1_edxy;   //[nTnP]
  float         TnP_l1_edz;   //[nTnP]
  float         TnP_l1_ip3d;   //[nTnP]
  float         TnP_l1_sip3d;   //[nTnP]
  float         TnP_l1_ptErr;   //[nTnP]
  int           TnP_l1_lostHits;   //[nTnP]
  int           TnP_l1_trackerLayers;   //[nTnP]
  int           TnP_l1_pixelLayers;   //[nTnP]
  float         TnP_l1_etaSc;   //[nTnP]
  int           TnP_l1_isGap;   //[nTnP]
  float         TnP_l1_r9;   //[nTnP]
  int           TnP_l1_convVeto;   //[nTnP]
  float         TnP_l1_mvaIdSpring15;   //[nTnP]
  float         TnP_l1_relIsoAfterFSR;   //[nTnP]
  float         TnP_l1_chargedHadIso03;   //[nTnP]
  int           TnP_l1_hasOwnFSR;   //[nTnP]
  int           TnP_l1_mcMatchId;   //[nTnP]
  int           TnP_l1_mcMatchAny;   //[nTnP]
  float         TnP_l1_mcPt;   //[nTnP]
  float         TnP_l1_mcPt1;   //[nTnP]
  int           TnP_l1_hlt1L;   //[nTnP]
  float         TnP_l1_p4WithFSR_pt;   //[nTnP]
  float         TnP_l1_p4WithFSR_eta;   //[nTnP]
  float         TnP_l1_p4WithFSR_phi;   //[nTnP]
  float         TnP_l1_p4WithFSR_mass;   //[nTnP]
  int           TnP_l2_pdgId;   //[nTnP]
  float         TnP_l2_pt;   //[nTnP]
  float         TnP_l2_eta;   //[nTnP]
  float         TnP_l2_phi;   //[nTnP]
  float         TnP_l2_mass;   //[nTnP]
  int           TnP_l2_charge;   //[nTnP]
  int           TnP_l2_tightId;   //[nTnP]
  int           TnP_l2_looseId;   //[nTnP]
  float         TnP_l2_dxy;   //[nTnP]
  float         TnP_l2_dz;   //[nTnP]
  float         TnP_l2_edxy;   //[nTnP]
  float         TnP_l2_edz;   //[nTnP]
  float         TnP_l2_ip3d;   //[nTnP]
  float         TnP_l2_sip3d;   //[nTnP]
  float         TnP_l2_ptErr;   //[nTnP]
  int           TnP_l2_lostHits;   //[nTnP]
  int           TnP_l2_trackerLayers;   //[nTnP]
  int           TnP_l2_pixelLayers;   //[nTnP]
  float         TnP_l2_etaSc;   //[nTnP]
  int           TnP_l2_isGap;   //[nTnP]
  float         TnP_l2_r9;   //[nTnP]
  int           TnP_l2_convVeto;   //[nTnP]
  float         TnP_l2_mvaIdSpring15;   //[nTnP]
  float         TnP_l2_relIsoAfterFSR;   //[nTnP]
  float         TnP_l2_chargedHadIso03;   //[nTnP]
  int           TnP_l2_hasOwnFSR;   //[nTnP]
  int           TnP_l2_mcMatchId;   //[nTnP]
  int           TnP_l2_mcMatchAny;   //[nTnP]
  float         TnP_l2_mcPt;   //[nTnP]
  float         TnP_l2_mcPt1;
  int           TnP_l2_hlt1L;   //[nTnP]
  float         TnP_l2_p4WithFSR_pt;   //[nTnP]
  float         TnP_l2_p4WithFSR_eta;   //[nTnP]
  float         TnP_l2_p4WithFSR_phi;   //[nTnP]
  float         TnP_l2_p4WithFSR_mass;   //[nTnP]
  ClassDef(ZtnP,1)
};
}
#endif
