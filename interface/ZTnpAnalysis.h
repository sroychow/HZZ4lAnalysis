#ifndef __ZTnpAnalysis__hh
#define __ZTnpAnalysis__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
//#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "PhysicsObjSelector.h"
#include "ZCandidate.h"

using std::ofstream;
using std::ifstream;
using std::ostream;
class ZTnpAnalysis: public PhysicsObjSelector {
    
public:
  enum EventType {
    mmmm = 0, eeee, eemm
  };
  enum MassRange {
    low = 0, high, fullblind, unblind 
  };
  ZTnpAnalysis();
  virtual ~ZTnpAnalysis();
  
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;

  void selectEvent() {};
  virtual void bookHistograms();
  void getZtnPpair();  
public:
  std::vector<vhtm::Vertex> vtxList_;
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >   tightMuon_;
  std::vector<vhtm::PackedPFCandidate>  fsrPVec_;
  std::vector<vhtm::TriggerObject> tObj_;
  std::vector<vhtm::GenParticle> genObj_;
  std::string treeFname_;
  unsigned long int run;
  unsigned long int lumi;
  unsigned long int event;
  int isData;
  
  TTree* outTree_;
  TFile* outTreeFile_;
  const static Int_t kMaxTnP = 8;
  //////////////////////////////////////////////////////////////
   int nTnP;
   float         TnP_pt[kMaxTnP];   
   float         TnP_eta[kMaxTnP];   
   float         TnP_phi[kMaxTnP];   
   float         TnP_mass[kMaxTnP];   
   int           TnP_hasFSR[kMaxTnP];   
   float         TnP_mll[kMaxTnP];   
   int           TnP_l1_pdgId[kMaxTnP];   
   float         TnP_l1_pt[kMaxTnP];   
   float         TnP_l1_eta[kMaxTnP];   
   float         TnP_l1_phi[kMaxTnP];   
   float         TnP_l1_mass[kMaxTnP];   
   int           TnP_l1_charge[kMaxTnP];   
   int           TnP_l1_tightId[kMaxTnP];   
   int           TnP_l1_looseId[kMaxTnP];   
   float         TnP_l1_dxy[kMaxTnP];   
   float         TnP_l1_dz[kMaxTnP];   
   float         TnP_l1_edxy[kMaxTnP];   
   float         TnP_l1_edz[kMaxTnP];   
   float         TnP_l1_ip3d[kMaxTnP];   
   float         TnP_l1_sip3d[kMaxTnP];   
   float         TnP_l1_ptErr[kMaxTnP];   
   int           TnP_l1_lostHits[kMaxTnP];   
   int           TnP_l1_trackerLayers[kMaxTnP];   
   int           TnP_l1_pixelLayers[kMaxTnP];   
   float         TnP_l1_etaSc[kMaxTnP];   
   int           TnP_l1_isGap[kMaxTnP];   
   float         TnP_l1_r9[kMaxTnP];   
   int           TnP_l1_convVeto[kMaxTnP];   
   float         TnP_l1_mvaIdSpring15[kMaxTnP];   
   float         TnP_l1_relIsoAfterFSR[kMaxTnP];   
   float         TnP_l1_chargedHadIso03[kMaxTnP];   
   int           TnP_l1_hasOwnFSR[kMaxTnP]; 

   int           TnP_l1_mcMatchId[kMaxTnP];   //[nTnP]
   int           TnP_l1_mcMatchAny[kMaxTnP];   //[nTnP]
   float         TnP_l1_mcPt[kMaxTnP];   //[nTnP]
   float         TnP_l1_mcPt1[kMaxTnP]; 
 
   int           TnP_l1_hlt1L[kMaxTnP];   
   float         TnP_l1_p4WithFSR_pt[kMaxTnP];   
   float         TnP_l1_p4WithFSR_eta[kMaxTnP];   
   float         TnP_l1_p4WithFSR_phi[kMaxTnP];   
   float         TnP_l1_p4WithFSR_mass[kMaxTnP];   
   int           TnP_l2_pdgId[kMaxTnP];   
   float         TnP_l2_pt[kMaxTnP];   
   float         TnP_l2_eta[kMaxTnP];   
   float         TnP_l2_phi[kMaxTnP];   
   float         TnP_l2_mass[kMaxTnP];   
   int           TnP_l2_charge[kMaxTnP];   
   int           TnP_l2_tightId[kMaxTnP];   
   int           TnP_l2_looseId[kMaxTnP];   
   float         TnP_l2_dxy[kMaxTnP];   
   float         TnP_l2_dz[kMaxTnP];   
   float         TnP_l2_edxy[kMaxTnP];   
   float         TnP_l2_edz[kMaxTnP];   
   float         TnP_l2_ip3d[kMaxTnP];   
   float         TnP_l2_sip3d[kMaxTnP];   
   float         TnP_l2_ptErr[kMaxTnP];   
   int           TnP_l2_lostHits[kMaxTnP];   
   int           TnP_l2_trackerLayers[kMaxTnP];   
   int           TnP_l2_pixelLayers[kMaxTnP];   
   float         TnP_l2_etaSc[kMaxTnP];   
   int           TnP_l2_isGap[kMaxTnP];   
   float         TnP_l2_r9[kMaxTnP];   
   int           TnP_l2_convVeto[kMaxTnP];   
   float         TnP_l2_mvaIdSpring15[kMaxTnP];   
   float         TnP_l2_relIsoAfterFSR[kMaxTnP];   
   float         TnP_l2_chargedHadIso03[kMaxTnP];   
   int           TnP_l2_hasOwnFSR[kMaxTnP];

   int           TnP_l2_mcMatchId[kMaxTnP];   //[nTnP]
   int           TnP_l2_mcMatchAny[kMaxTnP];   //[nTnP]
   float         TnP_l2_mcPt[kMaxTnP];   //[nTnP]
   float         TnP_l2_mcPt1[kMaxTnP];   

   int           TnP_l2_hlt1L[kMaxTnP];   
   float         TnP_l2_p4WithFSR_pt[kMaxTnP];   
   float         TnP_l2_p4WithFSR_eta[kMaxTnP];   
   float         TnP_l2_p4WithFSR_phi[kMaxTnP];   
   float         TnP_l2_p4WithFSR_mass[kMaxTnP];
  //////////////////////////////////////////////////////////////
};
#endif
