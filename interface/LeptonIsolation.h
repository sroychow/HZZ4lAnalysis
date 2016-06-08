#ifndef __LeptonIsolation__hh 
#define __LeptonIsolation__hh
#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

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
#include "PhysicsObjSelector.h"
#include "ZCandidate.h"
using std::string;
class LeptonIsolation : public PhysicsObjSelector {
    
public:

  LeptonIsolation();
  virtual ~LeptonIsolation();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  bool readJob(const string& jobFile, int& nFiles);  
  void clearLists();
  void selectEvent(){};
  virtual void bookHistograms();
  //void bookHistograms(TString folder);
  void bookDefaultIsoHistograms(TString dir);
  void readRhofitSlope();
  
  //void fillHistoforCone(std::string c,double mupt,double ch_had,double ch_lep,double nu_had,double nu_photon,double pu,
  //                      int ngoodVtx,double fGridRhoFastjetAll,const std::map<std::string, double>& lepCutMap,TString& leptype);
  template <typename T>
    void getLeptonIsolationInfo( std::vector<T> lepvec, double rho,int ngoodVtx,const std::map<std::string, double>& lepCutMap,
                                 TString leptype );

  void getMuonDefaultIso(std::vector<vhtm::Muon> muvec, int ngoodVtx, TString type, 
                         const double etaLow = 0.,const double etaHigh = 4.);
  void getElectronDefaultIso(std::vector<vhtm::Electron> elevec, int ngoodVtx, TString type, 
                             const double etaLow = 0., const double etaHigh = 4.);

  double pfisoEle( const vhtm::Electron& ele, double eventRho,double fsrPhotonEtSum, const string& cone);
  double pfisoMu(const vhtm::Muon& mu, double fsrPhotonEtSum, const string& cone);
  void getZllP4fromGen();
  template <typename T>
    void recoZSelector(const std::vector<T>& lepVec);
  template <typename T>
    void dogenMatching(std::vector<T> lepvec, std::vector<T>& sigvec, std::vector<T>& bkgvec);
  template <typename T>
    void dorecoZMatching(std::vector<T> lepvec, std::vector<T>& sigvec, std::vector<T>& bkgvec);
  template <typename T>
    void getClosestGenPartPdg( const T& lep );

  std::vector<vhtm::Vertex> vtxList;  
  bool _dumpEvent;
  bool dogenZmatching_;
  bool dorecoZmatching_;
  //bool useTightLepton_;
  std::map<std::string,std::string> cone_;
  std::vector<ZCandidate> Zllp4vec_,recoZllp4vec_;
  std::vector<vhtm::Muon> signalMuvec_,extraMuvec_;
  std::vector<vhtm::Electron> signalElevec_,extraElevec_;
};
#endif
