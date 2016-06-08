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
#include <typeinfo>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "HZZ4lUtil.h"
#include "AnaUtil.h"
#include "LeptonIsolation.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;
using std::setw;
using std::setiosflags;
using std::resetiosflags;
using std::ios;
using std::stoi;
using namespace vhtm;
// -----------
// Constructor
// -----------
LeptonIsolation::LeptonIsolation()
  : PhysicsObjSelector(),
    _dumpEvent(false), 
    dogenZmatching_(false),
    dorecoZmatching_(false)
{
  cone_["c15"]="015";
  cone_["c20"]="020";
  cone_["c25"]="025";
  cone_["c30"]="030";
  cone_["c35"]="035";
  cone_["c40"]="040";
  cone_["c45"]="045";
}
// ----------
// Destructor
// ----------
LeptonIsolation::~LeptonIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool LeptonIsolation::beginJob() 
{ 
  PhysicsObjSelector::beginJob();
  
  // Open the output ROOT file
  histf()->cd();
  PhysicsObjSelector::bookHistograms();
  bookHistograms();
  return true;
}
// ---------------
// Book Common histograms
// ---------------
void LeptonIsolation::bookHistograms() 
{
  histf()->cd();
  histf()->mkdir("Event");
  histf()->cd("Event");
  new TH1D("nGoodvtx","Number of Vertices",50,0,50);
  //Rho Values
  new TH1F("fGridRhoAll","Grid Rho for event",100,0.,100.);
  new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
  new TH1F("fGridRhoFastjetAllCalo","Grid Rho for event",100,0.,100.);
  new TH1F("fGridRhoFastjetCentralCalo","Grid Rho for event",100,0.,100.);
  new TH1F("fGridRhoFastjetCentralChargedPileUp","Grid Rho for event",100,0.,100.);
  new TH1F("fGridRhoFastjetCentralNeutral","Grid Rho for event",100,0.,100.);
  new TH1F("nLep","number of Tight Leptons(sig+bkg)",20,-0.5,19.5);
  new TH1F("mbkglZl1_mumu","Inv mass of bkg Leptons with leading lepton in reco Z",200,0.,200.);
  new TH1F("mbkglZl2_mumu","Inv mass of bkg Leptons with subleading lepton in reco Z",200,0.,200.);
  new TH1F("mbkglZl1_ee","Inv mass of bkg Leptons with leading lepton in reco Z",200,0.,200.);
  new TH1F("mbkglZl2_ee","Inv mass of bkg Leptons with subleading lepton in reco Z",200,0.,200.);
  bookDefaultIsoHistograms("SignalLepton");
  bookDefaultIsoHistograms("BkgLepton");
  bookDefaultIsoHistograms("SignalLepton_barrel");
  bookDefaultIsoHistograms("BkgLepton_barrel");
  bookDefaultIsoHistograms("SignalLepton_endcap");
  bookDefaultIsoHistograms("BkgLepton_endcap");
}
//------------------------------------
// book histograms for efficinecy study from default isolations
//------------------------------------
void LeptonIsolation::bookDefaultIsoHistograms(TString dir) 
{
  histf()->cd();
  histf()->mkdir(dir);
  histf()->cd(dir);
  new TH1F("nMu","number of Muons",20,-0.5,19.5);
  new TH1F("nEle","number of Electrons",20,-0.5,19.5);
  
  new TH1F("muPt", " Muon pt",200,-0.5,199.5);
  new TH1F("elePt"," Electron pt",200,-0.5,199.5);
  //default isolation from miniAOD
  new TProfile("IsoRelMuonOvsnVtx", "Delta-beta Corrected Relative Isolation Default vs nVertex", 80, 0, 80, 0.,1.);
  new TProfile("IsoRelElectronOvsnVtx", "Rho Corrected Relative Isolation Default vs nVertex", 80, 0, 80, 0.,1.);
  new TH1F("mupfrelIso_015","Muon rel Isolation for cone size 0.15",100,0.,1.);
  new TH1F("mupfrelIso_020","Muon rel Isolation for cone size 0.20",100,0.,1.);
  new TH1F("mupfrelIso_025","Muon rel Isolation for cone size 0.25",100,0.,1.);
  new TH1F("mupfrelIso_030","Muon rel Isolation for cone size 0.30",100,0.,1.);
  new TH1F("mupfrelIso_035","Muon rel Isolation for cone size 0.35",100,0.,1.);
  new TH1F("mupfrelIso_040","Muon rel Isolation for cone size 0.40",100,0.,1.);
  new TH1F("mupfrelIso_045","Muon rel Isolation for cone size 0.45",100,0.,1.);
  
  new TH1F("elepfrelIso_015","Electron rel Isolation for cone size 0.15",100,0.,1.);
  new TH1F("elepfrelIso_020","Electron rel Isolation for cone size 0.20",100,0.,1.);
  new TH1F("elepfrelIso_025","Electron rel Isolation for cone size 0.25",100,0.,1.);
  new TH1F("elepfrelIso_030","Electron rel Isolation for cone size 0.30",100,0.,1.);
  new TH1F("elepfrelIso_035","Electron rel Isolation for cone size 0.35",100,0.,1.);
  new TH1F("elepfrelIso_040","Electron rel Isolation for cone size 0.40",100,0.,1.);
  new TH1F("elepfrelIso_045","Electron rel Isolation for cone size 0.45",100,0.,1.);
  //just for comparison
  new TProfile("IsoRelMuonOvsnVtx_pf", "Delta-beta Corrected Relative Isolation PF vs nVertex", 80, 0, 80, 0.,1.);
  new TProfile("IsoRelElectronOvsnVtx_pf", "Rho Corrected Relative Isolation PF vs nVertex", 80, 0, 80, 0.,1.);
  //Delta-beta counter for muons
  new TH1F("mu_db_015", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_020", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_025", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_030", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_035", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_040", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  new TH1F("mu_db_045", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
  //rho counter for electrons
  new TH1F("ele_rho_015", "no of leptons passing rho cut for cone size 0.15", 7, -0.5, 6.5);
  new TH1F("ele_rho_020", "no of leptons passing rho cut for cone size 0.20", 7, -0.5, 6.5);
  new TH1F("ele_rho_025", "no of leptons passing rho cut for cone size 0.25", 7, -0.5, 6.5);
  new TH1F("ele_rho_030", "no of leptons passing rho cut for cone size 0.30", 7, -0.5, 6.5);
  new TH1F("ele_rho_035", "no of leptons passing rho cut for cone size 0.35", 7, -0.5, 6.5);
  new TH1F("ele_rho_040", "no of leptons passing rho cut for cone size 0.40", 7, -0.5, 6.5);
  new TH1F("ele_rho_045", "no of leptons passing rho cut for cone size 0.45", 7, -0.5, 6.5);
  if (dir.Contains("Signal")) {
    //Generator Level 
    new TH1F("genZMass","gen Z mass",200,0.,200.);
    new TH1F("genZl1Pt","leading lepton Pt from Z at gen ",200,0.,200.);
    new TH1F("genZl2Pt","leading lepton Pt from Z at gen ",200,0.,200.);
    //Reco Level
    new TH1F("recoZMass","reco Z mass",200,0.,200.);
    new TH1F("recoZl1Pt","leading lepton Pt from Z at reco ",200,0.,200.);
    new TH1F("recoZl2Pt","leading lepton Pt from Z at reco ",200,0.,200.);
    new TH1F("recoZflav","RecoZ flav",3,-0.5,2.5);
  }
  if (dir.Contains("BkgLepton")) {
    new TH1F("recoProbeGenPdg","PdgId of closest gen match particle of probe leptons ",1000,0.,1000.);
  }
}

//--------------------
double LeptonIsolation::pfisoMu(const vhtm::Muon& mu, double fsrPhotonEtSum, const string& cone) {
  //reference 
  //mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt 
  //- fsrPhotonEtSum - 0.5 * mu.sumPUPt));
  if (mu.isolationMap.find(cone) != mu.isolationMap.end()) 
    return (mu.isolationMap.at(cone).at(0) + std::max(0., mu.isolationMap.at(cone).at(2) + 
						      mu.isolationMap.at(cone).at(3) - fsrPhotonEtSum - 0.5 * mu.isolationMap.at(cone).at(4)));
  return 10000.;
}
//-----------------------------------------------------------------------------------
//Fill muon default isolation histograms
void LeptonIsolation::getMuonDefaultIso(std::vector<vhtm::Muon> muvec, int ngoodVtx, 
                                        TString type, const double etaLow, const double etaHigh) {
  histf()->cd();
  histf()->cd(type);
  AnaUtil::fillHist1D("nMu", muvec.size());
  for (const auto& mu : muvec ) {
    if( std::abs(mu.eta) < etaLow || std::abs(mu.eta) > etaHigh )   continue;
//    std::cout << "std::abs(mu.eta) <<std::endl; 
    AnaUtil::fillProfile("IsoRelMuonOvsnVtx",ngoodVtx,
                          HZZ4lUtil::pfiso(mu,0.)/mu.pt);
    AnaUtil::fillHist1D("muPt",mu.pt);

    for( auto& c : cone_ ) {
      double riso = pfisoMu(mu, 0., c.first)/mu.pt;
      AnaUtil::fillHist1D("mupfrelIso_" + c.second,riso);
      if( c.first == "c40") 
        AnaUtil::fillProfile("IsoRelMuonOvsnVtx_pf",ngoodVtx,riso);
      AnaUtil::fillHist1D("mu_db_" + c.second,0);
      ///Added parts for iso efficiency
      if( riso < 0.5 ) {
        AnaUtil::fillHist1D("mu_db_" + c.second,1);
        if( riso < 0.45 ){
          AnaUtil::fillHist1D("mu_db_" + c.second,2);
          if( riso < 0.4 ){
            AnaUtil::fillHist1D("mu_db_" + c.second,3);
            if( riso < 0.35 ){
              AnaUtil::fillHist1D("mu_db_" + c.second,4);
              if( riso < 0.3 ){
                AnaUtil::fillHist1D("mu_db_" + c.second,5);
                if( riso < 0.25 ){ 
                  AnaUtil::fillHist1D("mu_db_" + c.second,6);
                }
              }  
            }
          }
        }
      }
    }
  }
} 
//------------------------------------------------------------------------------------
double LeptonIsolation::pfisoEle( const vhtm::Electron& ele, double eventRho, 
                               double fsrPhotonEtSum, const string& cone) {
//reference 
//ele.chargedHadronIso + std::max(0., ele.neutralHadronIso + ele.photonIso  - fsrPhotonEtSum
//                                    - getEleRhoEffectiveArea(std::fabs(ele.eta)) * eventRho))
  double chhadiso = ele.isolationMap.at(cone).at(0);
  double nhadiso = ele.isolationMap.at(cone).at(2);
  double phiso = ele.isolationMap.at(cone).at(3);
  if( ele.isolationMap.find(cone) != ele.isolationMap.end() ) 
    return ( chhadiso + std::max(0., nhadiso + phiso - fsrPhotonEtSum
             - HZZ4lUtil::getEleRhoEffectiveArea(std::fabs(ele.eta)) * eventRho));
  return 10000.;
}

//Fill electron default isolation histograms
void LeptonIsolation::getElectronDefaultIso(std::vector<vhtm::Electron> elevec, 
                                            int ngoodVtx,TString type, const double etaLow, const double etaHigh) {
  histf()->cd();
  histf()->cd(type);
  AnaUtil::fillHist1D("nEle", elevec.size());
  for (const auto& ele : elevec ) {
    //std::cout << "etaLow=" << etaLow << "  etaHigh=" << etaHigh << "  " <<  std::abs(ele.eta) <<std::endl; 
    if( std::abs(ele.eta) < etaLow || std::abs(ele.eta) > etaHigh )   continue;
   
    AnaUtil::fillHist1D("elePt",ele.pt);

    AnaUtil::fillProfile("IsoRelElectronOvsnVtx",ngoodVtx,
                          HZZ4lUtil::pfiso(ele,getEventGridRho(),0.)/ele.pt);
    for( auto& c : cone_ ) {
      double riso = pfisoEle(ele, getEventGridRho(), 0., c.first)/ele.pt;
      AnaUtil::fillHist1D("elepfrelIso_" + c.second,riso);
      if( c.first == "c40") 
        AnaUtil::fillProfile("IsoRelElectronOvsnVtx_pf",ngoodVtx,riso);
      AnaUtil::fillHist1D("ele_rho_" + c.second,0);
      ///Added parts for iso efficiency
      if( riso < 0.5 ) {
        AnaUtil::fillHist1D("ele_rho_" + c.second,1);
        if( riso < 0.45 ){
          AnaUtil::fillHist1D("ele_rho_" + c.second,2);
          if( riso < 0.4 ){
            AnaUtil::fillHist1D("ele_rho_" + c.second,3);
            if( riso < 0.35 ){
              AnaUtil::fillHist1D("ele_rho_" + c.second,4);
              if( riso < 0.3 ){
                AnaUtil::fillHist1D("ele_rho_" + c.second,5);
                if( riso < 0.25 ){ 
                  AnaUtil::fillHist1D("ele_rho_" + c.second,6);
                }
              }  
            }
          }
        }
      }
    }  
  }
} 

// -------------------
// The main event loop
// -------------------
void LeptonIsolation::clearLists() {
 PhysicsObjSelector::clear();
 Zllp4vec_.clear();
 recoZllp4vec_.clear();
 extraMuvec_.clear();
 signalMuvec_.clear();
 extraElevec_.clear();
 signalElevec_.clear();
 vtxList.clear();
}
void LeptonIsolation::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = max(10000, nEvents()/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  //long int nRecoMuons=0,nGoodRecoMuons=0;
  std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
  for (int ev = 0; ev < nEvents(); ++ev) {

    clearEvent();
    clearLists();
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    //For data or for MC without pileup
    puevWt_ = 1;
    /*
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
    */
    // Show status of the run
    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << chain()->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

   op.verbose = (logOption() >> 1 & 0x1); 
   findVtxInfo(vtxList, op, fLog());
   
   if( vtxList.empty() )             continue;
  
   int ngoodVtx=vtxList.size();
   //to check triggered and non-triggered cases
   if (useTrigger() && !isTriggered(true, false))      continue;
   //dumpEvent(true);   
   //check if there is Z at gen level
   getZllP4fromGen();
   //no Z found at gen level
   if( !dorecoZmatching_ && dogenZmatching_ && Zllp4vec_.empty() )           continue;
   histf()->cd();
   histf()->cd("Event"); 
   AnaUtil::fillHist1D("nGoodvtx",ngoodVtx,1);
   AnaUtil::fillHist1D("fGridRhoAll",evt.fGridRhoAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAll",evt.fGridRhoFastjetAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAllCalo",evt.fGridRhoFastjetAllCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralCalo",evt.fGridRhoFastjetCentralCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralChargedPileUp",evt.fGridRhoFastjetCentralChargedPileUp,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralNeutral",evt.fGridRhoFastjetCentralNeutral,1);

   // main analysis object selection without iso
   findObjects(puevWt_);
   //!!!!!!!!!!!!!!Be very Careful!!!!!!!!!!!!!!!!!
   //when using loose leptons
   //const auto& loosemuVec_ = getLooseMuList();
   //when using tight leptons
   const auto& idmuVec_ = getTightMuList();

   //when using loose leptons
   //const auto& looseeleVec_ = getLooseEleList();
   //when using tight leptons
   const auto& ideleVec_ = getTightEleList();

   histf()->cd();
   histf()->cd("Event"); 
   AnaUtil::fillHist1D("nLep",ideleVec_.size() + idmuVec_.size(),1);

   TLorentzVector rZl1P4, rZl2P4;
   int rZl1Charge, rZl2Charge;
   if( dorecoZmatching_ ) {   
     //for reco Z with isolated leptons!!
     if( !idmuVec_.empty() )  {
       std::vector<vhtm::Muon> Isolep;
       for( auto& mu: idmuVec_ ) {
         if( pfisoMu(mu,0.,"c40")/mu.pt >= 0.4 )     continue;
         Isolep.push_back(mu);         
       }
       if( Isolep.size() >= 2)     recoZSelector<vhtm::Muon>(Isolep);
       //std::cout << "Iso mu pair found" << std::endl;
     }
     if( !ideleVec_.empty())  {
       std::vector<vhtm::Electron> Isolep;
       for( auto& ele: ideleVec_ ) {
         if( pfisoEle(ele,getEventGridRho(),0.,"c40")/ele.pt >= 0.5 )     continue;
         Isolep.push_back(ele);         
       }
       if( Isolep.size() >= 2)     recoZSelector<vhtm::Electron>(Isolep);
       //std::cout << "Iso ele pair found" << std::endl;
     }
     if(recoZllp4vec_.empty())   continue;

     //std::cout << "Reco Z found" << std::endl;
     histf()->cd();
     histf()->cd("SignalLepton");
     AnaUtil::fillHist1D("recoZMass",recoZllp4vec_[0].mass);
     AnaUtil::fillHist1D("recoZflav",recoZllp4vec_[0].flavour);
     if( recoZllp4vec_[0].l1P4.Pt() > recoZllp4vec_[0].l2P4.Pt() ) {
       rZl1P4 = recoZllp4vec_[0].l1P4;
       rZl2P4 = recoZllp4vec_[0].l2P4;
       rZl1Charge = recoZllp4vec_[0].l1Charge;
       rZl2Charge = recoZllp4vec_[0].l2Charge;
       AnaUtil::fillHist1D("recoZl1Pt",recoZllp4vec_[0].l1P4.Pt());
       AnaUtil::fillHist1D("recoZl2Pt",recoZllp4vec_[0].l2P4.Pt());
     }
     else {
       rZl2P4 = recoZllp4vec_[0].l1P4;
       rZl1P4 = recoZllp4vec_[0].l2P4;
       rZl2Charge = recoZllp4vec_[0].l1Charge;
       rZl1Charge = recoZllp4vec_[0].l2Charge;
       AnaUtil::fillHist1D("recoZl1Pt",recoZllp4vec_[0].l2P4.Pt());
       AnaUtil::fillHist1D("recoZl2Pt",recoZllp4vec_[0].l1P4.Pt());
     }
   } 

   if( !idmuVec_.empty() ) {
      if( dorecoZmatching_ ) {
        dorecoZMatching<vhtm::Muon>(idmuVec_,signalMuvec_,extraMuvec_);
        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton");
        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton");

        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton_barrel",1.2,2.4);

        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton_endcap",0,1.2);
        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton_endcap",1.2,2.4);

        for( auto& m: extraMuvec_) 
          getClosestGenPartPdg<vhtm::Muon>(m); 
        if( recoZllp4vec_[0].flavour == HZZ4lUtil::ZType::mumu ) {
          histf()->cd();
          histf()->cd("Event");
          for( auto& m: extraMuvec_) {
            if( rZl1Charge + m.charge == 0 ) 
              AnaUtil::fillHist1D("mbkglZl1_mumu", (rZl1P4 + HZZ4lUtil::getP4(m)).M());
            else
              AnaUtil::fillHist1D("mbkglZl2_mumu", (rZl2P4 + HZZ4lUtil::getP4(m)).M());
          }
        }
        
      }
      else if( dogenZmatching_ && !dorecoZmatching_ ) {
        dogenMatching<vhtm::Muon>(idmuVec_,signalMuvec_,extraMuvec_);

        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton");
        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton");
 
        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
        getMuonDefaultIso(signalMuvec_,ngoodVtx,"SignalLepton_endcap",1.2,2.4);

        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton_barrel",0,1.2);
        getMuonDefaultIso(extraMuvec_,ngoodVtx,"BkgLepton_endcap",1.2,2.4);
       } else {
         if (AnaUtil::cutValue(evselCutMap(), "Isbkg")) {
           getMuonDefaultIso(idmuVec_,ngoodVtx,"BkgLepton");
           getMuonDefaultIso(idmuVec_,ngoodVtx,"BkgLepton_barrel",0.,1.2);
           getMuonDefaultIso(idmuVec_,ngoodVtx,"BkgLepton_endcap",1.2,2.4);
         }
         else {
          getMuonDefaultIso(idmuVec_,ngoodVtx,"SignalLepton");
          getMuonDefaultIso(idmuVec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
          getMuonDefaultIso(idmuVec_,ngoodVtx,"SignalLepton_endcap",1.2,2.4);
         }
       }
      //getLeptonIsolationInfo<vhtm::Muon>( idmuVec_, evt.fGridRhoFastjetAll,
      //                                    ngoodVtx,muonCutMap(),"muon" );
      
   }
   if( !ideleVec_.empty() ) { 
     if( dorecoZmatching_ ) {//for DY case without gen matching
         dorecoZMatching<vhtm::Electron>(ideleVec_,signalElevec_,extraElevec_);
         getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton",0.,6.);
         getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton",0.,6.);

         getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
         getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton_barrel",0.,1.2);

         getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton_endcap",1.2,2.5);
         getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton_endcap",1.2,2.5);

         for( auto& e: extraElevec_) 
           getClosestGenPartPdg<vhtm::Electron>(e);
         if( recoZllp4vec_[0].flavour == HZZ4lUtil::ZType::ee ) {
           histf()->cd();
           histf()->cd("Event");
           for( auto& e: extraElevec_) { 
              if( rZl1Charge + e.charge == 0 ) 
                AnaUtil::fillHist1D("mbkglZl1_ee", (rZl1P4 + HZZ4lUtil::getP4(e)).M());
              else
                AnaUtil::fillHist1D("mbkglZl2_ee", (rZl2P4 + HZZ4lUtil::getP4(e)).M());
            }
         }
     }
     else if(dogenZmatching_ && !dorecoZmatching_ ) {//for Signal,DY case with genZ matching
       dogenMatching<vhtm::Electron>(ideleVec_,signalElevec_,extraElevec_);
       getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton");
       getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton");

       getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
       getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton_barrel",0.,1.2);

       getElectronDefaultIso(signalElevec_,ngoodVtx,"SignalLepton_endcap",1.2,2.5);
       getElectronDefaultIso(extraElevec_,ngoodVtx,"BkgLepton_endcap",1.2,2.5);
     } else {//for qcd
       if (AnaUtil::cutValue(evselCutMap(), "Isbkg")) {
         getElectronDefaultIso(ideleVec_,ngoodVtx,"BkgLepton");
         getElectronDefaultIso(ideleVec_,ngoodVtx,"BkgLepton_barrel",0.,1.2);
         getElectronDefaultIso(ideleVec_,ngoodVtx,"BkgLepton_endcap",1.2,2.5);
       }
       else {
         getElectronDefaultIso(ideleVec_,ngoodVtx,"SignalLepton");
         getElectronDefaultIso(ideleVec_,ngoodVtx,"SignalLepton_barrel",0.,1.2);
         getElectronDefaultIso(ideleVec_,ngoodVtx,"SignalLepton_endcap",1.2,2.5);
       }
     }
     //getLeptonIsolationInfo<vhtm::Electron>( ideleVec_, evt.fGridRhoFastjetAll,ngoodVtx,
     //                                         electronCutMap(),"electron" );
   }
  }
  // Analysis is over
  endJob();
}

void LeptonIsolation::getZllP4fromGen() {
//gen tau filter
   bool foundTau= false;
   for (const auto& gp: *genParticleColl()) {
     int pdgid = std::abs(gp.pdgId);
     if( pdgid == 15 ) {
       foundTau= true;
       break;
     }
   }
  if(!foundTau) {
  //std::cout << "New Event" << std::endl;
  histf()->cd();
  histf()->cd("SignalLepton");
  for (const auto& gp: *genParticleColl()) {
    int pdgid = std::abs(gp.pdgId);
    
    //int status = gp.status;
    if( pdgid != 23 )   continue;
    
    //std::cout << "Zfound at gen" << std::endl;
    if( gp.daughtIndices.size() < 2 )   continue;
    TLorentzVector ZP4 = HZZ4lUtil::getP4(gp);
    ZCandidate ztemp; 
    ztemp.mass = ZP4.M();
    AnaUtil::fillHist1D("genZMass", ztemp.mass);
    int nZl = 0;
    for (auto di: gp.daughtIndices ) {
      if (di >= ngenparticle()) continue;
      const GenParticle& dgp = genParticleColl()->at(di);
      int pid = std::abs(dgp.pdgId);
      if( pid == 23)      continue;
      if ( pid == 11)  ztemp.flavour = HZZ4lUtil::ZType::ee;
      else if ( pid == 13)  ztemp.flavour = HZZ4lUtil::ZType::mumu;
      else  ztemp.flavour = HZZ4lUtil::ZType::wrong;
      //std::cout << "Zdau pid=" << pid <<",Pt=" << HZZ4lUtil::getP4(dgp).Pt() << std::endl;
      //std::cout << nZl << std::endl;
      if( nZl == 0 ) {
        ztemp.l1P4 = HZZ4lUtil::getP4(dgp);
        ztemp.l1Charge = dgp.charge;
        nZl++;
      }
      else {
        ztemp.l2P4 = HZZ4lUtil::getP4(dgp);
        ztemp.l2Charge = dgp.charge;
      }
    }
    //std::cout << "Pdg" << pdgid << ",status=" << gp.status << ",dausize=" << gp.daughtIndices.size() 
    //          <<",GenZl1Pt=" << ztemp.l1P4.Pt() << "," << "GenZl2Pt=" << ztemp.l2P4.Pt() << std::endl;
    if( ztemp.l1P4.Pt() > ztemp.l2P4.Pt() ) {
      AnaUtil::fillHist1D("genZl1Pt", ztemp.l1P4.Pt());
      AnaUtil::fillHist1D("genZl2Pt", ztemp.l2P4.Pt());
    }
    else {
      AnaUtil::fillHist1D("genZl1Pt", ztemp.l2P4.Pt());
      AnaUtil::fillHist1D("genZl2Pt", ztemp.l1P4.Pt());
    }
   Zllp4vec_.push_back(ztemp);
  }
  }//found tau 
}

template <typename T>
void LeptonIsolation::getClosestGenPartPdg( const T& lep ) {
  TLorentzVector lepP4 = HZZ4lUtil::getP4(lep);
  int closestPdg=-1;
  double closestdR = 9999.;
  for (const auto& gp: *genParticleColl()) {
    TLorentzVector genP4 = HZZ4lUtil::getP4(gp);
    if( genP4.Pt() > 0. ) {
      if( lepP4.DeltaR(genP4) < closestdR ) 
        closestPdg = std::abs(gp.pdgId);
    }
  }
  histf()->cd();
  histf()->cd("BkgLepton");
  AnaUtil::fillHist1D("recoProbeGenPdg", closestPdg);
}

//do gen matching of the leptons from Z
template <typename T>
void LeptonIsolation::dogenMatching(std::vector<T> lepvec, std::vector<T>& sigvec, 
                                    std::vector<T>& bkgvec) {
  for(auto& lep : lepvec) {
    TLorentzVector lepP4 = HZZ4lUtil::getP4(lep);
    bool ismatched = false;
    for( auto& genZ: Zllp4vec_) {
      bool m1 = lepP4.DeltaR(genZ.l1P4) < 0.3 && std::fabs(genZ.l1P4.Pt() - lepP4.Pt()) < 25.;
      bool m2 = lepP4.DeltaR(genZ.l2P4) < 0.3 && std::fabs(genZ.l2P4.Pt() - lepP4.Pt()) < 25.;
      if( m1 || m2 ) {
        ismatched = true;
        break;
      }
    }//genz loop
    if(ismatched)  sigvec.push_back(lep);
    else  bkgvec.push_back(lep);
  }//lep loop
}

//find Z at reco level 
template <typename T>
void LeptonIsolation::recoZSelector(const std::vector<T>& lepVec) {
  for (unsigned int i = 0; i < lepVec.size(); ++i) {
    const auto& ip = lepVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip);
    for (unsigned int j = i+1; j < lepVec.size(); ++j) {
      const auto& jp = lepVec[j];
      // opposite charge
      if (ip.charge + jp.charge != 0) continue; 
      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp);
      if( lep1P4.DeltaR(lep2P4) < 0.2 )   continue;
      double Zmass = (lep1P4 + lep2P4).M();
      //std::cout << "Zmass=" << Zmass << std::endl;
      if( Zmass < AnaUtil::cutValue(evselCutMap(), "dyRecoZMlow") || 
          Zmass > AnaUtil::cutValue(evselCutMap(), "dyRecoZMhigh") )   continue;
      //std::cout << "Masscut Passed!!!" << std::endl;
      bool ptThA = lep1P4.Pt() > AnaUtil::cutValue(evselCutMap(), "l1Ptcut") && 
                   lep2P4.Pt() > AnaUtil::cutValue(evselCutMap(), "l2Ptcut");
      bool ptThB = lep2P4.Pt() > AnaUtil::cutValue(evselCutMap(), "l1Ptcut") && 
                   lep1P4.Pt() > AnaUtil::cutValue(evselCutMap(), "l2Ptcut");
      if( !ptThA && !ptThB )              continue;
      //std::cout << "Pt Cut Passed!!!" <<std::endl;
      ZCandidate ztmp;
      ztmp.l1P4 = lep1P4;
      ztmp.l1Charge = ip.charge;
      ztmp.l2P4 = lep2P4;
      ztmp.l2Charge = jp.charge;
      ztmp.mass = Zmass;
      ztmp.massDiff = std::fabs(Zmass - HZZ4lUtil::MZnominal);
      if (typeid(jp) == typeid(vhtm::Muon))
        ztmp.flavour = HZZ4lUtil::ZType::mumu;
      else if (typeid(jp) == typeid(vhtm::Electron))
        ztmp.flavour = HZZ4lUtil::ZType::ee;
      if( recoZllp4vec_.empty() )   recoZllp4vec_.push_back(ztmp);
      else {
        if( ztmp.massDiff < recoZllp4vec_[0].massDiff ) {
          recoZllp4vec_.clear();
          recoZllp4vec_.push_back(ztmp);
        }
      }
    }//lep j loop
  }//lep i loop  
}

template <typename T>
void LeptonIsolation::dorecoZMatching(std::vector<T> lepvec, std::vector<T>& sigvec, std::vector<T>& bkgvec) {
  for (auto& l: lepvec) {
   TLorentzVector lP4 = HZZ4lUtil::getP4(l);
   if( lP4 == recoZllp4vec_[0].l1P4 || lP4 == recoZllp4vec_[0].l2P4 )
     sigvec.push_back(l);
   else {
     if( lP4.DeltaR(recoZllp4vec_[0].l1P4) > 0.5 && lP4.DeltaR(recoZllp4vec_[0].l2P4) > 0.5 )
       bkgvec.push_back(l);
   }
  } 
}
bool LeptonIsolation::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    string key = tokens[0];
    string value = tokens[1];
    if (key == "dogenZmatching")
      dogenZmatching_ = (stoi(value.c_str()) > 0) ? true : false;
    if (key == "dorecoZmatching")
      dorecoZmatching_ = (stoi(value.c_str()) > 0) ? true : false;
    tokens.clear();
  }
  // Close the file
  fin.close();
  printJob();
  return true;
}
void LeptonIsolation::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
