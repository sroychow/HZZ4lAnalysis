#ifndef __ZCandidate__hh
#define __ZCandidate__hh

#include "TLorentzVector.h"
struct ZCandidate {
  int l1Index;
  TLorentzVector l1P4;
  int l1Charge;
  double l1Isolation;
  TLorentzVector l1FsrP4;
  
  int l2Index;
  TLorentzVector l2P4;
  int l2Charge;
  double l2Isolation;
  TLorentzVector l2FsrP4;
  
  int fsrWithLep;
  TLorentzVector fsrPhoP4;
  int flavour; // 0: mumu, 1: ee, 2: others
  double mass;
  double massDiff;
};
#endif
