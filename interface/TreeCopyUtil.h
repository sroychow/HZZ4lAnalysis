#ifndef __TreeCopyUtil__hh 
#define __TreeCopyUtil__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

class TTree;
class TChain;
class TFile;

#include "AnaBase.h"
#include "PhysicsObjects.h"
using std::string;
class TreeCopyUtil : public AnaBase {
    
 public:

  TreeCopyUtil();
  virtual ~TreeCopyUtil();
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();
  void selectEvent() {}
 private:
   std::string cloneFilename_;
  
};
#endif
