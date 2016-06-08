#include "configana.h"
#include "TreeCopyUtil.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include <string>
#include <iomanip>

using std::cout;
using std::setw;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace vhtm;

// -----------
// Constructor
// -----------
TreeCopyUtil::TreeCopyUtil()
  : AnaBase()
{
}
// ----------
// Destructor
// ----------
TreeCopyUtil::~TreeCopyUtil() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool TreeCopyUtil::beginJob() 
{ 
  AnaBase::beginJob();
  return true;
}
// -------------------
// The main event loop
// -------------------
void TreeCopyUtil::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = std::max(10000, nEvents()/1000);

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int fevt = (firstEvent() > -1) ? firstEvent() : 0;
  int levt = (lastEvent() > -1) ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent();
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

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
    if (doEventClone()) 
      cloner()->cloneEvent();
  }
  // Analysis is over
  endJob();
}

void TreeCopyUtil::endJob() {
  closeFiles();
  if (doEventClone())
    cloner()->writeTree();
}
