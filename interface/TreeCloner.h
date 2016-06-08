#ifndef __TreeCloner_h
#define __TreeCloner_h
#include <string>
class TTree;
class TChain;
class TFile;
class TreeCloner {
  public :
    TreeCloner( const std::string& cloneFilename);
    ~TreeCloner();
    void cloneTree(TTree* originalTree);
    void cloneEvent();
    void writeTree();
  private :
    std::string cloneFilename_;
    TFile* cloneFile_;
    TTree* newCloneTree_;
};
#endif
