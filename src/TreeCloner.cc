#include "TreeCloner.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

TreeCloner::TreeCloner( const std::string& cloneFilename ) {
  cloneFilename_ =cloneFilename;
}
TreeCloner::~TreeCloner() {
  // delete cloneFile_;
  // delete newCloneTree_;
}
void TreeCloner::cloneTree(TTree* originalTree) {
  cloneFile_ = TFile::Open(cloneFilename_.c_str(),"RECREATE");
  cloneFile_->mkdir("treeCreator");
  cloneFile_->cd("treeCreator");
  newCloneTree_ = originalTree->CloneTree(0);
}
void TreeCloner::cloneEvent() {
  newCloneTree_->Fill();
}

void TreeCloner::writeTree() {
  newCloneTree_->AutoSave();
  cloneFile_->Close();
}
