//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 26 10:44:49 2023 by ROOT version 6.26/11
// from TTree AnalysisxTree/Go4FileStore
// found on file: new.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TTamex_FullEvent.h"

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxOutputEvent = 1;

   // Declaration of leaf types
 //TTamex_FullEvent *OutputEvent_;
   TTamex_FullEvent *OutputEvent_TGo4EventElement;
   Double_t        OutputEvent_fTimeDiff[12];
   vector<unsigned int> OutputEvent_SSY;
   vector<unsigned int> OutputEvent_SFP;
   vector<unsigned int> OutputEvent_TAM;
   vector<int>     OutputEvent_CHA;
   vector<int>     OutputEvent_Edge_type;
   vector<int>     OutputEvent_Coarse_ct;
   vector<int>     OutputEvent_Fine_time;

   // List of branches
   TBranch        *b_OutputEvent_TGo4EventElement;   //!
   TBranch        *b_OutputEvent_fTimeDiff;   //!
   TBranch        *b_OutputEvent_SSY;   //!
   TBranch        *b_OutputEvent_SFP;   //!
   TBranch        *b_OutputEvent_TAM;   //!
   TBranch        *b_OutputEvent_CHA;   //!
   TBranch        *b_OutputEvent_Edge_type;   //!
   TBranch        *b_OutputEvent_Coarse_ct;   //!
   TBranch        *b_OutputEvent_Fine_time;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     ShowDetail(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("new.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("new.root");
      }
      f->GetObject("AnalysisxTree",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   OutputEvent_TGo4EventElement = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("OutputEvent.TGo4EventElement", &OutputEvent_TGo4EventElement, &b_OutputEvent_TGo4EventElement);
   fChain->SetBranchAddress("OutputEvent.fTimeDiff[12]", OutputEvent_fTimeDiff, &b_OutputEvent_fTimeDiff);
   fChain->SetBranchAddress("OutputEvent.SSY", &OutputEvent_SSY, &b_OutputEvent_SSY);
   fChain->SetBranchAddress("OutputEvent.SFP", &OutputEvent_SFP, &b_OutputEvent_SFP);
   fChain->SetBranchAddress("OutputEvent.TAM", &OutputEvent_TAM, &b_OutputEvent_TAM);
   fChain->SetBranchAddress("OutputEvent.CHA", &OutputEvent_CHA, &b_OutputEvent_CHA);
   fChain->SetBranchAddress("OutputEvent.Edge_type", &OutputEvent_Edge_type, &b_OutputEvent_Edge_type);
   fChain->SetBranchAddress("OutputEvent.Coarse_ct", &OutputEvent_Coarse_ct, &b_OutputEvent_Coarse_ct);
   fChain->SetBranchAddress("OutputEvent.Fine_time", &OutputEvent_Fine_time, &b_OutputEvent_Fine_time);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void test::ShowDetail(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);

	fprintf(stdout, "\nOutputEvent_TAM\t\t");
	for (UInt_t iTAM:OutputEvent_TAM) fprintf(stdout,"%u\t",iTAM);
	fprintf(stdout, "\nOutputEvent_CHA\t\t");
	for (Int_t iCHA:OutputEvent_CHA) fprintf(stdout,"%d\t",iCHA);
	fprintf(stdout, "\nOutputEvent_Edge_type\t");
	for (Int_t iEdge_type:OutputEvent_Edge_type) fprintf(stdout,"%d\t",iEdge_type);
	fprintf(stdout, "\nOutputEvent_Coarse_ct\t");
	for (Int_t iCoarse_ct:OutputEvent_Coarse_ct) fprintf(stdout,"%d\t",iCoarse_ct);
	fprintf(stdout, "\nOutputEvent_Fine_time\t");
	for (Int_t iFine_time:OutputEvent_Fine_time) fprintf(stdout,"%d\t",iFine_time);
	fprintf(stdout, "\n");

}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
