//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 13 22:23:53 2023 by ROOT version 6.26/11
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
   Double_t        OutputEvent_fTimeDiff[36];
   vector<unsigned int> OutputEvent_flip_SSY;
   vector<unsigned int> OutputEvent_flip_SFP;
   vector<unsigned int> OutputEvent_flip_TAM;
   vector<int>     OutputEvent_flip_TCHA;
   vector<int>     OutputEvent_flip_Edge_type;
   vector<int>     OutputEvent_flip_Coarse_ct;
   vector<int>     OutputEvent_flip_TDL;
   vector<unsigned int> OutputEvent_hit_SSY;
   vector<unsigned int> OutputEvent_hit_SFP;
   vector<unsigned int> OutputEvent_hit_TAM;
   vector<int>     OutputEvent_hit_PCHA;
   vector<double>  OutputEvent_hit_STOT;
   vector<double>  OutputEvent_hit_STle;
   vector<double>  OutputEvent_hit_FTOT;
   vector<double>  OutputEvent_hit_FTle;
   vector<double>  OutputEvent_hit_TTS;

   // List of branches
   TBranch        *b_OutputEvent_TGo4EventElement;   //!
   TBranch        *b_OutputEvent_fTimeDiff;   //!
   TBranch        *b_OutputEvent_flip_SSY;   //!
   TBranch        *b_OutputEvent_flip_SFP;   //!
   TBranch        *b_OutputEvent_flip_TAM;   //!
   TBranch        *b_OutputEvent_flip_TCHA;   //!
   TBranch        *b_OutputEvent_flip_Edge_type;   //!
   TBranch        *b_OutputEvent_flip_Coarse_ct;   //!
   TBranch        *b_OutputEvent_flip_TDL;   //!
   TBranch        *b_OutputEvent_hit_SSY;   //!
   TBranch        *b_OutputEvent_hit_SFP;   //!
   TBranch        *b_OutputEvent_hit_TAM;   //!
   TBranch        *b_OutputEvent_hit_PCHA;   //!
   TBranch        *b_OutputEvent_hit_STOT;   //!
   TBranch        *b_OutputEvent_hit_STle;   //!
   TBranch        *b_OutputEvent_hit_FTOT;   //!
   TBranch        *b_OutputEvent_hit_FTle;   //!
   TBranch        *b_OutputEvent_hit_TTS;   //!

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
   fChain->SetBranchAddress("OutputEvent.fTimeDiff[36]", OutputEvent_fTimeDiff, &b_OutputEvent_fTimeDiff);
   fChain->SetBranchAddress("OutputEvent.flip_SSY", &OutputEvent_flip_SSY, &b_OutputEvent_flip_SSY);
   fChain->SetBranchAddress("OutputEvent.flip_SFP", &OutputEvent_flip_SFP, &b_OutputEvent_flip_SFP);
   fChain->SetBranchAddress("OutputEvent.flip_TAM", &OutputEvent_flip_TAM, &b_OutputEvent_flip_TAM);
   fChain->SetBranchAddress("OutputEvent.flip_TCHA", &OutputEvent_flip_TCHA, &b_OutputEvent_flip_TCHA);
   fChain->SetBranchAddress("OutputEvent.flip_Edge_type", &OutputEvent_flip_Edge_type, &b_OutputEvent_flip_Edge_type);
   fChain->SetBranchAddress("OutputEvent.flip_Coarse_ct", &OutputEvent_flip_Coarse_ct, &b_OutputEvent_flip_Coarse_ct);
   fChain->SetBranchAddress("OutputEvent.flip_TDL", &OutputEvent_flip_TDL, &b_OutputEvent_flip_TDL);
   fChain->SetBranchAddress("OutputEvent.hit_SSY", &OutputEvent_hit_SSY, &b_OutputEvent_hit_SSY);
   fChain->SetBranchAddress("OutputEvent.hit_SFP", &OutputEvent_hit_SFP, &b_OutputEvent_hit_SFP);
   fChain->SetBranchAddress("OutputEvent.hit_TAM", &OutputEvent_hit_TAM, &b_OutputEvent_hit_TAM);
   fChain->SetBranchAddress("OutputEvent.hit_PCHA", &OutputEvent_hit_PCHA, &b_OutputEvent_hit_PCHA);
   fChain->SetBranchAddress("OutputEvent.hit_STOT", &OutputEvent_hit_STOT, &b_OutputEvent_hit_STOT);
   fChain->SetBranchAddress("OutputEvent.hit_STle", &OutputEvent_hit_STle, &b_OutputEvent_hit_STle);
   fChain->SetBranchAddress("OutputEvent.hit_FTOT", &OutputEvent_hit_FTOT, &b_OutputEvent_hit_FTOT);
   fChain->SetBranchAddress("OutputEvent.hit_FTle", &OutputEvent_hit_FTle, &b_OutputEvent_hit_FTle);
   fChain->SetBranchAddress("OutputEvent.hit_TTS", &OutputEvent_hit_TTS, &b_OutputEvent_hit_TTS);
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
   fChain->Show(entry,36);
}

void test::ShowDetail(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry,36);

    fprintf(stdout, "\nOutputEvent_flip_TAM\t\t");
    for (UInt_t iTAM:OutputEvent_flip_TAM) fprintf(stdout,"%u\t",iTAM);
    fprintf(stdout, "\nOutputEvent_flip_TCHA\t\t");
    for (Int_t iCHA:OutputEvent_flip_TCHA) fprintf(stdout,"%d\t",iCHA);
    fprintf(stdout, "\nOutputEvent_flip_Edge_type\t");
    for (Int_t iEdge_type:OutputEvent_flip_Edge_type) fprintf(stdout,"%d\t",iEdge_type);
    fprintf(stdout, "\nOutputEvent_flip_Coarse_ct\t");
    for (Int_t iCoarse_ct:OutputEvent_flip_Coarse_ct) fprintf(stdout,"%d\t",iCoarse_ct);
    fprintf(stdout, "\nOutputEvent_flip_TDL\t\t");
    for (Int_t iTDL:OutputEvent_flip_TDL) fprintf(stdout,"%d\t",iTDL);
    fprintf(stdout, "\n");


    fprintf(stdout, "\nOutputEvent_hit_SSY\t\t");
    for (UInt_t iSSY:OutputEvent_hit_SSY) fprintf(stdout,"%u\t",iSSY);
    fprintf(stdout, "\nOutputEvent_hit_SFP\t\t");
    for (UInt_t iSFP:OutputEvent_hit_SFP) fprintf(stdout,"%u\t",iSFP);
    fprintf(stdout, "\nOutputEvent_hit_TAM\t\t");
    for (UInt_t iTAM:OutputEvent_hit_TAM) fprintf(stdout,"%u\t",iTAM);
    fprintf(stdout, "\nOutputEvent_hit_PCHA\t\t");
    for (Int_t iPCHA:OutputEvent_hit_PCHA) fprintf(stdout,"%d\t",iPCHA);
    fprintf(stdout, "\nOutputEvent_hit_STOT\t\t");
    for (Double_t iSTOT:OutputEvent_hit_STOT) fprintf(stdout,"%.0f\t",iSTOT);
    fprintf(stdout, "\nOutputEvent_hit_STle\t\t");
    for (Double_t iSTle:OutputEvent_hit_STle) fprintf(stdout,"%.0f\t",iSTle);
    fprintf(stdout, "\nOutputEvent_hit_FTOT\t\t");
    for (Double_t iFTOT:OutputEvent_hit_FTOT) fprintf(stdout,"%.0f\t",iFTOT);
    fprintf(stdout, "\nOutputEvent_hit_FTle\t\t");
    for (Double_t iFTle:OutputEvent_hit_FTle) fprintf(stdout,"%.0f\t",iFTle);
    fprintf(stdout, "\nOutputEvent_hit_TTS\t\t");
    for (Double_t iTTS:OutputEvent_hit_TTS) fprintf(stdout,"%0.f\t",iTTS);
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
