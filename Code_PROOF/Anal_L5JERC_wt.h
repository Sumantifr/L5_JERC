//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 31 15:28:52 2018 by ROOT version 5.34/36
// from TTree T1/WPrimeNtuple
// found on file: root://se01.indiacms.res.in//store/user/chatterj/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_crab_L5JERC_TTBar_SemiLeptonic_JECtxt/180725_134430/0000/rootuple_jerc_l5_106.root
//////////////////////////////////////////////////////////

#ifndef Anal_L5JERC_wt_h
#define Anal_L5JERC_wt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include "TLorentzVector.h"
#include "TRandom.h"

#include <TProofOutputFile.h>


int getbinid(double val, int nbmx, float* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double eta_to_theta(double eta){
return(2*atan(exp(-2*eta)));
}

double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Anal_L5JERC_wt : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           irun;
   Int_t           ilumi;
   UInt_t          ievt;
   Int_t           nprim;
   Int_t           nvert;
   Int_t           ndofct;
   Int_t           nchict;
   Double_t        Rho;
   Double_t        event_weight;
   Float_t         qscale;
   Int_t           npu_vert;
  
   static const int njetmx =20;
  

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ilumi;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_nprim;   //!
   TBranch        *b_nvert;   //!
   TBranch        *b_ndofct;   //!
   TBranch        *b_nchict;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_npu_vert;   //!
 

   Anal_L5JERC_wt(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Anal_L5JERC_wt() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   static const int noptbins = 68;
   static const int netarange = 10 ;
   static const int nopbins = 3;
   static const int nomassbins = 3;
   static const int netaea = 7;
   static const int selneta = 3;
   static const int selnmass = 4;
   
   float min_wmass = 0; 
   float max_wmass = 200;
   static const int wmassbins = 100;
   float min_tmass = 0;
   float max_tmass = 500.;
   static const int tmassbins = 250;

   float ptcut = 25.;
   float ycut = 2.4;
   float bream_rad = 0.1;
   float btagvalue = 0.8385;
//   float weight;
   
   float ydiv = 1.4;
//float seleta[seleta+1] = {1.4,2.5,5.};
  float selmass[selnmass+1] = {220,270,320,370};
  
   
   const float pival = acos(-1.);
   
   TProofOutputFile *OutFile;
   TFile *fileOut; 
//  TProofOutputFile *OutTxt;
   ofstream fp ;

   TH1D *h_nvert;
   TH1D *h_ndofct;
   TH1D *h_nchict;
   TH1D *h_npuvert;

  float weight;
  
   double tot_weight = 0;  
   int nevent_total = 0;
   
   char name[100];
   char title[100];
   
   bool isMC;

   ClassDef(Anal_L5JERC_wt,0);
};

#endif

#ifdef Anal_L5JERC_wt_cxx
void Anal_L5JERC_wt::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("irun", &irun, &b_irun);
   fChain->SetBranchAddress("ilumi", &ilumi, &b_ilumi);
   fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
   fChain->SetBranchAddress("nprim", &nprim, &b_nprim);
   fChain->SetBranchAddress("nvert", &nvert, &b_nvert);
   fChain->SetBranchAddress("ndofct", &ndofct, &b_ndofct);
   fChain->SetBranchAddress("nchict", &nchict, &b_nchict);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("npu_vert", &npu_vert, &b_npu_vert);
}

Bool_t Anal_L5JERC_wt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Anal_L5JERC_wt_cxx
