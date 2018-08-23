#define Anal_L5JERC_wt_cxx
// The class definition in Anal_L5JERC_wt.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("Anal_L5JERC_wt.C")
// Root > T->Process("Anal_L5JERC_wt.C","some options")
// Root > T->Process("Anal_L5JERC_wt.C+")
//

#include "Anal_L5JERC_wt.h"
#include <TH2.h>
#include <TStyle.h>
#include "TProofServ.h"

void Anal_L5JERC_wt::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
 /*
    float ptcut = 25.;
	float ycut = 2.4;
	float bream_rad = 0.1;
	float btagvalue = 0.8385;
	float weight;
	
	float ydiv = 1.4;
	float selmass[selnmass+1] = {220,270,320,370};
*/
} 

void Anal_L5JERC_wt::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   const char* JetEta[netarange] = {"|eta|<0.5","0.5<|eta|<1.0","1.0<|eta|<1.5","1.5<|eta|<2.0","2.0<|eta|<2.5","2.5<|eta|<3.0","3.0<|eta|<3.2","3.2<|eta|<3.7","3.7<|eta|<4.2","4.2<|eta|<4.7"} ;
   float etarng[netarange+1] ={0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 3.7, 4.2, 4.7};

   float ptbins[noptbins+1] = {30, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
     4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000} ;

  float etaea[netaea+1] = {0,1.0,1.479,2.0,2.2,2.3,2.4,5};
  float EA_chhad[netaea] = {0.057,0.052,0.045,0.044,0.043,0.038,0.034};
  float EA_neuhad[netaea] = {0.062,0.112,0.077,0.033,0.007,0.007,0.012};
  float EA_pho[netaea] = {0.122,0.110,0.061,0.077,0.099,0.109,0.134};
  float EA_el[netaea] = {0.1566,0.1626,0.1073,0.0854,0.1051,0.1204,0.1524};
/*
  float min_wmass = 0; 
  float max_wmass = 200;
  int wmassbins = 100;
  float min_tmass = 0;
  float max_tmass = 500.;
  int tmassbins = 250;

  float ptcut = 25.;  
  float ycut = 2.4;
  float bream_rad = 0.1;
  float btagvalue = 0.8385;
  float weight;

  float ydiv = 1.4;
//float seleta[seleta+1] = {1.4,2.5,5.};
  float selmass[selnmass+1] = {220,270,320,370};
*/ 
    char infile[100];
	char outfile[100];
	char outfilx[100];
	sprintf(infile,"try1.log");
	int len = strlen(infile);
    strncpy(outfilx, infile, len-4);
    outfilx[len-4]='\0';
    sprintf (outfile,"%s.root",outfilx);

	OutFile = new TProofOutputFile(outfile);
 // fOutput->Add(OutFile);
  
	fileOut = OutFile->OpenFile("RECREATE");

  if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",
              OutFile->GetDir(), OutFile->GetFileName());
    }
   
   char text_file[100];
   sprintf(text_file,"%s.txt",outfilx);
 //  fp.open(text_file) ; 
 //  OutTxt = new TProofOutputFile(text_file);
 //  fp = OutTxt->OpenFile("");
    
  isMC = true;
  
   sprintf(name,"h_Nvert");
   h_nvert = new TH1D(name,name,100,-0.1,99.9);
   h_nvert->Sumw2();
   
   sprintf(name,"h_Ndofct");
   h_ndofct = new TH1D(name,name,100,-0.1,99.9);
   h_ndofct->Sumw2();
    
   sprintf(name,"h_Nchict");
   h_nchict = new TH1D(name,name,100,-0.1,99.9);
   h_nchict->Sumw2(); 
   
   sprintf(name,"h_NPUVert");
   h_npuvert = new TH1D(name,name,100,-0.1,99.9);
   h_npuvert->Sumw2(); 
	   
}

Bool_t Anal_L5JERC_wt::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Anal_L5JERC_wt::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   
   GetEntry(entry);
   
   nevent_total += 1;
   
  // tot_weight += event_weight;
   
   if(isMC){
   weight = event_weight;
   }else{
	   weight = 1;
	   }


h_nvert->Fill(nvert,weight);
h_ndofct->Fill(ndofct,weight);
h_nchict->Fill(nchict,weight);
if(isMC){
	h_npuvert->Fill(npu_vert,weight);
	}
  
return kTRUE;
}

void Anal_L5JERC_wt::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
fileOut->cd();
fileOut->Write();

fOutput->Add(OutFile); 

fileOut->Close();

}

void Anal_L5JERC_wt::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
//std::cout<<"Total Number of events in "<<fileOut->GetName()<<" is "<<nevent_total<<endl;
}
