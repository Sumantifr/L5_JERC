#define Anal_L5JERC_cxx
// The class definition in Anal_L5JERC.h has been generated automatically
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
// Root > T->Process("Anal_L5JERC.C")
// Root > T->Process("Anal_L5JERC.C","some options")
// Root > T->Process("Anal_L5JERC.C+")
//

#include "Anal_L5JERC.h"
#include <TH2.h>
#include <TStyle.h>
#include "TProofServ.h"

//#define TTBar

void Anal_L5JERC::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
 /*
    float ptcut = 25.;
	float ycut = 2.4;
	float brem_rad = 0.1;
	float btagvalue = 0.8385;
	float weight;
	
	float ydiv = 1.4;
	float selmass[selnmass+1] = {220,270,320,370};
*/
} 

void Anal_L5JERC::SlaveBegin(TTree * /*tree*/)
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
  float brem_rad = 0.1;
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
    
  Tout = new TTree("Tout", "PartonInfo");

  Tout->Branch("nevent_total", &nevent_total, "nevent_total/I");
  Tout->Branch("weightev", &weightev, "weightev/F");

  Tout->Branch("pq1_pt", &pq1_pt, "pq1_pt/F");  
  Tout->Branch("pq1_eta", &pq1_eta, "pq1_eta/F");  
  Tout->Branch("pq1_phi", &pq1_phi, "pq1_phi/F");  
  Tout->Branch("pq1_e", &pq1_e, "pq1_e/F");    

  Tout->Branch("pq2_pt", &pq2_pt, "pq2_pt/F");  
  Tout->Branch("pq2_eta", &pq2_eta, "pq2_eta/F");  
  Tout->Branch("pq2_phi", &pq2_phi, "pq2_phi/F");  
  Tout->Branch("pq2_e", &pq2_e, "pq2_e/F");   

  Tout->Branch("whad_pt", &whad_pt, "whad_pt/F");  
  Tout->Branch("whad_eta", &whad_eta, "whad_eta/F");  
  Tout->Branch("whad_phi", &whad_phi, "whad_phi/F");  
  Tout->Branch("whad_e", &whad_e, "whad_e/F");   
  Tout->Branch("whad_m", &whad_m, "whad_m/F");   
  
  Tout->Branch("bhad_pt", &bhad_pt, "bhad_pt/F");  
  Tout->Branch("bhad_eta", &bhad_eta, "bhad_eta/F");  
  Tout->Branch("bhad_phi", &bhad_phi, "bhad_phi/F");  
  Tout->Branch("bhad_e", &bhad_e, "bhad_e/F");   
  Tout->Branch("bhad_m", &bhad_m, "bhad_m/F");

  Tout->Branch("plep_pt", &plep_pt, "plep_pt/F");  
  Tout->Branch("plep_eta", &plep_eta, "plep_eta/F");  
  Tout->Branch("plep_phi", &plep_phi, "plep_phi/F");  
  Tout->Branch("plep_e", &plep_e, "plep_e/F");

  Tout->Branch("blep_pt", &blep_pt, "blep_pt/F");  
  Tout->Branch("blep_eta", &blep_eta, "blep_eta/F");  
  Tout->Branch("blep_phi", &blep_phi, "blep_phi/F");  
  Tout->Branch("blep_e", &blep_e, "blep_e/F");   
  Tout->Branch("blep_m", &blep_m, "blep_m/F"); 

  Tout->Branch("thad_pt", &thad_pt, "thad_pt/F");  
  Tout->Branch("thad_eta", &thad_eta, "thad_eta/F");  
  Tout->Branch("thad_phi", &thad_phi, "thad_phi/F");  
  Tout->Branch("thad_e", &thad_e, "thad_e/F");   
  Tout->Branch("thad_m", &thad_m, "thad_m/F");   
 
  Tout->Branch("tlep_pt", &tlep_pt, "tlep_pt/F");  
  Tout->Branch("tlep_eta", &tlep_eta, "tlep_eta/F");  
  Tout->Branch("tlep_phi", &tlep_phi, "tlep_phi/F");  
  Tout->Branch("tlep_e", &tlep_e, "tlep_e/F");   
  Tout->Branch("tlep_m", &tlep_m, "tlep_m/F");  

  Tout->Branch("delR_whadbhad", &delR_whadbhad, "delR_whadbhad/F"); 
  Tout->Branch("delPhi_whadbhad", &delPhi_whadbhad, "delPhi_whadbhad/F");  
  Tout->Branch("delEta_whadbhad", &delEta_whadbhad, "delEta_whadbhad/F"); 

  Tout->Branch("delR_whadblep", &delR_whadblep, "delR_whadblep/F"); 
  Tout->Branch("delPhi_whadblep", &delPhi_whadblep, "delPhi_whadblep/F");  
  Tout->Branch("delEta_whadblep", &delEta_whadblep, "delEta_whadblep/F");
  
  Tout->Branch("delR_whadplep", &delR_whadplep, "delR_whadplep/F"); 
  Tout->Branch("delPhi_whadplep", &delPhi_whadplep, "delPhi_whadplep/F");  
  Tout->Branch("delEta_whadplep", &delEta_whadplep, "delEta_whadplep/F");
  
  Tout->Branch("delR_bhadplep", &delR_bhadplep, "delR_bhadplep/F");
  Tout->Branch("delPhi_bhadplep", &delPhi_bhadplep, "delPhi_bhadplep/F");
  Tout->Branch("delEta_bhadplep", &delEta_bhadplep, "delEta_bhadplep/F");
  
  Tout->Branch("delR_blepplep", &delR_blepplep, "delR_blepplep/F");
  Tout->Branch("delPhi_blepplep", &delPhi_blepplep, "delPhi_blepplep/F");
  Tout->Branch("delEta_blepplep", &delEta_blepplep, "delEta_blepplep/F");
    
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
    
   sprintf(name,"NMuon_Pre"); 
   h_nmuon_pre = new TH1D(name,name,10,-0.1,9.9);
   h_nmuon_pre->Sumw2();
   
   sprintf(name,"Muon_TightID_Pre"); 
   h_idmuon_pre = new TH1D(name,name,2,-0.1,1.9); 
   h_idmuon_pre->Sumw2();
   
   sprintf(name,"NElectron_Pre");
   h_nelec_pre = new TH1D(name,name,10,-0.1,9.9);
   h_nelec_pre->Sumw2();
   
   sprintf(name,"Electron_TightID_Pre"); 
   h_idelec_pre = new TH1D(name,name,2,-0.1,1.9); 
   h_idelec_pre->Sumw2();
   
   sprintf(name,"NPhoton_Pre");
   h_npho_pre = new TH1D(name,name,10,-0.1,9.9);
   h_npho_pre->Sumw2();
   
   sprintf(name,"NAK4Jets_Pre"); 
   h_npfjetAK4_pre = new TH1D(name,name,10,-0.1,9.9);
   h_npfjetAK4_pre->Sumw2();
   
   sprintf(name,"NAK4BJets_Pre"); 
   h_nbjetAK4_pre = new TH1D(name,name,10,-0.1,9.9);
   h_nbjetAK4_pre->Sumw2();
   
   for(int jk=0;jk<netarange; jk++){
	  sprintf(name,"JetpT_AK8_EtaBin%i",jk+1);
	  h_jetptAK8[jk] = new TH1D(name,name,noptbins,ptbins);
	  h_jetptAK8[jk]->Sumw2();
	  
	  sprintf(name,"JetMass_AK8_EtaBin%i",jk+1);
	  h_jetmassAK8[jk] = new TH1D(name,name,150,0,300);
	  h_jetmassAK8[jk]->Sumw2();
	  
	  sprintf(name,"JetSDMass_AK8_EtaBin%i",jk+1);
	  h_jetsdmassAK8[jk] = new TH1D(name,name,150,0,300);
	  h_jetsdmassAK8[jk]->Sumw2();
	  
	  sprintf(name,"JetMass_Selected_AK8_EtaBin%i",jk+1);
	  h_jetmassAK8_sel[jk] = new TH1D(name,name,150,0,300);
	  h_jetmassAK8_sel[jk]->Sumw2();
	  
	  sprintf(name,"JetSDMass_Selected_AK8_EtaBin%i",jk+1);
	  h_jetsdmassAK8_sel[jk] = new TH1D(name,name,150,0,300);
	  h_jetsdmassAK8_sel[jk]->Sumw2();
	  
	  sprintf(name,"JetMass_SelectedPt_AK8_EtaBin%i",jk+1);
	  h_jetmassAK8_selpt[jk] = new TH1D(name,name,150,0,300);
	  h_jetmassAK8_selpt[jk]->Sumw2();
	  
	  sprintf(name,"JetSDMass_SelectedPt_AK8_EtaBin%i",jk+1);
	  h_jetsdmassAK8_selpt[jk] = new TH1D(name,name,150,0,300);
	  h_jetsdmassAK8_selpt[jk]->Sumw2();
	  
	  sprintf(name,"JetBTag_AK8_EtaBin%i",jk+1);
	  h_jetbtagAK8[jk] = new TH1D(name,name,20,0,2.);
	  h_jetbtagAK8[jk]->Sumw2();
	  
	  sprintf(name,"Jettau21_AK8_EtaBin%i",jk+1);
	  h_jettau21AK8[jk] = new TH1D(name,name,50,0.,1.);
	  h_jettau21AK8[jk]->Sumw2();
	  
	  sprintf(name,"Jettau32_AK8_EtaBin%i",jk+1);
	  h_jettau32AK8[jk] = new TH1D(name,name,50,0.,1.);
	  h_jettau32AK8[jk]->Sumw2();
    }
 
    
    sprintf(name,"Jety_AK8");
	h_jetyAK8 = new TH1D(name,name,100,-5.,5.);
	h_jetyAK8->Sumw2();
	
	sprintf(name,"Jetphi_AK8");
	h_jetphiAK8 = new TH1D(name,name,50,-M_PI,M_PI);
    h_jetphiAK8->Sumw2();
    
	for(int jk=0;jk<netarange; jk++){  
	  sprintf(name,"JetpT_AK4_EtaBin%i",jk+1);
	  h_jetptAK4[jk] = new TH1D(name,name,noptbins,ptbins);
	  h_jetptAK4[jk]->Sumw2();
	  
	  sprintf(name,"JetMass_AK4_EtaBin%i",jk+1);
	  h_jetmassAK4[jk] = new TH1D(name,name,noptbins,ptbins);
	  h_jetmassAK4[jk]->Sumw2();
	  
	  sprintf(name,"JetSDMass_AK4_EtaBin%i",jk+1);
	  h_jetsdmassAK4[jk] = new TH1D(name,name,noptbins,ptbins);
	  h_jetsdmassAK4[jk]->Sumw2();
	  
	  sprintf(name,"JetBTag_AK4_EtaBin%i",jk+1);
	  h_jetbtagAK4[jk] = new TH1D(name,name,20,0,2.);
	  h_jetbtagAK4[jk]->Sumw2();
	  
	  sprintf(name,"Jettau21_AK4_EtaBin%i",jk+1);
	  h_jettau21AK4[jk] = new TH1D(name,name,50,0.,1.);
	  h_jettau21AK4[jk]->Sumw2();
	  
	  sprintf(name,"Jettau32_AK4_EtaBin%i",jk+1);
	  h_jettau32AK4[jk] = new TH1D(name,name,50,0.,1.);
	  h_jettau32AK4[jk]->Sumw2();
	  
	  for(int kl=0; kl<noptbins; kl++){
	  sprintf(name,"JetPtReso_AK4_PtBin%i_EtaBin%i",kl+1,jk+1);
	  h_resopt_AK4[jk][kl] = new TH1D(name,name,100,-2.,2);
	  h_resopt_AK4[jk][kl]->Sumw2();
			}
	
	  for(int kl=0; kl<noptbins; kl++){
	  sprintf(name,"DiJetMass_AK4_PtBin%i_EtaBin%i",kl+1,jk+1);
	  h_massdijet_AK4[jk][kl] = new TH1D(name,name,100,-2.,2);
	  h_massdijet_AK4[jk][kl]->Sumw2();
			}
			
	  }

  sprintf(name,"Jety_AK4");
  h_jetyAK4 = new TH1D(name,name,100,-5.,5.);
  h_jetyAK4->Sumw2();
  
  sprintf(name,"Jetphi_AK4");
  h_jetphiAK4 = new TH1D(name,name,50,-M_PI,M_PI);
  h_jetphiAK4->Sumw2();

  sprintf(name,"NJet_AK4");
  h_jetnqbAK4 = new TH1D(name,name,10,-0.1,9.9);
  h_jetnqbAK4->Sumw2();
  
  sprintf(name,"NJet_B_AK4");
  h_jetnbAK4 = new TH1D(name,name,10,-0.1,9.9);
  h_jetnbAK4->Sumw2();
  
  sprintf(name,"NJet_Q_AK4");
  h_jetnqAK4 = new TH1D(name,name,10,-0.1,9.9);
  h_jetnqAK4->Sumw2();
  
  sprintf(name,"NMuon");
  h_nmuon = new TH1D(name,name,10,-0.1,9.9);
  h_nmuon->Sumw2();
  
  sprintf(name,"JetpT_AK4_Q1");
  h_jetpt_q1AK4 = new TH1D(name,name,noptbins,ptbins);
  h_jetpt_q1AK4->Sumw2();
  
  sprintf(name,"JetEta_AK4_Q1");
  h_jeteta_q1AK4 = new TH1D(name,name,100,-5.,5.);
  h_jeteta_q1AK4->Sumw2();
  
  sprintf(name,"JetPhi_AK4_Q1");
  h_jetphi_q1AK4 = new TH1D(name,name,50,-M_PI,M_PI);
  h_jetphi_q1AK4->Sumw2();
  
  sprintf(name,"JetpT_AK4_Q2");
  h_jetpt_q2AK4 = new TH1D(name,name,noptbins,ptbins);
  h_jetpt_q2AK4->Sumw2();
  
  sprintf(name,"JetEta_AK4_Q2");
  h_jeteta_q2AK4 = new TH1D(name,name,100,-5.,5.);
  h_jeteta_q2AK4->Sumw2();
  
  sprintf(name,"JetPhi_AK4_Q2");
  h_jetphi_q2AK4 = new TH1D(name,name,50,-M_PI,M_PI); 
  h_jetphi_q2AK4->Sumw2();
  
  sprintf(name,"JetpT_AK4_B1");
  h_jetpt_b1AK4 = new TH1D(name,name,noptbins,ptbins);
  h_jetpt_b1AK4->Sumw2();
  
  sprintf(name,"JetEta_AK4_B1");
  h_jeteta_b1AK4 = new TH1D(name,name,100,-5.,5.);
  h_jeteta_b1AK4->Sumw2();
  
  sprintf(name,"JetPhi_AK4_B1");
  h_jetphi_b1AK4 = new TH1D(name,name,50,-M_PI,M_PI);
  h_jetphi_b1AK4->Sumw2();
  
  sprintf(name,"JetpT_AK4_B2");
  h_jetpt_b2AK4 = new TH1D(name,name,noptbins,ptbins);
  h_jetpt_b2AK4->Sumw2();
  
  sprintf(name,"JetEta_AK4_B2");
  h_jeteta_b2AK4 = new TH1D(name,name,100,-5.,5.);
  h_jeteta_b2AK4->Sumw2();
  
  sprintf(name,"JetPhi_AK4_B2");
  h_jetphi_b2AK4 = new TH1D(name,name,50,-M_PI,M_PI); 
  h_jetphi_b2AK4->Sumw2();
  
  sprintf(name,"Passed_Muon_pT");
  h_mupt_pass  = new TH1D(name,name,200,0,2000);
  h_mupt_pass->Sumw2();
  
  sprintf(name,"Passed_Muon_Eta");
  h_mueta_pass = new TH1D(name,name,100,-5.,5.); 
  h_mueta_pass->Sumw2();
  
  sprintf(name,"Passed_Muon_Phi");
  h_muphi_pass = new TH1D(name,name,50,-M_PI,M_PI); 
  h_muphi_pass->Sumw2();
  
  sprintf(name,"PFMet_Et");
  h_metpt = new TH1D(name,name,60,0,300);
  h_metpt->Sumw2();
  
  sprintf(name,"PFMet_Phi");
  h_metphi = new TH1D(name,name,50,-M_PI,M_PI);
  h_metphi->Sumw2();
  
  sprintf(name,"PFMet_EtbySumEt");
  h_metbyEt = new TH1D(name,name,50,0,1.);
  h_metbyEt->Sumw2();
  
  sprintf(name,"PFMet_Significance");
  h_metSig = new TH1D(name,name,100,0,500);
  h_metSig->Sumw2();
  
  for(int jk=0;jk<netarange; jk++){

  sprintf(name,"Muon_number_EtaBin%i",jk+1);
  h_mun[jk] = new TH1D(name,name,20,0,20);
  h_mun[jk]->Sumw2();
  
  sprintf(name,"Muon_pT_EtaBin%i",jk+1);
  h_mupt[jk] = new TH1D(name,name,200,0,2000);
  h_mupt[jk]->Sumw2();

  }
  
  sprintf(name,"Muon_eta");
  h_mueta = new TH1D(name,name,100,-5.,5.);
  h_mueta->Sumw2();
  
  sprintf(name,"Muon_phi");
  h_muphi = new TH1D(name,name,50,-M_PI,M_PI);
  h_muphi->Sumw2();
  
  sprintf(name,"Muon_drbm");
  h_mudrbm = new TH1D(name,name,100,-5.,5.);
  h_mudrbm->Sumw2();
  
  sprintf(name,"Muon_trkvtx");
  h_mutrkvtx = new TH1D(name,name,100,-5.,5.);
  h_mutrkvtx->Sumw2();
  
  sprintf(name,"Muon_hbye");
  h_muhbye = new TH1D(name,name,100,0,1.);
  h_muhbye->Sumw2();
  
  sprintf(name,"Muon_emiso");
  h_muemiso = new TH1D(name,name,40,0,400);
  h_muemiso->Sumw2();
  
  sprintf(name,"Muon_hadiso");
  h_muhadiso = new TH1D(name,name,40,0,400);
  h_muhadiso->Sumw2();
  
  sprintf(name,"Muon_pfiso");
  h_mupfiso = new TH1D(name,name,20,0,200);
  h_mupfiso->Sumw2();
  
  sprintf(name,"Muon_kpt03");
  h_mukpt03 = new TH1D(name,name,20,0,2000);
  h_mukpt03->Sumw2();
  
  sprintf(name,"Muon_kpt05");
  h_mukpt05 = new TH1D(name,name,20,0,2000);
  h_mukpt05->Sumw2();
  
  sprintf(name,"Muon_hit");
  h_muhit = new TH1D(name,name,50,0,50);
  h_muhit->Sumw2();
  
  sprintf(name,"Muon_mst");
  h_mumst = new TH1D(name,name,10,0,10);
  h_mumst->Sumw2();
  
  sprintf(name,"Muon_trklay");
  h_mutrklay = new TH1D(name,name,20,0,20);
  h_mutrklay->Sumw2();
  
  sprintf(name,"WHad_Mass");
  h_whadmass = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmass->Sumw2();
  
  sprintf(name,"WHad_pt");
  h_whadpt = new TH1D(name,name,100,0,500);
  h_whadpt->Sumw2();
  
  sprintf(name,"WHad_y");
  h_whady = new TH1D(name,name,50,-5.,5.);
  h_whady->Sumw2();
  
  sprintf(name,"THad_Mass");
  h_thadmass = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmass->Sumw2();
  
  sprintf(name,"THad_pt");
  h_thadpt = new TH1D(name,name,100,0,500);
  h_thadpt->Sumw2();
  
  sprintf(name,"THad_y");
  h_thady = new TH1D(name,name,50,-5.,5.);
  h_thady->Sumw2();
  
  sprintf(name,"TW_Combined_2D_Mass");
  h_tw_2d = new TH2D(name,name,tmassbins,min_tmass,max_tmass,wmassbins,min_wmass,max_wmass);
  h_tw_2d->Sumw2();
  
  sprintf(name,"WHad_MScheme_Mass");
  h_whadmassA = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmassA->Sumw2();
  
  sprintf(name,"WHad_MScheme_pt");
  h_whadptA = new TH1D(name,name,500,0,500);
  h_whadptA->Sumw2();
  
  sprintf(name,"WHad_MScheme_y");
  h_whadyA = new TH1D(name,name,100,-5.,5.);
  h_whadyA->Sumw2();
  
  sprintf(name,"THad_MScheme_Mass");
  h_thadmassA = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmassA->Sumw2();
  
  sprintf(name,"THad_MScheme_pt");
  h_thadptA = new TH1D(name,name,100,0,500);
  h_thadptA->Sumw2();
  
  sprintf(name,"THad_MScheme_y");
  h_thadyA = new TH1D(name,name,50,-5.,5.);
  h_thadyA->Sumw2();
  
  sprintf(name,"TW_Combined_2D_MScheme_Mass");
  h_tw_2dA = new TH2D(name,name,tmassbins,min_tmass,max_tmass,wmassbins,min_wmass,max_wmass);
  h_tw_2dA->Sumw2();
  
  sprintf(name,"WHad_Truth_Mass");
  h_whadmass_truth = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmass_truth->Sumw2();
  
  sprintf(name,"WHad_Truth_pt");
  h_whadpt_truth = new TH1D(name,name,200,0,500);
  h_whadpt_truth->Sumw2();
  
  sprintf(name,"WHad_Truth_y");
  h_whady_truth = new TH1D(name,name,100,-5.,5.);
  h_whady_truth->Sumw2();
  
  sprintf(name,"THad_Truth_Mass");
  h_thadmass_truth = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmass_truth->Sumw2();
  
  sprintf(name,"THad_Truth_pt");
  h_thadpt_truth = new TH1D(name,name,200,0,500);
  h_thadpt_truth->Sumw2();
  
  sprintf(name,"THad_Truth_y");
  h_thady_truth = new TH1D(name,name,100,-5.,5.);
  h_thady_truth->Sumw2();
  
  sprintf(name,"TW_Combined_2D_Truth_Mass");
  h_tw_2d_truth = new TH2D(name,name,tmassbins,min_tmass,max_tmass,wmassbins,min_wmass,max_wmass);
  h_tw_2d_truth->Sumw2();
  
  sprintf(name,"WHad_Wrong_Mass");
  h_whadmass_wrong = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmass_wrong->Sumw2();
  
  sprintf(name,"WHad_Wrong_pt");
  h_whadpt_wrong = new TH1D(name,name,200,0,500);
  h_whadpt_wrong->Sumw2();
  
  sprintf(name,"WHad_Wrong_y");
  h_whady_wrong = new TH1D(name,name,100,-5.,5.);
  h_whady_wrong->Sumw2();
  
  sprintf(name,"THad_Wrong_Mass");
  h_thadmass_wrong = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmass_wrong->Sumw2();
  
  sprintf(name,"THad_Wrong_pt");
  h_thadpt_wrong = new TH1D(name,name,200,0,500);
  h_thadpt_wrong->Sumw2();
  
  sprintf(name,"THad_Wrong_y");
  h_thady_wrong = new TH1D(name,name,100,-5.,5.);
  h_thady_wrong->Sumw2();
  
  sprintf(name,"TW_Combined_2D_Wrong_Mass");
  h_tw_2d_wrong = new TH2D(name,name,tmassbins,min_tmass,max_tmass,wmassbins,min_wmass,max_wmass);
  h_tw_2d_wrong->Sumw2();
  
  sprintf(name,"WHad_UnMatched_Mass");
  h_whadmass_unmatch = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmass_unmatch->Sumw2();
  
  sprintf(name,"WHad_UnMatched_pt");
  h_whadpt_unmatch = new TH1D(name,name,200,0,500);
  h_whadpt_unmatch->Sumw2();
  
  sprintf(name,"WHad_UnMatched_y");
  h_whady_unmatch = new TH1D(name,name,100,-5.,5.);
  h_whady_unmatch->Sumw2();
  
  sprintf(name,"THad_UnMatched_Mass");
  h_thadmass_unmatch = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmass_unmatch->Sumw2();
  
  sprintf(name,"THad_UnMatched_pt");
  h_thadpt_unmatch = new TH1D(name,name,200,0,500);
  h_thadpt_unmatch->Sumw2();
  
  sprintf(name,"THad_UnMatched_y");
  h_thady_unmatch = new TH1D(name,name,100,-5.,5.);
  h_thady_unmatch->Sumw2();
  
  sprintf(name,"TW_Combined_2D_UnMatched_Mass");
  h_tw_2d_unmatch = new TH2D(name,name,tmassbins,min_tmass,max_tmass,wmassbins,min_wmass,max_wmass);
  h_tw_2d_unmatch->Sumw2();
  
  sprintf(name,"WHad_Parton_Mass");
  h_whadmass_parton = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
  h_whadmass_parton->Sumw2();
  
  sprintf(name,"THad_Parton_Mass");
  h_thadmass_parton = new TH1D(name,name,tmassbins,min_tmass,max_tmass);
  h_thadmass_parton->Sumw2();
  
  sprintf(name,"Likelihood_lambda");
  h_lambda = new TH1D(name,name,30,0,3);
  h_lambda->Sumw2();
  
  sprintf(name,"Likelihood_MScheme_lambda");
  h_lambdaA = new TH1D(name,name,30,0,3);
  h_lambdaA->Sumw2();
  
  sprintf(name,"Likelihood_truth_lambda");
  h_lambda_truth = new TH1D(name,name,30,0,3);
  h_lambda_truth->Sumw2();
  
  for(int eta=0; eta<selneta; eta++){
	  for(int ms=0; ms<selnmass; ms++){
		  for(int jec=0; jec<11; jec++){
			  sprintf(name,"Sel_WHad_Mass_JEC%i_MassBin%i_Eta%i",jec+1,ms+1,eta+1);
			  h_whadmass_sel[jec][ms][eta] = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
			  h_whadmass_sel[jec][ms][eta]->Sumw2();
			  
			  sprintf(name,"Sel_WHad_Truth_Mass_JEC%i_MassBin%i_Eta%i",jec+1,ms+1,eta+1);
			  h_whadmass_truth_sel[jec][ms][eta] = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
			  h_whadmass_truth_sel[jec][ms][eta]->Sumw2();
			  
			  sprintf(name,"Sel_WHad_Wrong_Mass_JEC%i_MassBin%i_Eta%i",jec+1,ms+1,eta+1);
			  h_whadmass_wrong_sel[jec][ms][eta] = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
			  h_whadmass_wrong_sel[jec][ms][eta]->Sumw2();
			  
			  sprintf(name,"Sel_WHad_UnMatched_Mass_JEC%i_MassBin%i_Eta%i",jec+1,ms+1,eta+1);
			  h_whadmass_unmatch_sel[jec][ms][eta] = new TH1D(name,name,wmassbins,min_wmass,max_wmass);
			  h_whadmass_unmatch_sel[jec][ms][eta]->Sumw2();
			  }
		  }
		  for(int ic=0; ic<2; ic++){
			sprintf(name,"Sel_Jety_%i_Eta%i",ic+1,eta+1);
			h_jety_2y[ic][eta] = new TH1D(name,name,100,-5.,5.);
			h_jety_2y[ic][eta]->Sumw2();
			
			sprintf(name,"Sel_JetPt_%i_Eta%i",ic+1,eta+1);
			h_jety_2pt[ic][eta] = new TH1D(name,name,200,0,2000);
			h_jety_2pt[ic][eta]->Sumw2();
		}
	  }
	   
}

Bool_t Anal_L5JERC::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Anal_L5JERC::GetEntry() or TBranch::GetEntry()
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
   
   tot_weight += event_weight;

	float etaea[netaea+1] = {0,1.0,1.479,2.0,2.2,2.3,2.4,5};
    float etarng[netarange+1] ={0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 3.7, 4.2, 4.7};
    float pbin_values[nopbins+1] = {100,300,500,10000};
    float mbin_values[nomassbins+1] = {0,140,200,3000};
    float ptbins[noptbins+1] = {30, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
     4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000} ;
  
  // float ydivs[selneta+1] = {0,0.5,1.0,1.5,2.0,2.5};
  
   float PV_weight[npvmax] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   
   float PV_weight_TT[npvmax] = {0.000220901,0.00693389,0.00745878,0.0475494,0.0773146,0.078206,0.0884815,0.0912987,0.0909217,0.103422,0.159877,0.316299,0.603483,0.986689,1.39064,1.72398,1.9567,2.08403,2.1358,2.13861,2.10244,2.04811,1.98871,1.95368,1.94097,1.95738,1.98176,1.98722,1.96425,1.90108,1.81445,1.69314,1.54288,1.37509,1.19599,1.01671,0.843972,0.68202,0.539962,0.417882,0.316269,0.234112,0.170225,0.121962,0.0855003,0.0589613,0.0401117,0.0269265,0.0178288,0.0116006,0.00744245,0.00471964,0.00294265,0.00180608,0.001089,0.000651485,0.000380997,0.000221757,0.000126294,7.14882e-05,3.99083e-05,2.20173e-05,1.19666e-05,6.49217e-06,3.45751e-06,1.83111e-06,9.56878e-07,4.98387e-07,2.55431e-07,1.31408e-07,6.61132e-08,3.30655e-08,1.63566e-08,8.00272e-09,3.86998e-09,1.84201e-09,8.6465e-10,4.0105e-10,1.83765e-10,8.08219e-11,3.55241e-11,1.52755e-11,6.40959e-12,2.64194e-12,1.06052e-12,4.17719e-13,1.57024e-13,5.7781e-14,2.1219e-14,7.52852e-15,2.59023e-15,8.91939e-16,3.03539e-16,9.64909e-17,1.69799e-17,0,0,0,0,0};
   
   float PV_weight_Zj[npvmax] = {0.0162932,0.00304803,0.0456651,0.116669,0.23593,0.415858,0.620988,0.819695,0.991391,1.1311,1.25117,1.34421,1.41921,1.47624,1.51854,1.5446,1.55676,1.55815,1.54714,1.52195,1.48079,1.43226,1.36623,1.29354,1.2103,1.12289,1.02907,0.933732,0.833599,0.736729,0.642214,0.552249,0.468871,0.393294,0.325187,0.26637,0.216105,0.173479,0.138931,0.109994,0.0876801,0.0695027,0.0545992,0.0437794,0.0346587,0.0273779,0.0217596,0.0175424,0.0142657,0.01143,0.00938838,0.00762231,0.00614745,0.00529929,0.00459658,0.00387351,0.00301634,0.0026131,0.00244303,0.00189471,0.00157309,0.00153251,0.00144049,0.00107836,0.0012727,0.000903082,0.00131467,0.000908631,0.00103471,0.000431562,0.000823078,0.000872884,0.000423908,0.000535174,0.00100365,0.00041925,0,0.000683158,0,0.00188225,0,0.00158966,0.00182679,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
   float PV_weight_Wj[npvmax] = {11.1533,0.00307276,0.0328439,0.115383,0.24584,0.435885,0.64987,0.856052,1.03057,1.16955,1.28619,1.37502,1.44646,1.49841,1.53449,1.56211,1.56804,1.56686,1.55455,1.52571,1.48458,1.43221,1.36504,1.29296,1.20965,1.12083,1.02561,0.930312,0.830668,0.73406,0.638972,0.55012,0.467473,0.391234,0.32338,0.265219,0.214419,0.172229,0.137771,0.108571,0.0860888,0.068222,0.0534077,0.0426865,0.0336501,0.0265847,0.0209903,0.0169463,0.0137073,0.0109752,0.00892725,0.00725307,0.00582303,0.00502794,0.00435825,0.00366336,0.00282718,0.00243509,0.00225749,0.00174102,0.00143779,0.00142844,0.00131558,0.000973996,0.00116354,0.000812439,0.00117203,0.000815989,0.000949941,0.000376483,0.000714085,0.00076455,0.000357915,0.000466578,0.000867344,0.000357342,0,0.000533584,0,0.00176488,0,0.00127968,0.00139087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   if(isMC){
   weight = event_weight;
   }else{
	   weight = 1;
	   }

weightev = weight;

if(isMC){
if(npu_vert>=0 && npu_vert<100){
	weight *= PV_weight[npu_vert];
    #ifdef TTBar
    weight *= PV_weight_TT[npu_vert];
    #endif
    #ifdef DYJets
    weight *= PV_weight_Zj[npu_vert];
    #endif
    #ifdef WJets
    weight *= PV_weight_Wj[npu_vert];
    #endif
	}
}

Tout->Fill();

h_nvert->Fill(nvert,weight);
h_ndofct->Fill(ndofct,weight);
h_nchict->Fill(nchict,weight);
if(isMC){
	h_npuvert->Fill(npu_vert,weight);
	}

for(int ht=0; ht<nHLTmx; ht++){

   ihlt[ht] = 0;
   prescl[ht] = 0;

   switch(ht){

   case 0:
   ihlt[ht] = ihlt01;
   prescl[ht] = prescl01;
   break;

   case 1:
   ihlt[ht] = ihlt02;
   prescl[ht] = prescl02;
   break;

   case 2:
   ihlt[ht] = ihlt03;
   prescl[ht] = prescl03;
   break;
 }
}//ht

int nmuon1 = 0;

for(int muon=0; muon<nmuons; muon++){

bool muon_id = Muon_TightID(muonisGL[muon],muonisPF[muon],muonchi[muon],muonhit[muon],muonmst[muon], muontrkvtx[muon], muondz[muon], muonpixhit[muon], muontrklay[muon]);
muonistight[nmuon1] = muon_id;
if(!muon_id) continue;
//if(!muonisMed[muon]) continue;
bool muon_iso_id = Muon_Iso_ID(muonpfiso[muon]);
if(!muon_iso_id) continue;
/*
if(!muonisLoose[muon]) continue;
*/
if(fabs(muonpt[0]) < ptcut) break;
if(fabs(muoneta[0]) > ycut) break;

muonpt[nmuon1] = muonpt[muon];
if(muonpt[nmuon1]<0) { muonq[nmuon1] = -1; }
else{
	if(muonpt[nmuon1]>0) { muonq[nmuon1] = +1; }
	else {muonq[nmuon1] = 0; }
	}
muonpt[nmuon1] = fabs(muonpt[nmuon1]);	
muone[nmuon1] = muone[muon];
muonp[nmuon1] = muonp[muon];
muoneta[nmuon1] = muoneta[muon];
muonphi[nmuon1] = muonphi[muon];

//if(isMC) { muonpt[nmuon1] *= Roch_Cor(muonq[nmuon1],muonpt[nmuon1],muoneta[nmuon1],muonphi[nmuon1]); }

nmuon1++;
if(nmuon1>=njetmx) break;
 }
 
nmuons = nmuon1;
 
int nelec1=0;

for(int elec=0; elec<nelecs; elec++){

bool elec_id = ELectron_ID(ele[elec],eleta[elec],elp[elec],elietaieta[elec],
						   eletain[elec],elphiin[elec],elmisshits[elec],elhovere[elec],elconvdist[elec]
						   ,Rho);
elecistight[elec] = elec_id;
if(!elec_id) continue;

if(fabs(elpt[0]) < ptcut) break;
if(abs(eleta[0])>ycut) break;

int eleaetatag = getbinid(abs(eleta[elec]),netaea,etaea);

bool elec_iso_id = Electron_Iso_ID(elsceta[elec],elpfiso[elec],Rho,eleaetatag);
if(!elec_iso_id) continue;

elpt[nelec1] = elpt[elec];
eleta[nelec1] = eleta[elec];
elphi[nelec1] = elphi[elec];
ele[nelec1] = ele[elec];
}  

nelecs = nelec1;

int npho1=0;
int nbrempho1=0;

for(int pho=0; pho<nphotons; pho++){

int phoeaetatag = getbinid(abs(phoeta[pho]),netaea,etaea);
bool pho_id = Photon_ID(phoe[pho],phoeta[pho],phohadbyem[pho],phoietaieta[pho],phochhadiso[pho],phoneuhadiso[pho],phophoiso[pho],Rho,phoeaetatag);
bool pho_iso_id = Photon_Iso_ID(phoe[pho],phoeta[pho],phohadiso[pho],phophoiso[pho]);

 float delphomu = brem_rad;
 float delphoe = brem_rad;
 int match_mu = -1;
 int match_e = -1;	
	
 for(int muon=0; muon<nmuons; muon++){
		if(delta2R(muoneta[muon],muonphi[muon],phoeta[pho],phophi[pho])< delphomu){
				delphomu = delta2R(muoneta[muon],muonphi[muon],phoeta[pho],phophi[pho]) ;
				match_mu = muon;
		}
	}
 for(int elec=0; elec<nelecs; elec++){
		if(delta2R(eleta[elec],elphi[elec],phoeta[pho],phophi[pho])< delphoe){
				delphoe = delta2R(eleta[elec],elphi[elec],phoeta[pho],phophi[pho]);
				match_e = elec;
		}
	 }
	 
	 if((match_mu<0)&&(match_e<0)){
		 
		 if(pho_id && pho_iso_id){
			phoe[npho1] = phoe[pho];
			phoeta[npho1] = phoeta[pho];
			phophi[npho1] = phophi[pho];
			npho1++;
			}
			else{continue;}
		 }
	 
	 if(((match_mu<0)&&(match_e>=0))||((match_mu>=0)&&(match_e>=0)&&(delphoe<delphomu)) ) { 
		 
			TLorentzVector tmp_e;
            tmp_e.SetPtEtaPhiE(elpt[match_e],eleta[match_e],elphi[match_e],ele[match_e]);
		    TLorentzVector tmp_pho;
		    tmp_pho.SetPtEtaPhiE(phoe[pho]*fabs(sin(eta_to_theta(phoeta[pho]))),phoeta[pho],phophi[pho],phoe[pho]);
		    tmp_e += tmp_pho;
			
			elpt[match_e] = tmp_e.Pt();
			eleta[match_e] = tmp_e.Eta();
			elphi[match_e] = tmp_e.Phi();
			ele[match_e] = tmp_e.E();
		  
		 }
		 
	 if(((match_mu>=0)&&(match_e<0))||((match_mu>=0)&&(match_e>=0)&&(delphomu<delphoe)) ) {
		 
			TLorentzVector tmp_mu;
            tmp_mu.SetPtEtaPhiE(muonpt[match_mu],muoneta[match_mu],muonphi[match_mu],muone[match_mu]);
		    TLorentzVector tmp_pho;
		    tmp_pho.SetPtEtaPhiE(phoe[pho]*fabs(sin(eta_to_theta(phoeta[pho]))),phoeta[pho],phophi[pho],phoe[pho]);
		    tmp_mu += tmp_pho;
			
			muonpt[match_mu] = tmp_mu.Pt();
			muoneta[match_mu] = tmp_mu.Eta();
			muonphi[match_mu] = tmp_mu.Phi();
			muone[match_mu] = tmp_mu.E();
		 
		 }	 
}

nphotons = npho1;

int npfjets1 = 0;

for(int jet=0; jet<npfjetAK8; jet++){

if(isMC && (pfjetAK8pt[0]> 3*qscale)) { break; }

if(pfjetAK8tightID[jet]!=1) continue;

bool clean_jet = true;

for(int muon=0; muon<nmuons; muon++){
	if(delta2R(pfjetAK8y[jet],pfjetAK8phi[jet],muoneta[muon],muonphi[muon])<0.8){clean_jet = false; break;}
	}
for(int elec=0; elec<nelecs; elec++){
	if(delta2R(pfjetAK8y[jet],pfjetAK8phi[jet],eleta[elec],elphi[elec])<0.8){clean_jet = false; break;}
	}
for(int pho=0; pho<nphotons; pho++){
	if(delta2R(pfjetAK8y[jet],pfjetAK8phi[jet],phoeta[pho],phophi[pho])<0.8){clean_jet = false; break;}
	}

if(!clean_jet) continue;

if(isMC){
pfjetAK8pt[jet] *= (1.+pfjetAK8reso[jet]) ;
}
pfjetAK8pt[jet] *= pfjetAK8JEC[jet];

if(pfjetAK8pt[jet]< ptcut) continue;
if(abs(pfjetAK8y[jet])> ycut) continue;

pfjetAK8pt[npfjets1] = pfjetAK8pt[jet] ;// * pfjetAK8JEC[jet];
pfjetAK8y[npfjets1] = pfjetAK8y[jet];
pfjetAK8phi[npfjets1] = pfjetAK8phi[jet];
pfjetAK8mass[npfjets1] = pfjetAK8mass[jet];
pfjetAK8btag[npfjets1] = pfjetAK8btag[jet];

pfjetAK8moment1[npfjets1] = pfjetAK8moment1[jet];
pfjetAK8chrad[npfjets1] = pfjetAK8chrad[jet];
pfjetAK8axis2[npfjets1] = pfjetAK8axis2[jet];
pfjetAK8pTD[npfjets1] = pfjetAK8pTD[jet];
pfjetAK8sdmass[npfjets1] = pfjetAK8sdmass[jet];
pfjetAK8tau1[npfjets1] = pfjetAK8tau1[jet];
pfjetAK8tau2[npfjets1] = pfjetAK8tau2[jet];
pfjetAK8tau3[npfjets1] = pfjetAK8tau3[jet];
pfjetAK8sub1pt[npfjets1] = pfjetAK8sub1pt[jet];
pfjetAK8sub1btag[npfjets1] = pfjetAK8sub1btag[jet];
pfjetAK8sub2pt[npfjets1] = pfjetAK8sub2pt[jet];
pfjetAK8sub2btag[npfjets1] = pfjetAK8sub2btag[jet];

npfjets1++;
 }
npfjetAK8 =  npfjets1;
 
npfjets1 = 0;
nbjetAK4 = 0;
nqjetAK4 = 0;

for(int jet=0; jet<npfjetAK4; jet++){
	pfjetAK4pt[jet] = (pfjetAK4pt[jet] * pfjetAK4JEC[jet]);
	pfjetAK4mass[jet] = pfjetAK4mass[jet] * pfjetAK4JEC[jet];
	if(isMC){
	pfjetAK4pt[jet] = pfjetAK4pt[jet] *(1+pfjetAK4reso[jet])  ;
	}
}

for(int jet=0; jet<npfjetAK4; jet++){
  for(int get=(jet+1); get<npfjetAK4; get++){
	 if(pfjetAK4pt[get]>pfjetAK4pt[jet]){
		 float tmppt = pfjetAK4pt[jet];
		 pfjetAK4pt[jet] = pfjetAK4pt[get];
		 pfjetAK4pt[get] = tmppt;
		 float tmpy = pfjetAK4eta[jet];
		 pfjetAK4eta[jet] = pfjetAK4eta[get];
		 pfjetAK4eta[get] = tmpy;
		 float tmpphi = pfjetAK4phi[jet];
		 pfjetAK4phi[jet] = pfjetAK4phi[get];
		 pfjetAK4phi[get] = tmpphi;
		 float tmpmass = pfjetAK4mass[jet];
		 pfjetAK4mass[jet] = pfjetAK4mass[get];
		 pfjetAK4mass[get] = tmpmass;
		 float tmpsdmass = pfjetAK4sdmass[jet];
		 pfjetAK4sdmass[jet] = pfjetAK4sdmass[get];
		 pfjetAK4sdmass[get] = tmpsdmass;
		 float tmpbtag = pfjetAK4btag[jet];
		 pfjetAK4btag[jet] = pfjetAK4btag[get];
		 pfjetAK4btag[get] = tmpbtag;
		 float tmppfjetAK4reso = pfjetAK4reso[jet];
		 pfjetAK4reso[jet] = pfjetAK4reso[get];
		 pfjetAK4reso[get] = tmppfjetAK4reso;
		 float tmppfjetAK4hadronflav = pfjetAK4hadronflav[jet];
		 pfjetAK4hadronflav[jet] = pfjetAK4hadronflav[get];
		 pfjetAK4hadronflav[get] = tmppfjetAK4hadronflav;
		 float tmppfjetAK4partonflav = pfjetAK4partonflav[jet];
		 pfjetAK4partonflav[jet] = pfjetAK4partonflav[get];
		 pfjetAK4partonflav[get] = tmppfjetAK4partonflav;
		 float tmppfjetAK4qgl = pfjetAK4qgl[jet];
		 pfjetAK4qgl[jet] = pfjetAK4qgl[get];
		 pfjetAK4qgl[get] = tmppfjetAK4qgl;
		 }
	 }
  }

for(int jet=0; jet<npfjetAK4; jet++){

if(pfjetAK4tightID[jet]!=1) continue;

bool clean_jet = true;

for(int muon=0; muon<nmuons; muon++){
	if(delta2R(pfjetAK4eta[jet],pfjetAK4phi[jet],muoneta[muon],muonphi[muon])<0.4){clean_jet = false; break;}
	}
	
for(int elec=0; elec<nelecs; elec++){
	if(delta2R(pfjetAK4eta[jet],pfjetAK4phi[jet],eleta[elec],elphi[elec])<0.4){clean_jet = false; break;}
	}

for(int pho=0; pho<nphotons; pho++){
	if(delta2R(pfjetAK4eta[jet],pfjetAK4phi[jet],phoeta[pho],phophi[pho])<0.4){clean_jet = false; break;}
	}

if(!clean_jet) continue;

if(isMC && (pfjetAK4pt[0] > 3*qscale)) { break; }

if(pfjetAK4pt[jet]< ptcut) continue;
if(abs(pfjetAK4eta[jet])> ycut) continue;
if((pfjetAK4btag[jet] < btagvalue)&&(pfjetAK4qgl[jet]<0.5)) continue; //targetting only quark jet

pfjetAK4pt[npfjets1] = pfjetAK4pt[jet] ;
pfjetAK4eta[npfjets1] = pfjetAK4eta[jet];
pfjetAK4phi[npfjets1] = pfjetAK4phi[jet];
pfjetAK4mass[npfjets1] = pfjetAK4mass[jet];
pfjetAK4btag[npfjets1] = pfjetAK4btag[jet];

if(pfjetAK4btag[jet] > btagvalue){
	bjetAK4pt[nbjetAK4] = pfjetAK4pt[jet] ;
	bjetAK4eta[nbjetAK4] = pfjetAK4eta[jet];
	bjetAK4phi[nbjetAK4] = pfjetAK4phi[jet];
	bjetAK4mass[nbjetAK4] = pfjetAK4mass[jet];
	bjetAK4btag[nbjetAK4] = pfjetAK4btag[jet];
	bjetAK4hadronflav[nbjetAK4] = pfjetAK4hadronflav[jet];
	bjetAK4partonflav[nbjetAK4] = pfjetAK4partonflav[jet];
	nbjetAK4++;
	}
	else{
			qjetAK4pt[nqjetAK4] = pfjetAK4pt[jet] ;
			qjetAK4eta[nqjetAK4] = pfjetAK4eta[jet];
			qjetAK4phi[nqjetAK4] = pfjetAK4phi[jet];
			qjetAK4mass[nqjetAK4] = pfjetAK4mass[jet];
			qjetAK4btag[nqjetAK4] = pfjetAK4btag[jet];
			qjetAK4hadronflav[nqjetAK4] = pfjetAK4hadronflav[jet];
			qjetAK4partonflav[nqjetAK4] = pfjetAK4partonflav[jet];
			nqjetAK4++;
		}

npfjets1++;
	}
npfjetAK4 =  npfjets1;
//nbjetAK4 = nbjetAK4;

nbgenjetAK4 = 0;
nqgenjetAK4 = 0;
int ngjets1 = 0;

for(int jet=0; jet<ngenjetAK4; jet++){
	
	if(genjetAK4pt[jet]< ptcut) continue;
	if(abs(genjetAK4y[jet])> ycut) continue;
	
	genjetAK4pt[ngjets1] = genjetAK4pt[jet];
	genjetAK4y[ngjets1] = genjetAK4y[jet];
	genjetAK4phi[ngjets1] = genjetAK4phi[jet];
	genjetAK4mass[ngjets1] = genjetAK4mass[jet];
	genjetAK4sdmass[ngjets1] = genjetAK4sdmass[jet];
	genjetAK4btag[ngjets1] =0;
	
	for(int igen=0; igen<ngenparticles; igen++){
	if(/*(abs(genpartstatus[igen])==23)&&genpartfromhard[igen]&&*/abs(genpartpdg[igen])==5){
         if(delta2R(genjetAK4y[ngjets1],genjetAK4phi[ngjets1],genparteta[igen],genpartphi[igen])<0.4){
			genjetAK4btag[ngjets1] = 1.;
			}
		}
	}
	
	if(genjetAK4btag[ngjets1] > btagvalue){
		
		bgenjetAK4pt[nbgenjetAK4] = genjetAK4pt[ngjets1];
		bgenjetAK4y[nbgenjetAK4] = genjetAK4y[ngjets1];
		bgenjetAK4phi[nbgenjetAK4] = genjetAK4phi[ngjets1];
		bgenjetAK4mass[nbgenjetAK4] = genjetAK4mass[ngjets1];
		bgenjetAK4sdmass[nbgenjetAK4] = genjetAK4sdmass[ngjets1];
		bgenjetAK4btag[nbgenjetAK4] = genjetAK4btag[ngjets1];
		nbgenjetAK4 ++;
		}else{
			qgenjetAK4pt[nqgenjetAK4] = genjetAK4pt[ngjets1];
			qgenjetAK4y[nqgenjetAK4] = genjetAK4y[ngjets1];
			qgenjetAK4phi[nqgenjetAK4] = genjetAK4phi[ngjets1];
			qgenjetAK4mass[nqgenjetAK4] = genjetAK4mass[ngjets1];
			qgenjetAK4sdmass[nqgenjetAK4] = genjetAK4sdmass[ngjets1];
			qgenjetAK4btag[nqgenjetAK4] = genjetAK4btag[ngjets1];
			nqgenjetAK4++;
			}
			
		ngjets1++;	
 }
ngenjetAK4 = ngjets1;

int ihlttrg[nHLTmx+1]= {0};

     if(trig_value>1){
     ihlttrg[nHLTmx] = 1;
     for (int ij=0; ij<nHLTmx; ij++) {
     ihlttrg[ij] = *(dec2bin(trig_value)+ij);
      }
     }      

    if(!isMC && (ihlttrg[nHLTmx] ==0) ) return kFALSE;
    float pre_wt = 1.e+10;
    if(!isMC){
    for(int ht=0; ht<nHLTmx; ht++){
    if(ihlt[ht] > 0){
		if(prescl[ht] < pre_wt) { pre_wt = prescl[ht]; }
			}
		}
	weight *= pre_wt;
	}

//     h_nmuon->Fill(nmuons,weight);
 
     for(int nmu=0; nmu<nmuons; nmu++){
 
	 int etamu = getbinid(abs(muoneta[nmu]),netarange,etarng);
	// if(isMC) { muonpt[nmu] *= Roch_Cor(muonq[nmu],muonpt[nmu],muoneta[nmu],muonphi[nmu]); }
     if(etamu>=0){
		 h_mupt[etamu]->Fill(muonpt[nmu],weight);
		 }
 
     h_mueta->Fill(muoneta[nmu],weight);
     h_muphi->Fill(muonphi[nmu],weight);
     h_mudrbm->Fill(muondrbm[nmu],weight);
     h_mutrkvtx->Fill(muontrkvtx[nmu],weight);
     h_muhbye->Fill(max(10.0,muonhcal[nmu]*1./(muonecal[nmu]+muonhcal[nmu])),weight);
     h_muemiso->Fill(muonemiso[nmu],weight);
     h_muhadiso->Fill(muonhadiso[nmu],weight);
     h_mupfiso->Fill(muonpfiso[nmu],weight);
	 h_mukpt03->Fill(muontkpt03[nmu],weight);
	 h_mukpt05->Fill(muontkpt05[nmu],weight);
	 h_muhit->Fill(muonhit[nmu],weight);
	 h_mumst->Fill(muonmst[nmu],weight);
	 h_mutrklay->Fill(muontrklay[nmu],weight);

	 }

     bool mu_pass = false;
	 if((nmuons==1)||(nmuons>1 && muonpt[1]<15.)) { mu_pass = true; }
	 
	 bool e_pass = false;
	 if((nelecs==1)||(nelecs>1 && elpt[1]<15.)) { e_pass = true; }
	
	 if(isMC){
	
	 for(int ijet=0; ijet<nqjetAK4; ijet++){
		 
		 int etatag = getbinid(abs(qjetAK4eta[ijet]),netarange,etarng);
		 
		 int match_gen = -1;
		 float min_delR = 0.2;
		 
		 for(int gjet=0; gjet<nqgenjetAK4; gjet++){
			 if(delta2R(qjetAK4eta[ijet],qjetAK4phi[ijet],qgenjetAK4y[gjet],qgenjetAK4phi[gjet]) < min_delR){
				 min_delR = delta2R(qjetAK4eta[ijet],qjetAK4phi[ijet],qgenjetAK4y[gjet],qgenjetAK4phi[gjet]);
				 match_gen = gjet;
				 }			 
			 }
			 
		 if(match_gen>=0){
			 int pttag = getbinid(genjetAK4pt[match_gen],noptbins,ptbins);
			 if(pttag>=0){
			 h_resopt_AK4[etatag][pttag]->Fill((qjetAK4pt[ijet]-qgenjetAK4pt[match_gen])*1./qgenjetAK4pt[match_gen],weight);
				}
			 }
		 
		 }
		 
	 }

	 h_nmuon_pre->Fill(nmuons,weight);
	 h_idmuon_pre->Fill(muonistight[0],weight);
	 h_idelec_pre->Fill(elecistight[0],weight);
	 h_nelec_pre->Fill(nelecs,weight);
	 h_npho_pre->Fill(nphotons,weight);
	 h_npfjetAK4_pre->Fill(npfjetAK4);
	 h_nbjetAK4_pre->Fill(nbjetAK4);

	 if(!mu_pass && !e_pass) return kFALSE;
	 if(mu_pass && e_pass) return kFALSE; 
	
//	 if(!mu_pass) return kFALSE;
	
//	 if(nphotons>0) return kFALSE;
	 if(npfjetAK4<4) return kFALSE;
	 if(nbjetAK4<2) return kFALSE;
	 if(nqjetAK4<2) return kFALSE;

     for(int ijet=0; ijet<npfjetAK8; ijet++){
		 
		 if(abs(pfjetAK8y[ijet]) > ycut) continue;
		 
		 h_jetyAK8->Fill(pfjetAK8y[ijet], weight);
		 h_jetphiAK8->Fill(pfjetAK8phi[ijet], weight);
		 
		 int etatag = getbinid(abs(pfjetAK8y[ijet]),netarange,etarng);
	
		 if(etatag>=0){
			 h_jetptAK8[etatag]->Fill(pfjetAK8pt[ijet],weight);
			 h_jetmassAK8[etatag]->Fill(pfjetAK8mass[ijet],weight);
			 h_jetsdmassAK8[etatag]->Fill(pfjetAK8sdmass[ijet],weight);
			 h_jetbtagAK8[etatag]->Fill(pfjetAK8btag[ijet],weight);
			 h_jettau21AK8[etatag]->Fill(pfjetAK8tau2[ijet]*1./pfjetAK8tau1[ijet],weight);
			 h_jettau32AK8[etatag]->Fill(pfjetAK8tau3[ijet]*1./pfjetAK8tau2[ijet],weight);
			 if((pfjetAK8tau3[ijet]*1./pfjetAK8tau2[ijet]<0.8)&&(pfjetAK8sub1btag[ijet]>btagvalue||pfjetAK8sub2btag[ijet]>btagvalue)){
			 h_jetmassAK8_sel[etatag]->Fill(pfjetAK8mass[ijet],weight);
			 h_jetsdmassAK8_sel[etatag]->Fill(pfjetAK8sdmass[ijet],weight);
			 if(pfjetAK8pt[ijet]>300){
				 h_jetmassAK8_selpt[etatag]->Fill(pfjetAK8mass[ijet],weight);
				 h_jetsdmassAK8_selpt[etatag]->Fill(pfjetAK8sdmass[ijet],weight);
					}
				}	
			 }
		 }
	 
	 TLorentzVector fourAK4_vec4[4];
	 float min_massdiff = 100;
	 int wtag[2] = {-1,-1};
	 
	 for(int ijet=0; ijet<nqjetAK4; ijet++){
		 for(int gjet=(ijet+1); gjet<nqjetAK4; gjet++){
		//	 if((pfjetAK4btag[ijet]<btagvalue)&&(pfjetAK4btag[gjet]<btagvalue)){
			 TLorentzVector tmp[2];
			 tmp[0].SetPtEtaPhiM(qjetAK4pt[ijet],qjetAK4eta[ijet],qjetAK4phi[ijet],qjetAK4mass[ijet]);
			 tmp[1].SetPtEtaPhiM(qjetAK4pt[gjet],qjetAK4eta[gjet],qjetAK4phi[gjet],qjetAK4mass[gjet]);
				if(abs((tmp[0]+tmp[1]).M()-81.)<min_massdiff) { 
					wtag[0] = ijet; wtag[1] = gjet;
					min_massdiff = abs((tmp[0]+tmp[1]).M()-81.);
					}
		//		}
			 }
		 }
	 
	 h_mupt_pass->Fill(muonpt[0],weight);
	 h_mueta_pass->Fill(muoneta[0],weight);
	 h_muphi_pass->Fill(muonphi[0],weight);

     h_nmuon->Fill(nmuons,weight);
     h_jetnqbAK4->Fill(npfjetAK4,weight);
     h_jetnqAK4->Fill(nqjetAK4,weight);
     h_jetnbAK4->Fill(nbjetAK4,weight);
	 
	 h_jetpt_q1AK4->Fill(qjetAK4pt[0],weight);
	 h_jeteta_q1AK4->Fill(qjetAK4eta[0],weight);
	 h_jetphi_q1AK4->Fill(qjetAK4phi[0],weight);
	 
	 h_jetpt_q2AK4->Fill(qjetAK4pt[1],weight);
	 h_jeteta_q2AK4->Fill(qjetAK4eta[1],weight);
	 h_jetphi_q2AK4->Fill(qjetAK4phi[1],weight);
	 
	 h_jetpt_b1AK4->Fill(bjetAK4pt[0],weight);
	 h_jeteta_b1AK4->Fill(bjetAK4eta[0],weight);
	 h_jetphi_b1AK4->Fill(bjetAK4phi[0],weight);
	 
	 h_jetpt_b2AK4->Fill(bjetAK4pt[1],weight);
	 h_jeteta_b2AK4->Fill(bjetAK4eta[1],weight);
	 h_jetphi_b2AK4->Fill(bjetAK4phi[1],weight);
	 
	 h_metpt->Fill(PFMET,weight);
	 h_metphi->Fill(PFMETPhi,weight);
	 h_metbyEt->Fill(PFMET*1./sumEt,weight);
	 h_metSig->Fill(MisEtSig,weight);
	 
	 for(int ijet=0; ijet<npfjetAK4; ijet++){
		 
		 h_jetyAK4->Fill(pfjetAK4eta[ijet], weight);
		 h_jetphiAK4->Fill(pfjetAK4phi[ijet], weight);
		 
		 int etatag = getbinid(abs(pfjetAK4eta[ijet]),netarange,etarng);
		 if(etatag>=0){
			 h_jetptAK4[etatag]->Fill(pfjetAK4pt[ijet],weight);
			 h_jetmassAK4[etatag]->Fill(pfjetAK4mass[ijet],weight);
			 h_jetsdmassAK4[etatag]->Fill(pfjetAK4sdmass[ijet],weight);
			 h_jetbtagAK4[etatag]->Fill(pfjetAK4btag[ijet],weight);
			 h_jettau21AK4[etatag]->Fill(pfjetAK4tau2[ijet]*1./pfjetAK4tau1[ijet],weight);
			 h_jettau32AK4[etatag]->Fill(pfjetAK4tau3[ijet]*1./pfjetAK4tau2[ijet],weight);
			 }
		 
		 }//ijet
		
	 int iw[2] ={0};
	 int iwc = 0;
	 int ib[2] = {0};
	 int iwb = 0;
	 
	 int iweta = -1;
	 int itm_bin[selnmass] = {-1,-1,-1,-1};
		
		  for(int ijet=0; ijet<nqjetAK4; ijet++){
	//		  if(abs(qjetAK4partonflav[ijet])<5){
				  if(iwc>1) break;
				  iw[iwc] = ijet;
				  fourAK4_vec4[iwc].SetPtEtaPhiM(qjetAK4pt[iw[iwc]],qjetAK4eta[iw[iwc]],qjetAK4phi[iw[iwc]],qjetAK4mass[iw[iwc]]);
				  iwc++;
	//			  }
		  }
	
	/*
		if((wtag[0]>-1)&&(wtag[1]>-1)){
		 fourAK4_vec4[0].SetPtEtaPhiM(pfjetAK4pt[wtag[0]],pfjetAK4eta[wtag[0]],pfjetAK4phi[wtag[0]],pfjetAK4mass[wtag[0]]);
		 fourAK4_vec4[1].SetPtEtaPhiM(pfjetAK4pt[wtag[1]],pfjetAK4eta[wtag[1]],pfjetAK4phi[wtag[1]],pfjetAK4mass[wtag[1]]);
		}
	*/ 
	// cout<<"wtags2 "<<iw[0]<<" "<<iw[1]<<endl;
	 
//	 fourAK4_vec4[0].SetPtEtaPhiM(qjetAK4pt[0],qjetAK4eta[0],qjetAK4phi[0],qjetAK4mass[0]);
//   fourAK4_vec4[1].SetPtEtaPhiM(qjetAK4pt[1],qjetAK4eta[1],qjetAK4phi[1],qjetAK4mass[1]);
     fourAK4_vec4[2].SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
     fourAK4_vec4[3].SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);

     TLorentzVector p4w = fourAK4_vec4[0]+fourAK4_vec4[1];
	 float mwj = p4w.M();
	 
     TLorentzVector p4t;

	 if(mwj>min_wmass && mwj<max_wmass){
	 
	 h_whadmass->Fill(p4w.M(),weight);
	 h_whadpt->Fill(p4w.Pt(),weight);
	 h_whady->Fill(p4w.Eta(),weight);
	 
	 }
	 
	 if(iwc>=2){
	 
	 TLorentzVector p4t_1 = p4w + fourAK4_vec4[2];
	 TLorentzVector p4t_2 = p4w + fourAK4_vec4[3];
	 
	 if(mu_pass){
		if(delta2R(muoneta[0],muonphi[0],fourAK4_vec4[2].Eta(),fourAK4_vec4[2].Phi())<delta2R(muoneta[0],muonphi[0],fourAK4_vec4[3].Eta(),fourAK4_vec4[3].Phi())){
		 p4t =  p4w + fourAK4_vec4[3];
			}
			else{
			   p4t =  p4w + fourAK4_vec4[2]; 
				}
		 }
		 
	 if(e_pass){
		if(delta2R(eleta[0],elphi[0],fourAK4_vec4[2].Eta(),fourAK4_vec4[2].Phi())<delta2R(eleta[0],elphi[0],fourAK4_vec4[3].Eta(),fourAK4_vec4[3].Phi())){
		 p4t =  p4w + fourAK4_vec4[3];
			}
		else{
			   p4t =  p4w + fourAK4_vec4[2]; 
				}
		 }  
			  
	 
//	 p4t = (abs(p4t_2.M() - 173.) < abs(p4t_1.M() - 173.)) ? p4t_2 : p4t_1;

     if(p4t.M()>min_tmass && p4t.M()<max_tmass){

	 h_thadmass->Fill(p4t.M(),weight);
	 h_thadpt->Fill(p4t.Pt(),weight);
	 h_thady->Fill(p4t.Eta(),weight);
	 h_tw_2d->Fill(p4t.M(),p4w.M(),weight);
	 /*
	 float prob;
	 if(p4w.M()>min_wmass && p4w.M()<max_wmass){ 
	 prob = h_topvs_w_true->GetBinContent(xaxis_true->FindBin(p4t.M()),yaxis_true->FindBin(p4w.M()));
	 prob = -log10(prob);
	 h_lambda->Fill(prob,weight);
		}
	*/   
	   }
   
//		FOR BB BE EE

	   if(abs(fourAK4_vec4[0].Eta())<ydiv && abs(fourAK4_vec4[1].Eta())<ydiv) {iweta = 0;}
		if((abs(fourAK4_vec4[0].Eta())<ydiv && abs(fourAK4_vec4[1].Eta())>ydiv)||(abs(fourAK4_vec4[0].Eta())>ydiv && abs(fourAK4_vec4[1].Eta())<ydiv)) {iweta = 1;}
		 if(abs(fourAK4_vec4[0].Eta())>ydiv && abs(fourAK4_vec4[1].Eta())>ydiv) {iweta = 2;} 
    
			if(iweta==0||iweta==2){
		    h_jety_2y[0][iweta]->Fill(fourAK4_vec4[0].Eta(),weight);
		    h_jety_2y[1][iweta]->Fill(fourAK4_vec4[1].Eta(),weight);
		    
		    h_jety_2pt[0][iweta]->Fill(fourAK4_vec4[0].Pt(),weight);
		    h_jety_2pt[1][iweta]->Fill(fourAK4_vec4[1].Pt(),weight);
			}
			if(iweta==1){
				if(abs(fourAK4_vec4[0].Eta()) < abs(fourAK4_vec4[1].Eta())) { 
					h_jety_2y[0][iweta]->Fill(fourAK4_vec4[0].Eta(),weight);
					h_jety_2y[1][iweta]->Fill(fourAK4_vec4[1].Eta(),weight);
					
					h_jety_2pt[0][iweta]->Fill(fourAK4_vec4[0].Pt(),weight);
					h_jety_2pt[1][iweta]->Fill(fourAK4_vec4[1].Pt(),weight);
					}
					else{
						h_jety_2y[0][iweta]->Fill(fourAK4_vec4[1].Eta(),weight);
						h_jety_2y[1][iweta]->Fill(fourAK4_vec4[0].Eta(),weight);
						
						h_jety_2pt[0][iweta]->Fill(fourAK4_vec4[1].Pt(),weight);
						h_jety_2pt[1][iweta]->Fill(fourAK4_vec4[0].Pt(),weight);
						}
				}
	
		//	iweta =  getbinid(max(fabs(fourAK4_vec4[0].Eta()),fabs(fourAK4_vec4[1].Eta())),selneta,ydivs);
			
			if(p4t.M()<selmass[0]) { itm_bin[0] = 0;}
			if(p4t.M()<selmass[1]) { itm_bin[1] = 1;}
			if(p4t.M()<selmass[2]) { itm_bin[2] = 2;} 
			if(p4t.M()<selmass[3]) { itm_bin[3] = 3;}
    
		for(int jc=0; jc<11; jc++){
			float jec = 1+(jc-5)*0.005;
			TLorentzVector tmp1 = fourAK4_vec4[0] * jec;
			TLorentzVector tmp2 = fourAK4_vec4[1] * jec;
			if(iweta>=0){
				if((tmp1+tmp2).M()>min_wmass && (tmp1+tmp2).M()<max_wmass){
				if(itm_bin[0]>=0){
				h_whadmass_sel[jc][itm_bin[0]][iweta]->Fill((tmp1+tmp2).M(),weight);
							 }
				if(itm_bin[1]>=0){
				h_whadmass_sel[jc][itm_bin[1]][iweta]->Fill((tmp1+tmp2).M(),weight);
							 }
				if(itm_bin[2]>=0){
				h_whadmass_sel[jc][itm_bin[2]][iweta]->Fill((tmp1+tmp2).M(),weight);
							 }
				if(itm_bin[3]>=0){
				h_whadmass_sel[jc][itm_bin[3]][iweta]->Fill((tmp1+tmp2).M(),weight);
							 }			 			 			 
						}
					}
				}
     }
	 
	 
	 for(int iv=0; iv<4; iv++){
		 fourAK4_vec4[iv].SetPxPyPzE(0,0,0,0);
		 }

	 p4w.SetPxPyPzE(0,0,0,0);
	 p4t.SetPxPyPzE(0,0,0,0);

	 min_massdiff = 150;
	 int toptag[3] = {-1,-1,-1};

	 for(int ijet=0; ijet<nqjetAK4; ijet++){
		 for(int gjet=(ijet+1); gjet<nqjetAK4; gjet++){
	//		 if((abs(qjetAK4partonflav[ijet])<5)&&(abs(qjetAK4partonflav[gjet])<5)) {
			 for(int bjet=0; bjet<nbjetAK4; bjet++){
		      TLorentzVector tmp[3];
		      tmp[0].SetPtEtaPhiM(qjetAK4pt[ijet],qjetAK4eta[ijet],qjetAK4phi[ijet],qjetAK4mass[ijet]);
		      tmp[1].SetPtEtaPhiM(qjetAK4pt[gjet],qjetAK4eta[gjet],qjetAK4phi[gjet],qjetAK4mass[gjet]);
		      tmp[2].SetPtEtaPhiM(bjetAK4pt[bjet],bjetAK4eta[bjet],bjetAK4phi[bjet],bjetAK4mass[bjet]);
		      float diff_tot = fabs((tmp[0]+tmp[1]).M()-81.)+fabs((tmp[0]+tmp[1]+tmp[2]).M()- 173);
		      if(diff_tot<min_massdiff) {
				  toptag[0] = ijet; toptag[1] = gjet; toptag[2] = bjet;
				  min_massdiff = diff_tot;
						}
					}
		//		}  //partonflav
			}
		 }	 
	//	 cout<<"toptags "<<toptag[0]<<'\t'<<toptag[1]<<'\t'<<toptag[2]<<" min_massdiff "<<min_massdiff<<endl;
		 
	 if((toptag[0]>=0)&&(toptag[1]>=0)&&(toptag[2]>=0)){
	//    cout<<"passed toptags "<<toptag[0]<<'\t'<<toptag[1]<<'\t'<<toptag[2]<<" min_massdiff "<<min_massdiff<<endl;
		 
		  TLorentzVector tmp[3];
		  tmp[0].SetPtEtaPhiM(qjetAK4pt[toptag[0]],qjetAK4eta[toptag[0]],qjetAK4phi[toptag[0]],qjetAK4mass[toptag[0]]);
		  tmp[1].SetPtEtaPhiM(qjetAK4pt[toptag[1]],qjetAK4eta[toptag[1]],qjetAK4phi[toptag[1]],qjetAK4mass[toptag[1]]);
		  tmp[2].SetPtEtaPhiM(bjetAK4pt[toptag[2]],bjetAK4eta[toptag[2]],bjetAK4phi[toptag[2]],bjetAK4mass[toptag[2]]);
		  
		  p4w = tmp[0]+tmp[1];
		  p4t = p4w+tmp[2];
		  
		  if(p4w.M()>min_wmass && p4w.M()<max_wmass){
		  h_whadmassA->Fill(p4w.M(),weight);
		  h_whadptA->Fill(p4w.Pt(),weight);
		  h_whadyA->Fill(p4w.Eta(),weight);
			}
	 
		  if(p4t.M()>min_tmass && p4t.M()<max_tmass){
			   h_thadmassA->Fill(p4t.M(),weight);
			   h_thadptA->Fill(p4t.Pt(),weight);
			   h_thadyA->Fill(p4t.Eta(),weight);
			   h_tw_2dA->Fill(p4t.M(),p4w.M(),weight);
			   
			  /* 
			  float prob;
			  if(p4w.M()>min_wmass && p4w.M()<max_wmass){ 
			  prob = h_topvs_w_true->GetBinContent(xaxis_true->FindBin(p4t.M()),yaxis_true->FindBin(p4w.M()));
			  prob = -log10(prob);
			  h_lambdaA->Fill(prob,weight);
					}
			 */  
			  }
		 }

//	 cout<<"w mass "<<mwj<<" top mass "<<p4t.M()<<endl;
	 
	 
	 if(isMC){
	 
	 ///// parton start /////
	 
	 TLorentzVector truthpart[4];
	 TLorentzVector wpartons[2];
	 TLorentzVector bpartons[2];
	 int wpartons_mom[2], bpartons_mom[2];
	 
	 TLorentzVector thad; TLorentzVector whad;  TLorentzVector bhad; TLorentzVector pq1; TLorentzVector pq2;
	 TLorentzVector tlep; TLorentzVector wlep;  TLorentzVector blep; TLorentzVector plep; TLorentzVector pnu;
	 int bhad_momid; int blep_momid;
	 
	 int itruthq = 0;
	 int iwq = 0;
	 int ibq = 0;
	 
	 for(int igen=0; igen<ngenparticles; igen++){
		 
		 if(/*(abs(genpartstatus[igen])!=22)||*/(abs(genpartstatus[igen])!=23)) continue;
		 if(!genpartfromhard[igen]) continue;
		 if(abs(genpartpdg[igen])>5) continue;
		 if((abs(genpartmompdg[igen])==24)||(abs(genpartmompdg[igen])==6)){
		 truthpart[itruthq].SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
		 itruthq++;
			}
		 if(abs(genpartmompdg[igen])==24 && abs(genpartpdg[igen])<5){
				wpartons[iwq].SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
				wpartons_mom[iwq] = genpartmompdg[igen];
				iwq++;
			 }
		 if(abs(genpartmompdg[igen])==6 && abs(genpartpdg[igen])==5){
				bpartons[ibq].SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
				bpartons_mom[ibq] = genpartmompdg[igen];
				ibq++;
			 }	 
		 }
	 if(iwq>1){
	 h_whadmass_parton->Fill((wpartons[0]+wpartons[1]).M(),weight);
	 
	 pq1 = wpartons[0];
	 pq2 = wpartons[1];
	 whad = (wpartons[0]+wpartons[1]);
	 
	 if((wpartons_mom[0] == wpartons_mom[1])&&(bpartons_mom[0]*bpartons_mom[1]<0)){
		 for(int ibp=0; ibp<2; ibp++){
			 if(bpartons_mom[ibp]*wpartons_mom[0] > 0){
				 bhad = bpartons[ibp];
				 bhad_momid = bpartons_mom[ibp];
				 }
				 else{
					 blep = bpartons[ibp];
					 blep_momid = bpartons_mom[ibp];
					 }
			 }
			 thad = whad + bhad;
		 }
	 }
	 
	 int plep_mom;
	 int nlep = 0;
	 
	 for(int igen=0; igen<ngenparticles; igen++){
		 
		 if(abs(genpartstatus[igen])!=23) continue;
		 if(!genpartfromhard[igen]) continue;
		 if(abs(genpartpdg[igen])<10 || abs(genpartpdg[igen])>20) continue;
	//	 cout<<"genpartid "<<genpartpdg[igen]<<endl;
		 if(((abs(genpartpdg[igen])==11)||(abs(genpartpdg[igen])==13)||(abs(genpartpdg[igen])==15)) && (abs(genpartmompdg[igen])==24)){
			 plep.SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
			 plep_mom = genpartmompdg[igen];
			 nlep++;
			 }
		 if(((abs(genpartpdg[igen])==12)||(abs(genpartpdg[igen])==14)||(abs(genpartpdg[igen])==16)) && (abs(genpartmompdg[igen])==24)){
			 pnu.SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
			 }	 
		 }
	 
		wlep = plep + pnu;
		if(plep_mom * blep_momid > 0){
		tlep = wlep + blep;
		}
	 
	 pq1_pt = pq1.Pt(); pq1_eta = pq1.Eta(); pq1_phi = pq1.Phi(); pq1_e = pq1.E();
	 pq2_pt = pq2.Pt(); pq2_eta = pq2.Eta(); pq2_phi = pq2.Phi(); pq2_e = pq1.E();
	 whad_pt = whad.Pt(); whad_eta = whad.Eta(); whad_phi = whad.Phi(); whad_e = whad.E(); whad_m = whad.M();
	 bhad_pt = bhad.Pt(); bhad_eta = bhad.Eta(); bhad_phi = bhad.Phi(); bhad_e = bhad.E(); bhad_m = bhad.M();
	 thad_pt = thad.Pt(); thad_eta = thad.Eta(); thad_phi = thad.Phi(); thad_e = thad.E(); thad_m = thad.M();
	 
	 if(nlep>0){
	 plep_pt = plep.Pt(); plep_eta = plep.Eta(); plep_phi = plep.Phi(); plep_e = plep.E();
	 wlep_pt = wlep.Pt(); wlep_eta = wlep.Eta(); wlep_phi = wlep.Phi(); wlep_e = wlep.E(); wlep_m = wlep.M();
	 blep_pt = blep.Pt(); blep_eta = blep.Eta(); blep_phi = blep.Phi(); blep_e = blep.E(); blep_m = blep.M();
	 tlep_pt = tlep.Pt(); tlep_eta = tlep.Eta(); tlep_phi = tlep.Phi(); tlep_e = tlep.E(); tlep_m = tlep.M();
	 
	 delR_whadbhad = delta2R(whad_eta,whad_phi,bhad_eta,bhad_phi); delPhi_whadbhad = PhiInRange(whad_phi-bhad_phi); delEta_whadbhad = abs(whad_eta - bhad_eta);
	 delR_whadblep = delta2R(whad_eta,whad_phi,blep_eta,blep_phi); delPhi_whadblep = PhiInRange(whad_phi-blep_phi); delEta_whadblep = abs(whad_eta - blep_eta);
	 delR_whadplep = delta2R(whad_eta,whad_phi,plep_eta,plep_phi); delPhi_whadplep = PhiInRange(whad_phi-plep_phi); delEta_whadplep = abs(whad_eta - plep_eta);
	 delR_bhadplep = delta2R(bhad_eta,bhad_phi,plep_eta,plep_phi); delPhi_bhadplep = PhiInRange(bhad_phi-plep_phi); delEta_bhadplep = abs(bhad_eta - plep_eta);
	 delR_blepplep = delta2R(blep_eta,blep_phi,plep_eta,plep_phi); delPhi_blepplep = PhiInRange(blep_phi-plep_phi); delEta_blepplep = abs(blep_eta - plep_eta);
	 }
	 
	 ///// parton end ////
	/* 
	 if(ibq>1){
	 if(delta2R((wpartons[0]+wpartons[1]).Eta(),(wpartons[0]+wpartons[1]).Phi(),bpartons[0].Eta(),bpartons[0].Phi()) < delta2R((wpartons[0]+wpartons[1]).Eta(),(wpartons[0]+wpartons[1]).Phi(),bpartons[1].Eta(),bpartons[1].Phi())){
		  h_thadmass_parton->Fill((wpartons[0]+wpartons[1]+bpartons[0]).M(),weight);
		 }else{
			   h_thadmass_parton->Fill((wpartons[0]+wpartons[1]+bpartons[1]).M(),weight);
			   TLorentzVector tmpvec4;
			   tmpvec4 = bpartons[0];
			   bpartons[0] = bpartons[1];
			   bpartons[1] = tmpvec4;
			  }
		  }
	 */
	 h_thadmass_parton->Fill(thad.M(),weight);
	 bpartons[0] = bhad;
	 bpartons[1] = blep;
	 
	// cout<<"ngenparticles "<<ngenparticles<<" iwq "<<iwq<<" ibq "<<ibq<<endl;
	 
	 TLorentzVector w_truthjet[3];
	 int match_w[2] = {-1,-1};
	 float mindelR = 0.4;
	 
	 if(iwq==2){
	 for(int iwc=0; iwc<iwq; iwc++){
		 if(iwc>1) break;
		 int iMatch = -1;
		 for(int ijet=0; ijet<nqjetAK4; ijet++){
		//	 if(pfjetAK4btag[ijet]<btagvalue){
			 if(delta2R(qjetAK4eta[ijet],qjetAK4phi[ijet],wpartons[iwc].Eta(),wpartons[iwc].Phi())<mindelR){
						iMatch = ijet;
						mindelR = delta2R(qjetAK4eta[ijet],qjetAK4phi[ijet],wpartons[iwc].Eta(),wpartons[iwc].Phi());
					}
		//		}
			}
		 match_w[iwc] = iMatch;
		 mindelR = 0.4;
		}
	  }
      
      int match_b[2]={-1,-1};
      TLorentzVector btruthjet[2];
	  TLorentzVector ttruthjet;
      
      if((match_w[0]>=0)&&(match_w[1]>=0)&&(match_w[0]!=match_w[1])){
		  if((match_w[0]<2)&&(match_w[1]<2)){
		  
		  w_truthjet[0].SetPtEtaPhiM(qjetAK4pt[match_w[0]],qjetAK4eta[match_w[0]],qjetAK4phi[match_w[0]],qjetAK4mass[match_w[0]]);
		  w_truthjet[1].SetPtEtaPhiM(qjetAK4pt[match_w[1]],qjetAK4eta[match_w[1]],qjetAK4phi[match_w[1]],qjetAK4mass[match_w[1]]);
		  w_truthjet[2] = w_truthjet[0] + w_truthjet[1];
		  
		  if(w_truthjet[2].M()>min_wmass && w_truthjet[2].M()<max_wmass){
		  
		  h_whadmass_truth->Fill(w_truthjet[2].M(),weight);
		  h_whadpt_truth->Fill(w_truthjet[2].Pt(),weight);
		  h_whady_truth->Fill(w_truthjet[2].Eta(),weight);
		  
			}
		  
		  for(int jc=0; jc<11; jc++){
						float jec = 1+(jc-5)*0.005;
							TLorentzVector tmp1 = w_truthjet[0] * jec;
							TLorentzVector tmp2 = w_truthjet[1] * jec;
							if(iweta>=0){
								if((tmp1+tmp2).M()>min_wmass && (tmp1+tmp2).M()<max_wmass){
								if(itm_bin[0]>=0){
								h_whadmass_truth_sel[jc][itm_bin[0]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}
								if(itm_bin[1]>=0){
								h_whadmass_truth_sel[jc][itm_bin[1]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}
								if(itm_bin[2]>=0){
								h_whadmass_truth_sel[jc][itm_bin[2]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}
								if(itm_bin[3]>=0){
								h_whadmass_truth_sel[jc][itm_bin[3]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}			 			 			 
									}
								}
							}
		  
		  if(ibq>0){
			  for(int ibc=0; ibc<ibq; ibc++){
				  if(ibc>1) break;
				  int iMatch = -1;
				  mindelR = 0.4;
				  for(int ibjet=0; ibjet<nbjetAK4; ibjet++){
					  if(ibjet==match_b[0]) continue;
						if(delta2R(bjetAK4eta[ibjet],bjetAK4phi[ibjet],bpartons[ibc].Eta(),bpartons[ibc].Phi())<mindelR){
							iMatch = ibjet;
							mindelR = delta2R(bjetAK4eta[ibjet],bjetAK4phi[ibjet],bpartons[ibc].Eta(),bpartons[ibc].Phi());
							}
						}
						match_b[ibc] = iMatch;
				  }
			  }
			  
			  bool had_b = false;
			  
			  if((match_b[0]>=0) && (match_b[1]>=0) && (match_b[0]!=match_b[1])){
			  if((match_b[0]<2)&&(match_b[1]<2)){
			  
			  btruthjet[0].SetPtEtaPhiM(bjetAK4pt[match_b[0]],bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],bjetAK4mass[match_b[0]]);
			  btruthjet[1].SetPtEtaPhiM(bjetAK4pt[match_b[1]],bjetAK4eta[match_b[1]],bjetAK4phi[match_b[1]],bjetAK4mass[match_b[1]]);
			  
			  if(mu_pass){
						if(delta2R(btruthjet[0].Eta(),btruthjet[0].Phi(),muoneta[0],muonphi[0]) < delta2R(btruthjet[1].Eta(),btruthjet[1].Phi(),muoneta[0],muonphi[0])){
								ttruthjet = (w_truthjet[2]+btruthjet[1]) ;
							}
							else {	ttruthjet = (w_truthjet[2]+btruthjet[0]) ; had_b = true;}
						}
			  if(e_pass){
						if(delta2R(btruthjet[0].Eta(),btruthjet[0].Phi(),eleta[0],elphi[0]) < delta2R(btruthjet[1].Eta(),btruthjet[1].Phi(),eleta[0],elphi[0])){
								ttruthjet = (w_truthjet[2]+btruthjet[1]) ;
							}
							else {	ttruthjet = (w_truthjet[2]+btruthjet[0]) ; had_b = true;}
						}	
			  
				if(had_b){
				    if(ttruthjet.M()>min_tmass && ttruthjet.M()<max_tmass){
				  
				    h_thadmass_truth->Fill(ttruthjet.M(),weight);
				    h_thadpt_truth->Fill(ttruthjet.Pt(),weight);
				    h_thady_truth->Fill(ttruthjet.Eta(),weight);
				    
				    h_tw_2d_truth->Fill(ttruthjet.M(),w_truthjet[2].M(),weight);
				    
						}
					}
					else{
						 if(ttruthjet.M()>min_tmass && ttruthjet.M()<max_tmass){
				  
							h_thadmass_wrong->Fill(ttruthjet.M(),weight);
							h_thadpt_wrong->Fill(ttruthjet.Pt(),weight);
							h_thady_wrong->Fill(ttruthjet.Eta(),weight);
				    
							h_tw_2d_wrong->Fill(ttruthjet.M(),w_truthjet[2].M(),weight);
				    
								}
						}
					}
						else{
							
							 TLorentzVector bjetmatch;
					
							 if(mu_pass){
							 if(delta2R(bjetAK4eta[0],bjetAK4phi[0],muoneta[0],muonphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],muoneta[0],muonphi[0])){
								bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
								ttruthjet = (w_truthjet[2]+bjetmatch) ;
								}
							 else {	
									bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
									ttruthjet = (w_truthjet[2]+bjetmatch);
								  }
								}
							if(e_pass){
							if(delta2R(bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],eleta[0],elphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],eleta[0],elphi[0])){
								bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
								ttruthjet = (w_truthjet[2]+bjetmatch) ;
								}
							else {	
								 bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
								 ttruthjet = (w_truthjet[2]+bjetmatch);
								 }
								}	
						
							h_thadmass_wrong->Fill(ttruthjet.M(),weight);
							h_thadpt_wrong->Fill(ttruthjet.Pt(),weight);
							h_thady_wrong->Fill(ttruthjet.Eta(),weight);
							
							h_tw_2d_wrong->Fill(ttruthjet.M(),w_truthjet[2].M(),weight);
							
							}
					}
					else{
					
					TLorentzVector bjetmatch;
					
					if(mu_pass){
						if(delta2R(bjetAK4eta[0],bjetAK4phi[0],muoneta[0],muonphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],muoneta[0],muonphi[0])){
								bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
								ttruthjet = (w_truthjet[2]+bjetmatch) ;
							}
							else {	
								 bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
								 ttruthjet = (w_truthjet[2]+bjetmatch);
								 }
						}
					if(e_pass){
						if(delta2R(bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],eleta[0],elphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],eleta[0],elphi[0])){
								bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
								ttruthjet = (w_truthjet[2]+bjetmatch) ;
							}
							else {	
								 bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
								 ttruthjet = (w_truthjet[2]+bjetmatch);
								 }
							 }	
						
							h_thadmass_unmatch->Fill(ttruthjet.M(),weight);
							h_thadpt_unmatch->Fill(ttruthjet.Pt(),weight);
							h_thady_unmatch->Fill(ttruthjet.Eta(),weight);
							
							h_tw_2d_unmatch->Fill(ttruthjet.M(),w_truthjet[2].M(),weight);
						
							}
						
				
			}
		  else{
			  TLorentzVector w_wrongjet[3];
			  w_wrongjet[0].SetPtEtaPhiM(qjetAK4pt[0],qjetAK4eta[0],qjetAK4phi[0],qjetAK4mass[0]);
			  w_wrongjet[1].SetPtEtaPhiM(qjetAK4pt[1],qjetAK4eta[1],qjetAK4phi[1],qjetAK4mass[1]);
			  w_wrongjet[2] = w_wrongjet[0] + w_wrongjet[1];
		  
			  if(w_wrongjet[2].M()>min_wmass && w_wrongjet[2].M()<max_wmass){
		  
				h_whadmass_wrong->Fill(w_wrongjet[2].M(),weight);
				h_whadpt_wrong->Fill(w_wrongjet[2].Pt(),weight);
				h_whady_wrong->Fill(w_wrongjet[2].Eta(),weight);
			  }
			  
			  for(int jc=0; jc<11; jc++){
						float jec = 1+(jc-5)*0.005;
							TLorentzVector tmp1 = w_wrongjet[0] * jec;
							TLorentzVector tmp2 = w_wrongjet[1] * jec;
							if(iweta>=0){
								if((tmp1+tmp2).M()>min_wmass && (tmp1+tmp2).M()<max_wmass){
								if(itm_bin[0]>=0){
								h_whadmass_wrong_sel[jc][itm_bin[0]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[1]>=0){
								h_whadmass_wrong_sel[jc][itm_bin[1]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[2]>=0){
								h_whadmass_wrong_sel[jc][itm_bin[2]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[3]>=0){
								h_whadmass_wrong_sel[jc][itm_bin[3]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}			 			 			 
									}
								}
							}
			  
			  if((match_b[0]>=0) && (match_b[1]>=0) && (match_b[0]!=match_b[1])){
				
				TLorentzVector bjetmatch;
					
				if(mu_pass){
					if(delta2R(bjetAK4eta[0],bjetAK4phi[0],muoneta[0],muonphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],muoneta[0],muonphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						   }
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							 }
						}
				if(e_pass){
					if(delta2R(bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],eleta[0],elphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],eleta[0],elphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						}
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							}
						}	
						
					h_thadmass_wrong->Fill(ttruthjet.M(),weight);
					h_thadpt_wrong->Fill(ttruthjet.Pt(),weight);
					h_thady_wrong->Fill(ttruthjet.Eta(),weight);
					
					h_tw_2d_wrong->Fill(ttruthjet.M(),w_wrongjet[2].M(),weight);
				
				}
				else{
					
					TLorentzVector bjetmatch;
					
					if(mu_pass){
					if(delta2R(bjetAK4eta[0],bjetAK4phi[0],muoneta[0],muonphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],muoneta[0],muonphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						   }
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							 }
						}
					if(e_pass){
					if(delta2R(bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],eleta[0],elphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],eleta[0],elphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						}
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							}
						}	
						
					h_thadmass_unmatch->Fill(ttruthjet.M(),weight);
					h_thadpt_unmatch->Fill(ttruthjet.Pt(),weight);
					h_thady_unmatch->Fill(ttruthjet.Eta(),weight);
					
					h_tw_2d_unmatch->Fill(ttruthjet.M(),w_wrongjet[2].M(),weight);
					
					}
			  
			}
		  }
		  else{
			  TLorentzVector w_wrongjet[3];
			  w_wrongjet[0].SetPtEtaPhiM(qjetAK4pt[0],qjetAK4eta[0],qjetAK4phi[0],qjetAK4mass[0]);
			  w_wrongjet[1].SetPtEtaPhiM(qjetAK4pt[1],qjetAK4eta[1],qjetAK4phi[1],qjetAK4mass[1]);
			  w_wrongjet[2] = w_wrongjet[0] + w_wrongjet[1];
		  
			  if(w_wrongjet[2].M()>min_wmass && w_wrongjet[2].M()<max_wmass){
		  
				h_whadmass_unmatch->Fill(w_wrongjet[2].M(),weight);
				h_whadpt_unmatch->Fill(w_wrongjet[2].Pt(),weight);
				h_whady_unmatch->Fill(w_wrongjet[2].Eta(),weight);
				
			  }
			  
			    for(int jc=0; jc<11; jc++){
						float jec = 1+(jc-5)*0.005;
							TLorentzVector tmp1 = w_wrongjet[0] * jec;
							TLorentzVector tmp2 = w_wrongjet[1] * jec;
							if(iweta>=0){
								if((tmp1+tmp2).M()>min_wmass && (tmp1+tmp2).M()<max_wmass){
								if(itm_bin[0]>=0){
								h_whadmass_unmatch_sel[jc][itm_bin[0]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[1]>=0){
								h_whadmass_unmatch_sel[jc][itm_bin[1]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[2]>=0){
								h_whadmass_unmatch_sel[jc][itm_bin[2]][iweta]->Fill((tmp1+tmp2).M(),weight);
									}
								if(itm_bin[3]>=0){
								h_whadmass_unmatch_sel[jc][itm_bin[3]][iweta]->Fill((tmp1+tmp2).M(),weight);
										}			 			 			 
									}
								}
							}
			
				TLorentzVector bjetmatch;
					
				if(mu_pass){
					if(delta2R(bjetAK4eta[0],bjetAK4phi[0],muoneta[0],muonphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],muoneta[0],muonphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						   }
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							 }
						}
				if(e_pass){
					if(delta2R(bjetAK4eta[match_b[0]],bjetAK4phi[match_b[0]],eleta[0],elphi[0]) < delta2R(bjetAK4eta[1],bjetAK4phi[1],eleta[0],elphi[0])){
						bjetmatch.SetPtEtaPhiM(bjetAK4pt[1],bjetAK4eta[1],bjetAK4phi[1],bjetAK4mass[1]);
						ttruthjet = (w_wrongjet[2]+bjetmatch) ;
						}
						else {	
							bjetmatch.SetPtEtaPhiM(bjetAK4pt[0],bjetAK4eta[0],bjetAK4phi[0],bjetAK4mass[0]);
							ttruthjet = (w_wrongjet[2]+bjetmatch);
							}
						}	
						
					h_thadmass_unmatch->Fill(ttruthjet.M(),weight);
					h_thadpt_unmatch->Fill(ttruthjet.Pt(),weight);
					h_thady_unmatch->Fill(ttruthjet.Eta(),weight);
					
					h_tw_2d_unmatch->Fill(ttruthjet.M(),w_wrongjet[2].M(),weight);
			}
		 
		 }//isMC
  
  // Tout->Fill();
  
   return kTRUE;
}

void Anal_L5JERC::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
fileOut->cd();
fileOut->Write();

fOutput->Add(OutFile); 

fileOut->Close();
/*
fp<<"Total Number of events in "<<fileOut->GetName()<<" is "<<nevent_total<<endl;
fp<<"Total Weight is "<<tot_weight<<endl;

fOutput->Add(OutTxt); 
fp.close();
*/
//TString msg = TString::Format("Total Weight is  %f",tot_weight);
//gProofServ->SendAsynMessage(msg);
/*
fp<<"Total Number of events in "<<fileOut->GetName()<<" is "<<nevent_total<<endl;
fp<<"Total Weight is "<<tot_weight<<endl;
*/
}

void Anal_L5JERC::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
//std::cout<<"Total Number of events in "<<fileOut->GetName()<<" is "<<nevent_total<<endl;
}
