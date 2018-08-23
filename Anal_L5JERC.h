//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 31 15:28:52 2018 by ROOT version 5.34/36
// from TTree T1/WPrimeNtuple
// found on file: root://se01.indiacms.res.in//store/user/chatterj/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_crab_L5JERC_TTBar_SemiLeptonic_JECtxt/180725_134430/0000/rootuple_jerc_l5_106.root
//////////////////////////////////////////////////////////

#ifndef Anal_L5JERC_h
#define Anal_L5JERC_h

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
/*
#include "/home/chatterj/t3store/L5_JERC/roccor.2017.v0/RoccoR.cc"
RoccoR  rc("/home/chatterj/t3store/L5_JERC/roccor.2017.v0/RoccoR2017v0.txt");

double Roch_Cor(float Q, float pt, float eta, float phi)
{
double dtSF = rc.kScaleDT(Q, pt, eta, phi,0,0);
return dtSF;
}
*/
int* dec2bin(int dec)
{
const int length = 10;//nHLTmx;
int input = dec;
int istep=0;
int div[length] = {0};
while(input>0){
  input = input/2;
  div[istep] = input%2;
  istep++;
 }

static int binary[10]; //nHLTmx

for(int ij=1; ij<(length+1); ij++){
 binary[ij-1] = div[length-ij];
 }
return binary;
}

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


bool Muon_TightID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst,
				  float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
bool tightid = false;
if(muonisGL && muonisPF){
	if(muonchi<10 && muonhit>0 && muonmst>1){
		if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
			if(muonpixhit>0 && muontrklay>5){
				tightid = true;
				}
			} 
		}
	}	
return tightid;
}

bool Muon_Iso_ID(float muonpfiso)
{
bool isoid = false;	
if(muonpfiso<0.15) { isoid = true; }
return isoid;
}

bool ELectron_ID(float ele, float eleta, float elp, float elietaieta, float eletain, float elphiin, float elmisshits, 
				 float elhovere, float elconvdist, float rho){
bool tightid = false;
if(abs(eleta)<1.479){
	if(elietaieta<0.0104 && fabs(eletain)<0.00353 && (fabs(elphiin)<0.0499) && elmisshits<=1){
		if(fabs(1./ele - 1./elp)<0.0258){
			if(elhovere < (0.026 + 1.12/ele + 0.0368*rho/ele)){
				if(elconvdist>=0){
					tightid = true;
					}
				}
			} //elpfiso<0.0361
		}
	}
	else{
		if(elietaieta<0.0305 && abs(eletain)<0.00567 && (elphiin<0.0165) && elmisshits<=1){
			if(abs(1./ele - 1./elp)<0.0158){
				if(elhovere < (0.026 + 0.5/ele + 0.201*rho/ele)){
					if(elconvdist>=0){
						tightid = true;
						}
					}
				} //elpfiso<0.0361
			}
		}	
return tightid;	
}

float etaea[8] = {0,1.0,1.479,2.0,2.2,2.3,2.4,5};

bool Electron_Iso_ID(float eleta, float elpfiso, float Rho, int ieta)
{
bool isoid = false;	
elpfiso -= Rho*etaea[ieta]; if(elpfiso<0) { elpfiso = 0; }
if((elpfiso<0.0361 && abs(eleta)<1.479)||(elpfiso<0.094 && abs(eleta)>1.479)) { isoid = true; }
return isoid;
}

bool Photon_ID(float phoe, float phoeta, float phohadbyem, float phoietaieta, 
			   float phochhadiso, float phoneuhadiso, float phophoiso, float Rho, int ieta){
bool tightid = false;

phochhadiso -= Rho*etaea[ieta]; if(phochhadiso<0) { phochhadiso = 0; }
phoneuhadiso -= Rho*etaea[ieta]; if(phoneuhadiso<0) { phoneuhadiso = 0; }
phophoiso -= Rho*etaea[ieta]; if(phophoiso<0) { phophoiso = 0; }

float phopt = phoe *1./ cosh(phoeta) ;
if(abs(phoeta)<1.479){
	if(phohadbyem<0.02 && phoietaieta<0.0103){
		if((phochhadiso<1.158)&&(phoneuhadiso<(1.267 + 0.0126*phopt + 0.000026*phopt*phopt))&&(phophoiso<(2.065 + 0.0035*phopt))){
		tightid = true;
			}
		}
	}else{
		if(phohadbyem<0.025 && phoietaieta<0.0217){
		tightid = true;
		}
	}
	return tightid;
}

bool Photon_Iso_ID(float phoe, float phoeta, float phohadiso, float phophoiso)
{
bool isoid = false;
float phopt = phoe * sin(eta_to_theta(phoeta));
if(abs(phoeta)<1.479){
	if(phohadiso<(1.267 + 0.0126*phopt + 0.000026*phopt*phopt) && phophoiso<(2.065 + 0.0035*phopt)){
		isoid = true;
		}
	}else{
		if(phohadiso<(8.916 + 0.0119*phopt + 0.000025*phopt*phopt ) && phophoiso<(23.272 + 0.0040*phopt)){
			isoid = true;
			}
		}
		return isoid;
}

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Anal_L5JERC : public TSelector {
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
   Int_t           trig_value;
   Bool_t          ihlt01;
   Bool_t          ihlt02;
   Bool_t          ihlt03;
   Float_t		   prescl01;
   Float_t		   prescl02;
   Float_t		   prescl03;
   Float_t         PFMET;
   Float_t         PFMETPhi;
   Float_t         MisEtSig;
   Float_t         PFCHMET;
   Float_t         PFCHMETPhi;
   Float_t         CHMisEtSig;
   Float_t         sumEt;
  
   static const int njetmx =20;
   
   Int_t           npfjetAK8;
   Bool_t          pfjetAK8looseID[njetmx];   //[npfjetAK8]
   Bool_t          pfjetAK8tightID[njetmx];   //[npfjetAK8]
   Bool_t          pfjetAK8tightLVID[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8pt[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8y[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8phi[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8mass[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8JEC[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8JECL1[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8JECL2[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8JECL3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8JECL2L3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8btag[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8reso[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8resoup[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8resodn[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8moment1[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8moment2[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8moment3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8chrad[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8axis2[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8pTD[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8beta[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8betaStar[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sdmass[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8tau1[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8tau2[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8tau3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub1pt[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub1y[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub1phi[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub1mass[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub1btag[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub2pt[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub2y[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub2phi[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub2mass[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK8sub2btag[njetmx];   //[npfjetAK8]
   Int_t           pfjetAK8hadronflav[njetmx];   //[npfjetAK8]
   Int_t           pfjetAK8partonflav[njetmx];   //[npfjetAK8]
   Int_t           pfjetAK8partonpdg[njetmx];   //[npfjetAK8]
   Int_t           pfjetAK8qgl[njetmx];   //[npfjetAK8]
   Int_t           npfjetAK4;
   Bool_t          pfjetAK4looseID[njetmx];   //[npfjetAK4]
   Bool_t          pfjetAK4tightID[njetmx];   //[npfjetAK4]
   Bool_t          pfjetAK4tightLVID[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4pt[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4eta[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4y[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4phi[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4mass[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4JEC[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4JECL1[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK4JECL2[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK4JECL3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK4JECL2L3[njetmx];   //[npfjetAK8]
   Float_t         pfjetAK4btag[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4reso[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4resoup[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4resodn[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4moment1[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4moment2[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4moment3[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4chrad[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4axis2[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4pTD[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4beta[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4betaStar[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4sdmass[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4tau1[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4tau2[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4tau3[njetmx];   //[npfjetAK4]
   Int_t           pfjetAK4hadronflav[njetmx];   //[npfjetAK4]
   Int_t           pfjetAK4partonflav[njetmx];   //[npfjetAK4]
   Int_t           pfjetAK4partonpdg[njetmx];   //[npfjetAK4]
   Float_t         pfjetAK4qgl[njetmx];   //[npfjetAK4]
   Float_t		   pfjetAK4PUID[njetmx];   //[npfjetAK4]
   Float_t         GENMET;
   Float_t         GENMETPhi;
   Float_t         GENMETSig;
   Int_t           ngenjetAK8;
   Float_t         genjetAK8pt[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8y[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8phi[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8mass[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sdmass[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8btag[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8tau1[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8tau2[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8tau3[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8moment1[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8moment2[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8moment3[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8chrad[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8axis2[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8pTD[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub1pt[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub1y[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub1phi[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub1mass[njetmx];   //[ngenjetAK8]
   Float_t		   genjetAK8sub1btag[njetmx];
   Float_t         genjetAK8sub2pt[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub2y[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub2phi[njetmx];   //[ngenjetAK8]
   Float_t         genjetAK8sub2mass[njetmx];   //[ngenjetAK8]
   Float_t		   genjetAK8sub2btag[njetmx];
   Int_t           ngenjetAK4;
   Float_t         genjetAK4pt[njetmx];   //[ngenjetAK4]
   Float_t         genjetAK4y[njetmx];   //[ngenjetAK4]
   Float_t         genjetAK4phi[njetmx];   //[ngenjetAK4]
   Float_t         genjetAK4mass[njetmx];   //[ngenjetAK4]
   Float_t         genjetAK4sdmass[njetmx];   //[ngenjetAK4]
   Float_t         genjetAK4btag[njetmx];   //[ngenjetAK4]
   Int_t           ngenparticles;
   Int_t           genpartstatus[njetmx];   //[ngenparticles]
   Int_t           genpartpdg[njetmx];   //[ngenparticles]
   Int_t           genpartmompdg[njetmx];   //[ngenparticles]
   Int_t           genpartdaugno[njetmx];   //[ngenparticles]
   Bool_t          genpartfromhard[njetmx];   //[ngenparticles]
   Bool_t          genpartfromhardbFSR[njetmx];   //[ngenparticles]
   Float_t         genpartpt[njetmx];   //[ngenparticles]
   Float_t         genparteta[njetmx];   //[ngenparticles]
   Float_t         genpartphi[njetmx];   //[ngenparticles]
   Float_t         genpartm[njetmx];   //[ngenparticles]
   Float_t         genpartq[njetmx];   //[ngenparticles]
   Int_t           nmuons;
   Bool_t          muonisPF[njetmx];   //[nmuons]
   Bool_t          muonisGL[njetmx];   //[nmuons]
   Bool_t          muonisTRK[njetmx];   //[nmuons]
   Bool_t          muonisLoose[njetmx];   //[nmuons]
   Bool_t          muonisGoodGL[njetmx];   //[nmuons]
   Bool_t          muonisMed[njetmx];   //[nmuons]
   Bool_t 		   muonistight[njetmx];  //[nmuons]
   Float_t         muonpt[njetmx];   //[nmuons]
   Float_t         muonp[njetmx];   //[nmuons]
   Float_t         muone[njetmx];   //[nmuons]
   Float_t         muoneta[njetmx];   //[nmuons]
   Float_t         muonphi[njetmx];   //[nmuons]
   Int_t		   muonq[njetmx];   //[nmuons]
   Float_t         muondrbm[njetmx];   //[nmuons]
   Float_t         muontrkvtx[njetmx];   //[nmuons]
   Float_t         muondz[njetmx];   //[nmuons]
   Float_t         muonpter[njetmx];   //[nmuons]
   Float_t         muonchi[njetmx];   //[nmuons]
   Int_t           muonndf[njetmx];   //[nmuons]
   Float_t         muonecal[njetmx];   //[nmuons]
   Float_t         muonhcal[njetmx];   //[nmuons]
   Float_t		   muonhit[njetmx]; 
   Float_t         muonemiso[njetmx];   //[nmuons]
   Float_t         muonhadiso[njetmx];   //[nmuons]
   Float_t         muonpfiso[njetmx];   //[nmuons]
   Float_t         muontkpt03[njetmx];   //[nmuons]
   Float_t         muontkpt05[njetmx];   //[nmuons]
   Float_t         muonposmatch[njetmx];   //[nmuons]
   Float_t         muontrkink[njetmx];   //[nmuons]
   Float_t         muonsegcom[njetmx];   //[nmuons]
   Float_t         muonpixhit[njetmx];   //[nmuons]
   Float_t         muonmst[njetmx];   //[nmuons]
   Float_t         muontrklay[njetmx];   //[nmuons]
   Float_t         muonvalfrac[njetmx];   //[nmuons]
   Int_t           nelecs;
   Bool_t 		   elecistight[njetmx];   //[nelecs]
   Float_t         elpt[njetmx];   //[nelecs]
   Float_t         eleta[njetmx];   //[nelecs]
   Float_t         elphi[njetmx];   //[nelecs]
   Float_t         elp[njetmx];   //[nelecs]
   Float_t         ele[njetmx];   //[nelecs]
   Float_t		   elsceta[njetmx];   //[nelecs]
   Bool_t          elmvaid[njetmx];   //[nelecs]
   Float_t         eldxy[njetmx];   //[nelecs]
   Float_t         eldz[njetmx];   //[nelecs]
   Float_t         elhovere[njetmx];   //[nelecs]
   Float_t         elchi[njetmx];   //[nelecs]
   Int_t           elndf[njetmx];   //[nelecs]
   Float_t         eltkpt03[njetmx];   //[nelecs]
   Float_t         eltkpt04[njetmx];   //[nelecs]
   Float_t         elemiso04[njetmx];   //[nelecs]
   Float_t         elhadiso04[njetmx];   //[nelecs]
   Float_t         eletain[njetmx];   //[nelecs]
   Float_t         elphiin[njetmx];   //[nelecs]
   Float_t         elfbrem[njetmx];   //[nelecs]
   Float_t         elhadisodepth03[njetmx];   //[nelecs]
   Float_t         eleoverp[njetmx];   //[nelecs]
   Float_t         elietaieta[njetmx];   //[nelecs]
   Float_t         elmisshits[njetmx];   //[nelecs]
   Float_t         elchhadiso[njetmx];   //[nelecs]
   Float_t         elneuhadiso[njetmx];   //[nelecs]
   Float_t         elphoiso[njetmx];   //[nelecs]
   Float_t         elpuchhadiso[njetmx];   //[nelecs]
   Float_t         elpfiso[njetmx];   //[nelecs]
   Float_t         elconvdist[njetmx];   //[nelecs]
   Float_t         elconvdoct[njetmx];   //[nelecs]
   Int_t           nphotons;
   Float_t         phoe[njetmx];   //[nphotons]
   Float_t         phoeta[njetmx];   //[nphotons]
   Float_t         phophi[njetmx];   //[nphotons]
   Bool_t          phomvaid[njetmx];   //[nphotons]
   Float_t         phoe1by9[njetmx];   //[nphotons]
   Float_t         phoe9by25[njetmx];   //[nphotons]
   Float_t         photrkiso[njetmx];   //[nphotons]
   Float_t         phoemiso[njetmx];   //[nphotons]
   Float_t         phohadiso[njetmx];   //[nphotons]
   Float_t         phochhadiso[njetmx];   //[nphotons]
   Float_t         phoneuhadiso[njetmx];   //[nphotons]
   Float_t         phophoiso[njetmx];   //[nphotons]
   Float_t         phoPUiso[njetmx];   //[nphotons]
   Float_t         phohadbyem[njetmx];   //[nphotons]
   Float_t         phoietaieta[njetmx];   //[nphotons]

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
   TBranch        *b_trig_value;   //!
   TBranch        *b_ihlt01;   //!
   TBranch        *b_ihlt02;   //!
   TBranch        *b_ihlt03;   //!
   TBranch        *b_prescl01;   //!
   TBranch        *b_prescl02;   //!
   TBranch        *b_prescl03;   //!
   TBranch        *b_miset;   //!
   TBranch        *b_misphi;   //!
   TBranch        *b_misetsig;   //!
   TBranch        *b_chmiset;   //!
   TBranch        *b_chmisphi;   //!
   TBranch        *b_chmisetsig;   //!
   TBranch        *b_sumEt;   //!
   TBranch        *b_npfjetAK8;   //!
   TBranch        *b_pfjetAK8looseID;   //!
   TBranch        *b_pfjetAK8tightID;   //!
   TBranch        *b_pfjetAK8tightLVID;   //!
   TBranch        *b_pfjetAK8pt;   //!
   TBranch        *b_pfjetAK8y;   //!
   TBranch        *b_pfjetAK8phi;   //!
   TBranch        *b_pfjetAK8mass;   //!
   TBranch        *b_pfjetAK8JEC;   //!
   TBranch        *b_pfjetAK8JECL1;   //!
   TBranch        *b_pfjetAK8JECL2;   //!
   TBranch        *b_pfjetAK8JECL3;   //!
   TBranch        *b_pfjetAK8JECL2L3;   //!
   TBranch        *b_pfjetAK8btag;   //!
   TBranch        *b_pfjetAK8reso;   //!
   TBranch        *b_pfjetAK8resoup;   //!
   TBranch        *b_pfjetAK8resodn;   //!
   TBranch        *b_pfjetAK8moment1;   //!
   TBranch        *b_pfjetAK8moment2;   //!
   TBranch        *b_pfjetAK8moment3;   //!
   TBranch        *b_pfjetAK8chrad;   //!
   TBranch        *b_pfjetAK8axis2;   //!
   TBranch        *b_pfjetAK8pTD;   //!
   TBranch        *b_pfjetAK8beta;   //!
   TBranch        *b_pfjetAK8betaStar;   //!
   TBranch        *b_pfjetAK8sdmass;   //!
   TBranch        *b_pfjetAK8tau1;   //!
   TBranch        *b_pfjetAK8tau2;   //!
   TBranch        *b_pfjetAK8tau3;   //!
   TBranch        *b_pfjetAK8sub1pt;   //!
   TBranch        *b_pfjetAK8sub1y;   //!
   TBranch        *b_pfjetAK8sub1phi;   //!
   TBranch        *b_pfjetAK8sub1mass;   //!
   TBranch        *b_pfjetAK8sub1btag;   //!
   TBranch        *b_pfjetAK8sub2pt;   //!
   TBranch        *b_pfjetAK8sub2y;   //!
   TBranch        *b_pfjetAK8sub2phi;   //!
   TBranch        *b_pfjetAK8sub2mass;   //!
   TBranch        *b_pfjetAK8sub2btag;   //!
   TBranch        *b_pfjetAK8hadronflav;   //!
   TBranch        *b_pfjetAK8partonflav;   //!
   TBranch        *b_pfjetAK8partonpdg;   //!
   TBranch        *b_npfjetAK4;   //!
   TBranch        *b_pfjetAK4looseID;   //!
   TBranch        *b_pfjetAK4tightID;   //!
   TBranch        *b_pfjetAK4tightLVID;   //!
   TBranch        *b_pfjetAK4pt;   //!
   TBranch        *b_pfjetAK4eta;   //!
   TBranch        *b_pfjetAK4y;   //!
   TBranch        *b_pfjetAK4phi;   //!
   TBranch        *b_pfjetAK4mass;   //!
   TBranch        *b_pfjetAK4JEC;   //!
   TBranch        *b_pfjetAK4JECL1;   //!
   TBranch        *b_pfjetAK4JECL2;   //!
   TBranch        *b_pfjetAK4JECL3;   //!
   TBranch        *b_pfjetAK4JECL2L3;   //!
   TBranch        *b_pfjetAK4btag;   //!
   TBranch        *b_pfjetAK4reso;   //!
   TBranch        *b_pfjetAK4resoup;   //!
   TBranch        *b_pfjetAK4resodn;   //!
   TBranch        *b_pfjetAK4moment1;   //!
   TBranch        *b_pfjetAK4moment2;   //!
   TBranch        *b_pfjetAK4moment3;   //!
   TBranch        *b_pfjetAK4chrad;   //!
   TBranch        *b_pfjetAK4axis2;   //!
   TBranch        *b_pfjetAK4pTD;   //!
   TBranch        *b_pfjetAK4beta;   //!
   TBranch        *b_pfjetAK4betaStar;   //!
   TBranch        *b_pfjetAK4sdmass;   //!
   TBranch        *b_pfjetAK4tau1;   //!
   TBranch        *b_pfjetAK4tau2;   //!
   TBranch        *b_pfjetAK4tau3;   //!
   TBranch        *b_pfjetAK4hadronflav;   //!
   TBranch        *b_pfjetAK4partonflav;   //!
   TBranch        *b_pfjetAK4partonpdg;   //!
   TBranch        *b_pfjetAK4qgl;   //!
   TBranch        *b_pfjetAK4PUID;   //!
   TBranch        *b_genmiset;   //!
   TBranch        *b_genmisphi;   //!
   TBranch        *b_genmisetsig;   //!
   TBranch        *b_ngenjetAK8;   //!
   TBranch        *b_genjetAK8pt;   //!
   TBranch        *b_genjetAK8y;   //!
   TBranch        *b_genjetAK8phi;   //!
   TBranch        *b_genjetAK8mass;   //!
   TBranch        *b_genjetAK8sdmass;   //!
   TBranch        *b_genjetAK8btag;   //!
   TBranch        *b_genjetAK8tau1;   //!
   TBranch        *b_genjetAK8tau2;   //!
   TBranch        *b_genjetAK8tau3;   //!
   TBranch        *b_genjetAK8moment1;   //!
   TBranch        *b_genjetAK8moment2;   //!
   TBranch        *b_genjetAK8moment3;   //!
   TBranch        *b_genjetAK8chrad;   //!
   TBranch        *b_genjetAK8axis2;   //!
   TBranch        *b_genjetAK8pTD;   //!
   TBranch        *b_genjetAK8sub1pt;   //!
   TBranch        *b_genjetAK8sub1y;   //!
   TBranch        *b_genjetAK8sub1phi;   //!
   TBranch        *b_genjetAK8sub1mass;   //!
   TBranch        *b_genjetAK8sub2pt;   //!
   TBranch        *b_genjetAK8sub2y;   //!
   TBranch        *b_genjetAK8sub2phi;   //!
   TBranch        *b_genjetAK8sub2mass;   //!
   TBranch        *b_ngenjetAK4;   //!
   TBranch        *b_genjetAK4pt;   //!
   TBranch        *b_genjetAK4y;   //!
   TBranch        *b_genjetAK4phi;   //!
   TBranch        *b_genjetAK4mass;   //!
   TBranch        *b_genjetAK4sdmass;   //!
   TBranch        *b_genjetAK4btag;   //!
   TBranch        *b_ngenparticles;   //!
   TBranch        *b_genpartstatus;   //!
   TBranch        *b_genpartpdg;   //!
   TBranch        *b_genpartmompdg;   //!
   TBranch        *b_genpartdaugno;   //!
   TBranch        *b_genpartfromhard;   //!
   TBranch        *b_genpartfromhardbFSR;   //!
   TBranch        *b_genpartpt;   //!
   TBranch        *b_genparteta;   //!
   TBranch        *b_genpartphi;   //!
   TBranch        *b_genpartm;   //!
   TBranch        *b_genpartq;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_muonisPF;   //!
   TBranch        *b_muonisGL;   //!
   TBranch        *b_muonisTRK;   //!
   TBranch        *b_muonisLoose;   //!
   TBranch        *b_muonisGoodGL;   //!
   TBranch        *b_muonisMed;   //!
   TBranch        *b_muonpt;   //!
   TBranch        *b_muonp;   //!
   TBranch        *b_muone;   //!
   TBranch        *b_muoneta;   //!
   TBranch        *b_muonphi;   //!
   TBranch        *b_muondrbm;   //!
   TBranch        *b_muontrkvtx;   //!
   TBranch        *b_muondz;   //!
   TBranch        *b_muonpter;   //!
   TBranch        *b_muonchi;   //!
   TBranch        *b_muonndf;   //!
   TBranch        *b_muonecal;   //!
   TBranch        *b_muonhcal;   //!
   TBranch        *b_muonemiso;   //!
   TBranch        *b_muonhadiso;   //!
   TBranch        *b_muonpfiso;   //!
   TBranch        *b_muontkpt03;   //!
   TBranch        *b_muontkpt05;   //!
   TBranch        *b_muonposmatch;   //!
   TBranch        *b_muontrkink;   //!
   TBranch        *b_muonsegcom;   //!
   TBranch        *b_muonhit;   //!
   TBranch        *b_muonpixhit;   //!
   TBranch        *b_muonmst;   //!
   TBranch        *b_muontrklay;   //!
   TBranch        *b_muonvalfrac;   //!
   TBranch        *b_nelecs;   //!
   TBranch        *b_elpt;   //!
   TBranch        *b_eleta;   //!
   TBranch        *b_elphi;   //!
   TBranch        *b_elp;   //!
   TBranch        *b_ele;   //!
   TBranch        *b_elmvaid;   //!
   TBranch        *b_eldxy;   //!
   TBranch        *b_eldz;   //!
   TBranch        *b_elhovere;   //!
   TBranch        *b_elchi;   //!
   TBranch        *b_elndf;   //!
   TBranch        *b_eltkpt03;   //!
   TBranch        *b_eltkpt04;   //!
   TBranch        *b_elemiso04;   //!
   TBranch        *b_elhadiso04;   //!
   TBranch        *b_eletain;   //!
   TBranch        *b_elphiin;   //!
   TBranch        *b_elfbrem;   //!
   TBranch        *b_elhadisodepth03;   //!
   TBranch        *b_eleoverp;   //!
   TBranch        *b_elietaieta;   //!
   TBranch        *b_elmisshits;   //!
   TBranch        *b_elchhadiso;   //!
   TBranch        *b_elneuhadiso;   //!
   TBranch        *b_elphoiso;   //!
   TBranch        *b_elpuchhadiso;   //!
   TBranch        *b_elpfiso;   //!
   TBranch        *b_elconvdist;   //!
   TBranch        *b_elconvdoct;   //!
   TBranch        *b_nphotons;   //!
   TBranch        *b_phoe;   //!
   TBranch        *b_phoeta;   //!
   TBranch        *b_phophi;   //!
   TBranch        *b_phomvaid;   //!
   TBranch        *b_phoe1by9;   //!
   TBranch        *b_phoe9by25;   //!
   TBranch        *b_photrkiso;   //!
   TBranch        *b_phoemiso;   //!
   TBranch        *b_phohadiso;   //!
   TBranch        *b_phochhadiso;   //!
   TBranch        *b_phoneuhadiso;   //!
   TBranch        *b_phophoiso;   //!
   TBranch        *b_phoPUiso;   //!
   TBranch        *b_phohadbyem;   //!
   TBranch        *b_phoietaieta;   //!

   Anal_L5JERC(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Anal_L5JERC() { }
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
   
   static const int npvmax = 100;
   
   float min_wmass = 0; 
   float max_wmass = 200;
   static const int wmassbins = 100;
   float min_tmass = 0;
   float max_tmass = 500.;
   static const int tmassbins = 250;

   float ptcut = 25.;
   float ycut = 2.4;
   float brem_rad = 0.1;
   float btagvalue = 0.8838;
//   float weight;
   
   float ydiv = 1.4;
 //  float ydivs[selneta+1] = {0,0.5,1.0,1.5,2.0,2.5};
   
//float seleta[seleta+1] = {1.4,2.5,5.};
  float selmass[selnmass+1] = {220,270,320,370};
  
   static const int nHLTmx= 3; 
   
   const float pival = acos(-1.);
   
   TProofOutputFile *OutFile;
   TFile *fileOut; 
//  TProofOutputFile *OutTxt;
   ofstream fp ;

   TTree *Tout ;

   TH1D *h_nvert;
   TH1D *h_ndofct;
   TH1D *h_nchict;
   TH1D *h_npuvert;

   TH1D *h_jetnAK8[netarange];
   TH1D *h_jetptAK8[netarange];
   TH1D *h_jetmassAK8[netarange];
   TH1D *h_jetsdmassAK8[netarange];  
   TH1D *h_jetbtagAK8[netarange]; 
   TH1D* h_jettau21AK8[netarange]; 
   TH1D* h_jettau32AK8[netarange]; 
   TH1D* h_jetmassAK8_sel[netarange];
   TH1D*	h_jetsdmassAK8_sel[netarange];  
   TH1D* h_jetmassAK8_selpt[netarange];
   TH1D*	h_jetsdmassAK8_selpt[netarange];  
   TH1D* h_jetyAK8;
   TH1D* h_jetphiAK8;
   
   TH1D *h_jetnAK4[netarange];
   TH1D *h_jetptAK4[netarange];
   TH1D *h_jetmassAK4[netarange];
   TH1D *h_jetsdmassAK4[netarange];  
   TH1D *h_jetbtagAK4[netarange]; 
   TH1D* h_jettau21AK4[netarange]; 
   TH1D* h_jettau32AK4[netarange]; 
   TH1D* h_jetyAK4;
   TH1D* h_jetphiAK4; 

   TH1D *h_nmuon_pre;
   TH1D *h_idmuon_pre;
   TH1D *h_nelec_pre;
   TH1D *h_idelec_pre;
   TH1D *h_npho_pre;
   TH1D *h_npfjetAK4_pre;
   TH1D *h_nbjetAK4_pre;

   TH1D *h_jetnqbAK4;
   TH1D *h_jetnqAK4;
   TH1D *h_jetnbAK4;
   TH1D *h_nmuon;
  
   TH1D *h_jetpt_q1AK4;
   TH1D *h_jeteta_q1AK4;
   TH1D *h_jetphi_q1AK4;
  
   TH1D *h_jetpt_q2AK4;
   TH1D *h_jeteta_q2AK4;
   TH1D *h_jetphi_q2AK4;

   TH1D *h_jetpt_b1AK4;
   TH1D *h_jeteta_b1AK4;
   TH1D *h_jetphi_b1AK4;
  
   TH1D *h_jetpt_b2AK4;
   TH1D *h_jeteta_b2AK4;
   TH1D *h_jetphi_b2AK4;

   TH1D *h_metpt;
   TH1D *h_metphi;
   TH1D *h_metbyEt;
   TH1D *h_metSig;

   TH1D *h_lambda;
   TH1D *h_lambdaA;
   TH1D *h_lambda_truth;

   TH1D *h_mun[netarange];
   TH1D *h_mupt[netarange];
   TH1D *h_mueta;
   TH1D *h_muphi;
   TH1D *h_mudrbm;
   TH1D *h_mutrkvtx;
   TH1D *h_muhbye;
   TH1D *h_muemiso;
   TH1D *h_muhadiso;
   TH1D *h_mupfiso;
   TH1D *h_mukpt03;
   TH1D *h_mukpt05;
   TH1D *h_muhit;
   TH1D *h_mumst;
   TH1D *h_mutrklay;

   TH1D *h_mupt_pass;
   TH1D *h_mueta_pass;
   TH1D *h_muphi_pass;

   TH1D *h_whadmass;
   TH1D *h_whadpt;
   TH1D *h_whady;
  
   TH1D *h_thadmass;
   TH1D *h_thadpt;
   TH1D *h_thady;
   TH2D *h_tw_2d;
  
   TH1D *h_whadmassA;
   TH1D *h_whadptA;
   TH1D *h_whadyA;
  
   TH1D *h_thadmassA;
   TH1D *h_thadptA;
   TH1D *h_thadyA;
   TH2D *h_tw_2dA;
  
   TH1D *h_whadmass_truth;
   TH1D *h_whadpt_truth;
   TH1D *h_whady_truth;
  
   TH1D *h_thadmass_truth;
   TH1D *h_thadpt_truth;
   TH1D *h_thady_truth;
   TH2D *h_tw_2d_truth;

   TH1D *h_whadmass_wrong;
   TH1D *h_whadpt_wrong;
   TH1D *h_whady_wrong;
  
   TH1D *h_thadmass_wrong;
   TH1D *h_thadpt_wrong;
   TH1D *h_thady_wrong;
   TH2D *h_tw_2d_wrong;
  
   TH1D *h_whadmass_unmatch;
   TH1D *h_whadpt_unmatch;
   TH1D *h_whady_unmatch;
  
   TH1D *h_thadmass_unmatch;
   TH1D *h_thadpt_unmatch;
   TH1D *h_thady_unmatch;
   TH2D *h_tw_2d_unmatch;
  
   TH1D *h_whadmass_parton;
   TH1D *h_thadmass_parton;

   TH1D *h_whadmass_sel[11][selnmass][selneta];
   TH1D *h_whadmass_truth_sel[11][selnmass][selneta];
   TH1D *h_whadmass_wrong_sel[11][selnmass][selneta];
   TH1D *h_whadmass_unmatch_sel[11][selnmass][selneta];
   
   TH1D *h_jety_2y[2][selneta];
   TH1D *h_jety_2pt[2][selneta];

   TH1D* h_resopt_AK4[netarange][noptbins];
   TH1D* h_massdijet_AK4[netarange][wmassbins];

   int nbjetAK4;
   float bjetAK4pt[njetmx], bjetAK4eta[njetmx], bjetAK4phi[njetmx], bjetAK4btag[njetmx], bjetAK4mass[njetmx], bjetAK4hadronflav[njetmx], bjetAK4partonflav[njetmx];

   int nqjetAK4;
   float qjetAK4pt[njetmx], qjetAK4eta[njetmx], qjetAK4phi[njetmx], qjetAK4btag[njetmx], qjetAK4mass[njetmx], qjetAK4hadronflav[njetmx], qjetAK4partonflav[njetmx];

   int nbgenjetAK4;
   float bgenjetAK4pt[njetmx], bgenjetAK4y[njetmx], bgenjetAK4phi[njetmx], bgenjetAK4btag[njetmx], bgenjetAK4mass[njetmx], bgenjetAK4sdmass[njetmx];

   int nqgenjetAK4;
   float qgenjetAK4pt[njetmx], qgenjetAK4y[njetmx], qgenjetAK4phi[njetmx], qgenjetAK4btag[njetmx], qgenjetAK4mass[njetmx], qgenjetAK4sdmass[njetmx];
 
   bool ihlt[nHLTmx];
   float prescl[nHLTmx];
    
  float thad_pt, thad_eta, thad_phi, thad_e, thad_m;
  float whad_pt, whad_eta, whad_phi, whad_e, whad_m;
  float tlep_pt, tlep_eta, tlep_phi, tlep_e, tlep_m;
  float wlep_pt, wlep_eta, wlep_phi, wlep_e, wlep_m;
  float bhad_pt, bhad_eta, bhad_phi, bhad_e, bhad_m;
  float blep_pt, blep_eta, blep_phi, blep_e, blep_m;
  float plep_pt, plep_eta, plep_phi, plep_e;
  float pq1_pt, pq1_eta, pq1_phi, pq1_e;
  float pq2_pt, pq2_eta, pq2_phi, pq2_e;

  float delR_whadbhad, delEta_whadbhad, delPhi_whadbhad;
  float delR_whadblep, delEta_whadblep, delPhi_whadblep;
  float delR_whadplep, delEta_whadplep, delPhi_whadplep;
  float delR_bhadplep, delEta_bhadplep, delPhi_bhadplep;
  float delR_blepplep, delEta_blepplep, delPhi_blepplep;
 
  float weight, weightev;
  
   double tot_weight = 0;  
   int nevent_total = 0;
   
   char name[100];
   char title[100];
   
   bool isMC;

   ClassDef(Anal_L5JERC,0);
};

#endif

#ifdef Anal_L5JERC_cxx
void Anal_L5JERC::Init(TTree *tree)
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
   fChain->SetBranchAddress("trig_value", &trig_value, &b_trig_value);
   fChain->SetBranchAddress("ihlt01", &ihlt01, &b_ihlt01);
   fChain->SetBranchAddress("ihlt02", &ihlt02, &b_ihlt02);
   fChain->SetBranchAddress("ihlt03", &ihlt03, &b_ihlt03);
   fChain->SetBranchAddress("prescl01", &prescl01, &b_prescl01);
   fChain->SetBranchAddress("prescl02", &prescl02, &b_prescl02);
   fChain->SetBranchAddress("prescl03", &prescl03, &b_prescl03);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_miset);
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_misphi);
   fChain->SetBranchAddress("MisEtSig", &MisEtSig, &b_misetsig);
   fChain->SetBranchAddress("PFCHMET", &PFCHMET, &b_chmiset);
   fChain->SetBranchAddress("PFCHMETPhi", &PFCHMETPhi, &b_chmisphi);
   fChain->SetBranchAddress("CHMisEtSig", &CHMisEtSig, &b_chmisetsig);
   fChain->SetBranchAddress("sumEt", &sumEt, &b_sumEt);
   fChain->SetBranchAddress("npfjetAK8", &npfjetAK8, &b_npfjetAK8);
   fChain->SetBranchAddress("pfjetAK8looseID", pfjetAK8looseID, &b_pfjetAK8looseID);
   fChain->SetBranchAddress("pfjetAK8tightID", pfjetAK8tightID, &b_pfjetAK8tightID);
   fChain->SetBranchAddress("pfjetAK8tightLVID", pfjetAK8tightLVID, &b_pfjetAK8tightLVID);
   fChain->SetBranchAddress("pfjetAK8pt", pfjetAK8pt, &b_pfjetAK8pt);
   fChain->SetBranchAddress("pfjetAK8y", pfjetAK8y, &b_pfjetAK8y);
   fChain->SetBranchAddress("pfjetAK8phi", pfjetAK8phi, &b_pfjetAK8phi);
   fChain->SetBranchAddress("pfjetAK8mass", pfjetAK8mass, &b_pfjetAK8mass);
   fChain->SetBranchAddress("pfjetAK8JEC", pfjetAK8JEC, &b_pfjetAK8JEC);
   fChain->SetBranchAddress("pfjetAK8JECL1", pfjetAK8JECL1, &b_pfjetAK8JECL1);
   fChain->SetBranchAddress("pfjetAK8JECL2", pfjetAK8JECL2, &b_pfjetAK8JECL2);
   fChain->SetBranchAddress("pfjetAK8JECL3", pfjetAK8JECL3, &b_pfjetAK8JECL3);
   fChain->SetBranchAddress("pfjetAK8JECL2L3", pfjetAK8JECL2L3, &b_pfjetAK8JECL2L3);
   fChain->SetBranchAddress("pfjetAK8btag", pfjetAK8btag, &b_pfjetAK8btag);
   fChain->SetBranchAddress("pfjetAK8reso", pfjetAK8reso, &b_pfjetAK8reso);
   fChain->SetBranchAddress("pfjetAK8resoup", pfjetAK8resoup, &b_pfjetAK8resoup);
   fChain->SetBranchAddress("pfjetAK8resodn", pfjetAK8resodn, &b_pfjetAK8resodn);
   fChain->SetBranchAddress("pfjetAK8moment1", pfjetAK8moment1, &b_pfjetAK8moment1);
   fChain->SetBranchAddress("pfjetAK8moment2", pfjetAK8moment2, &b_pfjetAK8moment2);
   fChain->SetBranchAddress("pfjetAK8moment3", pfjetAK8moment3, &b_pfjetAK8moment3);
   fChain->SetBranchAddress("pfjetAK8chrad", pfjetAK8chrad, &b_pfjetAK8chrad);
   fChain->SetBranchAddress("pfjetAK8axis2", pfjetAK8axis2, &b_pfjetAK8axis2);
   fChain->SetBranchAddress("pfjetAK8pTD", pfjetAK8pTD, &b_pfjetAK8pTD);
   fChain->SetBranchAddress("pfjetAK8beta", pfjetAK8beta, &b_pfjetAK8beta);
   fChain->SetBranchAddress("pfjetAK8betaStar", pfjetAK8betaStar, &b_pfjetAK8betaStar);
   fChain->SetBranchAddress("pfjetAK8sdmass", pfjetAK8sdmass, &b_pfjetAK8sdmass);
   fChain->SetBranchAddress("pfjetAK8tau1", pfjetAK8tau1, &b_pfjetAK8tau1);
   fChain->SetBranchAddress("pfjetAK8tau2", pfjetAK8tau2, &b_pfjetAK8tau2);
   fChain->SetBranchAddress("pfjetAK8tau3", pfjetAK8tau3, &b_pfjetAK8tau3);
   fChain->SetBranchAddress("pfjetAK8sub1pt", pfjetAK8sub1pt, &b_pfjetAK8sub1pt);
   fChain->SetBranchAddress("pfjetAK8sub1y", pfjetAK8sub1y, &b_pfjetAK8sub1y);
   fChain->SetBranchAddress("pfjetAK8sub1phi", pfjetAK8sub1phi, &b_pfjetAK8sub1phi);
   fChain->SetBranchAddress("pfjetAK8sub1mass", pfjetAK8sub1mass, &b_pfjetAK8sub1mass);
   fChain->SetBranchAddress("pfjetAK8sub1btag", pfjetAK8sub1btag, &b_pfjetAK8sub1btag);
   fChain->SetBranchAddress("pfjetAK8sub2pt", pfjetAK8sub2pt, &b_pfjetAK8sub2pt);
   fChain->SetBranchAddress("pfjetAK8sub2y", pfjetAK8sub2y, &b_pfjetAK8sub2y);
   fChain->SetBranchAddress("pfjetAK8sub2phi", pfjetAK8sub2phi, &b_pfjetAK8sub2phi);
   fChain->SetBranchAddress("pfjetAK8sub2mass", pfjetAK8sub2mass, &b_pfjetAK8sub2mass);
   fChain->SetBranchAddress("pfjetAK8sub2btag", pfjetAK8sub2btag, &b_pfjetAK8sub2btag);
   fChain->SetBranchAddress("pfjetAK8hadronflav", pfjetAK8hadronflav, &b_pfjetAK8hadronflav);
   fChain->SetBranchAddress("pfjetAK8partonflav", pfjetAK8partonflav, &b_pfjetAK8partonflav);
   fChain->SetBranchAddress("pfjetAK8partonpdg", pfjetAK8partonpdg, &b_pfjetAK8partonpdg);
   fChain->SetBranchAddress("npfjetAK4", &npfjetAK4, &b_npfjetAK4);
   fChain->SetBranchAddress("pfjetAK4looseID", pfjetAK4looseID, &b_pfjetAK4looseID);
   fChain->SetBranchAddress("pfjetAK4tightID", pfjetAK4tightID, &b_pfjetAK4tightID);
   fChain->SetBranchAddress("pfjetAK4tightLVID", pfjetAK4tightLVID, &b_pfjetAK4tightLVID);
   fChain->SetBranchAddress("pfjetAK4pt", pfjetAK4pt, &b_pfjetAK4pt);
   fChain->SetBranchAddress("pfjetAK4eta", pfjetAK4eta, &b_pfjetAK4eta);
   fChain->SetBranchAddress("pfjetAK4y", pfjetAK4y, &b_pfjetAK4y);
   fChain->SetBranchAddress("pfjetAK4phi", pfjetAK4phi, &b_pfjetAK4phi);
   fChain->SetBranchAddress("pfjetAK4mass", pfjetAK4mass, &b_pfjetAK4mass);
   fChain->SetBranchAddress("pfjetAK4JEC", pfjetAK4JEC, &b_pfjetAK4JEC);
   fChain->SetBranchAddress("pfjetAK4JECL1", pfjetAK4JECL1, &b_pfjetAK4JECL1);
   fChain->SetBranchAddress("pfjetAK4JECL2", pfjetAK4JECL2, &b_pfjetAK4JECL2);
   fChain->SetBranchAddress("pfjetAK4JECL3", pfjetAK4JECL3, &b_pfjetAK4JECL3);
   fChain->SetBranchAddress("pfjetAK4JECL2L3", pfjetAK4JECL2L3, &b_pfjetAK4JECL2L3);
   fChain->SetBranchAddress("pfjetAK4btag", pfjetAK4btag, &b_pfjetAK4btag);
   fChain->SetBranchAddress("pfjetAK4reso", pfjetAK4reso, &b_pfjetAK4reso);
   fChain->SetBranchAddress("pfjetAK4resoup", pfjetAK4resoup, &b_pfjetAK4resoup);
   fChain->SetBranchAddress("pfjetAK4resodn", pfjetAK4resodn, &b_pfjetAK4resodn);
   fChain->SetBranchAddress("pfjetAK4moment1", pfjetAK4moment1, &b_pfjetAK4moment1);
   fChain->SetBranchAddress("pfjetAK4moment2", pfjetAK4moment2, &b_pfjetAK4moment2);
   fChain->SetBranchAddress("pfjetAK4moment3", pfjetAK4moment3, &b_pfjetAK4moment3);
   fChain->SetBranchAddress("pfjetAK4chrad", pfjetAK4chrad, &b_pfjetAK4chrad);
   fChain->SetBranchAddress("pfjetAK4axis2", pfjetAK4axis2, &b_pfjetAK4axis2);
   fChain->SetBranchAddress("pfjetAK4pTD", pfjetAK4pTD, &b_pfjetAK4pTD);
   fChain->SetBranchAddress("pfjetAK4beta", pfjetAK4beta, &b_pfjetAK4beta);
   fChain->SetBranchAddress("pfjetAK4betaStar", pfjetAK4betaStar, &b_pfjetAK4betaStar);
   fChain->SetBranchAddress("pfjetAK4sdmass", pfjetAK4sdmass, &b_pfjetAK4sdmass);
   fChain->SetBranchAddress("pfjetAK4tau1", pfjetAK4tau1, &b_pfjetAK4tau1);
   fChain->SetBranchAddress("pfjetAK4tau2", pfjetAK4tau2, &b_pfjetAK4tau2);
   fChain->SetBranchAddress("pfjetAK4tau3", pfjetAK4tau3, &b_pfjetAK4tau3);
   fChain->SetBranchAddress("pfjetAK4hadronflav", pfjetAK4hadronflav, &b_pfjetAK4hadronflav);
   fChain->SetBranchAddress("pfjetAK4partonflav", pfjetAK4partonflav, &b_pfjetAK4partonflav);
   fChain->SetBranchAddress("pfjetAK4partonpdg", pfjetAK4partonpdg, &b_pfjetAK4partonpdg);
   fChain->SetBranchAddress("pfjetAK4qgl", pfjetAK4qgl, &b_pfjetAK4qgl);
   fChain->SetBranchAddress("pfjetAK4PUID", pfjetAK4PUID, &b_pfjetAK4PUID);
   if(isMC){
   fChain->SetBranchAddress("GENMET", &GENMET, &b_genmiset);
   fChain->SetBranchAddress("GENMETPhi", &GENMETPhi, &b_genmisphi);
   fChain->SetBranchAddress("GENMETSig", &GENMETSig, &b_genmisetsig);
   fChain->SetBranchAddress("ngenjetAK8", &ngenjetAK8, &b_ngenjetAK8);
   fChain->SetBranchAddress("genjetAK8pt", genjetAK8pt, &b_genjetAK8pt);
   fChain->SetBranchAddress("genjetAK8y", genjetAK8y, &b_genjetAK8y);
   fChain->SetBranchAddress("genjetAK8phi", genjetAK8phi, &b_genjetAK8phi);
   fChain->SetBranchAddress("genjetAK8mass", genjetAK8mass, &b_genjetAK8mass);
   fChain->SetBranchAddress("genjetAK8sdmass", genjetAK8sdmass, &b_genjetAK8sdmass);
   fChain->SetBranchAddress("genjetAK8btag", genjetAK8btag, &b_genjetAK8btag);
   fChain->SetBranchAddress("genjetAK8tau1", genjetAK8tau1, &b_genjetAK8tau1);
   fChain->SetBranchAddress("genjetAK8tau2", genjetAK8tau2, &b_genjetAK8tau2);
   fChain->SetBranchAddress("genjetAK8tau3", genjetAK8tau3, &b_genjetAK8tau3);
   fChain->SetBranchAddress("genjetAK8moment1", genjetAK8moment1, &b_genjetAK8moment1);
   fChain->SetBranchAddress("genjetAK8moment2", genjetAK8moment2, &b_genjetAK8moment2);
   fChain->SetBranchAddress("genjetAK8moment3", genjetAK8moment3, &b_genjetAK8moment3);
   fChain->SetBranchAddress("genjetAK8chrad", genjetAK8chrad, &b_genjetAK8chrad);
   fChain->SetBranchAddress("genjetAK8axis2", genjetAK8axis2, &b_genjetAK8axis2);
   fChain->SetBranchAddress("genjetAK8pTD", genjetAK8pTD, &b_genjetAK8pTD);
   fChain->SetBranchAddress("genjetAK8sub1pt", genjetAK8sub1pt, &b_genjetAK8sub1pt);
   fChain->SetBranchAddress("genjetAK8sub1y", genjetAK8sub1y, &b_genjetAK8sub1y);
   fChain->SetBranchAddress("genjetAK8sub1phi", genjetAK8sub1phi, &b_genjetAK8sub1phi);
   fChain->SetBranchAddress("genjetAK8sub1mass", genjetAK8sub1mass, &b_genjetAK8sub1mass);
   fChain->SetBranchAddress("genjetAK8sub2pt", genjetAK8sub2pt, &b_genjetAK8sub2pt);
   fChain->SetBranchAddress("genjetAK8sub2y", genjetAK8sub2y, &b_genjetAK8sub2y);
   fChain->SetBranchAddress("genjetAK8sub2phi", genjetAK8sub2phi, &b_genjetAK8sub2phi);
   fChain->SetBranchAddress("genjetAK8sub2mass", genjetAK8sub2mass, &b_genjetAK8sub2mass);
   fChain->SetBranchAddress("ngenjetAK4", &ngenjetAK4, &b_ngenjetAK4);
   fChain->SetBranchAddress("genjetAK4pt", genjetAK4pt, &b_genjetAK4pt);
   fChain->SetBranchAddress("genjetAK4y", genjetAK4y, &b_genjetAK4y);
   fChain->SetBranchAddress("genjetAK4phi", genjetAK4phi, &b_genjetAK4phi);
   fChain->SetBranchAddress("genjetAK4mass", genjetAK4mass, &b_genjetAK4mass);
   fChain->SetBranchAddress("genjetAK4sdmass", genjetAK4sdmass, &b_genjetAK4sdmass);
   fChain->SetBranchAddress("genjetAK4btag", genjetAK4btag, &b_genjetAK4btag);
   fChain->SetBranchAddress("ngenparticles", &ngenparticles, &b_ngenparticles);
   fChain->SetBranchAddress("genpartstatus", genpartstatus, &b_genpartstatus);
   fChain->SetBranchAddress("genpartpdg", genpartpdg, &b_genpartpdg);
   fChain->SetBranchAddress("genpartmompdg", genpartmompdg, &b_genpartmompdg);
   fChain->SetBranchAddress("genpartdaugno", genpartdaugno, &b_genpartdaugno);
   fChain->SetBranchAddress("genpartfromhard", genpartfromhard, &b_genpartfromhard);
   fChain->SetBranchAddress("genpartfromhardbFSR", genpartfromhardbFSR, &b_genpartfromhardbFSR);
   fChain->SetBranchAddress("genpartpt", genpartpt, &b_genpartpt);
   fChain->SetBranchAddress("genparteta", genparteta, &b_genparteta);
   fChain->SetBranchAddress("genpartphi", genpartphi, &b_genpartphi);
   fChain->SetBranchAddress("genpartm", genpartm, &b_genpartm);
   fChain->SetBranchAddress("genpartq", genpartq, &b_genpartq);
	}
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("muonisPF", muonisPF, &b_muonisPF);
   fChain->SetBranchAddress("muonisGL", muonisGL, &b_muonisGL);
   fChain->SetBranchAddress("muonisTRK", muonisTRK, &b_muonisTRK);
   fChain->SetBranchAddress("muonisLoose", muonisLoose, &b_muonisLoose);
   fChain->SetBranchAddress("muonisGoodGL", muonisGoodGL, &b_muonisGoodGL);
   fChain->SetBranchAddress("muonisMed", muonisMed, &b_muonisMed);
   fChain->SetBranchAddress("muonpt", muonpt, &b_muonpt);
   fChain->SetBranchAddress("muonp", muonp, &b_muonp);
   fChain->SetBranchAddress("muone", muone, &b_muone);
   fChain->SetBranchAddress("muoneta", muoneta, &b_muoneta);
   fChain->SetBranchAddress("muonphi", muonphi, &b_muonphi);
   fChain->SetBranchAddress("muondrbm", muondrbm, &b_muondrbm);
   fChain->SetBranchAddress("muontrkvtx", muontrkvtx, &b_muontrkvtx);
   fChain->SetBranchAddress("muondz", muondz, &b_muondz);
   fChain->SetBranchAddress("muonpter", muonpter, &b_muonpter);
   fChain->SetBranchAddress("muonchi", muonchi, &b_muonchi);
   fChain->SetBranchAddress("muonndf", muonndf, &b_muonndf);
   fChain->SetBranchAddress("muonecal", muonecal, &b_muonecal);
   fChain->SetBranchAddress("muonhcal", muonhcal, &b_muonhcal);
   fChain->SetBranchAddress("muonemiso", muonemiso, &b_muonemiso);
   fChain->SetBranchAddress("muonhadiso", muonhadiso, &b_muonhadiso);
   fChain->SetBranchAddress("muonpfiso", muonpfiso, &b_muonpfiso);
   fChain->SetBranchAddress("muontkpt03", muontkpt03, &b_muontkpt03);
   fChain->SetBranchAddress("muontkpt05", muontkpt05, &b_muontkpt05);
   fChain->SetBranchAddress("muonposmatch", muonposmatch, &b_muonposmatch);
   fChain->SetBranchAddress("muontrkink", muontrkink, &b_muontrkink);
   fChain->SetBranchAddress("muonsegcom", muonsegcom, &b_muonsegcom);
   fChain->SetBranchAddress("muonthit", muonhit, &b_muonhit);
   fChain->SetBranchAddress("muonpixhit", muonpixhit, &b_muonpixhit);
   fChain->SetBranchAddress("muonmst", muonmst, &b_muonmst);
   fChain->SetBranchAddress("muontrklay", muontrklay, &b_muontrklay);
   fChain->SetBranchAddress("muonvalfrac", muonvalfrac, &b_muonvalfrac);
   fChain->SetBranchAddress("nelecs", &nelecs, &b_nelecs);
   fChain->SetBranchAddress("elpt", elpt, &b_elpt);
   fChain->SetBranchAddress("eleta", eleta, &b_eleta);
   fChain->SetBranchAddress("elphi", elphi, &b_elphi);
   fChain->SetBranchAddress("elp", elp, &b_elp);
   fChain->SetBranchAddress("ele", ele, &b_ele);
   fChain->SetBranchAddress("elmvaid", elmvaid, &b_elmvaid);
   fChain->SetBranchAddress("eldxy", eldxy, &b_eldxy);
   fChain->SetBranchAddress("eldz", eldz, &b_eldz);
   fChain->SetBranchAddress("elhovere", elhovere, &b_elhovere);
   fChain->SetBranchAddress("elchi", elchi, &b_elchi);
   fChain->SetBranchAddress("elndf", elndf, &b_elndf);
   fChain->SetBranchAddress("eltkpt03", eltkpt03, &b_eltkpt03);
   fChain->SetBranchAddress("eltkpt04", eltkpt04, &b_eltkpt04);
   fChain->SetBranchAddress("elemiso04", elemiso04, &b_elemiso04);
   fChain->SetBranchAddress("elhadiso04", elhadiso04, &b_elhadiso04);
   fChain->SetBranchAddress("eletain", eletain, &b_eletain);
   fChain->SetBranchAddress("elphiin", elphiin, &b_elphiin);
   fChain->SetBranchAddress("elfbrem", elfbrem, &b_elfbrem);
   fChain->SetBranchAddress("elhadisodepth03", elhadisodepth03, &b_elhadisodepth03);
   fChain->SetBranchAddress("eleoverp", eleoverp, &b_eleoverp);
   fChain->SetBranchAddress("elietaieta", elietaieta, &b_elietaieta);
   fChain->SetBranchAddress("elmisshits", elmisshits, &b_elmisshits);
   fChain->SetBranchAddress("elchhadiso", elchhadiso, &b_elchhadiso);
   fChain->SetBranchAddress("elneuhadiso", elneuhadiso, &b_elneuhadiso);
   fChain->SetBranchAddress("elphoiso", elphoiso, &b_elphoiso);
   fChain->SetBranchAddress("elpuchhadiso", elpuchhadiso, &b_elpuchhadiso);
   fChain->SetBranchAddress("elpfiso", elpfiso, &b_elpfiso);
   fChain->SetBranchAddress("elconvdist", elconvdist, &b_elconvdist);
   fChain->SetBranchAddress("elconvdoct", elconvdoct, &b_elconvdoct);
   fChain->SetBranchAddress("nphotons", &nphotons, &b_nphotons);
   fChain->SetBranchAddress("phoe", phoe, &b_phoe);
   fChain->SetBranchAddress("phoeta", phoeta, &b_phoeta);
   fChain->SetBranchAddress("phophi", phophi, &b_phophi);
   fChain->SetBranchAddress("phomvaid", phomvaid, &b_phomvaid);
   fChain->SetBranchAddress("phoe1by9", phoe1by9, &b_phoe1by9);
   fChain->SetBranchAddress("phoe9by25", phoe9by25, &b_phoe9by25);
   fChain->SetBranchAddress("photrkiso", photrkiso, &b_photrkiso);
   fChain->SetBranchAddress("phoemiso", phoemiso, &b_phoemiso);
   fChain->SetBranchAddress("phohadiso", phohadiso, &b_phohadiso);
   fChain->SetBranchAddress("phochhadiso", phochhadiso, &b_phochhadiso);
   fChain->SetBranchAddress("phoneuhadiso", phoneuhadiso, &b_phoneuhadiso);
   fChain->SetBranchAddress("phophoiso", phophoiso, &b_phophoiso);
   fChain->SetBranchAddress("phoPUiso", phoPUiso, &b_phoPUiso);
   fChain->SetBranchAddress("phohadbyem", phohadbyem, &b_phohadbyem);
   fChain->SetBranchAddress("phoietaieta", phoietaieta, &b_phoietaieta);
}

Bool_t Anal_L5JERC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Anal_L5JERC_cxx
