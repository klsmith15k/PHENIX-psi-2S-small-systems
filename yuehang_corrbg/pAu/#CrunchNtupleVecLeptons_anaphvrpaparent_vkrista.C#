//separate PYTHIA processes//for bb
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TText.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TChain.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "TGraphErrors.h"

#define PI  3.141592653589793238
#define me  0.000510998
#define mm  0.1056584
#define mB  5.279

using namespace std;
int _cc_or_bb;
TH1D *h_00_mass_mm_cc_FG12;
TH1D *h_01_mass_mm_cc_FG12;
TH1D *h_11_mass_mm_cc_FG12;
TH1D *h_00_mass_mm_cc_FG11;
TH1D *h_01_mass_mm_cc_FG11;
TH1D *h_11_mass_mm_cc_FG11;
TH1D *h_00_mass_mm_cc_FG22;
TH1D *h_01_mass_mm_cc_FG22;
TH1D *h_11_mass_mm_cc_FG22;
TH1D *h_00_mass_mm_bb_FG12;
TH1D *h_01_mass_mm_bb_FG12;
TH1D *h_11_mass_mm_bb_FG12;
TH1D *h_00_mass_mm_bb_FG11;
TH1D *h_01_mass_mm_bb_FG11;
TH1D *h_11_mass_mm_bb_FG11;
TH1D *h_00_mass_mm_bb_FG22;
TH1D *h_01_mass_mm_bb_FG22;
TH1D *h_11_mass_mm_bb_FG22;

TH2D *h_00_mass_pt_mm_cc_FG12;
TH2D *h_01_mass_pt_mm_cc_FG12;
TH2D *h_11_mass_pt_mm_cc_FG12;
TH2D *h_00_mass_pt_mm_cc_FG11;
TH2D *h_01_mass_pt_mm_cc_FG11;
TH2D *h_11_mass_pt_mm_cc_FG11;
TH2D *h_00_mass_pt_mm_cc_FG22;
TH2D *h_01_mass_pt_mm_cc_FG22;
TH2D *h_11_mass_pt_mm_cc_FG22;
TH2D *h_00_mass_pt_mm_bb_FG12;
TH2D *h_01_mass_pt_mm_bb_FG12;
TH2D *h_11_mass_pt_mm_bb_FG12;
TH2D *h_00_mass_pt_mm_bb_FG11;
TH2D *h_01_mass_pt_mm_bb_FG11;
TH2D *h_11_mass_pt_mm_bb_FG11;
TH2D *h_00_mass_pt_mm_bb_FG22;
TH2D *h_01_mass_pt_mm_bb_FG22;
TH2D *h_11_mass_pt_mm_bb_FG22;

TH1D *ffit_ccrpA;

TH1D *h_11_single_pt_m_bb_FGLS;
TH1D *h_11_single_pt_B_bb_FGLS;

TH1D *w_11_single_pt_m_bb_FGLS;
TH1D *w_11_single_pt_B_bb_FGLS;

TH1D *r_11_single_pt_m_bb_FGLS;
TH1D *r_11_single_pt_B_bb_FGLS;

TH2D *h_00_mass_dphi_mm_cc_FG12;
TH2D *h_01_mass_dphi_mm_cc_FG12;
TH2D *h_11_mass_dphi_mm_cc_FG12;

TH2D *h_00_mass_dphi_mm_bb_FG12;
TH2D *h_01_mass_dphi_mm_bb_FG12;
TH2D *h_11_mass_dphi_mm_bb_FG12;

TH2D *h_00_mass_dphi_mm_bb_FGLS;
TH2D *h_01_mass_dphi_mm_bb_FGLS;
TH2D *h_11_mass_dphi_mm_bb_FGLS;

TH2D *h_11_d00_mass_dphi_mm_bb_FG12;
TH2D *h_11_d01_mass_dphi_mm_bb_FG12;
TH2D *h_11_d11_mass_dphi_mm_bb_FG12;

TH2D *h_11_d00_mass_dphi_mm_bb_FGLS;
TH2D *h_11_d01_mass_dphi_mm_bb_FGLS;
TH2D *h_11_d11_mass_dphi_mm_bb_FGLS;

TH2D *h_11_mass_dphi_em_cc_FG12;
TH2D *h_11_mass_dphi_em_bb_FG12;
TH2D *h_11_mass_dphi_em_bb_FGLS;

TH2D *h_11_mass_dphi_ee_cc_FG12;
TH2D *h_11_mass_dphi_ee_bb_FG12;
TH2D *h_11_mass_dphi_ee_bb_FGLS;

TH1D *h_11_dphi_em_cc_FG12_varbin;
TH1D *h_11_dphi_em_bb_FG12_varbin;
TH1D *h_11_dphi_em_bb_FGLS_varbin;

TH1D *h_11_dphi_ee_cc_FG12_varbin;
TH1D *h_11_dphi_ee_bb_FG12_varbin;
TH1D *h_11_dphi_ee_bb_FGLS_varbin;

TH1D *h_11_pair_pt_m_bb_FGLS;
TH1D *w_11_pair_pt_m_bb_FGLS;
TH1D *r_11_pair_pt_m_bb_FGLS;


double weight1;
double weight2;
double highptcut;

TCanvas *c1;
TDatime dt;
TText *tdate;
TText *tdate1;

Char_t ofile[50];

void CrunchNtupleVecLeptons_anaphvrpaparent_vkrista(int a);
void book_histo();
void write_histo();

int CheckAccS(double rap);
int CheckAccN(double rap);
int CheckP(double pt, double pz);
int CheckAccEMe(double eta, double pt);
int CheckAccEMm(double eta, double pt);
int CheckButsykAcc_plpl(double phi, int q, double pT);


void CrunchNtupleVecLeptons_anaphvrpaparent_vkrista(int _crpamode)
{
    // Create a timer object to benchmark this loop
    TStopwatch timer;
    timer.Start();

    tdate = new TText(0.91,0.91,dt.AsString());
    tdate->SetNDC();
    tdate->SetTextSize(0.03);
    tdate->SetTextAlign(31);

    tdate1 = new TText(0.99,0.97,dt.AsString());
    tdate1->SetNDC();
    tdate1->SetTextSize(0.03);
    tdate1->SetTextAlign(31);

    gStyle->SetOptStat(1111111);
    gStyle->SetStatFormat("6.8g");
    gStyle->SetOptFit(0);
    //gStyle->SetPalette(1);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);

    int nChains = 0;

    // _cc_or_bb = 0; //0 = cc, 1 = bb
    _cc_or_bb = 1; //0 = cc, 1 = bb

    int _select_process = -1;//all processes

    // int _select_process = 0;
    // int _select_process = 1;
    // int _select_process = 2;

    int _applyccrpA = 0;
    _applyccrpA = _crpamode;

    TChain etree("leptons");
    Char_t path_ch[100];
    ifstream fin;


          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/msel_hf/filelist.txt");
          // strcpy (ofile,"output/pythia_mb_leptons_gps0049.root");
          // nChains = 50;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/msel_hf/filelist.txt");
          // strcpy (ofile,"output/pythia_mb_leptons_gps0099.root");
          // nChains = 100;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/pythiamsel5_bb/filelist.txt");
          // strcpy (ofile,"output/pythia_msel5_leptons_gps0004.root");
          // nChains = 5;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/ana_powheg_bb_singles_anav5/filelist.txt");
          // strcpy (ofile,"output/powheg_bb_leptons_gps0004.root");
          // nChains = 5;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/ana_powheg_cc_singles_anav5/filelist.txt");
          // strcpy (ofile,"output/powheg_cc_vf_leptons_gps0004.root");
          // nChains = 5;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/pythiamsel4_cc/filelist.txt");
          // strcpy (ofile,"output/pythia_msel4_leptons_gps0004.root");
          // nChains = 5;

          // fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/msel_hf/filelist.txt");
          // strcpy (ofile,"output/pythia_mb_leptons_test.root");
          // nChains = 1;


          //fvtx disk
          //CHARM
          // fin.open("/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/yhleung/msel1cc/filelist.txt");
          // strcpy (ofile,Form("output/pythia_mb_cc_leptons_vf_process_%d.root",_select_process));
          // nChains = 40;
          // nChains = 10;


          //BOTTOM
          // fin.open("/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/yhleung/msel1bb/filelist.txt");
          // strcpy (ofile,Form("output/pythia_mb_bb_process_%d.root",_select_process));
          // nChains = 200;

          // //POWHEG CHARM
          // fin.open("/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/yhleung/powheg/pp/run6/filelist.txt");
          // fin.open("/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/yhleung/powheg/pp/filelistcc.txt");
          // strcpy (ofile,"output/powheg_cc_vf_leptons_14664.root");
          // strcpy (ofile,"output/powheg_cc_vf_leptons_1000.root");
          // nChains = 14664;
          // nChains = 100;

          //POWHEG BOTTOM
          // fin.open("/gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/yhleung/powheg/pp/run6/filelist.txt");
          fin.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj3/yhleung/pythia/dst/powhegbbparent/run6/filelist.txt");
          strcpy (ofile,Form("output/powheg_bbparent_vrpa_leptons_5000_ccrpa%d_vkrista.root",_applyccrpA));
          // nChains = 7973;
          // nChains = 100;
          // nChains = 5000;
          // nChains = 1000;
          // nChains = 500;
          nChains = 5000;

    double nTrueEvents = 0.0;
    for(int i=0; i<nChains; i++){

      fin >> path_ch;
        etree.Add(path_ch);
        cout<<path_ch<<"  "<<etree.GetEntries() <<endl;
        nTrueEvents+=100000.0;

    }

    cout << " done TChain-ing.. " << endl;
    book_histo();


    Double_t  E=0;
    Double_t  m=0;
    Double_t  pt=0;

    Double_t  pE=0;
    Double_t  pm=0;
    Double_t  ppt=0;


    //Declaration of leaves types
    Float_t         ID;
    Float_t         px;
    Float_t         py;
    Float_t         pz;
    Float_t         pID;
    Float_t         ppx;
    Float_t         ppy;
    Float_t         ppz;
    Float_t         evt;
    // Float_t         decay;
    // Float_t         process;



    // Set branch addresses.
    etree.SetBranchAddress("ID",&ID);
    etree.SetBranchAddress("px",&px);
    etree.SetBranchAddress("py",&py);
    etree.SetBranchAddress("pz",&pz);
    etree.SetBranchAddress("pID",&pID);
    etree.SetBranchAddress("ppx",&ppx);
    etree.SetBranchAddress("ppy",&ppy);
    etree.SetBranchAddress("ppz",&ppz);
    etree.SetBranchAddress("evt",&evt);
    // etree.SetBranchAddress("decay",&decay);
    // etree.SetBranchAddress("process",&process);



    Long64_t nEntries = etree.GetEntries();//changed, why use int_t()??
    cout<<endl<<"Start processing "<<endl<<endl;
    cout<<"--------------------------------------------"<<endl;
    cout<<"Total entries "<<nEntries<<endl;

    TLorentzVector Parent, myParticle;
    vector<TLorentzVector> electrons;
    vector<TLorentzVector> positrons;
    vector<TLorentzVector> muplus;
    vector<TLorentzVector> muminus;

    vector<TLorentzVector> muplusparent;
    vector<TLorentzVector> muminusparent;



    // TF1 *ffit_ccrpA = new TF1("ffit_ccrpA","pol3",0,10);
    ffit_ccrpA = new TH1D("ffit_ccrpA","ffit_ccrpA",50,0,10);
    highptcut = 10.0;

    if(_applyccrpA>30){
    // ffit_ccrpA = new TF1("ffit_ccrpA","([0]+[1]*x+[2]*x*x+[3]*x*x*x)+[4]*exp(-(x-[5])*(x-[5])/2/[6]/[6])",0,9);
    }

    // if(_applyccrpA==0){ffit_ccrpA->SetParameters(1.,0.,0.,0.,0.);}//NO WEIGHT
    // if(_applyccrpA==1){ffit_ccrpA->SetParameters(-0.0126872,0.860661,-0.224038,0.016317);}
    // if(_applyccrpA==2){ffit_ccrpA->SetParameters(0.272555,1.23352,-0.333023,0.0243945);}
    // if(_applyccrpA==3){ffit_ccrpA->SetParameters(-0.382825,1.4936,-0.429175,0.0355012);}
    // if(_applyccrpA==4){ffit_ccrpA->SetParameters(-0.0832401,1.97126,-0.563057,0.0450929);}
    // if(_applyccrpA==5){ffit_ccrpA->SetParameters(-0.286107,1.08048,-0.331414,0.0287583);}
    // if(_applyccrpA==6){ffit_ccrpA->SetParameters(-0.0547509,1.39146,-0.42842,0.0364409);}
    // if(_applyccrpA==7){ffit_ccrpA->SetParameters(-0.318365,1.06362,-0.28865,0.0236927);}
    // if(_applyccrpA==8){ffit_ccrpA->SetParameters(-0.129494,1.39015,-0.365864,0.0285548);}
    // if(_applyccrpA==9){ffit_ccrpA->SetParameters(-0.350567,1.51046,-0.471939,0.0405668);}
    // if(_applyccrpA==10){ffit_ccrpA->SetParameters(-0.00849729,1.97257,-0.625613,0.052979);}

    // if(_applyccrpA==11){
    //   ffit_ccrpA->SetParameters(-1.1807,2.13109,-0.602156,0.0484389);
    //   highptcut = 6.0;
    // }  //MOD5, CHARM_RDAU MACRO
    // if(_applyccrpA==12){
    //   ffit_ccrpA->SetParameters(-1.55339,3.16207,-0.889124,0.0703654);
    //   highptcut = 6.0;
    // }  //MOD5, CHARM_RDAU MACRO

    //too flat
    // if(_applyccrpA==13){ffit_ccrpA->SetParameters(-1.93313,2.94119,-0.792159,0.0605706);}  //MOD6, CHARM_RDAU MACRO
    // if(_applyccrpA==14){ffit_ccrpA->SetParameters(-2.76851,4.47554,-1.19857,0.0905468);}  //MOD6, CHARM_RDAU MACRO

    // if(_applyccrpA==13){
    //   ffit_ccrpA->SetParameters(-2.98739,4.34731,-1.29765,0.11103);
    //   highptcut = 5.4;
    // }  //MOD6, CHARM_RDAU MACRO
    // if(_applyccrpA==14){
    //   ffit_ccrpA->SetParameters(-4.42997,6.68073,-1.98296,0.168061);
    //   highptcut = 5.4;
    // }  //MOD6, CHARM_RDAU MACRO



    // if(_applyccrpA==15){
    //   ffit_ccrpA->SetParameters(-4.36175,6.32752,-2.04938,0.188445);
    //   highptcut = 4.0;
    // }  //MOD7, CHARM_RDAU MACRO
    // if(_applyccrpA==16){
    //   ffit_ccrpA->SetParameters(-6.72109,9.96389,-3.22126,0.294969);
    //   highptcut = 4.0;
    // }  //MOD7, CHARM_RDAU MACRO


    // if(_applyccrpA==17){
    //   ffit_ccrpA->SetParameters(-6.16783,8.91388,-2.99214,0.282811);
    //   highptcut = 3.6;
    // }  //MOD8, CHARM_RDAU MACRO
    // if(_applyccrpA==18){
    //   ffit_ccrpA->SetParameters(-9.84593,14.4214,-4.83837,0.456287);
    //   highptcut = 3.6;
    // }  //MOD8, CHARM_RDAU MACRO

    // if(_applyccrpA==19){
    //   ffit_ccrpA->SetParameters(-8.0605,11.4569,-3.87331,0.368229);
    //   highptcut = 3.6;
    // }  //MOD9, CHARM_RDAU MACRO
    // if(_applyccrpA==20){
    //   ffit_ccrpA->SetParameters(-12.9894,18.6404,-6.29432,0.596965);
    //   highptcut = 3.6;
    // }  //MOD9, CHARM_RDAU MACRO

    if(_applyccrpA==21){
      ffit_ccrpA->SetBinContent(1,0.5);
      ffit_ccrpA->SetBinContent(2,0.55);
      ffit_ccrpA->SetBinContent(3,0.6);
      ffit_ccrpA->SetBinContent(4,0.7);
      ffit_ccrpA->SetBinContent(5,0.8);
      ffit_ccrpA->SetBinContent(6,0.9);
      ffit_ccrpA->SetBinContent(7,1.1);
      ffit_ccrpA->SetBinContent(8,1.3);
      ffit_ccrpA->SetBinContent(9,1.5);
      ffit_ccrpA->SetBinContent(10,1.7);

      ffit_ccrpA->SetBinContent(11,1.8);
      ffit_ccrpA->SetBinContent(12,1.7);
      ffit_ccrpA->SetBinContent(13,1.5);
      ffit_ccrpA->SetBinContent(14,1.3);
      ffit_ccrpA->SetBinContent(15,1.1);
      ffit_ccrpA->SetBinContent(16,1.0);
      ffit_ccrpA->SetBinContent(17,1.0);
      ffit_ccrpA->SetBinContent(18,1.0);
      ffit_ccrpA->SetBinContent(19,1.0);
      ffit_ccrpA->SetBinContent(20,1.0);

      ffit_ccrpA->SetBinContent(21,1.);
      ffit_ccrpA->SetBinContent(22,1.);
      ffit_ccrpA->SetBinContent(23,1.);
      ffit_ccrpA->SetBinContent(24,1.);
      ffit_ccrpA->SetBinContent(25,1.);
      ffit_ccrpA->SetBinContent(26,1.);
      ffit_ccrpA->SetBinContent(27,1.);
      ffit_ccrpA->SetBinContent(28,1.);
      ffit_ccrpA->SetBinContent(29,1.);
      ffit_ccrpA->SetBinContent(30,1.);

      ffit_ccrpA->SetBinContent(31,1.);
      ffit_ccrpA->SetBinContent(32,1.);
      ffit_ccrpA->SetBinContent(33,1.);
      ffit_ccrpA->SetBinContent(34,1.);
      ffit_ccrpA->SetBinContent(35,1.);
      ffit_ccrpA->SetBinContent(36,1.);
      ffit_ccrpA->SetBinContent(37,1.);
      ffit_ccrpA->SetBinContent(38,1.);
      ffit_ccrpA->SetBinContent(39,1.);
      ffit_ccrpA->SetBinContent(40,1.);

      ffit_ccrpA->SetBinContent(41,1.);
      ffit_ccrpA->SetBinContent(42,1.);
      ffit_ccrpA->SetBinContent(43,1.);
      ffit_ccrpA->SetBinContent(44,1.);
      ffit_ccrpA->SetBinContent(45,1.);
      ffit_ccrpA->SetBinContent(46,1.);
      ffit_ccrpA->SetBinContent(47,1.);
      ffit_ccrpA->SetBinContent(48,1.);
      ffit_ccrpA->SetBinContent(49,1.);
      ffit_ccrpA->SetBinContent(50,1.);

      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO
    if(_applyccrpA==22){
      ffit_ccrpA->SetBinContent(1,0.25);
      ffit_ccrpA->SetBinContent(2,0.3);
      ffit_ccrpA->SetBinContent(3,0.35);
      ffit_ccrpA->SetBinContent(4,0.4);
      ffit_ccrpA->SetBinContent(5,0.45);
      ffit_ccrpA->SetBinContent(6,0.5);
      ffit_ccrpA->SetBinContent(7,0.55);
      ffit_ccrpA->SetBinContent(8,0.6);
      ffit_ccrpA->SetBinContent(9,0.7);
      ffit_ccrpA->SetBinContent(10,0.8);

      ffit_ccrpA->SetBinContent(11,0.9);
      ffit_ccrpA->SetBinContent(12,1.1);
      ffit_ccrpA->SetBinContent(13,1.3);
      ffit_ccrpA->SetBinContent(14,1.5);
      ffit_ccrpA->SetBinContent(15,1.7);
      ffit_ccrpA->SetBinContent(16,1.8);
      ffit_ccrpA->SetBinContent(17,1.7);
      ffit_ccrpA->SetBinContent(18,1.5);
      ffit_ccrpA->SetBinContent(19,1.3);
      ffit_ccrpA->SetBinContent(20,1.1);

      ffit_ccrpA->SetBinContent(21,1.);
      ffit_ccrpA->SetBinContent(22,1.);
      ffit_ccrpA->SetBinContent(23,1.);
      ffit_ccrpA->SetBinContent(24,1.);
      ffit_ccrpA->SetBinContent(25,1.);
      ffit_ccrpA->SetBinContent(26,1.);
      ffit_ccrpA->SetBinContent(27,1.);
      ffit_ccrpA->SetBinContent(28,1.);
      ffit_ccrpA->SetBinContent(29,1.);
      ffit_ccrpA->SetBinContent(30,1.);

      ffit_ccrpA->SetBinContent(31,1.);
      ffit_ccrpA->SetBinContent(32,1.);
      ffit_ccrpA->SetBinContent(33,1.);
      ffit_ccrpA->SetBinContent(34,1.);
      ffit_ccrpA->SetBinContent(35,1.);
      ffit_ccrpA->SetBinContent(36,1.);
      ffit_ccrpA->SetBinContent(37,1.);
      ffit_ccrpA->SetBinContent(38,1.);
      ffit_ccrpA->SetBinContent(39,1.);
      ffit_ccrpA->SetBinContent(40,1.);

      ffit_ccrpA->SetBinContent(41,1.);
      ffit_ccrpA->SetBinContent(42,1.);
      ffit_ccrpA->SetBinContent(43,1.);
      ffit_ccrpA->SetBinContent(44,1.);
      ffit_ccrpA->SetBinContent(45,1.);
      ffit_ccrpA->SetBinContent(46,1.);
      ffit_ccrpA->SetBinContent(47,1.);
      ffit_ccrpA->SetBinContent(48,1.);
      ffit_ccrpA->SetBinContent(49,1.);
      ffit_ccrpA->SetBinContent(50,1.);

      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO



    if(_applyccrpA==27){

      ffit_ccrpA->SetBinContent(1,1.6);
      ffit_ccrpA->SetBinContent(2,1.6);
      ffit_ccrpA->SetBinContent(3,1.7);
      ffit_ccrpA->SetBinContent(4,1.7);
      ffit_ccrpA->SetBinContent(5,1.6);
      ffit_ccrpA->SetBinContent(6,1.6);
      ffit_ccrpA->SetBinContent(7,1.4);
      ffit_ccrpA->SetBinContent(8,1.2);
      ffit_ccrpA->SetBinContent(9,1.1);
      ffit_ccrpA->SetBinContent(10,1.0);

      ffit_ccrpA->SetBinContent(11,0.9);
      ffit_ccrpA->SetBinContent(12,0.85);
      ffit_ccrpA->SetBinContent(13,0.8);
      ffit_ccrpA->SetBinContent(14,0.75);
      ffit_ccrpA->SetBinContent(15,0.75);
      ffit_ccrpA->SetBinContent(16,0.75);
      ffit_ccrpA->SetBinContent(17,0.75);
      ffit_ccrpA->SetBinContent(18,0.75);
      ffit_ccrpA->SetBinContent(19,0.75);
      ffit_ccrpA->SetBinContent(20,0.75);

      ffit_ccrpA->SetBinContent(21,0.75);
      ffit_ccrpA->SetBinContent(22,0.75);
      ffit_ccrpA->SetBinContent(23,0.75);
      ffit_ccrpA->SetBinContent(24,0.75);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);

      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);

      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==28){

      ffit_ccrpA->SetBinContent(1,1.0);
      ffit_ccrpA->SetBinContent(2,1.2);
      ffit_ccrpA->SetBinContent(3,1.4);
      ffit_ccrpA->SetBinContent(4,1.6);
      ffit_ccrpA->SetBinContent(5,1.6);
      ffit_ccrpA->SetBinContent(6,1.7);
      ffit_ccrpA->SetBinContent(7,1.7);
      ffit_ccrpA->SetBinContent(8,1.6);
      ffit_ccrpA->SetBinContent(9,1.6);
      ffit_ccrpA->SetBinContent(10,1.4);

      ffit_ccrpA->SetBinContent(11,1.2);
      ffit_ccrpA->SetBinContent(12,1.1);
      ffit_ccrpA->SetBinContent(13,1.0);
      ffit_ccrpA->SetBinContent(14,0.9);
      ffit_ccrpA->SetBinContent(15,0.85);
      ffit_ccrpA->SetBinContent(16,0.8);
      ffit_ccrpA->SetBinContent(17,0.75);
      ffit_ccrpA->SetBinContent(18,0.75);
      ffit_ccrpA->SetBinContent(19,0.75);
      ffit_ccrpA->SetBinContent(20,0.75);

      ffit_ccrpA->SetBinContent(21,0.75);
      ffit_ccrpA->SetBinContent(22,0.75);
      ffit_ccrpA->SetBinContent(23,0.75);
      ffit_ccrpA->SetBinContent(24,0.75);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);

      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);

      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==29){

      ffit_ccrpA->SetBinContent(1,0.5);
      ffit_ccrpA->SetBinContent(2,0.6);
      ffit_ccrpA->SetBinContent(3,0.8);
      ffit_ccrpA->SetBinContent(4,1.0);
      ffit_ccrpA->SetBinContent(5,1.2);
      ffit_ccrpA->SetBinContent(6,1.4);
      ffit_ccrpA->SetBinContent(7,1.55);
      ffit_ccrpA->SetBinContent(8,1.65);
      ffit_ccrpA->SetBinContent(9,1.7);
      ffit_ccrpA->SetBinContent(10,1.7);

      ffit_ccrpA->SetBinContent(11,1.65);
      ffit_ccrpA->SetBinContent(12,1.55);
      ffit_ccrpA->SetBinContent(13,1.4);
      ffit_ccrpA->SetBinContent(14,1.25);
      ffit_ccrpA->SetBinContent(15,1.1);
      ffit_ccrpA->SetBinContent(16,1.0);
      ffit_ccrpA->SetBinContent(17,0.9);
      ffit_ccrpA->SetBinContent(18,0.83);
      ffit_ccrpA->SetBinContent(19,0.8);
      ffit_ccrpA->SetBinContent(20,0.775);

      ffit_ccrpA->SetBinContent(21,0.763);
      ffit_ccrpA->SetBinContent(22,0.75);
      ffit_ccrpA->SetBinContent(23,0.75);
      ffit_ccrpA->SetBinContent(24,0.75);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);

      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);

      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==30){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.42);
      ffit_ccrpA->SetBinContent(3,0.45);
      ffit_ccrpA->SetBinContent(4,0.5);
      ffit_ccrpA->SetBinContent(5,0.6);
      ffit_ccrpA->SetBinContent(6,0.8);
      ffit_ccrpA->SetBinContent(7,1.0);
      ffit_ccrpA->SetBinContent(8,1.2);
      ffit_ccrpA->SetBinContent(9,1.4);
      ffit_ccrpA->SetBinContent(10,1.55);
      ffit_ccrpA->SetBinContent(11,1.65);
      ffit_ccrpA->SetBinContent(12,1.7);
      ffit_ccrpA->SetBinContent(13,1.7);
      ffit_ccrpA->SetBinContent(14,1.65);
      ffit_ccrpA->SetBinContent(15,1.55);
      ffit_ccrpA->SetBinContent(16,1.4);
      ffit_ccrpA->SetBinContent(17,1.25);
      ffit_ccrpA->SetBinContent(18,1.1);
      ffit_ccrpA->SetBinContent(19,1.0);
      ffit_ccrpA->SetBinContent(20,0.9);
      ffit_ccrpA->SetBinContent(21,0.83);
      ffit_ccrpA->SetBinContent(22,0.8);
      ffit_ccrpA->SetBinContent(23,0.775);
      ffit_ccrpA->SetBinContent(24,0.763);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);



      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==31){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.4);
      ffit_ccrpA->SetBinContent(3,0.4);
      ffit_ccrpA->SetBinContent(4,0.4);
      ffit_ccrpA->SetBinContent(5,0.42);
      ffit_ccrpA->SetBinContent(6,0.45);
      ffit_ccrpA->SetBinContent(7,0.5);
      ffit_ccrpA->SetBinContent(8,0.6);
      ffit_ccrpA->SetBinContent(9,0.8);
      ffit_ccrpA->SetBinContent(10,1.0);
      ffit_ccrpA->SetBinContent(11,1.2);
      ffit_ccrpA->SetBinContent(12,1.4);
      ffit_ccrpA->SetBinContent(13,1.55);
      ffit_ccrpA->SetBinContent(14,1.65);
      ffit_ccrpA->SetBinContent(15,1.7);
      ffit_ccrpA->SetBinContent(16,1.7);
      ffit_ccrpA->SetBinContent(17,1.65);
      ffit_ccrpA->SetBinContent(18,1.55);
      ffit_ccrpA->SetBinContent(19,1.4);
      ffit_ccrpA->SetBinContent(20,1.25);
      ffit_ccrpA->SetBinContent(21,1.1);
      ffit_ccrpA->SetBinContent(22,1.0);
      ffit_ccrpA->SetBinContent(23,0.9);
      ffit_ccrpA->SetBinContent(24,0.83);
      ffit_ccrpA->SetBinContent(25,0.8);
      ffit_ccrpA->SetBinContent(26,0.775);
      ffit_ccrpA->SetBinContent(27,0.763);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==32){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.4);
      ffit_ccrpA->SetBinContent(3,0.4);
      ffit_ccrpA->SetBinContent(4,0.4);
      ffit_ccrpA->SetBinContent(5,0.4);
      ffit_ccrpA->SetBinContent(6,0.4);
      ffit_ccrpA->SetBinContent(7,0.4);
      ffit_ccrpA->SetBinContent(8,0.42);
      ffit_ccrpA->SetBinContent(9,0.45);
      ffit_ccrpA->SetBinContent(10,0.5);
      ffit_ccrpA->SetBinContent(11,0.6);
      ffit_ccrpA->SetBinContent(12,0.8);
      ffit_ccrpA->SetBinContent(13,1.0);
      ffit_ccrpA->SetBinContent(14,1.2);
      ffit_ccrpA->SetBinContent(15,1.4);
      ffit_ccrpA->SetBinContent(16,1.55);
      ffit_ccrpA->SetBinContent(17,1.65);
      ffit_ccrpA->SetBinContent(18,1.7);
      ffit_ccrpA->SetBinContent(19,1.7);

      ffit_ccrpA->SetBinContent(20,1.65);
      ffit_ccrpA->SetBinContent(21,1.55);
      ffit_ccrpA->SetBinContent(22,1.4);
      ffit_ccrpA->SetBinContent(23,1.25);
      ffit_ccrpA->SetBinContent(24,1.1);
      ffit_ccrpA->SetBinContent(25,1.0);
      ffit_ccrpA->SetBinContent(26,0.9);
      ffit_ccrpA->SetBinContent(27,0.83);
      ffit_ccrpA->SetBinContent(28,0.8);
      ffit_ccrpA->SetBinContent(29,0.775);

      ffit_ccrpA->SetBinContent(30,0.763);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);
      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);

      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);
      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);


      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO



    if(_applyccrpA==40){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.42);
      ffit_ccrpA->SetBinContent(3,0.45);
      ffit_ccrpA->SetBinContent(4,0.5);
      ffit_ccrpA->SetBinContent(5,0.6);
      ffit_ccrpA->SetBinContent(6,0.8);
      ffit_ccrpA->SetBinContent(7,1.0);
      ffit_ccrpA->SetBinContent(8,1.2);
      ffit_ccrpA->SetBinContent(9,1.43);
      ffit_ccrpA->SetBinContent(10,1.72);
      ffit_ccrpA->SetBinContent(11,1.92);
      ffit_ccrpA->SetBinContent(12,2.0);
      ffit_ccrpA->SetBinContent(13,2.0);

      ffit_ccrpA->SetBinContent(14,1.92);
      ffit_ccrpA->SetBinContent(15,1.72);
      ffit_ccrpA->SetBinContent(16,1.43);
      ffit_ccrpA->SetBinContent(17,1.25);
      ffit_ccrpA->SetBinContent(18,1.1);
      ffit_ccrpA->SetBinContent(19,1.0);
      ffit_ccrpA->SetBinContent(20,0.9);
      ffit_ccrpA->SetBinContent(21,0.83);
      ffit_ccrpA->SetBinContent(22,0.8);
      ffit_ccrpA->SetBinContent(23,0.775);

      ffit_ccrpA->SetBinContent(24,0.763);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);

      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);

      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==41){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.42);
      ffit_ccrpA->SetBinContent(3,0.45);

      ffit_ccrpA->SetBinContent(4,0.5);
      ffit_ccrpA->SetBinContent(5,0.6);
      ffit_ccrpA->SetBinContent(6,0.8);
      ffit_ccrpA->SetBinContent(7,1.0);
      ffit_ccrpA->SetBinContent(8,1.2);
      ffit_ccrpA->SetBinContent(9,1.5);
      ffit_ccrpA->SetBinContent(10,1.9);
      ffit_ccrpA->SetBinContent(11,2.2);
      ffit_ccrpA->SetBinContent(12,2.3);
      ffit_ccrpA->SetBinContent(13,2.3);

      ffit_ccrpA->SetBinContent(14,2.2);
      ffit_ccrpA->SetBinContent(15,1.9);
      ffit_ccrpA->SetBinContent(16,1.5);
      ffit_ccrpA->SetBinContent(17,1.25);
      ffit_ccrpA->SetBinContent(18,1.1);
      ffit_ccrpA->SetBinContent(19,1.0);
      ffit_ccrpA->SetBinContent(20,0.9);
      ffit_ccrpA->SetBinContent(21,0.83);
      ffit_ccrpA->SetBinContent(22,0.8);
      ffit_ccrpA->SetBinContent(23,0.775);

      ffit_ccrpA->SetBinContent(24,0.763);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);

      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);

      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==42){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.42);
      ffit_ccrpA->SetBinContent(3,0.45);

      ffit_ccrpA->SetBinContent(4,0.5);
      ffit_ccrpA->SetBinContent(5,0.6);
      ffit_ccrpA->SetBinContent(6,0.8);
      ffit_ccrpA->SetBinContent(7,1.0);
      ffit_ccrpA->SetBinContent(8,1.15);
      ffit_ccrpA->SetBinContent(9,1.29);
      ffit_ccrpA->SetBinContent(10,1.4);
      ffit_ccrpA->SetBinContent(11,1.5);
      ffit_ccrpA->SetBinContent(12,1.55);
      ffit_ccrpA->SetBinContent(13,1.55);

      ffit_ccrpA->SetBinContent(14,1.5);
      ffit_ccrpA->SetBinContent(15,1.4);
      ffit_ccrpA->SetBinContent(16,1.3);
      ffit_ccrpA->SetBinContent(17,1.2);
      ffit_ccrpA->SetBinContent(18,1.1);
      ffit_ccrpA->SetBinContent(19,1.0);
      ffit_ccrpA->SetBinContent(20,0.9);
      ffit_ccrpA->SetBinContent(21,0.83);
      ffit_ccrpA->SetBinContent(22,0.8);
      ffit_ccrpA->SetBinContent(23,0.775);

      ffit_ccrpA->SetBinContent(24,0.763);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);

      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);

      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO

    if(_applyccrpA==43){

      ffit_ccrpA->SetBinContent(1,0.4);
      ffit_ccrpA->SetBinContent(2,0.42);
      ffit_ccrpA->SetBinContent(3,0.45);

      ffit_ccrpA->SetBinContent(4,0.5);
      ffit_ccrpA->SetBinContent(5,0.6);
      ffit_ccrpA->SetBinContent(6,0.8);
      ffit_ccrpA->SetBinContent(7,1.0);
      ffit_ccrpA->SetBinContent(8,1.12);
      ffit_ccrpA->SetBinContent(9,1.22);
      ffit_ccrpA->SetBinContent(10,1.3);
      ffit_ccrpA->SetBinContent(11,1.37);
      ffit_ccrpA->SetBinContent(12,1.4);
      ffit_ccrpA->SetBinContent(13,1.4);

      ffit_ccrpA->SetBinContent(14,1.37);
      ffit_ccrpA->SetBinContent(15,1.31);
      ffit_ccrpA->SetBinContent(16,1.24);
      ffit_ccrpA->SetBinContent(17,1.16);
      ffit_ccrpA->SetBinContent(18,1.07);
      ffit_ccrpA->SetBinContent(19,0.99);
      ffit_ccrpA->SetBinContent(20,0.9);
      ffit_ccrpA->SetBinContent(21,0.83);
      ffit_ccrpA->SetBinContent(22,0.8);
      ffit_ccrpA->SetBinContent(23,0.775);

      ffit_ccrpA->SetBinContent(24,0.763);
      ffit_ccrpA->SetBinContent(25,0.75);
      ffit_ccrpA->SetBinContent(26,0.75);
      ffit_ccrpA->SetBinContent(27,0.75);
      ffit_ccrpA->SetBinContent(28,0.75);
      ffit_ccrpA->SetBinContent(29,0.75);
      ffit_ccrpA->SetBinContent(30,0.75);
      ffit_ccrpA->SetBinContent(31,0.75);
      ffit_ccrpA->SetBinContent(32,0.75);
      ffit_ccrpA->SetBinContent(33,0.75);

      ffit_ccrpA->SetBinContent(34,0.75);
      ffit_ccrpA->SetBinContent(35,0.75);
      ffit_ccrpA->SetBinContent(36,0.75);
      ffit_ccrpA->SetBinContent(37,0.75);
      ffit_ccrpA->SetBinContent(38,0.75);
      ffit_ccrpA->SetBinContent(39,0.75);
      ffit_ccrpA->SetBinContent(40,0.75);
      ffit_ccrpA->SetBinContent(41,0.75);
      ffit_ccrpA->SetBinContent(42,0.75);
      ffit_ccrpA->SetBinContent(43,0.75);

      ffit_ccrpA->SetBinContent(44,0.75);
      ffit_ccrpA->SetBinContent(45,0.75);
      ffit_ccrpA->SetBinContent(46,0.75);
      ffit_ccrpA->SetBinContent(47,0.75);
      ffit_ccrpA->SetBinContent(48,0.75);
      ffit_ccrpA->SetBinContent(49,0.75);
      ffit_ccrpA->SetBinContent(50,0.75);

      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO


    if(_applyccrpA==50){

      ffit_ccrpA->SetBinContent(1,0.65);
      ffit_ccrpA->SetBinContent(2,0.66);
      ffit_ccrpA->SetBinContent(3,0.68);
      ffit_ccrpA->SetBinContent(4,0.72);
      ffit_ccrpA->SetBinContent(5,0.8);
      ffit_ccrpA->SetBinContent(6,0.92);
      ffit_ccrpA->SetBinContent(7,1.05);
      ffit_ccrpA->SetBinContent(8,1.18);
      ffit_ccrpA->SetBinContent(9,1.3);
      ffit_ccrpA->SetBinContent(10,1.4);
      ffit_ccrpA->SetBinContent(11,1.47);
      ffit_ccrpA->SetBinContent(12,1.5);
      ffit_ccrpA->SetBinContent(13,1.5);

      ffit_ccrpA->SetBinContent(14,1.47);
      ffit_ccrpA->SetBinContent(15,1.4);
      ffit_ccrpA->SetBinContent(16,1.3);
      ffit_ccrpA->SetBinContent(17,1.18);
      ffit_ccrpA->SetBinContent(18,1.05);
      ffit_ccrpA->SetBinContent(19,0.96);
      ffit_ccrpA->SetBinContent(20,0.91);
      ffit_ccrpA->SetBinContent(21,0.88);
      ffit_ccrpA->SetBinContent(22,0.87);
      ffit_ccrpA->SetBinContent(23,0.86);

      ffit_ccrpA->SetBinContent(24,0.85);
      ffit_ccrpA->SetBinContent(25,0.85);
      ffit_ccrpA->SetBinContent(26,0.85);
      ffit_ccrpA->SetBinContent(27,0.85);
      ffit_ccrpA->SetBinContent(28,0.85);
      ffit_ccrpA->SetBinContent(29,0.85);
      ffit_ccrpA->SetBinContent(30,0.85);
      ffit_ccrpA->SetBinContent(31,0.85);
      ffit_ccrpA->SetBinContent(32,0.85);
      ffit_ccrpA->SetBinContent(33,0.85);

      ffit_ccrpA->SetBinContent(34,0.85);
      ffit_ccrpA->SetBinContent(35,0.85);
      ffit_ccrpA->SetBinContent(36,0.85);
      ffit_ccrpA->SetBinContent(37,0.85);
      ffit_ccrpA->SetBinContent(38,0.85);
      ffit_ccrpA->SetBinContent(39,0.85);
      ffit_ccrpA->SetBinContent(40,0.85);
      ffit_ccrpA->SetBinContent(41,0.85);
      ffit_ccrpA->SetBinContent(42,0.85);
      ffit_ccrpA->SetBinContent(43,0.85);

      ffit_ccrpA->SetBinContent(44,0.85);
      ffit_ccrpA->SetBinContent(45,0.85);
      ffit_ccrpA->SetBinContent(46,0.85);
      ffit_ccrpA->SetBinContent(47,0.85);
      ffit_ccrpA->SetBinContent(48,0.85);
      ffit_ccrpA->SetBinContent(49,0.85);
      ffit_ccrpA->SetBinContent(50,0.85);



      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO

    if(_applyccrpA==51){

      ffit_ccrpA->SetBinContent(1,0.75);
      ffit_ccrpA->SetBinContent(2,0.76);
      ffit_ccrpA->SetBinContent(3,0.78);
      ffit_ccrpA->SetBinContent(4,0.82);
      ffit_ccrpA->SetBinContent(5,0.9);
      ffit_ccrpA->SetBinContent(6,1.0);
      ffit_ccrpA->SetBinContent(7,1.08);
      ffit_ccrpA->SetBinContent(8,1.15);
      ffit_ccrpA->SetBinContent(9,1.21);
      ffit_ccrpA->SetBinContent(10,1.25);
      ffit_ccrpA->SetBinContent(11,1.28);
      ffit_ccrpA->SetBinContent(12,1.3);
      ffit_ccrpA->SetBinContent(13,1.3);

      ffit_ccrpA->SetBinContent(14,1.28);
      ffit_ccrpA->SetBinContent(15,1.25);
      ffit_ccrpA->SetBinContent(16,1.21);
      ffit_ccrpA->SetBinContent(17,1.15);
      ffit_ccrpA->SetBinContent(18,1.09);
      ffit_ccrpA->SetBinContent(19,1.04);
      ffit_ccrpA->SetBinContent(20,1.0);
      ffit_ccrpA->SetBinContent(21,0.97);
      ffit_ccrpA->SetBinContent(22,0.96);
      ffit_ccrpA->SetBinContent(23,0.95);

      ffit_ccrpA->SetBinContent(24,0.95);
      ffit_ccrpA->SetBinContent(25,0.95);
      ffit_ccrpA->SetBinContent(26,0.95);
      ffit_ccrpA->SetBinContent(27,0.95);
      ffit_ccrpA->SetBinContent(28,0.95);
      ffit_ccrpA->SetBinContent(29,0.95);
      ffit_ccrpA->SetBinContent(30,0.95);
      ffit_ccrpA->SetBinContent(31,0.95);
      ffit_ccrpA->SetBinContent(32,0.95);
      ffit_ccrpA->SetBinContent(33,0.95);

      ffit_ccrpA->SetBinContent(34,0.95);
      ffit_ccrpA->SetBinContent(35,0.95);
      ffit_ccrpA->SetBinContent(36,0.95);
      ffit_ccrpA->SetBinContent(37,0.95);
      ffit_ccrpA->SetBinContent(38,0.95);
      ffit_ccrpA->SetBinContent(39,0.95);
      ffit_ccrpA->SetBinContent(40,0.95);
      ffit_ccrpA->SetBinContent(41,0.95);
      ffit_ccrpA->SetBinContent(42,0.95);
      ffit_ccrpA->SetBinContent(43,0.95);

      ffit_ccrpA->SetBinContent(44,0.95);
      ffit_ccrpA->SetBinContent(45,0.95);
      ffit_ccrpA->SetBinContent(46,0.95);
      ffit_ccrpA->SetBinContent(47,0.95);
      ffit_ccrpA->SetBinContent(48,0.95);
      ffit_ccrpA->SetBinContent(49,0.95);
      ffit_ccrpA->SetBinContent(50,0.95);



      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO

    if(_applyccrpA==52){

      ffit_ccrpA->SetBinContent(1,0.3);
      ffit_ccrpA->SetBinContent(2,0.32);
      ffit_ccrpA->SetBinContent(3,0.35);
      ffit_ccrpA->SetBinContent(4,0.4);
      ffit_ccrpA->SetBinContent(5,0.5);
      ffit_ccrpA->SetBinContent(6,0.65);
      ffit_ccrpA->SetBinContent(7,0.85);
      ffit_ccrpA->SetBinContent(8,1.1);
      ffit_ccrpA->SetBinContent(9,1.35);
      ffit_ccrpA->SetBinContent(10,1.6);
      ffit_ccrpA->SetBinContent(11,1.78);
      ffit_ccrpA->SetBinContent(12,1.85);
      ffit_ccrpA->SetBinContent(13,1.85);

      ffit_ccrpA->SetBinContent(14,1.78);
      ffit_ccrpA->SetBinContent(15,1.6);
      ffit_ccrpA->SetBinContent(16,1.35);
      ffit_ccrpA->SetBinContent(17,1.1);
      ffit_ccrpA->SetBinContent(18,0.9);
      ffit_ccrpA->SetBinContent(19,0.8);
      ffit_ccrpA->SetBinContent(20,0.75);
      ffit_ccrpA->SetBinContent(21,0.705);
      ffit_ccrpA->SetBinContent(22,0.68);
      ffit_ccrpA->SetBinContent(23,0.66);

      ffit_ccrpA->SetBinContent(24,0.65);
      ffit_ccrpA->SetBinContent(25,0.65);
      ffit_ccrpA->SetBinContent(26,0.65);
      ffit_ccrpA->SetBinContent(27,0.65);
      ffit_ccrpA->SetBinContent(28,0.65);
      ffit_ccrpA->SetBinContent(29,0.65);
      ffit_ccrpA->SetBinContent(30,0.65);
      ffit_ccrpA->SetBinContent(31,0.65);
      ffit_ccrpA->SetBinContent(32,0.65);
      ffit_ccrpA->SetBinContent(33,0.65);

      ffit_ccrpA->SetBinContent(34,0.65);
      ffit_ccrpA->SetBinContent(35,0.65);
      ffit_ccrpA->SetBinContent(36,0.65);
      ffit_ccrpA->SetBinContent(37,0.65);
      ffit_ccrpA->SetBinContent(38,0.65);
      ffit_ccrpA->SetBinContent(39,0.65);
      ffit_ccrpA->SetBinContent(40,0.65);
      ffit_ccrpA->SetBinContent(41,0.65);
      ffit_ccrpA->SetBinContent(42,0.65);
      ffit_ccrpA->SetBinContent(43,0.65);

      ffit_ccrpA->SetBinContent(44,0.65);
      ffit_ccrpA->SetBinContent(45,0.65);
      ffit_ccrpA->SetBinContent(46,0.65);
      ffit_ccrpA->SetBinContent(47,0.65);
      ffit_ccrpA->SetBinContent(48,0.65);
      ffit_ccrpA->SetBinContent(49,0.65);
      ffit_ccrpA->SetBinContent(50,0.65);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO

    if(_applyccrpA==53){

      ffit_ccrpA->SetBinContent(1,0.2);
      ffit_ccrpA->SetBinContent(2,0.22);
      ffit_ccrpA->SetBinContent(3,0.25);
      ffit_ccrpA->SetBinContent(4,0.3);
      ffit_ccrpA->SetBinContent(5,0.38);
      ffit_ccrpA->SetBinContent(6,0.5);
      ffit_ccrpA->SetBinContent(7,0.7);
      ffit_ccrpA->SetBinContent(8,1.0);
      ffit_ccrpA->SetBinContent(9,1.3);
      ffit_ccrpA->SetBinContent(10,1.6);
      ffit_ccrpA->SetBinContent(11,1.9);
      ffit_ccrpA->SetBinContent(12,2.0);
      ffit_ccrpA->SetBinContent(13,2.0);

      ffit_ccrpA->SetBinContent(14,1.9);
      ffit_ccrpA->SetBinContent(15,1.6);
      ffit_ccrpA->SetBinContent(16,1.3);
      ffit_ccrpA->SetBinContent(17,1.0);
      ffit_ccrpA->SetBinContent(18,0.84);
      ffit_ccrpA->SetBinContent(19,0.74);
      ffit_ccrpA->SetBinContent(20,0.67);
      ffit_ccrpA->SetBinContent(21,0.62);
      ffit_ccrpA->SetBinContent(22,0.58);
      ffit_ccrpA->SetBinContent(23,0.56);

      ffit_ccrpA->SetBinContent(24,0.55);
      ffit_ccrpA->SetBinContent(25,0.55);
      ffit_ccrpA->SetBinContent(26,0.55);
      ffit_ccrpA->SetBinContent(27,0.55);
      ffit_ccrpA->SetBinContent(28,0.55);
      ffit_ccrpA->SetBinContent(29,0.55);
      ffit_ccrpA->SetBinContent(30,0.55);
      ffit_ccrpA->SetBinContent(31,0.55);
      ffit_ccrpA->SetBinContent(32,0.55);
      ffit_ccrpA->SetBinContent(33,0.55);

      ffit_ccrpA->SetBinContent(34,0.55);
      ffit_ccrpA->SetBinContent(35,0.55);
      ffit_ccrpA->SetBinContent(36,0.55);
      ffit_ccrpA->SetBinContent(37,0.55);
      ffit_ccrpA->SetBinContent(38,0.55);
      ffit_ccrpA->SetBinContent(39,0.55);
      ffit_ccrpA->SetBinContent(40,0.55);
      ffit_ccrpA->SetBinContent(41,0.55);
      ffit_ccrpA->SetBinContent(42,0.55);
      ffit_ccrpA->SetBinContent(43,0.55);

      ffit_ccrpA->SetBinContent(44,0.55);
      ffit_ccrpA->SetBinContent(45,0.55);
      ffit_ccrpA->SetBinContent(46,0.55);
      ffit_ccrpA->SetBinContent(47,0.55);
      ffit_ccrpA->SetBinContent(48,0.55);
      ffit_ccrpA->SetBinContent(49,0.55);
      ffit_ccrpA->SetBinContent(50,0.55);


      highptcut = 10.0;
    }  //MODXX, CHARM_RDAU MACRO



    if(_applyccrpA==0)
    {
      for(int ibin=0;ibin<=ffit_ccrpA->GetNbinsX();ibin++){
        ffit_ccrpA->SetBinContent(ibin,1.);
      }
    }
    // if(_applyccrpA==23){
    //   ffit_ccrpA->SetParameters(-2.89284,3.70122,-1.07711,0.0932336);
    //   highptcut = 4.4;
    // }  //MODXX, CHARM_RDAU MACRO
    // if(_applyccrpA==24){
    //   ffit_ccrpA->SetParameters(-4.78068,6.14392,-1.78562,0.153931);
    //   highptcut = 4.4;
    // }  //MODXX, CHARM_RDAU MACRO
    //
    //
    // if(_applyccrpA==31){
    //   ffit_ccrpA->SetParameters(0.272555,1.23352,-0.333023,0.0243945,1.5*1.5,2,0.1);
    //   highptcut = 6.0;
    // }  //MODXX, CHARM_RDAU MACRO
    // if(_applyccrpA==32){
    //   ffit_ccrpA->SetParameters(-0.0126872,0.860661,-0.224038,0.016317,0.9*1.5,2,0.1);
    //   highptcut = 6.0;
    // }  //MODXX, CHARM_RDAU MACRO


    int prev_evt = -1;
    bool e_PHENIX_N=0;
    bool p_PHENIX_N=0;
    bool e_PHENIX_S=0;
    bool p_PHENIX_S=0;

    bool e_P=0;
    bool p_P=0;

    double _dphi=-999;
    double _dphi_varbin=-999;

    Long64_t nentries =0; //what is wrong?

    for (Long64_t iEntry = 0; iEntry<=nEntries; iEntry++)
    {
        // nbytes = etree.GetEntry(iEntry);//original
        etree.GetEntry(iEntry);// what is wrong?????


        if(evt!=prev_evt)
        {
            nentries++;
            if(nentries%100000 ==0) cout<<"Processed FG12 "<<nentries<<endl;
            // if(nentries%100 ==0) cout<<"Processed FG12 "<<nentries<<endl;

            unsigned int ie, ip;


            //_________________________UNLIKE-SIGN________________________________________
            for(ie=0; ie<muminus.size(); ie++){
                for(ip=0; ip<muplus.size(); ip++){

                    Parent = muminus[ie]+muplus[ip];

                    e_PHENIX_S = CheckAccS(muminus[ie].PseudoRapidity());
                    p_PHENIX_S = CheckAccS(muplus[ip].PseudoRapidity());
                    e_PHENIX_N = CheckAccN(muminus[ie].PseudoRapidity());
                    p_PHENIX_N = CheckAccN(muplus[ip].PseudoRapidity());

                    e_P = CheckP(muminus[ie].Pt(),muminus[ie].Pz());
                    p_P = CheckP(muplus[ip].Pt(),muplus[ip].Pz());


                    if(muminusparent[ie].Pt()<highptcut){weight1 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muminusparent[ie].Pt()));}
                    else{weight1 = 1;}

                    if(muplusparent[ip].Pt()<highptcut){weight2 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muplusparent[ip].Pt()));}
                    else{weight2 = 1;}


                    // cout<<"parent:"<<muminusparent[ie].Pt()<<" muon:"<<muminus[ie].Pt()<<endl;

                    _dphi = fabs(atan2(muminus[ie].Py(),muminus[ie].Px())-atan2(muplus[ip].Py(),muplus[ip].Px()));
                    if(_dphi>TMath::Pi()) _dphi = 2*TMath::Pi()-_dphi;

                    //replace dphi with opening angle
                    // _dphi = muminus[ie].Angle(muplus[ip].Vect())*1.;



                    if(_cc_or_bb==0){//cc
                      h_00_mass_mm_cc_FG12->Fill(Parent.M() ,weight1*weight2);
                      h_00_mass_pt_mm_cc_FG12->Fill(Parent.M(), Parent.Pt() ,weight1*weight2);
                      h_00_mass_dphi_mm_cc_FG12->Fill(Parent.M(), _dphi ,weight1*weight2);

                      if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
                        h_01_mass_pt_mm_cc_FG12    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
                        h_01_mass_mm_cc_FG12       ->Fill(Parent.M(), weight1*weight2);
                        h_01_mass_dphi_mm_cc_FG12  ->Fill(Parent.M(), _dphi, weight1*weight2);

                        if(e_P && p_P){
                        h_11_mass_pt_mm_cc_FG12    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
                        h_11_mass_mm_cc_FG12       ->Fill(Parent.M(),weight1*weight2);
                        h_11_mass_dphi_mm_cc_FG12  ->Fill(Parent.M(), _dphi,weight1*weight2);

                        // cout<<"13 "<<muminus[ie].Px()<<" "<<muminus[ie].Py()<<" "<<muminus[ie].Pz()<<endl;
                        // cout<<"-13 "<<muplus[ip].Px()<<" "<<muplus[ip].Py()<<" "<<muplus[ip].Pz()<<endl;

                        }//momentum cut

                      }// IS in phenix if_statement
                      }//cc


          if(_cc_or_bb==1){//bb
            h_00_mass_mm_bb_FG12->Fill(Parent.M(),weight1*weight2);
            h_00_mass_pt_mm_bb_FG12->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
            h_00_mass_dphi_mm_bb_FG12->Fill(Parent.M(), _dphi,weight1*weight2);

            if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
              h_01_mass_pt_mm_bb_FG12    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
              h_01_mass_mm_bb_FG12       ->Fill(Parent.M(),weight1*weight2);
              h_01_mass_dphi_mm_bb_FG12  ->Fill(Parent.M(), _dphi,weight1*weight2);

              if(e_P && p_P){
              h_11_mass_pt_mm_bb_FG12    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
              h_11_mass_mm_bb_FG12       ->Fill(Parent.M(),weight1*weight2);
              h_11_mass_dphi_mm_bb_FG12  ->Fill(Parent.M(), _dphi,weight1*weight2);



              // if(muminus_decay[ie]==0 && muplus_decay[ip]==0){h_11_d00_mass_dphi_mm_bb_FG12->Fill(Parent.M(), _dphi);}
              // if(muminus_decay[ie]==1 && muplus_decay[ip]==1){h_11_d11_mass_dphi_mm_bb_FG12->Fill(Parent.M(), _dphi);}
              // if(muminus_decay[ie]==0 && muplus_decay[ip]==1){h_11_d01_mass_dphi_mm_bb_FG12->Fill(Parent.M(), _dphi);}
              // if(muminus_decay[ie]==1 && muplus_decay[ip]==0){h_11_d01_mass_dphi_mm_bb_FG12->Fill(Parent.M(), _dphi);}

              }//momentum cut

            }// IS in phenix if_statement
            }//bb
    }//ip loop
  }//ie loop




            //____________________________END UNLIKE-SIGN________________________________________

            //__________________________________LIKE-SIGN 11-muminus_____________________________________
for(ie=0; ie<muminus.size(); ie++){
  for(ip=ie+1; ip<muminus.size(); ip++){

    Parent = muminus[ie]+muminus[ip];

    e_PHENIX_S = CheckAccS(muminus[ie].PseudoRapidity());
    p_PHENIX_S = CheckAccS(muminus[ip].PseudoRapidity());
    e_PHENIX_N = CheckAccN(muminus[ie].PseudoRapidity());
    p_PHENIX_N = CheckAccN(muminus[ip].PseudoRapidity());

    e_P = CheckP(muminus[ie].Pt(),muminus[ie].Pz());
    p_P = CheckP(muminus[ip].Pt(),muminus[ip].Pz());


    if(muminusparent[ie].Pt()<highptcut){weight1 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muminusparent[ie].Pt()));}
    else{weight1 = 1;}

    if(muminusparent[ip].Pt()<highptcut){weight2 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muminusparent[ip].Pt()));}
    else{weight2 = 1;}

    // if(muminus[ie].Pt()<2.3 && muminus[ie].Pt()>1.7){weight1 = weight1*2.2;}
    // if(muminus[ip].Pt()<2.3 && muminus[ip].Pt()>1.7){weight2 = weight2*2.2;}



    _dphi = fabs(atan2(muminus[ie].Py(),muminus[ie].Px())-atan2(muminus[ip].Py(),muminus[ip].Px()));
    if(_dphi>TMath::Pi()) _dphi = 2*TMath::Pi()-_dphi;

    //replace dphi with opening angle
    // _dphi = muminus[ie].Angle(muminus[ip].Vect())*1.;


    if(_cc_or_bb==0){//cc

      // cout<<"parents:"<<eparentID[0]<<","<<eparentID[1]<<" "<<evt<<endl;
      // cout<<"pz:"<<muminus[0].Pz()<<","<<muminus[1].Pz()<<endl;

      h_00_mass_mm_cc_FG11->Fill(Parent.M(),weight1*weight2);
      h_00_mass_pt_mm_cc_FG11->Fill(Parent.M(), Parent.Pt(),weight1*weight2);

      if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
        h_01_mass_pt_mm_cc_FG11    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_01_mass_mm_cc_FG11       ->Fill(Parent.M(),weight1*weight2);

        if(e_P && p_P){
        h_11_mass_pt_mm_cc_FG11    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_11_mass_mm_cc_FG11       ->Fill(Parent.M(),weight1*weight2);
        }//momentum cut

      }// IS in phenix if_statement
      }//cc

    if(_cc_or_bb==1){//bb
      h_00_mass_mm_bb_FG11->Fill(Parent.M(),weight1*weight2);
      h_00_mass_pt_mm_bb_FG11->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
      h_00_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);

      if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
        h_01_mass_pt_mm_bb_FG11    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_01_mass_mm_bb_FG11       ->Fill(Parent.M(),weight1*weight2);
        h_01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);

        if(e_P && p_P){
        h_11_mass_pt_mm_bb_FG11    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_11_mass_mm_bb_FG11       ->Fill(Parent.M(),weight1*weight2);
        h_11_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);



        if(Parent.M()>3.5){
        // cout<<"1:"<<muminus[ie].Pt()<<" "<<weight1<<" 2:"<<muminus[ip].Pt()<<" "<<weight2<<endl;

        // double counting may occur, need to be fixed!
        // h_11_single_pt_m_bb_FGLS ->Fill(muminus[ie].Pt());
        // h_11_single_pt_B_bb_FGLS ->Fill(muminusparent[ie].Pt());
        //
        // h_11_single_pt_m_bb_FGLS ->Fill(muminus[ip].Pt());
        // h_11_single_pt_B_bb_FGLS ->Fill(muminusparent[ip].Pt());
        //
        // w_11_single_pt_m_bb_FGLS ->Fill(muminus[ie].Pt(), weight1);
        // w_11_single_pt_B_bb_FGLS ->Fill(muminusparent[ie].Pt(), weight1);
        //
        // w_11_single_pt_m_bb_FGLS ->Fill(muminus[ip].Pt(), weight2);
        // w_11_single_pt_B_bb_FGLS ->Fill(muminusparent[ip].Pt(), weight2);

        }

        // if(muminus_decay[ie]==0 && muminus_decay[ip]==0){h_11_d00_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muminus_decay[ie]==1 && muminus_decay[ip]==1){h_11_d11_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muminus_decay[ie]==0 && muminus_decay[ip]==1){h_11_d01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muminus_decay[ie]==1 && muminus_decay[ip]==0){h_11_d01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}

        }//momentum cut

      }// IS in phenix if_statement
      }//bb

  }//ie loop
}//ip loop
            //___________________________________END LIKE-SIGN 11-muminus_____________________________________

            //__________________________________LIKE-SIGN 22-muplus_____________________________________
for(ie=0; ie<muplus.size(); ie++){
  for(ip=ie+1; ip<muplus.size(); ip++){

    Parent = muplus[ie]+muplus[ip];
    e_PHENIX_S = CheckAccS(muplus[ie].PseudoRapidity());
    p_PHENIX_S = CheckAccS(muplus[ip].PseudoRapidity());
    e_PHENIX_N = CheckAccN(muplus[ie].PseudoRapidity());
    p_PHENIX_N = CheckAccN(muplus[ip].PseudoRapidity());

    e_P = CheckP(muplus[ie].Pt(),muplus[ie].Pz());
    p_P = CheckP(muplus[ip].Pt(),muplus[ip].Pz());


    // if(muplus[ie].Pt()<highptcut){weight1 = ffit_ccrpA->Eval(muplus[ie].Pt());}
    // else{weight1 = ffit_ccrpA->Eval(highptcut);}
    //
    // if(muplus[ip].Pt()<highptcut){weight2 = ffit_ccrpA->Eval(muplus[ip].Pt());}
    // else{weight2 = ffit_ccrpA->Eval(highptcut);}

    if(muplusparent[ie].Pt()<highptcut){weight1 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muplusparent[ie].Pt()));}
    else{weight1 = 1;}

    if(muplusparent[ip].Pt()<highptcut){weight2 = ffit_ccrpA->GetBinContent(ffit_ccrpA->FindBin(muplusparent[ip].Pt()));}
    else{weight2 = 1;}

    // if(muplus[ie].Pt()<2.3 && muplus[ie].Pt()>1.7){weight1 = weight1*2.2;}
    // if(muplus[ip].Pt()<2.3 && muplus[ip].Pt()>1.7){weight2 = weight2*2.2;}


    // cout<<"1:"<<muplusparent[ie].Pt()<<" "<<muplus[ie].Pt()<<" "<<weight1<<endl;
    // cout<<"2:"<<muplusparent[ip].Pt()<<" "<<muplus[ip].Pt()<<" "<<weight2<<endl;

    _dphi = fabs(atan2(muplus[ie].Py(),muplus[ie].Px())-atan2(muplus[ip].Py(),muplus[ip].Px()));
    if(_dphi>TMath::Pi()) _dphi = 2*TMath::Pi()-_dphi;

    //replace dphi with opening angle
    // _dphi = muplus[ie].Angle(muplus[ip].Vect())*1.;


    if(_cc_or_bb==0){//cc
      h_00_mass_mm_cc_FG22->Fill(Parent.M(),weight1*weight2);
      h_00_mass_pt_mm_cc_FG22->Fill(Parent.M(), Parent.Pt(),weight1*weight2);

      if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
        h_01_mass_pt_mm_cc_FG22    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_01_mass_mm_cc_FG22       ->Fill(Parent.M(),weight1*weight2);

        if(e_P && p_P){
        h_11_mass_pt_mm_cc_FG22    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_11_mass_mm_cc_FG22       ->Fill(Parent.M(),weight1*weight2);
        }//momentum cut

      }// IS in phenix if_statement
    }//cc

    if(_cc_or_bb==1){//bb
      h_00_mass_mm_bb_FG22->Fill(Parent.M(),weight1*weight2);
      h_00_mass_pt_mm_bb_FG22->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
      h_00_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);

      if((e_PHENIX_S && p_PHENIX_S) || (e_PHENIX_N && p_PHENIX_N)){
        h_01_mass_pt_mm_bb_FG22    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_01_mass_mm_bb_FG22       ->Fill(Parent.M(),weight1*weight2);
        h_01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);

        if(e_P && p_P){
        h_11_mass_pt_mm_bb_FG22    ->Fill(Parent.M(), Parent.Pt(),weight1*weight2);
        h_11_mass_mm_bb_FG22       ->Fill(Parent.M(),weight1*weight2);
        h_11_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi,weight1*weight2);

        // cout<<"1:"<<muplus[ie].Pt()<<" "<<weight1<<" 2:"<<muplus[ip].Pt()<<" "<<weight2<<endl;

        // if(muplus_decay[ie]==0 && muplus_decay[ip]==0){h_11_d00_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muplus_decay[ie]==1 && muplus_decay[ip]==1){h_11_d11_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muplus_decay[ie]==0 && muplus_decay[ip]==1){h_11_d01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}
        // if(muplus_decay[ie]==1 && muplus_decay[ip]==0){h_11_d01_mass_dphi_mm_bb_FGLS->Fill(Parent.M(), _dphi);}

        if(Parent.M()>3.5){


        // double counting may occur, need to be fixed!
        // h_11_single_pt_m_bb_FGLS ->Fill(muplus[ie].Pt());
        // h_11_single_pt_B_bb_FGLS ->Fill(muplusparent[ie].Pt());
        //
        // h_11_single_pt_m_bb_FGLS ->Fill(muplus[ip].Pt());
        // h_11_single_pt_B_bb_FGLS ->Fill(muplusparent[ip].Pt());
        //
        // w_11_single_pt_m_bb_FGLS ->Fill(muplus[ie].Pt(), weight1);
        // w_11_single_pt_B_bb_FGLS ->Fill(muplusparent[ie].Pt(), weight1);
        //
        // w_11_single_pt_m_bb_FGLS ->Fill(muplus[ip].Pt(), weight2);
        // w_11_single_pt_B_bb_FGLS ->Fill(muplusparent[ip].Pt(), weight2);


        h_11_pair_pt_m_bb_FGLS ->Fill(Parent.Pt());
        w_11_pair_pt_m_bb_FGLS ->Fill(Parent.Pt(),weight1*weight2);

        // cout<<"1:"<<muminus[ie].Pt()<<" "<<weight1<<" 2:"<<muminus[ip].Pt()<<" "<<weight2<<endl;

        }



        }//momentum cut

      }// IS in phenix if_statement
    }//bb

  }//ie loop
}//ip loop
            //_________________________________END LIKE-SIGN 22-muplus______________________________________


//electron-muons

//end electron-muons

//unlike sign dielectrons

            electrons.clear();
            positrons.clear();
            muplus.clear();
            muminus.clear();

            muplusparent.clear();
            muminusparent.clear();


        }
        prev_evt = evt;

        if(TMath::Abs(ID)==11) m = me;
        else if(TMath::Abs(ID)==13) m = mm;
        pm = mB;
        //else cout << "wtf is the mass?!?!??" << endl;


        //evt, ID, pID, px, py, pz, in_STAR, in_PHENIXcd
        pt = TMath::Sqrt(px*px + py*py);
        E  = TMath::Sqrt(pt*pt + pz*pz + m*m);

        ppt = TMath::Sqrt(ppx*ppx + ppy*ppy);
        pE  = TMath::Sqrt(ppt*ppt + ppz*ppz + pm*pm);

        myParticle.SetPxPyPzE(px,py,pz,E);
        Parent.SetPxPyPzE(ppx,ppy,ppz,pE);


        if(_select_process==-1){
  	if(ID==11)
        {
            electrons.push_back(TLorentzVector(px,py,pz,E));
            // electrons_decay.push_back(decay);
            // if(cc==1){electrons_process.push_back(0);}
            // if(bb==1){electrons_process.push_back(1);}

	      }
        else if(ID==-11){
            positrons.push_back(TLorentzVector(px,py,pz,E));
            // positrons_decay.push_back(decay);
            // if(cc==1){positrons_process.push_back(0);}
            // if(bb==1){positrons_process.push_back(1);}

        }
        else if(ID==13){
            muminus.push_back(TLorentzVector(px,py,pz,E));
            muminusparent.push_back(TLorentzVector(ppx,ppy,ppz,pE));
            // muminus_decay.push_back(decay);
            // eparentID.push_back(pID);

            // if(cc==1){muminus_process.push_back(0);}
            // if(bb==1){muminus_process.push_back(1);}

            if(fabs(myParticle.Rapidity())>1.2 &&  fabs(myParticle.Rapidity())<2.2){

            h_11_single_pt_m_bb_FGLS ->Fill(myParticle.Pt());
            h_11_single_pt_B_bb_FGLS ->Fill(Parent.Pt());

            double wtt =0;
            if(Parent.Pt()<10){wtt = ffit_ccrpA->FindBin(Parent.Pt());}
            if(Parent.Pt()>=10){wtt = ffit_ccrpA->FindBin(10-0.01);}

            w_11_single_pt_m_bb_FGLS ->Fill(myParticle.Pt(), wtt );
            w_11_single_pt_B_bb_FGLS ->Fill(Parent.Pt(),     wtt );

            }



	            }
        else if(ID==-13){
            muplus.push_back(TLorentzVector(px,py,pz,E));
            muplusparent.push_back(TLorentzVector(ppx,ppy,ppz,pE));
            // muplus_decay.push_back(decay);
            // pparentID.push_back(pID);

            // if(cc==1){muplus_process.push_back(0);}
            // if(bb==1){muplus_process.push_back(1);}

            if(fabs(myParticle.Rapidity())>1.2 &&  fabs(myParticle.Rapidity())<2.2){

            h_11_single_pt_m_bb_FGLS ->Fill(myParticle.Pt());
            h_11_single_pt_B_bb_FGLS ->Fill(Parent.Pt());

            double wtt =0;
            if(Parent.Pt()<10){wtt = ffit_ccrpA->FindBin(Parent.Pt());}
            if(Parent.Pt()>=10){wtt = ffit_ccrpA->FindBin(10-0.01);}

            w_11_single_pt_m_bb_FGLS ->Fill(myParticle.Pt(),  wtt );
            w_11_single_pt_B_bb_FGLS ->Fill(Parent.Pt(),      wtt );

            }

        }
      }//select process


    }//entries




    h_00_mass_mm_cc_FG12   ->Sumw2();
    h_00_mass_pt_mm_cc_FG12   ->Sumw2();
    h_01_mass_mm_cc_FG12   ->Sumw2();
    h_01_mass_pt_mm_cc_FG12   ->Sumw2();
    h_11_mass_mm_cc_FG12   ->Sumw2();
    h_11_mass_pt_mm_cc_FG12   ->Sumw2();

    h_00_mass_mm_cc_FG11   ->Sumw2();
    h_00_mass_pt_mm_cc_FG11   ->Sumw2();
    h_01_mass_mm_cc_FG11   ->Sumw2();
    h_01_mass_pt_mm_cc_FG11   ->Sumw2();
    h_11_mass_mm_cc_FG11   ->Sumw2();
    h_11_mass_pt_mm_cc_FG11   ->Sumw2();

    h_00_mass_mm_cc_FG22   ->Sumw2();
    h_00_mass_pt_mm_cc_FG22   ->Sumw2();
    h_01_mass_mm_cc_FG22   ->Sumw2();
    h_01_mass_pt_mm_cc_FG22   ->Sumw2();
    h_11_mass_mm_cc_FG22   ->Sumw2();
    h_11_mass_pt_mm_cc_FG22   ->Sumw2();


    h_00_mass_mm_bb_FG12   ->Sumw2();
    h_00_mass_pt_mm_bb_FG12   ->Sumw2();
    h_01_mass_mm_bb_FG12   ->Sumw2();
    h_01_mass_pt_mm_bb_FG12   ->Sumw2();
    h_11_mass_mm_bb_FG12   ->Sumw2();
    h_11_mass_pt_mm_bb_FG12   ->Sumw2();

    h_00_mass_mm_bb_FG11   ->Sumw2();
    h_00_mass_pt_mm_bb_FG11   ->Sumw2();
    h_01_mass_mm_bb_FG11   ->Sumw2();
    h_01_mass_pt_mm_bb_FG11   ->Sumw2();
    h_11_mass_mm_bb_FG11   ->Sumw2();
    h_11_mass_pt_mm_bb_FG11   ->Sumw2();

    h_00_mass_mm_bb_FG22   ->Sumw2();
    h_00_mass_pt_mm_bb_FG22   ->Sumw2();
    h_01_mass_mm_bb_FG22   ->Sumw2();
    h_01_mass_pt_mm_bb_FG22   ->Sumw2();
    h_11_mass_mm_bb_FG22   ->Sumw2();
    h_11_mass_pt_mm_bb_FG22   ->Sumw2();

    h_00_mass_dphi_mm_bb_FGLS->Sumw2();
    h_01_mass_dphi_mm_bb_FGLS->Sumw2();
    h_11_mass_dphi_mm_bb_FGLS->Sumw2();

    h_11_d00_mass_dphi_mm_bb_FG12->Sumw2();
    h_11_d01_mass_dphi_mm_bb_FG12->Sumw2();
    h_11_d11_mass_dphi_mm_bb_FG12->Sumw2();

    h_11_d00_mass_dphi_mm_bb_FGLS->Sumw2();
    h_11_d01_mass_dphi_mm_bb_FGLS->Sumw2();
    h_11_d11_mass_dphi_mm_bb_FGLS->Sumw2();

    h_11_mass_dphi_em_cc_FG12->Sumw2();
    h_11_mass_dphi_em_bb_FG12->Sumw2();
    h_11_mass_dphi_em_bb_FGLS->Sumw2();

    h_11_mass_dphi_ee_cc_FG12->Sumw2();
    h_11_mass_dphi_ee_bb_FG12->Sumw2();
    h_11_mass_dphi_ee_bb_FGLS->Sumw2();

    h_11_dphi_em_cc_FG12_varbin->Sumw2();
    h_11_dphi_em_bb_FG12_varbin->Sumw2();
    h_11_dphi_em_bb_FGLS_varbin->Sumw2();

    h_11_dphi_ee_cc_FG12_varbin->Sumw2();
    h_11_dphi_ee_bb_FG12_varbin->Sumw2();
    h_11_dphi_ee_bb_FGLS_varbin->Sumw2();


    for(int ibin=1;ibin<=r_11_single_pt_m_bb_FGLS->GetNbinsX();ibin++){
      if(h_11_single_pt_m_bb_FGLS->GetBinContent(ibin)>0){
      r_11_single_pt_m_bb_FGLS->SetBinContent(ibin, w_11_single_pt_m_bb_FGLS->GetBinContent(ibin)/h_11_single_pt_m_bb_FGLS->GetBinContent(ibin) );
      }

      if(h_11_single_pt_B_bb_FGLS->GetBinContent(ibin)>0){
      r_11_single_pt_B_bb_FGLS->SetBinContent(ibin, w_11_single_pt_B_bb_FGLS->GetBinContent(ibin)/h_11_single_pt_B_bb_FGLS->GetBinContent(ibin) );
      }

    }

    for(int ibin=1;ibin<=r_11_pair_pt_m_bb_FGLS->GetNbinsX();ibin++){
      if(h_11_pair_pt_m_bb_FGLS->GetBinContent(ibin)>0){
      r_11_pair_pt_m_bb_FGLS->SetBinContent(ibin, w_11_pair_pt_m_bb_FGLS->GetBinContent(ibin)/h_11_pair_pt_m_bb_FGLS->GetBinContent(ibin) );
      }
    }



    write_histo();

    cout << "Total number of generated events in MCNLO: " << nTrueEvents << endl;
    cout << "histograms have been normalized by this number! " << endl << endl;

    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    printf("RealTime=%f minutes, CpuTime=%f minutes\n",rtime/60,ctime/60);


}

void book_histo()
{

  h_00_mass_pt_mm_cc_FG12     = new TH2D("h_00_mass_pt_mm_cc_FG12" ,"h_00_mass_pt_mm_cc_FG12"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_cc_FG12     = new TH2D("h_01_mass_pt_mm_cc_FG12" ,"h_01_mass_pt_mm_cc_FG12"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_cc_FG12     = new TH2D("h_11_mass_pt_mm_cc_FG12" ,"h_11_mass_pt_mm_cc_FG12"  ,1500,0,15, 800,0,8);

  h_00_mass_pt_mm_cc_FG11     = new TH2D("h_00_mass_pt_mm_cc_FG11" ,"h_00_mass_pt_mm_cc_FG11"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_cc_FG11     = new TH2D("h_01_mass_pt_mm_cc_FG11" ,"h_01_mass_pt_mm_cc_FG11"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_cc_FG11     = new TH2D("h_11_mass_pt_mm_cc_FG11" ,"h_11_mass_pt_mm_cc_FG11"  ,1500,0,15, 800,0,8);

  h_00_mass_pt_mm_cc_FG22     = new TH2D("h_00_mass_pt_mm_cc_FG22" ,"h_00_mass_pt_mm_cc_FG22"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_cc_FG22     = new TH2D("h_01_mass_pt_mm_cc_FG22" ,"h_01_mass_pt_mm_cc_FG22"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_cc_FG22     = new TH2D("h_11_mass_pt_mm_cc_FG22" ,"h_11_mass_pt_mm_cc_FG22"  ,1500,0,15, 800,0,8);


  h_00_mass_pt_mm_bb_FG12     = new TH2D("h_00_mass_pt_mm_bb_FG12" ,"h_00_mass_pt_mm_bb_FG12"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_bb_FG12     = new TH2D("h_01_mass_pt_mm_bb_FG12" ,"h_01_mass_pt_mm_bb_FG12"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_bb_FG12     = new TH2D("h_11_mass_pt_mm_bb_FG12" ,"h_11_mass_pt_mm_bb_FG12"  ,1500,0,15, 800,0,8);

  h_00_mass_pt_mm_bb_FG11     = new TH2D("h_00_mass_pt_mm_bb_FG11" ,"h_00_mass_pt_mm_bb_FG11"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_bb_FG11     = new TH2D("h_01_mass_pt_mm_bb_FG11" ,"h_01_mass_pt_mm_bb_FG11"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_bb_FG11     = new TH2D("h_11_mass_pt_mm_bb_FG11" ,"h_11_mass_pt_mm_bb_FG11"  ,1500,0,15, 800,0,8);

  h_00_mass_pt_mm_bb_FG22     = new TH2D("h_00_mass_pt_mm_bb_FG22" ,"h_00_mass_pt_mm_bb_FG22"  ,1500,0,15, 800,0,8);
  h_01_mass_pt_mm_bb_FG22     = new TH2D("h_01_mass_pt_mm_bb_FG22" ,"h_01_mass_pt_mm_bb_FG22"  ,1500,0,15, 800,0,8);
  h_11_mass_pt_mm_bb_FG22     = new TH2D("h_11_mass_pt_mm_bb_FG22" ,"h_11_mass_pt_mm_bb_FG22"  ,1500,0,15, 800,0,8);



  h_00_mass_mm_cc_FG12     = new TH1D("h_00_mass_mm_cc_FG12" ,"h_00_mass_mm_cc_FG12"  ,1500,0,15);
  h_01_mass_mm_cc_FG12     = new TH1D("h_01_mass_mm_cc_FG12" ,"h_01_mass_mm_cc_FG12"  ,1500,0,15);
  h_11_mass_mm_cc_FG12     = new TH1D("h_11_mass_mm_cc_FG12" ,"h_11_mass_mm_cc_FG12"  ,1500,0,15);

  h_00_mass_mm_cc_FG11     = new TH1D("h_00_mass_mm_cc_FG11" ,"h_00_mass_mm_cc_FG11"  ,1500,0,15);
  h_01_mass_mm_cc_FG11     = new TH1D("h_01_mass_mm_cc_FG11" ,"h_01_mass_mm_cc_FG11"  ,1500,0,15);
  h_11_mass_mm_cc_FG11     = new TH1D("h_11_mass_mm_cc_FG11" ,"h_11_mass_mm_cc_FG11"  ,1500,0,15);

  h_00_mass_mm_cc_FG22     = new TH1D("h_00_mass_mm_cc_FG22" ,"h_00_mass_mm_cc_FG22"  ,1500,0,15);
  h_01_mass_mm_cc_FG22     = new TH1D("h_01_mass_mm_cc_FG22" ,"h_01_mass_mm_cc_FG22"  ,1500,0,15);
  h_11_mass_mm_cc_FG22     = new TH1D("h_11_mass_mm_cc_FG22" ,"h_11_mass_mm_cc_FG22"  ,1500,0,15);


  h_00_mass_mm_bb_FG12     = new TH1D("h_00_mass_mm_bb_FG12" ,"h_00_mass_mm_bb_FG12"  ,1500,0,15);
  h_01_mass_mm_bb_FG12     = new TH1D("h_01_mass_mm_bb_FG12" ,"h_01_mass_mm_bb_FG12"  ,1500,0,15);
  h_11_mass_mm_bb_FG12     = new TH1D("h_11_mass_mm_bb_FG12" ,"h_11_mass_mm_bb_FG12"  ,1500,0,15);

  h_00_mass_mm_bb_FG11     = new TH1D("h_00_mass_mm_bb_FG11" ,"h_00_mass_mm_bb_FG11"  ,1500,0,15);
  h_01_mass_mm_bb_FG11     = new TH1D("h_01_mass_mm_bb_FG11" ,"h_01_mass_mm_bb_FG11"  ,1500,0,15);
  h_11_mass_mm_bb_FG11     = new TH1D("h_11_mass_mm_bb_FG11" ,"h_11_mass_mm_bb_FG11"  ,1500,0,15);

  h_00_mass_mm_bb_FG22     = new TH1D("h_00_mass_mm_bb_FG22" ,"h_00_mass_mm_bb_FG22"  ,1500,0,15);
  h_01_mass_mm_bb_FG22     = new TH1D("h_01_mass_mm_bb_FG22" ,"h_01_mass_mm_bb_FG22"  ,1500,0,15);
  h_11_mass_mm_bb_FG22     = new TH1D("h_11_mass_mm_bb_FG22" ,"h_11_mass_mm_bb_FG22"  ,1500,0,15);


  h_00_mass_dphi_mm_cc_FG12     = new TH2D("h_00_mass_dphi_mm_cc_FG12" ,"h_00_mass_dphi_mm_cc_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_01_mass_dphi_mm_cc_FG12     = new TH2D("h_01_mass_dphi_mm_cc_FG12" ,"h_01_mass_dphi_mm_cc_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_mm_cc_FG12     = new TH2D("h_11_mass_dphi_mm_cc_FG12" ,"h_11_mass_dphi_mm_cc_FG12"  ,1500,0,15, 180,0,TMath::Pi());

  h_00_mass_dphi_mm_bb_FG12     = new TH2D("h_00_mass_dphi_mm_bb_FG12" ,"h_00_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_01_mass_dphi_mm_bb_FG12     = new TH2D("h_01_mass_dphi_mm_bb_FG12" ,"h_01_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_mm_bb_FG12     = new TH2D("h_11_mass_dphi_mm_bb_FG12" ,"h_11_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());

  h_00_mass_dphi_mm_bb_FGLS = new TH2D("h_00_mass_dphi_mm_bb_FGLS" ,"h_00_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());
  h_01_mass_dphi_mm_bb_FGLS = new TH2D("h_01_mass_dphi_mm_bb_FGLS" ,"h_01_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_mm_bb_FGLS = new TH2D("h_11_mass_dphi_mm_bb_FGLS" ,"h_11_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());

  h_11_d00_mass_dphi_mm_bb_FG12 = new TH2D("h_11_d00_mass_dphi_mm_bb_FG12" ,"h_11_d00_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_d01_mass_dphi_mm_bb_FG12 = new TH2D("h_11_d01_mass_dphi_mm_bb_FG12" ,"h_11_d01_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_d11_mass_dphi_mm_bb_FG12 = new TH2D("h_11_d11_mass_dphi_mm_bb_FG12" ,"h_11_d11_mass_dphi_mm_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());

  h_11_d00_mass_dphi_mm_bb_FGLS = new TH2D("h_11_d00_mass_dphi_mm_bb_FGLS" ,"h_11_d00_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_d01_mass_dphi_mm_bb_FGLS = new TH2D("h_11_d01_mass_dphi_mm_bb_FGLS" ,"h_11_d01_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_d11_mass_dphi_mm_bb_FGLS = new TH2D("h_11_d11_mass_dphi_mm_bb_FGLS" ,"h_11_d11_mass_dphi_mm_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());


  h_11_mass_dphi_em_cc_FG12     = new TH2D("h_11_mass_dphi_em_cc_FG12" ,"h_11_mass_dphi_em_cc_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_em_bb_FG12     = new TH2D("h_11_mass_dphi_em_bb_FG12" ,"h_11_mass_dphi_em_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_em_bb_FGLS     = new TH2D("h_11_mass_dphi_em_bb_FGLS" ,"h_11_mass_dphi_em_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());

  h_11_mass_dphi_ee_cc_FG12   = new TH2D("h_11_mass_dphi_ee_cc_FG12" ,"h_11_mass_dphi_ee_cc_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_ee_bb_FG12   = new TH2D("h_11_mass_dphi_ee_bb_FG12" ,"h_11_mass_dphi_ee_bb_FG12"  ,1500,0,15, 180,0,TMath::Pi());
  h_11_mass_dphi_ee_bb_FGLS   = new TH2D("h_11_mass_dphi_ee_bb_FGLS" ,"h_11_mass_dphi_ee_bb_FGLS"  ,1500,0,15, 180,0,TMath::Pi());

  h_11_dphi_em_cc_FG12_varbin     = new TH1D("h_11_dphi_em_cc_FG12_varbin" ,"h_11_dphi_em_cc_FG12_varbin"  ,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  h_11_dphi_em_bb_FG12_varbin     = new TH1D("h_11_dphi_em_bb_FG12_varbin" ,"h_11_dphi_em_bb_FG12_varbin"  ,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  h_11_dphi_em_bb_FGLS_varbin     = new TH1D("h_11_dphi_em_bb_FGLS_varbin" ,"h_11_dphi_em_bb_FGLS_varbin"  ,100,-0.5*TMath::Pi(),1.5*TMath::Pi());

  h_11_dphi_ee_cc_FG12_varbin     = new TH1D("h_11_dphi_ee_cc_FG12_varbin" ,"h_11_dphi_ee_cc_FG12_varbin"  ,28,0.3,3.1);//strange binning
  h_11_dphi_ee_bb_FG12_varbin     = new TH1D("h_11_dphi_ee_bb_FG12_varbin" ,"h_11_dphi_ee_bb_FG12_varbin"  ,28,0.3,3.1);
  h_11_dphi_ee_bb_FGLS_varbin     = new TH1D("h_11_dphi_ee_bb_FGLS_varbin" ,"h_11_dphi_ee_bb_FGLS_varbin"  ,28,0.3,3.1);


  h_11_single_pt_m_bb_FGLS = new TH1D("h_11_single_pt_m_bb_FGLS" ,"h_11_single_pt_m_bb_FGLS"  ,50,0,10);//strange binning
  h_11_single_pt_B_bb_FGLS = new TH1D("h_11_single_pt_B_bb_FGLS" ,"h_11_single_pt_B_bb_FGLS"  ,50,0,10);//strange binning

  w_11_single_pt_m_bb_FGLS = new TH1D("w_11_single_pt_m_bb_FGLS" ,"w_11_single_pt_m_bb_FGLS"  ,50,0,10);//strange binning
  w_11_single_pt_B_bb_FGLS = new TH1D("w_11_single_pt_B_bb_FGLS" ,"w_11_single_pt_B_bb_FGLS"  ,50,0,10);//strange binning

  r_11_single_pt_m_bb_FGLS = new TH1D("r_11_single_pt_m_bb_FGLS" ,"r_11_single_pt_m_bb_FGLS"  ,50,0,10);//strange binning
  r_11_single_pt_B_bb_FGLS = new TH1D("r_11_single_pt_B_bb_FGLS" ,"r_11_single_pt_B_bb_FGLS"  ,50,0,10);//strange binning

  h_11_pair_pt_m_bb_FGLS = new TH1D("h_11_pair_pt_m_bb_FGLS" , "h_11_pair_pt_m_bb_FGLS" , 5,0,5);
  w_11_pair_pt_m_bb_FGLS = new TH1D("w_11_pair_pt_m_bb_FGLS" , "w_11_pair_pt_m_bb_FGLS" , 5,0,5);
  r_11_pair_pt_m_bb_FGLS = new TH1D("r_11_pair_pt_m_bb_FGLS" , "r_11_pair_pt_m_bb_FGLS" , 5,0,5);


  w_11_single_pt_m_bb_FGLS->SetLineColor(kRed);
  w_11_single_pt_B_bb_FGLS->SetLineColor(kRed);
}

void write_histo()
{
  TFile *outhistfile = new TFile (ofile, "RECREATE");
  outhistfile->cd();

  h_00_mass_mm_cc_FG12->Write();
  h_00_mass_pt_mm_cc_FG12   ->Write();
  h_01_mass_mm_cc_FG12->Write();
  h_01_mass_pt_mm_cc_FG12   ->Write();
  h_11_mass_mm_cc_FG12->Write();
  h_11_mass_pt_mm_cc_FG12   ->Write();

  h_00_mass_mm_cc_FG11->Write();
  h_00_mass_pt_mm_cc_FG11   ->Write();
  h_01_mass_mm_cc_FG11->Write();
  h_01_mass_pt_mm_cc_FG11   ->Write();
  h_11_mass_mm_cc_FG11->Write();
  h_11_mass_pt_mm_cc_FG11   ->Write();

  h_00_mass_mm_cc_FG22->Write();
  h_00_mass_pt_mm_cc_FG22   ->Write();
  h_01_mass_mm_cc_FG22->Write();
  h_01_mass_pt_mm_cc_FG22   ->Write();
  h_11_mass_mm_cc_FG22->Write();
  h_11_mass_pt_mm_cc_FG22   ->Write();

  h_00_mass_mm_bb_FG12->Write();
  h_00_mass_pt_mm_bb_FG12   ->Write();
  h_01_mass_mm_bb_FG12->Write();
  h_01_mass_pt_mm_bb_FG12   ->Write();
  h_11_mass_mm_bb_FG12->Write();
  h_11_mass_pt_mm_bb_FG12   ->Write();

  h_00_mass_mm_bb_FG11->Write();
  h_00_mass_pt_mm_bb_FG11   ->Write();
  h_01_mass_mm_bb_FG11->Write();
  h_01_mass_pt_mm_bb_FG11   ->Write();
  h_11_mass_mm_bb_FG11->Write();
  h_11_mass_pt_mm_bb_FG11   ->Write();

  h_00_mass_mm_bb_FG22->Write();
  h_00_mass_pt_mm_bb_FG22   ->Write();
  h_01_mass_mm_bb_FG22->Write();
  h_01_mass_pt_mm_bb_FG22   ->Write();
  h_11_mass_mm_bb_FG22->Write();
  h_11_mass_pt_mm_bb_FG22   ->Write();

  h_00_mass_dphi_mm_cc_FG12->Write();
  h_01_mass_dphi_mm_cc_FG12->Write();
  h_11_mass_dphi_mm_cc_FG12->Write();

  h_00_mass_dphi_mm_bb_FG12->Write();
  h_01_mass_dphi_mm_bb_FG12->Write();
  h_11_mass_dphi_mm_bb_FG12->Write();

  h_00_mass_dphi_mm_bb_FGLS->Write();
  h_01_mass_dphi_mm_bb_FGLS->Write();
  h_11_mass_dphi_mm_bb_FGLS->Write();

  h_11_d00_mass_dphi_mm_bb_FG12->Write();
  h_11_d01_mass_dphi_mm_bb_FG12->Write();
  h_11_d11_mass_dphi_mm_bb_FG12->Write();
  h_11_d00_mass_dphi_mm_bb_FGLS->Write();
  h_11_d01_mass_dphi_mm_bb_FGLS->Write();
  h_11_d11_mass_dphi_mm_bb_FGLS->Write();


  h_11_mass_dphi_em_cc_FG12->Write();
  h_11_mass_dphi_em_bb_FG12->Write();
  h_11_mass_dphi_em_bb_FGLS->Write();

  h_11_mass_dphi_ee_cc_FG12->Write();
  h_11_mass_dphi_ee_bb_FG12->Write();
  h_11_mass_dphi_ee_bb_FGLS->Write();

  h_11_dphi_em_cc_FG12_varbin->Write();
  h_11_dphi_em_bb_FG12_varbin->Write();
  h_11_dphi_em_bb_FGLS_varbin->Write();

  h_11_dphi_ee_cc_FG12_varbin->Write();
  h_11_dphi_ee_bb_FG12_varbin->Write();
  h_11_dphi_ee_bb_FGLS_varbin->Write();

  h_11_single_pt_m_bb_FGLS ->Write();
  h_11_single_pt_B_bb_FGLS ->Write();

  w_11_single_pt_m_bb_FGLS ->Write();
  w_11_single_pt_B_bb_FGLS ->Write();

  r_11_single_pt_m_bb_FGLS ->Write();
  r_11_single_pt_B_bb_FGLS ->Write();


  h_11_pair_pt_m_bb_FGLS  ->Write();
  w_11_pair_pt_m_bb_FGLS  ->Write();
  r_11_pair_pt_m_bb_FGLS  ->Write();


  ffit_ccrpA->Write();

  outhistfile->Close();
  delete outhistfile;
}

int CheckAccS(double rap)
{
  int pass = 0;
	if(rap<-1.2 && rap>-2.2)
	  {
	    pass= 1;
	  }
	else
	  {
	  pass = 0;
	  }
	return pass;
}

int CheckAccN(double rap)
{
  int pass = 0;
	if(rap>1.2 && rap<2.2)
	  {
	    pass= 1;
	  }
	else
	  {
	  pass = 0;
	  }
	return pass;
}
int CheckP(double pt, double pz)
{
  int pass = 0;
  if(sqrt(pt*pt+pz*pz)>3)
  {
    pass=1;
  }
  else
  {
    pass=0;
  }
  return pass;
}

int CheckAccEMe(double eta, double pt)
{
  int pass = 0;
	if((eta<0.5 && eta>-0.5) && pt>0.5)
	  {
	    pass= 1;
	  }
	else
	  {
	  pass = 0;
	  }
	return pass;
}

int CheckAccEMm(double eta, double pt)
{
  int pass = 0;
	if((fabs(eta)<2.1 && fabs(eta)>1.4) && pt>1.0)
	  {
	    pass= 1;
	  }
	else
	  {
	  pass = 0;
	  }
	return pass;
}

int CheckButsykAcc_plpl(double phi, int q, double pT)
{


    int pass = 0;


    if(pT<0.2)
        return 0;

    if(phi<=-3.14159/2.0)        phi = phi + 3.14159*2.0;
    if( (phi>=4.2)  && (q==11)  ) phi = phi - 3.14159*2.0;
    else if( (phi<=-1.0) && (q==-11) ) phi = phi + 3.14159*2.0;


    q = q/11;

    //WEST ARM
    if( (double(q)/pT >= -4.854*phi-2.859) &&
        (double(q)/pT >= -3.236*phi-1.906) &&
        (double(q)/pT <= -3.236*phi+3.178) &&
        (double(q)/pT <= -4.854*phi+4.767)  )
        pass = 2;

    //EAST ARM
    else if( (double(q)/pT >= -4.854*phi+10.49) &&
             (double(q)/pT >= -3.236*phi+6.990) &&
             (double(q)/pT <= -3.236*phi+12.07) &&
             (double(q)/pT <= -4.854*phi+18.11)  )
        pass = 1;

    return pass;


}
