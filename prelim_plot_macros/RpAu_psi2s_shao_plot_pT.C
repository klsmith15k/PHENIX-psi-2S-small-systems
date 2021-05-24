
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <Math/MinimizerOptions.h>
#include <fstream>
#include <iostream>
#include <TAttText.h>
#include <TFrame.h>
#include <TBox.h>


using namespace std;


void RpAu_psi2s_shao_plot_pT()
{
    gROOT->ForceStyle();

  bool north_arm = true;
  int arm = 0;
  double FACTOR = 0.0005;

  static const int nbins = 6;
  static const int narm = 2;

  if(north_arm)
    arm = 1;
    
  double eff_bbc = 0.55;
  double eff_jpsi = 0.79;

  double abs1 = 0;
  double abs2 = 0;
  double pi = 3.141593;


  const Double_t pT_array1[nbins] = {0.25, 1.0, 2.0, 3.0, 4.5, 6.05};
  const Double_t pT_array2[nbins] = {0.25, 1.0, 2.0, 3.0, 4.5, 3.5};
  const Double_t pT_array3[nbins] = {0.25, 1.0, 2.0, 3.0, 4.5, 3.5};
  const Double_t pT_array4[nbins] = {0.25, 1.0, 2.0, 3.0, 4.5, 3.5};

  double delta_pt[nbins] = {0.5,1.0,1.0,1.0,2.0, 2.0};
  
  double pp_jpsi_binshift[nbins] = {1,1,1,1,1, 1};
  double pAu_jpsi_binshift[nbins] = {1,1,1,1,1, 1};

  double pp_psi2s_binshift[nbins] = {1,1,1,1,1, 1};
  double pAu_psi2s_binshift[nbins] = {1,1,1,1,1, 1};

  double bias[nbins] = {0.86, 0.86, 0.86, 0.86, 0.86, 0.86};
  double ncoll[nbins] = {4.7, 4.7, 4.7, 4.7, 4.7, 4.7};

  double pAu_MB[narm][nbins] = {
    {2.09769926640*pow(10,11), 2.09769926640*pow(10,11), 2.09769926640*pow(10,11), 2.09769926640*pow(10,11), 2.09769926640*pow(10,11), 2.09769926640*pow(10,11)},
    {1.93066407424*pow(10,11), 1.93066407424*pow(10,11), 1.93066407424*pow(10,11), 1.93066407424*pow(10,11), 1.93066407424*pow(10,11), 1.93066407424*pow(10,11)}};
 
  double pp_MB_228[narm] = {1.092422549504*pow(10,12) ,1.073082815488*pow(10,12)};    // PPG228 good runlist

  double pp_MB[narm] = {9.16550312960000000*pow(10,11), 8.16541736960000000*pow(10,11)};   // PPG201 good runlist

 TFile *file_data1, *file_data2;

  TH1D *t1, *t2;

  TGraphAsymmErrors *f1;
  TGraphAsymmErrors *f2;

  TEfficiency *e1;
  TEfficiency *e2;

  // 1a  for fvtx Jpsi
    //====================================================================================
  // fid 188N cut

  // fid 188N cut
  double jpsi_fvtx[narm][nbins] = {
    {3076,1977,1253,806},
    {2317,1995,1544,1143}};

  double jpsi_fvtx_err[narm][nbins] = {
    {64,50,36,29},
    {50,46,39,35}};
 
 double jpsi_fvtx_trig_eff[2][4] = {0};
 double jpsi_fvtx_acc_eff[2][4] = {0};
 
  {
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/fvtx_tracks/sanghoon_simulations/Run15pAu_fvtx_jpsi_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int cent =0; cent < 4; cent++)
    	 jpsi_fvtx_acc_eff[0][cent] = e1->GetEfficiency(cent+1);
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/fvtx_tracks/sanghoon_simulations/Run15pAu_fvtx_jpsi_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int cent =0; cent < 4; cent++)
    	 jpsi_fvtx_acc_eff[1][cent] = e2->GetEfficiency(cent+1);
	 
       file_data1->Close();
       file_data2->Close();
    }
       
    //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_fvtx/sanghoon_simulations/Run15pAu_fvtx_jpsi_trigeff_S.root");
       file_data1->GetObject("h1",f1);
       
       jpsi_fvtx_trig_eff[0][0] = f1->Eval(10);
       jpsi_fvtx_trig_eff[0][1] = f1->Eval(30);
       jpsi_fvtx_trig_eff[0][2] = f1->Eval(50);
       jpsi_fvtx_trig_eff[0][3] = f1->Eval(72);
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_fvtx/sanghoon_simulations/Run15pAu_fvtx_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       jpsi_fvtx_trig_eff[1][0] = f2->Eval(10);
       jpsi_fvtx_trig_eff[1][1] = f2->Eval(30);
       jpsi_fvtx_trig_eff[1][2] = f2->Eval(50);
       jpsi_fvtx_trig_eff[1][3] = f2->Eval(72);
       
    	file_data1->Close();
    	file_data2->Close();
       //////////////////////////////////////
	
   }

  //=================
     // for 1b
    double jpsi_fvtx_pAu_inv_yield[narm][nbins] = {0};
    double jpsi_fvtx_pAu_inv_yield_err[narm][nbins] = {0};

    // calculcate pAu inv. yield for jpsi fvtx
    // calculcate pAu inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {
 	    jpsi_fvtx_pAu_inv_yield[arm][bin] =  bias[bin] * jpsi_fvtx[arm][bin] / ( ( FACTOR + jpsi_fvtx_trig_eff[arm][bin]) *(  FACTOR + jpsi_fvtx_acc_eff[arm][bin]) * pAu_MB[arm][bin] );
 	    jpsi_fvtx_pAu_inv_yield_err[arm][bin] =  bias[bin] * jpsi_fvtx_err[arm][bin] / ( ( FACTOR + jpsi_fvtx_trig_eff[arm][bin] ) * (  FACTOR + jpsi_fvtx_acc_eff[arm][bin] )* pAu_MB[arm][bin]);
 	    // cout << "jpsi_fvtx_pAu inv yield: " << jpsi_fvtx_pAu_inv_yield[0][bin] << endl;
 	  }
      }
  
     
  
  
       // for pp invariant yield for jpsi fvtx
    double pp_jpsi_fvtx_inv_yield[narm] = {0}; 
    double pp_jpsi_fvtx_inv_yield_err[narm] = {0}; 
    double pp_jpsi_fvtx[narm] = {13392,8192};  
    double pp_jpsi_fvtx_err[narm] = {215,198};

    double pp_jpsi_fvtx_acc_eff[narm] = {0}; 
    double pp_jpsi_fvtx_trig_eff[narm] = {0};
    
   

 // RUN15pp
     {
          
       /////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/sanghoon_simulations_final/Run15pp_fvtx_jpsi_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       //  pp_jpsi_fvtx_acc_eff[0] = t1->GetBinContent(1);
       pp_jpsi_fvtx_acc_eff[0] = 0.0265;
       //pp_jpsi_fvtx_acc_eff[0] = 0.0158;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/sanghoon_simulations_final/Run15pp_fvtx_jpsi_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       // pp_jpsi_fvtx_acc_eff[1] = t2->GetBinContent(1);
       pp_jpsi_fvtx_acc_eff[1] = 0.0197;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////
     }
     // RUN15pp
     {
     
       /////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_fvtx/sanghoon_simulations_final/Run15pp_fvtx_jpsi_trigeff_S.root");
       file_data1->GetObject("geff_y_int",t1);
       //pp_jpsi_fvtx_trig_eff[0] = t1->GetBinContent(1);
       pp_jpsi_fvtx_trig_eff[0] = 0.7526;
       // pp_jpsi_fvtx_trig_eff[0] = 0.7343;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_fvtx/sanghoon_simulations_final/Run15pp_fvtx_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       // pp_jpsi_fvtx_trig_eff[1] = t2->GetBinContent(1);
       pp_jpsi_fvtx_trig_eff[1] = 0.7255;
	        
       file_data1->Close();
       file_data2->Close();
       //////////////////////////////
     }


     // calculate pp inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	pp_jpsi_fvtx_inv_yield[arm] = pp_jpsi_fvtx[arm]*eff_bbc / ( ( FACTOR + pp_jpsi_fvtx_trig_eff[arm] ) * ( FACTOR + pp_jpsi_fvtx_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
 	pp_jpsi_fvtx_inv_yield_err[arm] = pp_jpsi_fvtx_err[arm]*eff_bbc / ( ( FACTOR + pp_jpsi_fvtx_trig_eff[arm] ) * ( FACTOR + pp_jpsi_fvtx_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi);
       cout << "pp jpsi_fvtx fvtx inv yield: " << pp_jpsi_fvtx_inv_yield[arm] << endl;
      }

   

    double RpAu_jpsi_fvtx_cent[narm][nbins] = {0};
    double RpAu_jpsi_fvtx_cent_err[narm][nbins] = {0};
    
    // jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_jpsi_fvtx_cent[arm][bin] =( jpsi_fvtx_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_jpsi_fvtx_inv_yield[arm] ) ;
 	    abs1 = jpsi_fvtx_pAu_inv_yield_err[arm][bin] / jpsi_fvtx_pAu_inv_yield[arm][bin];
 	    abs2 = pp_jpsi_fvtx_inv_yield_err[arm] / pp_jpsi_fvtx_inv_yield[arm];
 	    RpAu_jpsi_fvtx_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_jpsi_fvtx_cent[arm][bin];
 	  }
      }

    
    // for ppg228 jpsi 
    //======================

    // double jpsi_ppg228[narm][nbins] = {
    //   {4916,3208,2128,1386},
    //   {6078,5174,3968,3039}}; 

    // double jpsi_ppg228_err[narm][nbins] = {
    //   {91,111,75,43},
    //   {104,90,80,69}};

    // double jpsi_ppg228_trig_eff[narm][nbins] = {
    //   {0.7573,0.6802,0.6241,0.6156},
    //   {0.7397,0.7394,0.7370,0.7397}};  

    // double jpsi_ppg228_acc_eff[narm][nbins] = {
    //   {0.0261,0.0270,0.0278,0.0276},
    //   {0.0502,0.0504,0.0507,0.0502}};
    
    // double jpsi_ppg228_pAu_inv_yield[narm][nbins] = {0};
    // double jpsi_ppg228_pAu_inv_yield_err[narm][nbins] = {0};

    // double percent_err = 0;

    // for(int arm = 0; arm < 2; arm++)
    //   {
    // 	for(int bin = 0; bin < 4; bin++)
    // 	  {
    // 	    jpsi_ppg228_pAu_inv_yield[arm][bin] = bias[bin] * jpsi_ppg228[arm][bin] / ( jpsi_ppg228_trig_eff[arm][bin]*jpsi_ppg228_acc_eff[arm][bin] * pAu_MB[arm][bin] );
    // 	    jpsi_ppg228_pAu_inv_yield_err[arm][bin] = jpsi_ppg228_err[arm][bin] / ( jpsi_ppg228_trig_eff[arm][bin]*jpsi_ppg228_acc_eff[arm][bin] * pAu_MB[arm][bin] );

    // 	    // percent_err = ( jpsi_ppg228_pAu_inv_yield_err[arm][bin]/  jpsi_ppg228_pAu_inv_yield[arm][bin] ) * 100;
    // 	    // cout <<  "J/psi inv yield ppg228: " << jpsi_ppg228_pAu_inv_yield[arm][bin] << ", and error: " <<  jpsi_ppg228_pAu_inv_yield_err[arm][bin] << ", and % err: " << percent_err << endl;
    // 	  }
    //   }

    // double pp_jpsi_ppg228_inv_yield[narm] = {0}; 
    // double pp_jpsi_ppg228_inv_yield_err[narm] = {0};  
    // double pp_jpsi_ppg228[narm] = {28511,31452};
    // double pp_jpsi_ppg228_err[narm] = {203,215};
    // double pp_jpsi_ppg228_trig_eff[narm] = {0.7499,0.7139};
    // double pp_jpsi_ppg228_acc_eff[narm] = {0.0383,0.0453};

    // for(int arm = 0; arm < 2; arm++)
    //   {
    // 	pp_jpsi_ppg228_inv_yield[arm] = pp_jpsi_ppg228[arm]*eff_bbc / ( pp_jpsi_ppg228_trig_eff[arm]*pp_jpsi_ppg228_acc_eff[arm] * pp_MB_228[arm] * eff_jpsi );
    // 	pp_jpsi_ppg228_inv_yield_err[arm] = pp_jpsi_ppg228_err[arm]*eff_bbc / ( pp_jpsi_ppg228_trig_eff[arm]*pp_jpsi_ppg228_acc_eff[arm] * pp_MB_228[arm] * eff_jpsi );
    //   }

    double RpAu_jpsi_ppg228_cent[narm][18] = {
      {0.69, 0.78, 0.68, 0.78, 0.88, 0.78 , 0.92, 0.92, 1.02, 1.02, 1.12, 1.06, 1.06, 1.22, 1.16, 1.13, 1.31 , 1.71},
      {0.65, 0.61, 0.63, 0.64, 0.66, 0.71, 0.78, 0.73, 0.83, 0.83, 0.88, 0.92, 0.97, 1.05, 0.98, 1.07, 1.20, 0.96}};

    double RpAu_jpsi_ppg228_cent_err[narm][18] = {
      {0.07, 0.04, 0.04, 0.04, 0.03, 0.04, 0.06, 0.05, 0.05, 0.06, 0.07, 0.07, 0.08, 0.10, 0.11, 0.12, 0.13, 0.24},
      {0.05, 0.03, 0.03, 0.02, 0.02, 0.03, 0.04, 0.03, 0.03, 0.04, 0.05, 0.06, 0.06, 0.07, 0.07, 0.10, 0.12, 0.11}};

 // double RpAu_jpsi_ppg228_cent[narm][18] = {
 //      {10.69, 10.78, 0.68, 0.78, 0.88, 0.78 , 0.92, 0.92, 1.02, 1.02, 1.12, 1.06, 1.06, 1.22, 1.16, 1.13, 1.31 , 1.71},
 //      {10.65, 10.61, 0.63, 0.64, 0.66, 0.71, 0.78, 0.73, 0.83, 0.83, 0.88, 0.92, 0.97, 1.05, 0.98, 1.07, 1.20, 0.96}};

 //    double RpAu_jpsi_ppg228_cent_err[narm][18] = {
 //      {10.07, 10.04, 0.04, 0.04, 0.03, 0.04, 0.06, 0.05, 0.05, 0.06, 0.07, 0.07, 0.08, 0.10, 0.11, 0.12, 0.13, 0.24},
 //      {10.05, 10.03, 0.03, 0.02, 0.02, 0.03, 0.04, 0.03, 0.03, 0.04, 0.05, 0.06, 0.06, 0.07, 0.07, 0.10, 0.12, 0.11}};

    // for(int arm = 0; arm < 2; arm++)
    //   {
    // 	for(int bin = 0; bin < 4; bin++)
    // 	  {
    // 	    abs1 = 0;
    // 	    abs2 = 0;

    // 	    RpAu_jpsi_ppg228_cent[arm][bin] =( jpsi_ppg228_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_jpsi_ppg228_inv_yield[arm] ) ;
    // 	    abs1 = jpsi_ppg228_pAu_inv_yield_err[arm][bin] / jpsi_ppg228_pAu_inv_yield[arm][bin];
    // 	    abs2 = pp_jpsi_ppg228_inv_yield_err[arm] / pp_jpsi_ppg228_inv_yield[arm];
    // 	    RpAu_jpsi_ppg228_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_jpsi_ppg228_cent[arm][bin];  
    // 	  }
    //   }
     
    //======================================
  
  // with badprob cut N and S
     double jpsi_sng[narm][nbins] = {
      {1634,1139,791,547},
      {1624,1279,1010,840}};

    double jpsi_sng_err[narm][nbins] = {
      {74,63,45,37},
      {53,49,42,38}};
   
    //correct
    double jpsi_sng_acc_eff[2][4] = {0};
    double jpsi_sng_trig_eff[2][4] = {0};
  
    {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_tracks/sanghoon_simulations/Run15pAu_sng_jpsi_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int cent =0; cent < 4; cent++)
	 jpsi_sng_acc_eff[0][cent] = e1->GetEfficiency(cent+1);
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_tracks/sanghoon_simulations/Run15pAu_sng_jpsi_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int cent =0; cent < 4; cent++)
	 jpsi_sng_acc_eff[1][cent] = e2->GetEfficiency(cent+1);
	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sng/sanghoon_simulations/Run15pAu_sng_jpsi_trigeff_S.root");
       file_data1->GetObject("h1",f1);
       
       jpsi_sng_trig_eff[0][0] = f1->Eval(10);
       jpsi_sng_trig_eff[0][1] = f1->Eval(30);
       jpsi_sng_trig_eff[0][2] = f1->Eval(50);
       jpsi_sng_trig_eff[0][3] = f1->Eval(72);
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sng/sanghoon_simulations/Run15pAu_sng_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       jpsi_sng_trig_eff[1][0] = f2->Eval(10);
       jpsi_sng_trig_eff[1][1] = f2->Eval(30);
       jpsi_sng_trig_eff[1][2] = f2->Eval(50);
       jpsi_sng_trig_eff[1][3] = f2->Eval(72);
	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


 //=================
    
     double jpsi_sng_pAu_inv_yield[narm][nbins] = {0};
     double jpsi_sng_pAu_inv_yield_err[narm][nbins] = {0};
     
     // calculcate pAu inv. yield for jpsi sng
     for(int arm = 0; arm < 2; arm++)
       {
	 for(int bin = 0; bin < 4; bin++)
	   {
	     FACTOR = 0;
	     jpsi_sng_pAu_inv_yield[arm][bin] =  bias[bin] * jpsi_sng[arm][bin] / ( ( FACTOR + jpsi_sng_trig_eff[arm][bin] ) * ( FACTOR + jpsi_sng_acc_eff[arm][bin] ) * pAu_MB[arm][bin]  );
	     jpsi_sng_pAu_inv_yield_err[arm][bin] =  bias[bin] * jpsi_sng_err[arm][bin] / ( ( FACTOR + jpsi_sng_trig_eff[arm][bin] ) * ( FACTOR + jpsi_sng_acc_eff[arm][bin] ) * pAu_MB[arm][bin] );
	     // cout << "jpsi_sng_pAu inv yield: " << jpsi_sng_pAu_inv_yield[arm][bin] << "+\-" << jpsi_sng_pAu_inv_yield_err[arm][bin] << endl;
	   }
       }
     
     // jan 20
     double pp_jpsi_sng_inv_yield[narm] = {0}; 
     double pp_jpsi_sng_inv_yield_err[narm] = {0}; 
     double pp_jpsi_sng[narm] = {5326,5044}; 
     double pp_jpsi_sng_err[narm] = {130,115};   


     double pp_jpsi_sng_acc_eff[narm] = {0};  //0.0102,0.0086}; // *1.98/2.12 = 0.0202,0.0182
     double pp_jpsi_sng_trig_eff[narm] = {0};   // 0.7218,0.6758}; 

 // RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_tracks/sanghoon_simulations_final/Run15pp_sng_jpsi_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       //  pp_jpsi_sng_acc_eff[0] = t1->GetBinContent(1);
       pp_jpsi_sng_acc_eff[0] = 0.0071;
       // pp_jpsi_sng_acc_eff[0] = 0.0118;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_tracks/sanghoon_simulations_final/Run15pp_sng_jpsi_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       //pp_jpsi_sng_acc_eff[1] = t2->GetBinContent(1);
       pp_jpsi_sng_acc_eff[1] = 0.0083;
       //pp_jpsi_sng_acc_eff[1] = 0.0116;
       
       file_data1->Close();
       file_data2->Close();
       //////////////////////////////////////////
     }
     // RUN15pp
     {
     
       //////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sng/sanghoon_simulations_final/Run15pp_sng_jpsi_trigeff_S.root");
       file_data1->GetObject("geff_y_int",t1);
       //pp_jpsi_sng_trig_eff[0] = t1->GetBinContent(1);
       pp_jpsi_sng_trig_eff[0] = 0.7237;
       // pp_jpsi_sng_trig_eff[0] = 0.7281;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sng/sanghoon_simulations_final/Run15pp_sng_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //pp_jpsi_sng_trig_eff[1] = t2->GetBinContent(1);
       pp_jpsi_sng_trig_eff[1] = 0.6972;
       // pp_jpsi_sng_trig_eff[1] = 0.6961;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }
  
     // calculate pp inv. yield for jpsi sng
    for(int arm = 0; arm < 2; arm++)
      {
	FACTOR = 0;
 	pp_jpsi_sng_inv_yield[arm] = pp_jpsi_sng[arm]*eff_bbc / ( ( FACTOR +  pp_jpsi_sng_trig_eff[arm] ) * ( FACTOR + pp_jpsi_sng_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
 	pp_jpsi_sng_inv_yield_err[arm] = pp_jpsi_sng_err[arm]*eff_bbc / ( ( FACTOR + pp_jpsi_sng_trig_eff[arm] ) * ( FACTOR + pp_jpsi_sng_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi);
 	// cout << "pp jpsi_sng sng inv yield: " << pp_jpsi_sng_inv_yield[arm] << "+\-" << pp_jpsi_sng_inv_yield_err[arm] << endl;
      }

    double RpAu_jpsi_sng_cent[narm][nbins] = {0};
    double RpAu_jpsi_sng_cent_err[narm][nbins] = {0};
    
    // jpsi sng
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_jpsi_sng_cent[arm][bin] =( jpsi_sng_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_jpsi_sng_inv_yield[arm] ) ;
 	    abs1 = jpsi_sng_pAu_inv_yield_err[arm][bin] / jpsi_sng_pAu_inv_yield[arm][bin];
 	    abs2 = pp_jpsi_sng_inv_yield_err[arm] / pp_jpsi_sng_inv_yield[arm];
 	    RpAu_jpsi_sng_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_jpsi_sng_cent[arm][bin];

	    //cout <<  RpAu_jpsi_sng_cent_err[arm][bin] / RpAu_jpsi_sng_cent[arm][bin] << endl;
 	  }
      }

 // 4a for sng+dbl Jpsi
    //====================================================================================
//====================================================================================

//
    //double sngdbl_jpsi[narm][nbins] = {
    // {900, 4385, 3396, 1513, 869, 2166},
    // {1000, 4606, 3281, 1714, 957, 2410}};
  
    //double sngdbl_jpsi_err[narm][nbins] = {
      // {35, 77, 64, 44, 34, 52},
      // {35, 75, 62, 46, 34, 55}};

    double sngdbl_jpsi[narm][nbins] = {
      {10000, 4385, 3396, 100000,10000, 2166},
      {100000, 4606, 3281, 1000000,1000000, 2410}};

    double sngdbl_jpsi_err[narm][nbins] = {
      {35, 77, 64, 44, 34, 52},
      {35, 75, 62, 46, 34, 55}};
    
    // correct
    double sngdbl_jpsi_trig_eff[narm][nbins] = {0};
    double sngdbl_jpsi_acc_eff[narm][nbins] = {0};
      
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_jpsi_acc_eff[0][pt] = e1->GetEfficiency(pt+1);

       sngdbl_jpsi_acc_eff[0][0] = 0.0297;
       sngdbl_jpsi_acc_eff[0][1] = 0.0267;
       sngdbl_jpsi_acc_eff[0][2] = 0.0244;
       sngdbl_jpsi_acc_eff[0][3] = 0.0255;
       sngdbl_jpsi_acc_eff[0][4] = 0.0307;
       sngdbl_jpsi_acc_eff[0][5] = 0.0264;
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_jpsi_acc_eff[1][pt] = e2->GetEfficiency(pt+1);

       sngdbl_jpsi_acc_eff[1][0] = 0.0393;
       sngdbl_jpsi_acc_eff[1][1] = 0.0347;
       sngdbl_jpsi_acc_eff[1][2] = 0.0305;
       sngdbl_jpsi_acc_eff[1][3] = 0.0343;
       sngdbl_jpsi_acc_eff[1][4] = 0.0420;
       sngdbl_jpsi_acc_eff[1][5] = 0.0356;
	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_S.root");
       file_data1->GetObject("h1",f1);
       
       sngdbl_jpsi_trig_eff[0][0] = 0.6357;
       sngdbl_jpsi_trig_eff[0][1] = 0.6602;
       sngdbl_jpsi_trig_eff[0][2] = 0.6748;
       sngdbl_jpsi_trig_eff[0][3] = 0.6592;
       sngdbl_jpsi_trig_eff[0][4] = 0.6601;
       sngdbl_jpsi_trig_eff[0][5] = 0.6585;
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       sngdbl_jpsi_trig_eff[1][0] = 0.6971;
       sngdbl_jpsi_trig_eff[1][1] = 0.7299;
       sngdbl_jpsi_trig_eff[1][2] = 0.7639;
       sngdbl_jpsi_trig_eff[1][3] = 0.7732;
       sngdbl_jpsi_trig_eff[1][4] = 0.7880;
       sngdbl_jpsi_trig_eff[1][5] = 0.7776;
	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


    double sngdbl_jpsi_pAu_inv_yield[narm][nbins] = {0};
    double sngdbl_jpsi_pAu_inv_yield_err[narm][nbins] = {0};

    // calculcate pAu inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < nbins; bin++)
	  {
	   
	    sngdbl_jpsi_pAu_inv_yield[arm][bin] = bias[bin] * sngdbl_jpsi[arm][bin] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm][bin] ) * ( FACTOR + sngdbl_jpsi_acc_eff[arm][bin] ) * pAu_MB[arm][bin] * 2 * pi * delta_pt[bin] * pAu_jpsi_binshift[bin] );
	    sngdbl_jpsi_pAu_inv_yield_err[arm][bin] = bias[bin] * sngdbl_jpsi_err[arm][bin] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm][bin]) * ( FACTOR + sngdbl_jpsi_acc_eff[arm][bin]) * pAu_MB[arm][bin] * 2 * pi * delta_pt[bin] * pAu_jpsi_binshift[bin] );
	    cout << "sngdbl_jpsi_pAu inv yield: " << sngdbl_jpsi_pAu_inv_yield[arm][bin] << "+\-" << sngdbl_jpsi_pAu_inv_yield_err[arm][bin] << endl;
	  }
      }
    
    double pp_sngdbl_jpsi_inv_yield[narm][nbins] = {0}; 
    double pp_sngdbl_jpsi_inv_yield_err[narm][nbins] = {0}; 

    double pp_sngdbl_jpsi[narm][nbins] = {
      {1864, 8471, 5228, 2102, 971, 2863},
      {1327, 5863, 3691, 1590, 840, 2151}};

    double pp_sngdbl_jpsi_err[narm][nbins] = {
      {48, 103, 80, 109, 45, 150},
      {40, 104, 68, 43, 32, 72}};

    double pp_sngdbl_jpsi_acc_eff[narm][nbins] = {0};
    double pp_sngdbl_jpsi_trig_eff[narm][nbins] = {0};

    //RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       //pp_sngdbl_jpsi_acc_eff[0] = t1->GetBinContent(1);
       //pp_sngdbl_jpsi_acc_eff[0] = 0.0337;
       pp_sngdbl_jpsi_acc_eff[0][0] = 0.0386;
       pp_sngdbl_jpsi_acc_eff[0][1] = 0.0348;
       pp_sngdbl_jpsi_acc_eff[0][2] = 0.0310;
       pp_sngdbl_jpsi_acc_eff[0][3] = 0.0334;
       pp_sngdbl_jpsi_acc_eff[0][4] = 0.0388; 
       pp_sngdbl_jpsi_acc_eff[0][5] = 0.0343;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       //pp_sngdbl_jpsi_acc_eff[1] = t2->GetBinContent(1);
       //pp_sngdbl_jpsi_acc_eff[1] = 0.0281;
       pp_sngdbl_jpsi_acc_eff[1][0] = 0.0322;
       pp_sngdbl_jpsi_acc_eff[1][1] = 0.0286;
       pp_sngdbl_jpsi_acc_eff[1][2] = 0.0258;
       pp_sngdbl_jpsi_acc_eff[1][3] = 0.0278;
       pp_sngdbl_jpsi_acc_eff[1][4] = 0.0356; // double check
       pp_sngdbl_jpsi_acc_eff[1][5] = 0.0292;
       
       file_data1->Close();
       file_data2->Close();
       //////////////////////////////////////////
     }
     // RUN15pp
     {
     
       //////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_S.root");
       file_data1->GetObject("geff_y_int",t1);
       //pp_sngdbl_jpsi_trig_eff[0] = t1->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[0][0] = 0.7407;
       pp_sngdbl_jpsi_trig_eff[0][1] = 0.7421;
       pp_sngdbl_jpsi_trig_eff[0][2] = 0.7522;
       pp_sngdbl_jpsi_trig_eff[0][3] = 0.7413;
       pp_sngdbl_jpsi_trig_eff[0][4] = 0.7475;
       pp_sngdbl_jpsi_trig_eff[0][5] = 0.7425;
 
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //ypp_sngdbl_jpsi_trig_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[1][0] = 0.6626;
       pp_sngdbl_jpsi_trig_eff[1][1] = 0.7002;
       pp_sngdbl_jpsi_trig_eff[1][2] = 0.7258;
       pp_sngdbl_jpsi_trig_eff[1][3] = 0.7384;
       pp_sngdbl_jpsi_trig_eff[1][4] = 0.7517;
       pp_sngdbl_jpsi_trig_eff[1][5] = 0.7421;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }


    // calculate pp inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bins = 0; bins < nbins; bins++)
	  {
	    pp_sngdbl_jpsi_inv_yield[arm][bins] = pp_sngdbl_jpsi[arm][bins]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm][bins] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm][bins] ) * pp_MB[arm] * eff_jpsi  * 2 * pi * delta_pt[bins] * pp_jpsi_binshift[bins] );
	    pp_sngdbl_jpsi_inv_yield_err[arm][bins] = pp_sngdbl_jpsi_err[arm][bins]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm][bins] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm][bins] ) * pp_MB[arm] * eff_jpsi  * 2 * pi * delta_pt[bins] * pp_jpsi_binshift[bins] );
	    cout << "pp sngdbl_jpsi fvtx inv yield: " << pp_sngdbl_jpsi_inv_yield[arm][bins] << "+\-" << pp_sngdbl_jpsi_inv_yield_err[arm][bins] << endl;
	  }
      }


    double RpAu_sngdbl_jpsi_cent[narm][nbins] = {0};
    double RpAu_sngdbl_jpsi_cent_err[narm][nbins] = {0};
    
    // jpsi sngdbl
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < nbins; bin++)
	  {
	    abs1 = 0;
	    abs2 = 0;

	    RpAu_sngdbl_jpsi_cent[arm][bin] =( sngdbl_jpsi_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_sngdbl_jpsi_inv_yield[arm][bin] ) ;
	    abs1 = sngdbl_jpsi_pAu_inv_yield_err[arm][bin] / sngdbl_jpsi_pAu_inv_yield[arm][bin];
	    abs2 = pp_sngdbl_jpsi_inv_yield_err[arm][bin] / pp_sngdbl_jpsi_inv_yield[arm][bin];
	    RpAu_sngdbl_jpsi_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_jpsi_cent[arm][bin];
	  }
      }

// 4a for sng+dbl psi2s
    //====================================================================================

    //
    // double sngdbl_psi2s[narm][nbins] = {
    //   {12, 65, 92, 32, 36, 52},
    //   {18, 101, 100, 64, 49, 91}};
 
    // double sngdbl_psi2s_err[narm][nbins] = {
    //   {9, 23, 22, 14, 14, 18},
    //   {9, 20, 18, 14, 11, 17}};

    // changed binning
  double sngdbl_psi2s[narm][nbins] = {
    {12000, 57, 89, 10000, 16000, 52},  // 0-0.5, 0.5-1.5, 1.5-2.5, N.A, 4.5 - 6.5 , 2.5 - 4.5
    //   {18, 101, 100, 100064, 25, 91}};
  {18000, 101, 107, 100064, 25000, 91}};
     
    double sngdbl_psi2s_err[narm][nbins] = {
      {9, 23, 22, 14, 9, 18},
      {9, 20, 18, 14, 6, 17}};
 
    double sngdbl_psi2s_acc_eff[narm][nbins] = {0};
    double sngdbl_psi2s_trig_eff[narm][nbins] = {0};
    
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_psi2s_acc_eff[0][pt] = e1->GetEfficiency(pt+1);

       sngdbl_psi2s_acc_eff[0][0] = 0.0367;
       sngdbl_psi2s_acc_eff[0][1] = 0.0343;
       sngdbl_psi2s_acc_eff[0][2] = 0.0315;
       sngdbl_psi2s_acc_eff[0][3] = 0.0295;
       sngdbl_psi2s_acc_eff[0][4] = 0.0373; // 4.5 - 6.5
	 sngdbl_psi2s_acc_eff[0][5] = 0.0301;  // 2.5 - 4.5
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_psi2s_acc_eff[1][pt] = e2->GetEfficiency(pt+1);

       sngdbl_psi2s_acc_eff[1][0] = 0.0523;
       sngdbl_psi2s_acc_eff[1][1] = 0.0462;
       sngdbl_psi2s_acc_eff[1][2] = 0.0428;
       sngdbl_psi2s_acc_eff[1][3] = 0.0399;
       sngdbl_psi2s_acc_eff[1][4] = 0.0477;
       sngdbl_psi2s_acc_eff[1][5] = 0.0410;
	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_S.root");
       file_data1->GetObject("h2",f1);
       
       sngdbl_psi2s_trig_eff[0][0] = 0.6693;
       sngdbl_psi2s_trig_eff[0][1] = 0.6681;
       sngdbl_psi2s_trig_eff[0][2] = 0.6718;
       sngdbl_psi2s_trig_eff[0][3] = 0.6760;
       sngdbl_psi2s_trig_eff[0][4] = 0.6850;  // 4.5 - 6.5
       sngdbl_psi2s_trig_eff[0][5] = 0.6742; // 2.5 - 4.5
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("h2",f2);
       
       sngdbl_psi2s_trig_eff[1][0] = 0.7390;
       sngdbl_psi2s_trig_eff[1][1] = 0.7515;
       sngdbl_psi2s_trig_eff[1][2] = 0.7751;
       sngdbl_psi2s_trig_eff[1][3] = 0.7833;
       sngdbl_psi2s_trig_eff[1][4] = 0.7990;  // 4.5 - 6.5
       sngdbl_psi2s_trig_eff[1][5] = 0.7849;
	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }

    double sngdbl_psi2s_pAu_inv_yield[narm][nbins] = {0};
    double sngdbl_psi2s_pAu_inv_yield_err[narm][nbins] = {0};

    // calculcate pAu inv. yield for psi2s fvtx
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < nbins; bin++)
 	  {
 	    sngdbl_psi2s_pAu_inv_yield[arm][bin] = bias[bin] * sngdbl_psi2s[arm][bin] / ( ( FACTOR + sngdbl_psi2s_trig_eff[arm][bin] ) * ( FACTOR + sngdbl_psi2s_acc_eff[arm][bin] ) * pAu_MB[arm][bin] * 2 * pi * delta_pt[bin] * pAu_psi2s_binshift[bin] );
 	    sngdbl_psi2s_pAu_inv_yield_err[arm][bin] = bias[bin] * sngdbl_psi2s_err[arm][bin] / ( ( FACTOR + sngdbl_psi2s_trig_eff[arm][bin] ) * ( FACTOR + sngdbl_psi2s_acc_eff[arm][bin] ) * pAu_MB[arm][bin] * 2 * pi * delta_pt[bin] * pAu_psi2s_binshift[bin] );
 	    // cout << "sngdbl_psi2s_pAu inv yield: " << sngdbl_psi2s_pAu_inv_yield[arm][bin] << "+\-" << sngdbl_psi2s_pAu_inv_yield_err[arm][bin] << endl;
 	  }
      }

  
    double pp_sngdbl_psi2s_inv_yield[narm][nbins] = {0}; 
    double pp_sngdbl_psi2s_inv_yield_err[narm][nbins] = {0}; 
    
    // averaged 2nd gaussian
    // double pp_sngdbl_psi2s[narm][nbins] = {
    //   {43, 213, 201, 85, 19, 139},  // 0-0.5, 0.5-1.5, 1.5-2.5, N.A, 4.5 - 6.5 , 2.5 - 4.5
    //   {13, 175, 175, 48, 27, 96}};

    // double pp_sngdbl_psi2s_err[narm][nbins] = {
    //   {11, 24, 21, 15, 8, 19},.q

    //   {9, 22, 20, 12, 7, 16}};

  // individual 2nd gausian fits
    double pp_sngdbl_psi2s[narm][nbins] = {
      {43, 218, 197, 85, 50, 138},  // 0.0 - 0.5, 0.5-1.5, 1.5-2.5, -, -, 2.5-4.5
      {13, 173, 166, 48, 59, 96}};

    double pp_sngdbl_psi2s_err[narm][nbins] = {
      {11, 24, 22, 15, 12, 19},
      {9, 22, 20, 12, 13, 16}};

    double pp_sngdbl_psi2s_acc_eff[narm][nbins] = {0}; 
    double pp_sngdbl_psi2s_trig_eff[narm][nbins] = {0}; 
    
    //RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       // pp_sngdbl_psi2s_acc_eff[0] = t1->GetBinContent(1);
       // pp_sngdbl_psi2s_acc_eff[0] = 0.0417;
       pp_sngdbl_psi2s_acc_eff[0][0] = 0.0468;
       pp_sngdbl_psi2s_acc_eff[0][1] = 0.0436;
       pp_sngdbl_psi2s_acc_eff[0][2] = 0.0404;
       pp_sngdbl_psi2s_acc_eff[0][3] = 0.0384;
       pp_sngdbl_psi2s_acc_eff[0][4] = 0.0443; // 4.5 - 6.5
	 pp_sngdbl_psi2s_acc_eff[0][5] = 0.0394;  // 2.5 - 4.5
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       // pp_sngdbl_psi2s_acc_eff[1] = t2->GetBinContent(1);
       //pp_sngdbl_psi2s_acc_eff[1] = 0.0367;
       pp_sngdbl_psi2s_acc_eff[1][0] = 0.0434;
       pp_sngdbl_psi2s_acc_eff[1][1] = 0.0390;
       pp_sngdbl_psi2s_acc_eff[1][2] = 0.0354;
       pp_sngdbl_psi2s_acc_eff[1][3] = 0.0332;
       pp_sngdbl_psi2s_acc_eff[1][4] = 0.0395; 
       pp_sngdbl_psi2s_acc_eff[1][5] = 0.0338;
              
       file_data1->Close();
       file_data2->Close();
       //////////////////////////////////////////
     }
     // RUN15pp
     {
     
       //////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_trigeff_S.root");
       file_data1->GetObject("geff_y_int",t1);
       // pp_sngdbl_psi2s_trig_eff[0] = t1->GetBinContent(1);
       // pp_sngdbl_psi2s_trig_eff[0] = 0.7620;
       pp_sngdbl_psi2s_trig_eff[0][0] = 0.7506;
       pp_sngdbl_psi2s_trig_eff[0][1] = 0.7562;
       pp_sngdbl_psi2s_trig_eff[0][2] = 0.7703;
       pp_sngdbl_psi2s_trig_eff[0][3] = 0.7633;
       pp_sngdbl_psi2s_trig_eff[0][4] = 0.7455;  // 4.5 - 6.5
       pp_sngdbl_psi2s_trig_eff[0][5] = 0.7586; // 2.5 - 4.5
          
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //pp_sngdbl_psi2s_trig_eff[1] = t2->GetBinContent(1);
       // pp_sngdbl_psi2s_trig_eff[1] = 0.7352;
       pp_sngdbl_psi2s_trig_eff[1][0] = 0.7056;
       pp_sngdbl_psi2s_trig_eff[1][1] = 0.7149;
       pp_sngdbl_psi2s_trig_eff[1][2] = 0.7371;
       pp_sngdbl_psi2s_trig_eff[1][3] = 0.7576;
       pp_sngdbl_psi2s_trig_eff[1][4] = 0.7915;
       pp_sngdbl_psi2s_trig_eff[1][5] = 0.7634;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }


    // calculate pp inv. yield for psi2s fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bins = 0; bins < nbins; bins++)
	  {
	    pp_sngdbl_psi2s_inv_yield[arm][bins] = pp_sngdbl_psi2s[arm][bins] * eff_bbc / ( (FACTOR + pp_sngdbl_psi2s_trig_eff[arm][bins] ) * ( FACTOR + pp_sngdbl_psi2s_acc_eff[arm][bins] ) * pp_MB[arm] * eff_jpsi * 2 * pi * delta_pt[bins] * pp_psi2s_binshift[bins] );
	    pp_sngdbl_psi2s_inv_yield_err[arm][bins] = pp_sngdbl_psi2s_err[arm][bins] * eff_bbc / ( ( FACTOR + pp_sngdbl_psi2s_trig_eff[arm][bins] ) * ( FACTOR + pp_sngdbl_psi2s_acc_eff[arm][bins] ) * pp_MB[arm] * eff_jpsi * 2 * pi * delta_pt[bins] * pp_psi2s_binshift[bins] );
	    cout << "pp sngdbl_psi2s sngdbl inv yield: " << pp_sngdbl_psi2s_inv_yield[arm][bins] << endl; //<< "+\-" << pp_sngdbl_psi2s_inv_yield_err[arm][nbins] << endl;
	  }
      }
     double RpAu_sngdbl_psi2s_cent[narm][nbins] = {0};
     double RpAu_sngdbl_psi2s_cent_err[narm][nbins] = {0};
    
    // psi2s sngdbl
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < nbins; bin++)
 	  {
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_sngdbl_psi2s_cent[arm][bin] =( sngdbl_psi2s_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_sngdbl_psi2s_inv_yield[arm][bin] ) ;
 	    abs1 = sngdbl_psi2s_pAu_inv_yield_err[arm][bin] / sngdbl_psi2s_pAu_inv_yield[arm][bin];
 	    abs2 = pp_sngdbl_psi2s_inv_yield_err[arm][bin] / pp_sngdbl_psi2s_inv_yield[arm][bin];
 	    RpAu_sngdbl_psi2s_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_psi2s_cent[arm][bin];
	    cout << "RpAu psi2s[" << arm << "][" << bin << "]: " <<  RpAu_sngdbl_psi2s_cent[arm][bin] << ", ncoll: " << ncoll[bin] << ", raw psi2s yield: " <<  sngdbl_psi2s[arm][bin] << ", psi2s pAu inv yield: " <<  sngdbl_psi2s_pAu_inv_yield[arm][bin] <<  endl;
 	  }
      }

    ////////////////////////////// calcualte weighted MB sum
    double w_ave[2] = {0};
    double num[2] = {0};
    double denom[2] = {0};
    double weight =(4.7)/(0.86) ;

    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < nbins; bin++)
 	  {
	    if(bin !=3)
	      {
		num[arm] += RpAu_sngdbl_psi2s_cent[arm][bin] * pp_sngdbl_psi2s_inv_yield[arm][bin] * weight ;
		denom[arm] +=  pp_sngdbl_psi2s_inv_yield[arm][bin] * weight;
	      }
	    w_ave[arm] = num[arm]/denom[arm];
	 
	  }
      }
	////////////////////////////// calcualte weighted MB sum
   

  
  //====================================================================================
    TGraphAsymmErrors *p1; 
    TGraphAsymmErrors *p2;  // J[psi bkwd
    TGraphAsymmErrors *p3;
    TGraphAsymmErrors *p4;
    TGraphAsymmErrors *p5;  // CNM fwd
    TGraphAsymmErrors *p6;  // CNM bkwd

    double jpsi_transport[2][6] = {0};
    double psi2s_transport[2][6] = {0};
    double cnm_effects[2][6] = {0};

    // Ralph Rapp predciitons
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
      file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/tony_bestfit_parameters/sanghoon_corrbg/shao_psi2s_predictions/psi2S_RpAu_MB_predictions.root");
      file_data1->GetObject("RpAu_epps_MB_predictions_bkwd",p2);
      file_data1->GetObject("RpAu_epps_MB_predictions_fwd",p1);
   file_data1->GetObject("RpAu_ncteq_MB_predictions_bkwd",p4);
      file_data1->GetObject("RpAu_ncteq_MB_predictions_fwd",p3);
    
      file_data1->Close();
	
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


    const Double_t err_x[nbins] = {0,0,0,0,0,0};
    const Double_t err_x_228[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    const Double_t width = 0.145;
    const Double_t jpsi_ppg228_sys[narm][nbins] = {
      {0.04,0.04,0.05,0.04, 0.06, 0.06}, //  ->
      {0.05,0.07,0.07,0.08, 0.06, 0.06}};

    
    double pT_array_228[18] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375};
       
    double sys_err_228[2][18] = {
    {0.06, 0.07, 0.06, 0.07, 0.08, 0.07, 0.08, 0.09, 0.10, 0.10, 0.11, 0.10, 0.10, 0.12, 0.11, 0.10, 0.13, 0.17},
    {0.04,0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06, 0.07, 0.07, 0.06}};


   TGraphErrors *gr1 = new TGraphErrors(6,pT_array1,RpAu_sngdbl_psi2s_cent[1],err_x, RpAu_sngdbl_psi2s_cent_err[1]);
   TGraphErrors *gr3 = new TGraphErrors(6,pT_array2,RpAu_sngdbl_jpsi_cent[1],err_x, RpAu_sngdbl_jpsi_cent_err[1]);
   TGraphErrors *gr5 = new TGraphErrors(18,pT_array_228,RpAu_jpsi_ppg228_cent[1],err_x_228, RpAu_jpsi_ppg228_cent_err[1] );
   TGraphErrors *gr7 = new TGraphErrors(4,pT_array4,RpAu_jpsi_sng_cent[1],err_x,RpAu_jpsi_sng_cent_err[1]);
   TGraphErrors *gr9 = new TGraphErrors(6,pT_array4,RpAu_sngdbl_psi2s_cent[1],err_x, RpAu_sngdbl_psi2s_cent_err[1]);
   TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(6);
   TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(6);
  
   
   TGraphErrors *gr2 = new TGraphErrors(6,pT_array1,RpAu_sngdbl_psi2s_cent[1],err_x, RpAu_sngdbl_psi2s_cent_err[1]);
   TGraphErrors *gr4 = new TGraphErrors(6,pT_array2,RpAu_sngdbl_jpsi_cent[0],err_x, RpAu_sngdbl_jpsi_cent_err[0]);
   TGraphErrors *gr6 = new TGraphErrors(18,pT_array_228,RpAu_jpsi_ppg228_cent[0],err_x_228, RpAu_jpsi_ppg228_cent_err[0] );;
   TGraphErrors *gr8 = new TGraphErrors(4,pT_array4,RpAu_jpsi_sng_cent[0],err_x,RpAu_jpsi_sng_cent_err[0]);
    TGraphErrors *gr10 = new TGraphErrors(6,pT_array4,RpAu_sngdbl_psi2s_cent[0],err_x, RpAu_sngdbl_psi2s_cent_err[0]);
   TGraphAsymmErrors *gr4_sys = new TGraphAsymmErrors(6);
   TGraphAsymmErrors *gr6_sys = new TGraphAsymmErrors(6);
    

   for(int i = 0; i < 18; i++)
     {
       // gr3_sys->SetPoint(i, pT_array2[i],RpAu_sngdbl_jpsi_cent[1][i]);   // need to also include array index here as well as arm index
       // gr3_sys->SetPointError(i, width, width, jpsi_ppg228_sys[1][i], jpsi_ppg228_sys[1][i]);

       gr5_sys->SetPoint(i, pT_array_228[i],RpAu_jpsi_ppg228_cent[1][i]);   // need to also include array index here as well as arm index
       gr5_sys->SetPointError(i, width, width, sys_err_228[1][i], sys_err_228[1][i]);

       // gr4_sys->SetPoint(i, pT_array2[i],RpAu_sngdbl_jpsi_cent[0][i]);
       // gr4_sys->SetPointError(i, width, width, jpsi_ppg228_sys[0][i], jpsi_ppg228_sys[0][i]);

       gr6_sys->SetPoint(i, pT_array_228[i],RpAu_jpsi_ppg228_cent[0][i]);
       gr6_sys->SetPointError(i, width, width, sys_err_228[0][i], sys_err_228[0][i]);
	  
     }

      
    gr3_sys->SetMarkerColor(kCyan+1);                     
    gr3_sys->SetLineColor(kCyan+1);                         
    gr3_sys->SetFillStyle(0); 
    gr3_sys->SetLineWidth(1);

    gr4_sys->SetMarkerColor(kCyan+1);                     
    gr4_sys->SetLineColor(kCyan+1);                         
    gr4_sys->SetFillStyle(0); 
    gr4_sys->SetLineWidth(1);

    gr5_sys->SetMarkerColor(kBlack);                     
    gr5_sys->SetLineColor(kBlack);                         
    gr5_sys->SetFillStyle(0); 
    gr5_sys->SetLineWidth(1);

    gr6_sys->SetMarkerColor(kBlack);                     
    gr6_sys->SetLineColor(kBlack);                         
    gr6_sys->SetFillStyle(0); 
    gr6_sys->SetLineWidth(1);


    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,545);
  
    c1->SetTickx(1);
    c1->SetTicky(1);
   
    c1->cd();
    c1->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.15);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetLineWidth(2);
   
 
    gr1->SetMarkerColor(kBlue+1);
    gr1->SetTitle("");
    gr1->SetLineWidth(3);
    gr1->SetLineColor(kBlue+1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(2);
  
  
    gr1->GetXaxis()->SetTitleOffset(0.95);
    gr1->GetXaxis()->SetTitleFont(102);
    gr1->GetXaxis()->SetLabelFont(102);
    gr1->GetXaxis()->SetTitle("p_{T} [GeV/c^{2}]");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.06);
    
    gr1->GetXaxis()->SetRangeUser(0,6);
     gr1->GetXaxis()->SetNdivisions(5);
    // gr1->GetYaxis()->SetNdivisions(8);
   

    gr1->GetYaxis()->SetTitleOffset(0.9);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.075);
    gr1->GetYaxis()->SetLabelSize(0.065);
 
    gr1->GetYaxis()->SetTitle("R_{pAu}");
 
    gr1->GetYaxis()->SetRangeUser(0,2.5);
    gr1->GetYaxis()->SetNdivisions(7);
 
    gr2->SetMarkerColor(kBlue+1);
    gr2->SetLineWidth(3);
    gr2->SetLineColor(kBlue+1);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(2);                       

    gr3->SetMarkerColor(kCyan+1);
    gr3->SetLineWidth(3);
    gr3->SetLineColor(kCyan+1);
    gr3->SetMarkerStyle(20); 
    gr3->SetMarkerSize(2);                       

    gr4->SetMarkerColor(kCyan+1);
    gr4->SetLineWidth(3);
    gr4->SetLineColor(kCyan+1);
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerSize(2);                       

    gr5->SetMarkerColor(kBlack);                     
    gr5->SetLineColor(kBlack);                         
    gr5->SetLineWidth(3);                           
    gr5->SetMarkerStyle(21);                       
    // gr5->SetMarkerSize(1.75); 
    gr5->SetMarkerSize(2); 

    gr6->SetMarkerColor(kBlack);                     
    gr6->SetLineColor(kBlack);                         
    gr6->SetLineWidth(3);                           
    gr6->SetMarkerStyle(21);                       
    //  gr6->SetMarkerSize(1.75);                       
    gr6->SetMarkerSize(2);                       
    
    gr7->SetMarkerColor(kBlack);
    gr7->SetLineWidth(3);
    gr7->SetLineColor(kBlack);
    gr7->SetMarkerStyle(20);
    gr7->SetMarkerSize(2); 
  
    gr8->SetMarkerColor(kBlack);
    gr8->SetLineWidth(3);
    gr8->SetLineColor(kBlack);
    gr8->SetMarkerStyle(20);
    gr8->SetMarkerSize(2);

    gr9->SetMarkerColor(kBlue+1);
    gr9->SetLineWidth(3);
    gr9->SetLineColor(kBlue+1);
    gr9->SetMarkerStyle(20);
    gr9->SetMarkerSize(2); 
  
    gr10->SetMarkerColor(kBlue+1);
    gr10->SetLineWidth(3);
    gr10->SetLineColor(kBlue+1);
    gr10->SetMarkerStyle(20);
    gr10->SetMarkerSize(2);   

    /////////////////// predicitons
    p1->SetFillStyle(3344);
    p1->SetLineColor(kBlack);
    p1->SetFillColor(kBlack);
    p1->SetLineWidth(3);

    p2->SetFillStyle(3344);
    p2->SetLineColor(kBlack);
    p2->SetFillColor(kBlack);
    p2->SetLineWidth(3);
    
    // p5->SetFillStyle(3244);
    // p5->SetLineColor(kBlack);
    // p5->SetFillColor(kBlack);
    // p5->SetLineWidth(3);

    // p6->SetFillStyle(3244);
    // p6->SetLineColor(kBlack);
    // p6->SetFillColor(kBlack);
    // p6->SetLineWidth(3);

    p3->SetFillStyle(3001);
    p3->SetLineColor(kBlue+1);
    //  p3->SetFillColor(kBlue+1);
    p3->SetFillColorAlpha(4,0.3);
    p3->SetLineWidth(3);

    p4->SetFillStyle(3001);
    p4->SetLineColor(kBlue+1);
    //  p4->SetFillColor(kBlue+1);
    p4->SetFillColorAlpha(4,0.3);
    p4->SetLineWidth(3);
    
      
 
     if(north_arm)
      {
	gr1->SetMarkerColor(kWhite);
	gr1->SetLineColor(kWhite);
	gr1->Draw("AP");
	//	gr3->Draw("P");
       	gr5->Draw("P");
	//	gr7->Draw("P");
		gr9->Draw("P");
		//gr3_sys->Draw("e2same");
	       	gr5_sys->Draw("e2same");
		p1->Draw("e3same");
		p3->Draw("e3same");
		//	p5->Draw("e3same");
      }
     else
       {
	 gr1->SetMarkerColor(kWhite);
	 gr1->SetLineColor(kWhite);
	 gr1->Draw("AP");
	 // gr2->Draw("P");
	 //	 gr4->Draw("P");
		 gr6->Draw("P");
	 // gr8->Draw("P");
	 gr10->Draw("P");
	 //	 gr4_sys->Draw("e2same");
	  gr6_sys->Draw("e2same");
	 p2->Draw("e3same");
	 p4->Draw("e3same");
	 //p4->Draw("P");
       }

     TLatex l2;
    l2.SetTextSize(0.045);
    l2.SetTextAlign(12);
    l2.SetTextColor(kBlack);
  
    char text3[100];
    char text4[100];

    if(north_arm)
      sprintf(text3,"1.2 < y < 2.2");
    else
      sprintf(text3,"-2.2 < y < -1.2");

    sprintf(text4,"Mixed Events");
    l2.SetTextFont(102);  
    l2.DrawLatexNDC(0.64, 0.21, text3); // x0, y0
    // l2.SetTextColor(kBlue);
    // l2.DrawLatexNDC(0.70, 0.91, text4); // x0, y0

    TLine *line = new TLine(0,1,6,1);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");

    double global_sys = 0.121;  // RpAu J/psi PPG228 in both Norht and South

    TBox *box = new TBox(5.85,1-global_sys,6.0,1+global_sys);
    box->SetFillStyle(1000);
    box->SetFillColorAlpha(1,0.5);
    box->Draw();
    
    TLegend *leg2 = new TLegend(0.20, 0.62, 0.44, 0.93);  //(start x, start y, end x, end y)
    leg2->SetFillColor(19); 
    leg2->SetFillStyle(3003); 
    leg2->SetLineWidth(1);
    leg2->SetLineColor(0);
    leg2->SetTextSize(0.0425); 
    leg2->SetTextFont(102);

    if(north_arm)
      {

	//	leg2->AddEntry(gr1,"R_{pAu} J/#psi Fvtx Tracks", "p");
	leg2->AddEntry(gr5,"R_{pAu} J/#psi, p+Au #sqrt{s_{NN}}=200 GeV", "p");
	//	leg2->AddEntry(gr3,"R_{pAu} J/#psi Sng+Dbl", "p");
	//	leg2->AddEntry(gr7,"R_{pAu} J/#psi Sng Tracks", "p");
       	leg2->AddEntry(gr9,"R_{pAu} #psi(2S), p+Au #sqrt{s_{NN}}=200 GeV", "p");
	leg2->AddEntry(p1,"#psi(2S) nCTEQ15, Shao et al", "f");
	leg2->AddEntry(p3,"#psi(2S) EPPS16, Shao et al", "f");
	//	leg2->AddEntry(p1,"J/#psi Transport Model (Du&Rapp)", "f");
	//	leg2->AddEntry(p5,"CNM Effects (Du&Rapp)", "f");


      }
    else
      {

	//	leg2->AddEntry(gr2,"R_{pAu} J/#psi Fvtx Tracks", "p");
	leg2->AddEntry(gr6,"R_{pAu} J/#psi, p+Au #sqrt{s_{NN}}=200 GeV", "p");
	//	leg2->AddEntry(gr4,"R_{pAu} J/#psi Sng+Dbl", "p");
	//	leg2->AddEntry(gr8,"R_{pAu} J/#psi Sng Tracks", "p");
	leg2->AddEntry(gr10,"R_{pAu} #psi(2S), p+Au #sqrt{s_{NN}}=200 GeV", "p");
	leg2->AddEntry(p2,"#psi(2S) nCTEQ15, Shao et al", "f");
	leg2->AddEntry(p3,"#psi(2S) EPPS16, Shao et al", "f");
	//	leg2->AddEntry(p6,"CNM Effects (Du&Rapp)", "f");


      }
    leg2->Draw();


    cout << "North RpAu MB weighted average: " << w_ave[1] << endl;
    cout << "South RpAu MB weighted average: " << w_ave[0] << endl;

    if(north_arm)
      c1->SaveAs(Form("../psi2S_pdf/RpAu_psi2s_shao_plot_pT_fwd.pdf"));
    else
      c1->SaveAs(Form("../psi2S_pdf/RpAu_psi2s_shao_plot_pT_bkwd.pdf"));
}








