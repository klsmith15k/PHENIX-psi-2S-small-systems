

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


void pAu_psi2S_ratio_ncoll()
{
    gROOT->ForceStyle();

  bool north_arm = true;
  int arm = 0;
  double FACTOR = 0.0005;

  static const int nbins = 4;
  static const int narm = 2;

  if(north_arm)
    arm = 1;
    
  double eff_bbc = 0.55;
  double eff_jpsi = 0.79;

  double abs1 = 0;
  double abs2 = 0;

  double cent_bin[4] = {0.2, 0.2, 0.2, 0.24};

  const Double_t ncoll_array1[nbins] = {1.9,6.0,4.3,13};  // need all 4 points for ppg228 array
  const Double_t ncoll_array2[nbins] = {8.3,6.2,4.5,2.7};
  const Double_t ncoll_array3[nbins] = {8.2,6.1,4.4,2.6};
  const Double_t ncoll_array4[nbins] = {8.2,5.9,3.3,0};

  const Double_t ncoll_array6[6] = {12.7, 11.5, 9.81, 7.09, 4.28, 2.53};

  double bias_228[nbins] = {0.90,0.98,1.02,1.0};
  double ncoll_228[nbins] = {8.2,6.1,4.4,2.6};

  double bias[nbins] = {0.90,0.98,1.01,0};
  double ncoll[nbins] = {8.2,6.1,3.4,0};

  double pAu_MB_228[narm][nbins] = {
    {5.047069934*pow(10,10),4.99512054720000000e+10,4.97313174400000000e+10,5.96167043840000000e+10},
   {4.6346380480*pow(10,10),4.59845716160000000e+10,4.58064919040000000e+10,5.49289634240000000e+10} };

double pAu_MB[narm][nbins] = {
  {5.047069934*pow(10,10),4.99512054720000000e+10,1.09348021824*pow(10,11),0},
  {4.6346380480*pow(10,10),4.59845716160000000e+10,1.00735455328*pow(10,11), 0} };

 
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

    double jpsi_ppg228[narm][nbins] = {
      {4916,3208,2128,1386},
      {6078,5174,3968,3039}}; 

    double jpsi_ppg228_err[narm][nbins] = {
      {91,111,75,43},
      {104,90,80,69}};

    double jpsi_ppg228_trig_eff[narm][nbins] = {
      {0.7573,0.6802,0.6241,0.6156},
      {0.7397,0.7394,0.7370,0.7397}};  

    double jpsi_ppg228_acc_eff[narm][nbins] = {
      {0.0261,0.0270,0.0278,0.0276},
      {0.0502,0.0504,0.0507,0.0502}};
    
    double jpsi_ppg228_pAu_inv_yield[narm][nbins] = {0};
    double jpsi_ppg228_pAu_inv_yield_err[narm][nbins] = {0};

    double percent_err = 0;

    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < 4; bin++)
	  {
	    jpsi_ppg228_pAu_inv_yield[arm][bin] = bias_228[bin] * jpsi_ppg228[arm][bin] / ( jpsi_ppg228_trig_eff[arm][bin]*jpsi_ppg228_acc_eff[arm][bin] * pAu_MB_228[arm][bin] );
	    jpsi_ppg228_pAu_inv_yield_err[arm][bin] = jpsi_ppg228_err[arm][bin] / ( jpsi_ppg228_trig_eff[arm][bin]*jpsi_ppg228_acc_eff[arm][bin] * pAu_MB_228[arm][bin] );

	    // percent_err = ( jpsi_ppg228_pAu_inv_yield_err[arm][bin]/  jpsi_ppg228_pAu_inv_yield[arm][bin] ) * 100;
	    // cout <<  "J/psi inv yield ppg228: " << jpsi_ppg228_pAu_inv_yield[arm][bin] << ", and error: " <<  jpsi_ppg228_pAu_inv_yield_err[arm][bin] << ", and % err: " << percent_err << endl;
 	  }
      }

    ///////// weighted average of 40-60 and 60-84  ppg228 pAU inv yields
    double array_a = 0;
    double array_b = 0;
    double weight2 = 0;  // 0.85437;  // this is cent binwidth * Ncoll / bias correction for 40-60% bin
    double weight3 = 0;  // 0.6240;  // this is cent binwidth * Ncoll / bias correction for 60-84% bin

    double jpsi_ppg228_pAu_inv_yield_3bins[2][4] = {0};
    double jpsi_ppg228_pAu_inv_yield_3bins_err[2][4] = {0};

    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < 2; bin++)
	  {
	     jpsi_ppg228_pAu_inv_yield_3bins[arm][bin] = jpsi_ppg228_pAu_inv_yield[arm][bin] ;
	     jpsi_ppg228_pAu_inv_yield_3bins_err[arm][bin] = jpsi_ppg228_pAu_inv_yield_err[arm][bin] ;
	  }
	
	weight2 = cent_bin[2] * ncoll_228[2] / bias_228[2];
	weight3 = cent_bin[3] * ncoll_228[3] / bias_228[3];

	jpsi_ppg228_pAu_inv_yield_3bins[arm][2] = (  (jpsi_ppg228_pAu_inv_yield[arm][2] * weight2) + (jpsi_ppg228_pAu_inv_yield[arm][3] * weight3 ) ) / (  weight2 + weight3);
	jpsi_ppg228_pAu_inv_yield_3bins_err[arm][2] = (  (jpsi_ppg228_pAu_inv_yield_err[arm][2] * weight2) + (jpsi_ppg228_pAu_inv_yield_err[arm][3] * weight3 ) ) / (  weight2 + weight3 );

      }

    double pp_jpsi_ppg228_inv_yield[narm] = {0}; 
    double pp_jpsi_ppg228_inv_yield_err[narm] = {0};  
    double pp_jpsi_ppg228[narm] = {28511,31452};
    double pp_jpsi_ppg228_err[narm] = {203,215};
    double pp_jpsi_ppg228_trig_eff[narm] = {0.7499,0.7139};
    double pp_jpsi_ppg228_acc_eff[narm] = {0.0383,0.0453};

    for(int arm = 0; arm < 2; arm++)
      {
	pp_jpsi_ppg228_inv_yield[arm] = pp_jpsi_ppg228[arm]*eff_bbc / ( pp_jpsi_ppg228_trig_eff[arm]*pp_jpsi_ppg228_acc_eff[arm] * pp_MB_228[arm] * eff_jpsi );
	pp_jpsi_ppg228_inv_yield_err[arm] = pp_jpsi_ppg228_err[arm]*eff_bbc / ( pp_jpsi_ppg228_trig_eff[arm]*pp_jpsi_ppg228_acc_eff[arm] * pp_MB_228[arm] * eff_jpsi );
      }

    double RpAu_jpsi_ppg228_cent[narm][nbins] = {0};
    double RpAu_jpsi_ppg228_cent_err[narm][nbins] = {0};

    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_jpsi_ppg228_cent[arm][bin] =( jpsi_ppg228_pAu_inv_yield[arm][bin] ) / ( ncoll_228[bin] * pp_jpsi_ppg228_inv_yield[arm] ) ;
 	    abs1 = jpsi_ppg228_pAu_inv_yield_err[arm][bin] / jpsi_ppg228_pAu_inv_yield[arm][bin];
 	    abs2 = pp_jpsi_ppg228_inv_yield_err[arm] / pp_jpsi_ppg228_inv_yield[arm];
 	    RpAu_jpsi_ppg228_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_jpsi_ppg228_cent[arm][bin];  
 	  }
      }
     
    //======================================
  
      // NO badprob cut
    // double jpsi_sng[narm][nbins] = {
    //{2267,1652,1105,815},
    //{2281,1877,1444,1193}};

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

  
  
    // double sngdbl_jpsi[narm][nbins] = {
    // {5324,3587,2354,1615},
      //  {4606,3844,2975,2320}}; no sng probcut
    // {3932,3250,2535,1976}};  // with sng badprobcut

  double sngdbl_jpsi[narm][nbins] = {
    {4670,3115,2029,1345},
    {3937,3245,2539,1967}};
  
    double sngdbl_jpsi_err[narm][nbins] = {
      {98,77,56,45},
      {77,300,79,46}};
     
    // correct
    double sngdbl_jpsi_trig_eff[2][4] = {0};
    double sngdbl_jpsi_acc_eff[2][4] = {0};
      
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int cent =0; cent < 4; cent++)
	 sngdbl_jpsi_acc_eff[0][cent] = e1->GetEfficiency(cent+1);
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int cent =0; cent < 4; cent++)
	 sngdbl_jpsi_acc_eff[1][cent] = e2->GetEfficiency(cent+1);
	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_S.root");
       file_data1->GetObject("h1",f1);
       
       sngdbl_jpsi_trig_eff[0][0] = f1->Eval(10);
       sngdbl_jpsi_trig_eff[0][1] = f1->Eval(30);
       sngdbl_jpsi_trig_eff[0][2] = f1->Eval(50);
       sngdbl_jpsi_trig_eff[0][3] = f1->Eval(72);
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       sngdbl_jpsi_trig_eff[1][0] = f2->Eval(10);
       sngdbl_jpsi_trig_eff[1][1] = f2->Eval(30);
       sngdbl_jpsi_trig_eff[1][2] = f2->Eval(50);
       sngdbl_jpsi_trig_eff[1][3] = f2->Eval(72);
	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


    double sngdbl_jpsi_pAu_inv_yield[narm][nbins] = {0};
    double sngdbl_jpsi_pAu_inv_yield_err[narm][nbins] = {0};

    // calculcate pAu inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < 4; bin++)
	  {
	    sngdbl_jpsi_pAu_inv_yield[arm][bin] = bias[bin] * sngdbl_jpsi[arm][bin] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm][bin] ) * ( FACTOR + sngdbl_jpsi_acc_eff[arm][bin] ) * pAu_MB[arm][bin] );
	    sngdbl_jpsi_pAu_inv_yield_err[arm][bin] = bias[bin] * sngdbl_jpsi_err[arm][bin] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm][bin]) * ( FACTOR + sngdbl_jpsi_acc_eff[arm][bin]) * pAu_MB[arm][bin] );
	    cout << "sngdbl_jpsi_pAu inv yield: " << sngdbl_jpsi_pAu_inv_yield[arm][bin] << "+\-" << sngdbl_jpsi_pAu_inv_yield_err[arm][bin] << endl;
	  }
      }
    
    double pp_sngdbl_jpsi_inv_yield[narm] = {0}; 
    double pp_sngdbl_jpsi_inv_yield_err[narm] = {0}; 

    double pp_sngdbl_jpsi[narm] = {18692,13235};
    double pp_sngdbl_jpsi_err[narm] = {162,124};

    double pp_sngdbl_jpsi_acc_eff[narm] = {0};
    double pp_sngdbl_jpsi_trig_eff[narm] = {0};

    //RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       //pp_sngdbl_jpsi_acc_eff[0] = t1->GetBinContent(1);
       pp_sngdbl_jpsi_acc_eff[0] = 0.0337;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       //pp_sngdbl_jpsi_acc_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_acc_eff[1] = 0.0281;
       
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
       pp_sngdbl_jpsi_trig_eff[0] = 0.7467;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //ypp_sngdbl_jpsi_trig_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[1] = 0.7172;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }


    // calculate pp inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	pp_sngdbl_jpsi_inv_yield[arm] = pp_sngdbl_jpsi[arm]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm] * pp_MB[arm] * eff_jpsi ) );
	pp_sngdbl_jpsi_inv_yield_err[arm] = pp_sngdbl_jpsi_err[arm]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi);
	cout << "pp sngdbl_jpsi fvtx inv yield: " << pp_sngdbl_jpsi_inv_yield[arm] << "+\-" << pp_sngdbl_jpsi_inv_yield_err[arm] << endl;
      }


    double RpAu_sngdbl_jpsi_cent[narm][nbins] = {0};
    double RpAu_sngdbl_jpsi_cent_err[narm][nbins] = {0};
    
    // jpsi sngdbl
    for(int arm = 0; arm < 2; arm++)
      {
	for(int bin = 0; bin < 4; bin++)
	  {
	    abs1 = 0;
	    abs2 = 0;

	    RpAu_sngdbl_jpsi_cent[arm][bin] =( sngdbl_jpsi_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_sngdbl_jpsi_inv_yield[arm] ) ;
	    abs1 = sngdbl_jpsi_pAu_inv_yield_err[arm][bin] / sngdbl_jpsi_pAu_inv_yield[arm][bin];
	    abs2 = pp_sngdbl_jpsi_inv_yield_err[arm] / pp_sngdbl_jpsi_inv_yield[arm];
	    RpAu_sngdbl_jpsi_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_jpsi_cent[arm][bin];
	  }
      }

// 4a for sng+dbl psi2s
    //====================================================================================

 // // mixed events
      double sngdbl_psi2s[narm][nbins] = {
	{62,151,0,0},
	{92,79, 147,0}};  // three points North
       
    double sngdbl_psi2s_err[narm][nbins] = {
      {29,28,0,0},
      {20,17,21,0}};
    /////////////////////////////////////////////////////////////
  
    
    double sngdbl_psi2s_acc_eff[2][4] = {0};
    double sngdbl_psi2s_trig_eff[2][4] = {0};
    
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       // for(int cent =0; cent < 4; cent++)
       // sngdbl_psi2s_acc_eff[0][cent] = e1->GetEfficiency(cent+1);

       sngdbl_psi2s_acc_eff[0][0] = 0.0313;
       sngdbl_psi2s_acc_eff[0][1] = 0.0341;
       sngdbl_psi2s_acc_eff[0][2] = 0;
       sngdbl_psi2s_acc_eff[0][3] = 0;

	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       // for(int cent =0; cent < 4; cent++)
       // sngdbl_psi2s_acc_eff[1][cent] = e2->GetEfficiency(cent+1);

       sngdbl_psi2s_acc_eff[1][0] = 0.0434;
       sngdbl_psi2s_acc_eff[1][1] = 0.0434;
       sngdbl_psi2s_acc_eff[1][2] = 0.0460;
       sngdbl_psi2s_acc_eff[1][3] = 0;
	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_S.root");
       file_data1->GetObject("h2",f1);
       
       sngdbl_psi2s_trig_eff[0][0] = 0.7840;  // bkgd hits, prob00
       sngdbl_psi2s_trig_eff[0][1] = 0.6477;  // bkgd hits, prob00
	 sngdbl_psi2s_trig_eff[0][2] = 0;
       sngdbl_psi2s_trig_eff[0][3] = 0;
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("h2",f2);
       
       sngdbl_psi2s_trig_eff[1][0] = 0.7642;  // bkgd hits, prob00
       sngdbl_psi2s_trig_eff[1][1] = 0.7735;  // bkgd hits, prob00
       sngdbl_psi2s_trig_eff[1][2] = 0.7718;  // bkgd hits, prob00
       sngdbl_psi2s_trig_eff[1][3] = 0;
	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }

    double sngdbl_psi2s_pAu_inv_yield[narm][nbins] = {0};
    double sngdbl_psi2s_pAu_inv_yield_err[narm][nbins] = {0};

    // calculcate pAu inv. yield for psi2s fvtx
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {

	    FACTOR = 0;  // no prob cut in fvtx sim and prob cut used in sng data - so no adjustment	    

	    sngdbl_psi2s_pAu_inv_yield[arm][bin] = bias[bin] * sngdbl_psi2s[arm][bin] / ( (sngdbl_psi2s_trig_eff[arm][bin] ) * ( sngdbl_psi2s_acc_eff[arm][bin] ) * pAu_MB[arm][bin] );
 	    sngdbl_psi2s_pAu_inv_yield_err[arm][bin] = bias[bin] * sngdbl_psi2s_err[arm][bin] / ( sngdbl_psi2s_trig_eff[arm][bin]*sngdbl_psi2s_acc_eff[arm][bin] * pAu_MB[arm][bin] );
 	    // cout << "sngdbl_psi2s_pAu inv yield: " << sngdbl_psi2s_pAu_inv_yield[arm][bin] << "+\-" << sngdbl_psi2s_pAu_inv_yield_err[arm][bin] << endl;
 	  }
      }

  
    double pp_sngdbl_psi2s_inv_yield[narm] = {0}; 
    double pp_sngdbl_psi2s_inv_yield_err[narm] = {0}; 
    
   
     // mixed events
     double pp_sngdbl_psi2s[narm] = {586,464};
    double pp_sngdbl_psi2s_err[narm] = {38,36};

    // like-sign
    //   double pp_sngdbl_psi2s[narm] = {593,471};
    //   double pp_sngdbl_psi2s_err[narm] = {38,37};
    
   
    double pp_sngdbl_psi2s_acc_eff[narm] = {0}; 
    double pp_sngdbl_psi2s_trig_eff[narm] = {0}; 
    
    //RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       // pp_sngdbl_psi2s_acc_eff[0] = t1->GetBinContent(1);
    
       pp_sngdbl_psi2s_acc_eff[0] = 0.0423;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       // pp_sngdbl_psi2s_acc_eff[1] = t2->GetBinContent(1);
     
       pp_sngdbl_psi2s_acc_eff[1] = 0.0372;
       
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
     
       pp_sngdbl_psi2s_trig_eff[0] = 0.7637;  // bkgd hits, prob00
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //pp_sngdbl_psi2s_trig_eff[1] = t2->GetBinContent(1);
     
       pp_sngdbl_psi2s_trig_eff[1] = 0.7350;  // bkgd hits, prob00
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }



    // calculate pp inv. yield for psi2s fvtx
    for(int arm = 0; arm < 2; arm++)
      {

	FACTOR = 0;

 	pp_sngdbl_psi2s_inv_yield[arm] = pp_sngdbl_psi2s[arm]*eff_bbc / ( ( pp_sngdbl_psi2s_trig_eff[arm] ) * ( pp_sngdbl_psi2s_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
 	pp_sngdbl_psi2s_inv_yield_err[arm] = pp_sngdbl_psi2s_err[arm]*eff_bbc / ( ( pp_sngdbl_psi2s_trig_eff[arm] ) * ( pp_sngdbl_psi2s_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi);
 	// cout << "pp sngdbl_psi2s fvtx inv yield: " << pp_sngdbl_psi2s_inv_yield[arm] << "+\-" << pp_sngdbl_psi2s_inv_yield_err[arm] << endl;
      }

     double RpAu_sngdbl_psi2s_cent[narm][nbins] = {0};
     double RpAu_sngdbl_psi2s_cent_err[narm][nbins] = {0};
    
    // psi2s sngdbl
    for(int arm = 0; arm < 2; arm++)
      {
 	for(int bin = 0; bin < 4; bin++)
 	  {
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_sngdbl_psi2s_cent[arm][bin] =( sngdbl_psi2s_pAu_inv_yield[arm][bin] ) / ( ncoll[bin] * pp_sngdbl_psi2s_inv_yield[arm] ) ;
 	    abs1 = sngdbl_psi2s_pAu_inv_yield_err[arm][bin] / sngdbl_psi2s_pAu_inv_yield[arm][bin];
 	    abs2 = pp_sngdbl_psi2s_inv_yield_err[arm] / pp_sngdbl_psi2s_inv_yield[arm];
 	    RpAu_sngdbl_psi2s_cent_err[arm][bin] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_psi2s_cent[arm][bin];
	    cout << "RpAu psi2s[" << arm << "][" << bin << "]: " <<  RpAu_sngdbl_psi2s_cent[arm][bin] << ", ncoll: " << ncoll[bin] << ", raw psi2s yield: " <<  sngdbl_psi2s[arm][bin] << ", psi2s pAu inv yield: " <<  sngdbl_psi2s_pAu_inv_yield[arm][bin] <<  endl;
 	  }
      }



  ///////////////////////////////////////////////////////////
  

    double ratio_pAu[2][4] = {0};
    double ratio_pAu_err[2][4] = {0};

    // using PPG228 J/psi
    for( int arm = 0; arm < 2; arm++)
      {
    	for(int cent = 0; cent < 4; cent++)  // last bin is zero
    	  {
    	    double temp0 = 0;
    	    double tempN = 0;
    	    double tempD = 0;
    	    double tempRN = 0;
    	    double tempRD = 0;

    	    ratio_pAu[arm][cent] = sngdbl_psi2s_pAu_inv_yield[arm][cent] /  jpsi_ppg228_pAu_inv_yield_3bins[arm][cent];
	  
    	    tempN = sngdbl_psi2s_pAu_inv_yield_err[arm][cent] /  sngdbl_psi2s_pAu_inv_yield[arm][cent] ;
    	    tempD = jpsi_ppg228_pAu_inv_yield_3bins_err[arm][cent] /  jpsi_ppg228_pAu_inv_yield_3bins[arm][cent] ;

    	    temp0 = tempN * tempN + tempD * tempD;
    	    ratio_pAu_err[arm][cent] = sqrt( temp0 ) * ratio_pAu[arm][cent];

	    cout << "RpAu 228: " << ratio_pAu[arm][cent] << "+/-" << ratio_pAu[arm][cent]*ratio_pAu_err[arm][cent] << endl;
    	  }
      }



 
    // ALICE data from hepdata pub JHEP 02 (20121) 002 - inspire id 100166
      // ALICE data from hepdata pub JHEP 02 (20121) 002 - inspire id 100166
    ///////////////////////////////////////////////////////////////////////////////
    double ratio_alice[2][6] = {
      {0.0176, 0.0103, 0.0140, 0.0110, 0.0105, 0.0162},
      // {0.0251, 0.0182, 0.0173, 0.0130, 0.0148, 0.0246}};
      {0.0246, 0.0148, 0.0130, 0.0173, 0.0182, 0.0251}};

    double ratio_alice_err[2][6] = {
      {0.0034, 0.0031, 0.0023, 0.0027, 0.0031, 0.0056},
      // {0.0064, 0.0038, 0.0033, 0.0029, 0.0042, 0.0046}};
      {0.0046, 0.0042, 0.0029, 0.0033, 0.0038, 0.0064}};

    double ratio_alice_sys[2][6] = {
      {0.0014, 0.00086, 0.0014, 0.0010, 0.0012, 0.0016},
      // {0.0026, 0.0016, 0.0013, 0.0016, 0.0021, 0.0023}};
      {0.0023, 0.0021, 0.0016, 0.0013, 0.0016, 0.0026}};
    ///////////////////////////////////////////////////////////////////////////////


    const Double_t err_x[4] = {0,0,0, 0};
    const Double_t err_x6[6] = {0,0,0,0,0,0};
    const Double_t width = 0.25;
    const Double_t jpsi_ppg228_sys[narm][nbins] = {
      {0.04,0.04,0.05,0.04}, //  ->
      {0.05,0.07,0.07,0.08}};

   
    const Double_t ratio_percent[narm][nbins] = {
      {0.12, 0.12, 0, 0}, //  ->
      {0.09, 0.09, 0.09, 0}};

    double ratio_sys[2][4] = {0};

    for(int arm = 0; arm < 2; arm++)
     {
       for(int i = 0; i < 4; i++)
	 {
	   ratio_sys[arm][i] = ratio_percent[arm][i]*ratio_pAu[arm][i];
	 }
     }
   

    TGraphErrors *gr1 = new TGraphErrors(4,ncoll_array1,RpAu_jpsi_sng_cent[1],err_x,RpAu_jpsi_sng_cent_err[1]);

    TGraphErrors *gr3 = new TGraphErrors(6,ncoll_array6,ratio_alice[1],err_x6, ratio_alice_err[1]);

    TGraphErrors *gr5 = new TGraphErrors(3,ncoll_array4,ratio_pAu[1],err_x, ratio_pAu_err[1] );


    TGraphErrors *gr7 = new TGraphErrors(4,ncoll_array4,RpAu_jpsi_sng_cent[1],err_x,RpAu_jpsi_sng_cent_err[1]);
    TGraphErrors *gr9 = new TGraphErrors(4,ncoll_array4,RpAu_sngdbl_psi2s_cent[1],err_x, RpAu_sngdbl_psi2s_cent_err[1]);
    TGraphAsymmErrors *gr1_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr7_sys = new TGraphAsymmErrors(4);
  
   
    TGraphErrors *gr2 = new TGraphErrors(6,ncoll_array6,ratio_alice[0],err_x6, ratio_alice_err[0]);

    TGraphErrors *gr4 = new TGraphErrors(6,ncoll_array6,ratio_alice[0],err_x6, ratio_alice_err[0]);

    TGraphErrors *gr6 = new TGraphErrors(3,ncoll_array4, ratio_pAu[0],err_x, ratio_pAu_err[0]);


    TGraphErrors *gr8 = new TGraphErrors(4,ncoll_array4,RpAu_jpsi_sng_cent[0],err_x,RpAu_jpsi_sng_cent_err[0]);
    TGraphErrors *gr10 = new TGraphErrors(4,ncoll_array4,RpAu_sngdbl_psi2s_cent[0],err_x, RpAu_sngdbl_psi2s_cent_err[0]);
    TGraphAsymmErrors *gr2_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr4_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr6_sys = new TGraphAsymmErrors(4);
    TGraphAsymmErrors *gr8_sys = new TGraphAsymmErrors(4);

  

 for(int i = 0; i < 6; i++)
     {

       gr3_sys->SetPoint(i, ncoll_array6[i],ratio_alice[1][i]);  
       gr3_sys->SetPointError(i, width, width, ratio_alice_sys[1][i], ratio_alice_sys[1][i]);
  
       gr4_sys->SetPoint(i, ncoll_array6[i],ratio_alice[0][i]);  
       gr4_sys->SetPointError(i, width, width, ratio_alice_sys[0][i], ratio_alice_sys[0][i]);

     }

   for(int i = 0; i < 4; i++)
     {
       gr5_sys->SetPoint(i, ncoll_array4[i],ratio_pAu[1][i]);
       gr5_sys->SetPointError(i, width, width, ratio_sys[1][i], ratio_sys[1][i]);

       gr6_sys->SetPoint(i, ncoll_array4[i],ratio_pAu[0][i]);
       gr6_sys->SetPointError(i, width, width,ratio_sys[0][i], ratio_sys[0][i]);
	  
     }

   gr1_sys->SetMarkerColor(kRed+1);                     
    gr1_sys->SetLineColor(kRed+1);                         
    gr1_sys->SetFillStyle(0); 
   
    gr2_sys->SetMarkerColor(kRed+1);                     
    gr2_sys->SetLineColor(kRed+1);                         
    gr2_sys->SetFillStyle(0); 

    /////////////////////////////////////////////////////

    gr3_sys->SetMarkerColorAlpha(kBlack, 0.8);                     
    gr3_sys->SetLineColorAlpha(kBlack, 0.8);                         
    gr3_sys->SetFillStyle(0); 
    gr3_sys->SetLineWidth(1);

    gr4_sys->SetMarkerColorAlpha(kBlack, 0.8);                     
    gr4_sys->SetLineColorAlpha(kBlack, 0.8);                         
    gr4_sys->SetFillStyle(0); 
    gr4_sys->SetLineWidth(1);

    gr5_sys->SetMarkerColorAlpha(kRed+1, 0.8);                     
    gr5_sys->SetLineColorAlpha(kRed+1, 0.8);                         
    gr5_sys->SetFillStyle(0); 
    gr5_sys->SetLineWidth(1);

    gr6_sys->SetMarkerColorAlpha(kRed+1, 0.8);                     
    gr6_sys->SetLineColorAlpha(kRed+1, 0.8);                         
    gr6_sys->SetFillStyle(0); 
    gr6_sys->SetLineWidth(1);

    /////////////////////////////////////////////////////

    gr7_sys->SetMarkerColor(kBlack);                     
    gr7_sys->SetLineColor(kBlack);                         
    gr7_sys->SetFillStyle(0); 

    gr8_sys->SetMarkerColor(kBlack);                     
    gr8_sys->SetLineColor(kBlack);                         
    gr8_sys->SetFillStyle(0); 
   

    
    
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
  
  
    gr1->GetXaxis()->SetTitleOffset(1.05);
    gr1->GetXaxis()->SetTitleFont(102);
    gr1->GetXaxis()->SetLabelFont(102);
    gr1->GetXaxis()->SetTitle("#LTN_{coll}#GT");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.05);

    gr1->GetXaxis()->SetRangeUser(0.9,13.9);
    gr1->GetXaxis()->SetNdivisions(5);

    gr1->GetYaxis()->SetTitleOffset(0.95);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.065);
    gr1->GetYaxis()->SetLabelSize(0.05);


    gr1->GetYaxis()->SetTitle("#psi(2S)/#psi(1S) Ratio");
 
    gr1->GetYaxis()->SetRangeUser(0,0.055);
    gr1->GetYaxis()->SetNdivisions(7);
 
    gr2->SetMarkerColor(kBlue+1);
    gr2->SetLineWidth(3);
    gr2->SetLineColor(kBlue+1);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(2);                                 

    //////////////////////////////////////////////////
    gr3->SetMarkerColorAlpha(kBlack, 0.8);
    gr3->SetLineWidth(3);
    gr3->SetLineColorAlpha(kBlack, 0.8);
    gr3->SetMarkerStyle(20); 
    gr3->SetMarkerSize(2);                       

    gr4->SetMarkerColorAlpha(kBlack, 0.8);
    gr4->SetLineWidth(3);
    gr4->SetLineColorAlpha(kBlack, 0.8);
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerSize(2);                       
    //////////////////////////////////////////////////

    gr5->SetMarkerColorAlpha(kRed+1, 0.8);
    gr5->SetLineColorAlpha(kRed+1, 0.8);
    gr5->SetLineWidth(3);                           
    gr5->SetMarkerStyle(21);                       
    gr5->SetMarkerSize(2); 

    gr6->SetMarkerColorAlpha(kRed+1, 0.8);                     
    gr6->SetLineColorAlpha(kRed+1, 0.8);
    gr6->SetLineWidth(3);                           
    gr6->SetMarkerStyle(21);                       
    gr6->SetMarkerSize(2);                       
    //////////////////////////////////////////////////

      
 
     if(north_arm)
      {
	gr1->SetMarkerColor(kWhite);
	gr1->SetLineColor(kWhite);
	gr1->Draw("AP");
	gr3->Draw("P");
	gr3_sys->Draw("e2same");
       	gr5->Draw("P");
	gr5_sys->Draw("e2same");

      }
     else
       {
	 gr1->SetMarkerColor(kWhite);
	 gr1->SetLineColor(kWhite);
	 gr1->Draw("AP");
	 // gr2->Draw("P");
	 gr4->Draw("P");
	 gr4_sys->Draw("e2same");
	 gr6->Draw("P");
	 gr6_sys->Draw("e2same");
       }

   
     TLatex l2;
    l2.SetTextSize(0.0425);
    l2.SetTextAlign(12);
    l2.SetTextColor(kBlack);
  
    char text3[100];
    char text4[100];
    char text5[100];
    char text6[100];
    char text7[100];

    if(north_arm)
      {
	sprintf(text3,"1.2 < y < 2.2");
	sprintf(text4,"2.03 < y < 3.53");
	sprintf(text5,"Inclusive");
	sprintf(text6, "(JHEP02 (2021) 002)");
	sprintf(text7, "(PRC 102, 014902)");
      }
    else
      {
	sprintf(text3,"-2.2 < y < -1.2");
	sprintf(text4,"-4.46 < y < -2.96");
	sprintf(text5,"Inclusive");
	sprintf(text6, "(JHEP02 (2021) 002)");
	sprintf(text7, "(PRC 102, 014902)");
      }
     
    l2.SetTextFont(102);  
    l2.DrawLatexNDC(0.245, 0.64, text3); // x0, y0
    l2.DrawLatexNDC(0.245, 0.905, text5); // x0, y0
    l2.DrawLatexNDC(0.245, 0.796, text4); // x0, y0
    l2.SetTextSize(0.04);
    l2.DrawLatexNDC(0.24, 0.75, text6); // x0, y0
    // l2.DrawLatexNDC(0.24, 0.595, text7); // x0, y0


    TLine *line = new TLine(0.8,0.0247,14, 0.0247);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");



    double global_sys = 0.101*0.0247; 

    TBox *box = new TBox(13.75, 0.0247-global_sys, 13.95, 0.0247+global_sys);
 
    box->SetFillStyle(1000);
    box->SetFillColorAlpha(kRed+1,0.5);
    //  box->Draw();

    double global_sys0 = 0.084*0.0223; 

    TBox *box0 = new TBox(13.55, 0.0247-global_sys0, 13.75, 0.0247+global_sys0);
    box0->SetFillStyle(1000);
    box0->SetFillColorAlpha(13,0.5);
    // box0->Draw();


    double global_sys2 = 0.020*.0223;

    TBox *box2 = new TBox(13.35, 0.0247-global_sys2, 13.55, 0.0247+global_sys2);
    box2->SetFillStyle(1000);
    box2->SetFillColorAlpha(1,0.5);
    // box2->Draw();


    TLegend *leg = new TLegend(0.18, 0.835, 0.44, 0.86);  //(start x, start y, end x, end y)

    leg->SetLineColor(0);
    leg->SetTextSize(0.0425); 
    leg->SetTextFont(102);

    if(north_arm)
      {

	leg->AddEntry(gr3,"ALICE p+Pb #sqrt{s_{NN}}=8.16 TeV", "p");
      }
    else
      {

	leg->AddEntry(gr4,"ALICE p+Pb #sqrt{s_{NN}}=8.16 TeV", "p");

      }
    leg->Draw();
    
    TLegend *leg1 = new TLegend(0.18, 0.675, 0.44, 0.705);  //(start x, start y, end x, end y)

    leg1->SetLineColor(0);
    leg1->SetTextSize(0.0425); 
    leg1->SetTextFont(102);

    if(north_arm)
      {
	
  	leg1->AddEntry(gr5,"PHENIX p+Au #sqrt{s_{NN}}=200 GeV", "p");
      }
    else
      {
  	leg1->AddEntry(gr6,"PHENIX p+Au #sqrt{s_{NN}}=200 GeV", "p");

      }
    leg1->Draw();



    if(north_arm)
      c1->SaveAs(Form("../psi2S_pdf/RpAu_ratio_alice_plot_fwd.pdf"));
    else
      c1->SaveAs(Form("../psi2S_pdf/RpAu_ratio_alice_plot_bkwd.pdf"));




}











