
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


void RpAu_psi2s_transport_plot_2pts()
{
    gROOT->ForceStyle();

  bool north_arm = false;
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

  double cent_binwidth[4] = {0.20,0.64, 0, 0};

  const Double_t ncoll_array1[nbins] = {8.1,6.0,4.3,2.7};  // need all 4 points for ppg228 array
  const Double_t ncoll_array2[nbins] = {8.3,6.2,4.6,2.7};
  const Double_t ncoll_array3[nbins] = {8.1,6.0,4.3,2.5};  // need all 4 points for ppg228 array
  const Double_t ncoll_array4[nbins] = {8.2,4.2, 9,0};

  double bias_228[nbins] = {0.90,0.98,1.02,1.0};
  double ncoll_228[nbins] = {8.2,6.1,4.4,2.6};

  double bias[nbins] = {0.90,1.0,0,0};
  double ncoll[nbins] = {8.2,4.3,0,0};

  double pAu_MB_228[narm][nbins] = {
    {5.047069934*pow(10,10),4.99512054720000000e+10,4.97313174400000000e+10,5.96167043840000000e+10},
   {4.6346380480*pow(10,10),4.59845716160000000e+10,4.58064919040000000e+10,5.49289634240000000e+10} };

double pAu_MB[narm][nbins] = {
  {5.047069934*pow(10,10),1.59299227296e+11,0,0},
  {4.6346380480*pow(10,10),1.46720026944e+11,0,0} };

 
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
	    FACTOR = 0;  // sngdbl pAu has prob cut on both arms and both tracks - April 25th
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
       pp_sngdbl_jpsi_trig_eff[0] = 0.7458;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //ypp_sngdbl_jpsi_trig_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[1] = 0.7162;
       
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

////////////////////////////// calcualte weighted MB sum
    double w_ave2[2] = {0};
    double num2[2] = {0};
    double denom2[2] = {0};
    double weight =(4.7)/(0.86) ;

    for(int arm = 0; arm < 1; arm++)
      {
 	for(int bin = 0; bin < 2; bin++)
 	  {
	     if(bin < 3)
	      {
		num2[arm] += RpAu_sngdbl_psi2s_cent[arm][bin] * ncoll[bin] * weight * cent_binwidth[bin] ;
		denom2[arm] += ncoll[bin] * weight * cent_binwidth[bin] ;
	      }
	    w_ave2[arm] = num2[arm]/denom2[arm];
	 
	  }
      }
	////////////////////////////// calcualte weighted MB sum
   


  //====================================================================================
    TGraphAsymmErrors *p1;  // J/psi fwd = { 8.1, 6.7, 5.1, 3.3, 1.6,  1.2}
    TGraphAsymmErrors *p2;  // J[psi bkwd
    TGraphAsymmErrors *p3;  // psi2s fwd
    TGraphAsymmErrors *p4;  // psi2s bkwd
    TGraphAsymmErrors *p5;  // CNM fwd
    TGraphAsymmErrors *p6;  // CNM bkwd

    double jpsi_transport[2][6] = {0};
    double psi2s_transport[2][6] = {0};
    double cnm_effects[2][6] = {0};

    // Ralph Rapp predciitons
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/psi2S_root_files/RpAu_transport_ncoll.root");
       file_data1->GetObject("RpAu_transport_bkwd",p2);
       file_data1->GetObject("RpAu_transport_fwd",p1);
       file_data1->GetObject("RpAu_2s_transport_bkwd",p4);
       file_data1->GetObject("RpAu_2s_transport_fwd",p3);
       file_data1->GetObject("RpAu_CNM_bkwd",p6);
       file_data1->GetObject("RpAu_CNM_fwd",p5);

	file_data1->Close();
	
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


    const Double_t err_x[nbins] = {0,0,0,0};
    const Double_t width = 0.145;
    const Double_t jpsi_ppg228_sys[narm][nbins] = {
      {0.04,0.04,0.05,0.04}, //  ->
      {0.05,0.07,0.07,0.08}};

     const Double_t psi2s_percent[narm][nbins] = {
      {0.16, 0.16, 0.0, 0}, //  ->
      {0.13, 0.13, 0.14, 0}};

    double psi2s_sys[2][4] = {0};

    for(int arm = 0; arm < 2; arm++)
     {
       for(int i = 0; i < 4; i++)
	 {
	   psi2s_sys[arm][i] = psi2s_percent[arm][i]*RpAu_sngdbl_psi2s_cent[arm][i];
	 }
     }
   

   TGraphErrors *gr1 = new TGraphErrors(4,ncoll_array1,RpAu_jpsi_fvtx_cent[1],err_x,RpAu_jpsi_fvtx_cent_err[1]);
   TGraphErrors *gr3 = new TGraphErrors(4,ncoll_array2,RpAu_sngdbl_jpsi_cent[1],err_x, RpAu_sngdbl_jpsi_cent_err[1]);
   TGraphErrors *gr5 = new TGraphErrors(4,ncoll_array3,RpAu_jpsi_ppg228_cent[1],err_x, RpAu_jpsi_ppg228_cent_err[1] );
   TGraphErrors *gr7 = new TGraphErrors(4,ncoll_array4,RpAu_jpsi_sng_cent[1],err_x,RpAu_jpsi_sng_cent_err[1]);
   TGraphErrors *gr9 = new TGraphErrors(4,ncoll_array4,RpAu_sngdbl_psi2s_cent[1],err_x, RpAu_sngdbl_psi2s_cent_err[1]);
   TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(4);
   TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(4);
  
   
   TGraphErrors *gr2 = new TGraphErrors(4,ncoll_array1,RpAu_jpsi_fvtx_cent[0],err_x,RpAu_jpsi_fvtx_cent_err[0]);
   TGraphErrors *gr4 = new TGraphErrors(4,ncoll_array2,RpAu_sngdbl_jpsi_cent[0],err_x, RpAu_sngdbl_jpsi_cent_err[0]);
   TGraphErrors *gr6 = new TGraphErrors(4,ncoll_array3,RpAu_jpsi_ppg228_cent[0],err_x, jpsi_ppg228_sys[0]);
   TGraphErrors *gr8 = new TGraphErrors(4,ncoll_array4,RpAu_jpsi_sng_cent[0],err_x,RpAu_jpsi_sng_cent_err[0]);
    TGraphErrors *gr10 = new TGraphErrors(4,ncoll_array4,RpAu_sngdbl_psi2s_cent[0],err_x, RpAu_sngdbl_psi2s_cent_err[0]);
   TGraphAsymmErrors *gr10_sys = new TGraphAsymmErrors(4);
   TGraphAsymmErrors *gr6_sys = new TGraphAsymmErrors(4);
    

   
   for(int i = 0; i < 4; i++)
     {
       gr3_sys->SetPoint(i, ncoll_array4[i],RpAu_sngdbl_psi2s_cent[1][i]);   // need to also include array index here as well as arm index
       gr3_sys->SetPointError(i, width, width, psi2s_sys[1][i], psi2s_sys[1][i]);

       gr5_sys->SetPoint(i, ncoll_array3[i],RpAu_jpsi_ppg228_cent[1][i]);   // need to also include array index here as well as arm index
       gr5_sys->SetPointError(i, width, width, jpsi_ppg228_sys[1][i], jpsi_ppg228_sys[1][i]);

       gr10_sys->SetPoint(i, ncoll_array4[i],RpAu_sngdbl_psi2s_cent[0][i]);
       gr10_sys->SetPointError(i, width, width, psi2s_sys[0][i], psi2s_sys[0][i]);

       gr6_sys->SetPoint(i, ncoll_array3[i],RpAu_jpsi_ppg228_cent[0][i]);
       gr6_sys->SetPointError(i, width, width, jpsi_ppg228_sys[0][i], jpsi_ppg228_sys[0][i]);
	  
     }

      
      
    gr3_sys->SetMarkerColor(kRed+1);                     
    gr3_sys->SetLineColor(kRed+1);                         
    gr3_sys->SetFillStyle(0); 
    gr3_sys->SetLineWidth(1);

    gr10_sys->SetMarkerColor(kRed+1);                     
    gr10_sys->SetLineColor(kRed+1);                         
    gr10_sys->SetFillStyle(0); 
    gr10_sys->SetLineWidth(1);

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
   
 
    gr1->SetMarkerColor(kRed+1);
    gr1->SetTitle("");
    gr1->SetLineWidth(3);
    gr1->SetLineColor(kRed+1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(2);
  
  
    gr1->GetXaxis()->SetTitleOffset(1.05);
    gr1->GetXaxis()->SetTitleFont(102);
    gr1->GetXaxis()->SetLabelFont(102);
    gr1->GetXaxis()->SetTitle("#LTN_{coll}#GT");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.05);
    
    gr1->GetXaxis()->SetRangeUser(0.1,10);
     gr1->GetXaxis()->SetNdivisions(5);
    // gr1->GetYaxis()->SetNdivisions(8);
   

    gr1->GetYaxis()->SetTitleOffset(0.8);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.075);
    gr1->GetYaxis()->SetLabelSize(0.05);
 
    gr1->GetYaxis()->SetTitle("R_{pAu}");
 
    gr1->GetYaxis()->SetRangeUser(0,2.5);
    gr1->GetYaxis()->SetNdivisions(7);

    if(north_arm)
      gr1->GetXaxis()->SetRangeUser(0, 100);
    else
      gr1->GetXaxis()->SetRangeUser(0, 100);


    gr2->SetMarkerColor(kRed+1);
    gr2->SetLineWidth(3);
    gr2->SetLineColor(kRed+1);
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
    gr5->SetMarkerStyle(20);                       
    gr5->SetMarkerSize(2); 

    gr6->SetMarkerColorAlpha(kBlack, 0.8);                     
    gr6->SetLineColorAlpha(kBlack, 0.8);                         
    gr6->SetLineWidth(3);                           
    gr6->SetMarkerStyle(20);                       
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

    gr9->SetMarkerColor(kRed+1);
    gr9->SetLineWidth(3);
    gr9->SetLineColor(kRed+1);
    gr9->SetMarkerStyle(21);
    gr9->SetMarkerSize(2); 
  
    gr10->SetMarkerColorAlpha(kRed+1, 0.8);
    gr10->SetLineWidth(3);
    gr10->SetLineColorAlpha(kRed+1, 0.8);
    gr10->SetMarkerStyle(21);
    gr10->SetMarkerSize(2);   

    /////////////////// predicitons
 
    // predicitons - psi(2S) South
    p4->SetFillStyle(3001);
    p4->SetLineColor(kRed+1);
    p4->SetFillColorAlpha(kRed+1,0.3);
    p4->SetLineWidth(1);
   
    // predicitons - J/psi South
    p2->SetFillStyle(3224);
    p2->SetLineColor(kBlack);
    p2->SetFillColorAlpha(kBlack, 0.5);
    p2->SetLineWidth(1);
   
 
     if(north_arm)
      {
	gr1->SetMarkerColor(kWhite);
	gr1->SetLineColor(kWhite);
	gr1->Draw("AP");
	//	gr3->Draw("P");
	gr5->Draw("P");
	//	gr7->Draw("P");
		gr9->Draw("P");
		gr3_sys->Draw("e2same");
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
	 // gr4->Draw("P");
	 gr6->Draw("P");
	 // gr8->Draw("P");
	 gr10->Draw("P");
	 gr10_sys->Draw("e2same");
	 gr6_sys->Draw("e2same");
      	 p2->Draw("e3same");
	 p4->Draw("e3same");
	 //	 p6->Draw("e3same");
       }

   TLatex l2;
    l2.SetTextSize(0.0425);
    l2.SetTextAlign(12);
    //   l2.SetTextColor(kBlue+4);
  
    char text3[100];
    char text4[100];
    char text5[100];

    if(north_arm)
      sprintf(text3,"1.2 < y < 2.2, Inclusive");
    else
      sprintf(text3,"-2.2 < y < -1.2, Inclusive");

    sprintf(text4,"Inclusive");
    sprintf(text5, "(Phys. Rev. C 102, 014902)");
    l2.SetTextFont(102);  
    l2.DrawLatexNDC(0.24, 0.91, text3); // x0, y0
    //  l2.DrawLatexNDC(0.76, 0.91, text4); // x0, y0
    l2.SetTextSize(0.04);
    l2.DrawLatexNDC(0.24, 0.81, text5); // x0, y0



    double global_sys = 0.101;  // RpAu J/psi PPG228 in both Norht and South

    TBox *box = new TBox(8.5, 1-global_sys, 8.7, 1+global_sys);
    box->SetFillStyle(1000);
    box->SetFillColorAlpha(1,0.5);
    box->Draw();


    TLine *line = new TLine(2.2,1,8.6,1);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");

     TLegend *leg = new TLegend(0.18, 0.84, 0.44, 0.89);  //(start x, start y, end x, end y)
    leg->SetFillColor(19); 
    leg->SetFillStyle(3003); 
    leg->SetLineWidth(1);
    leg->SetLineColor(0);
    leg->SetTextSize(0.0425); 
    leg->SetTextFont(102);

    if(north_arm)
      {

	leg->AddEntry(gr5,"R_{pAu} J/#psi PPG228 No Fvtx", "p");

      }
    else
      {
	leg->AddEntry(gr6,"R_{pAu} J/#psi, p+Au #sqrt{s_{NN}}=200 GeV", "p");
      }

    leg->Draw();



    TLegend *leg2 = new TLegend(0.18, 0.61, 0.44, 0.79);  //(start x, start y, end x, end y)
    leg2->SetFillColor(19); 
    leg2->SetFillStyle(3003); 
    leg2->SetLineWidth(1);
    leg2->SetLineColor(0);
    leg2->SetTextSize(0.0425); 
    leg2->SetTextFont(102);

    if(north_arm)
      {
	//	leg2->AddEntry(gr1,"R_{pAu} J/#psi Fvtx Tracks", "p");
	//	leg2->AddEntry(gr3,"R_{pAu} J/#psi Sng+Dbl", "p");
	//	leg2->AddEntry(gr5,"R_{pAu} J/#psi PPG228 No Fvtx", "p");
	//	leg2->AddEntry(gr7,"R_{pAu} J/#psi Sng Tracks", "p");
       	//leg2->AddEntry(gr9,"R_{pAu} #psi(2S) Sng+Dbl Tracks", "p");
	leg2->AddEntry(p1,"J/#psi Transport Model (Du&Rapp)", "f");
	leg2->AddEntry(p5,"CNM Effects (Du&Rapp)", "f");
	leg2->AddEntry(p3,"#psi(2S) Transport Model (Du&Rapp)", "f");

      }
    else
      {
	//	leg2->AddEntry(gr2,"R_{pAu} J/#psi Fvtx Tracks", "p");
	//	leg2->AddEntry(gr4,"R_{pAu} J/#psi Sng+Dbl", "p");
	//	leg2->AddEntry(gr6,"R_{pAu} J/#psi, p+Au #sqrt{s_{NN}}=200 GeV", "p");
	//	leg2->AddEntry(gr8,"R_{pAu} J/#psi Sng Tracks", "p");
	leg2->AddEntry(gr10,"R_{pAu} #psi(2S), p+Au #sqrt{s_{NN}}=200 GeV", "p");
	leg2->AddEntry(p2,"J/#psi Transport Model (Du&Rapp)", "f");
	//leg2->AddEntry(p6,"CNM Effects (Du&Rapp)", "f");
	leg2->AddEntry(p4,"#psi(2S) Transport Model (Du&Rapp)", "f");

      }
    leg2->Draw();


    cout << "North psi(2S) sngdbl RpAu MB weighted average: " << w_ave2[1] << endl;
    cout << "South psi(2S) sngdbl RpAu MB weighted average: " << w_ave2[0] << endl;

 if(north_arm)
      c1->SaveAs(Form("../psi2S_pdf/RpAu_psi2s_rapp_plot_fwd.pdf"));
    else
      c1->SaveAs(Form("../psi2S_pdf/RpAu_psi2s_rapp_plot_bkwd.pdf"));



  TFile *h;

  h = new TFile("RpAu_psi2s_transport_ncoll_S.root", "RECREATE");

  h->cd();

  gr6->SetName("jpsi_data");
  gr6->SetMarkerColor(kBlack);
  gr6->SetLineColor(kBlack);
  gr6->Write();

  gr6_sys->SetName("jpsi_sys");
  gr6_sys->SetMarkerColor(kBlack);
  gr6_sys->SetLineColor(kBlack);
  gr6_sys->Write();

  gr10->SetName("psi2s_data");
  gr10->SetMarkerColor(kBlack);
  gr10->SetLineColor(kBlack);
  gr10->Write();

  gr10_sys->SetName("psi2s_sys");
  gr10_sys->SetMarkerColor(kBlack);
  gr10_sys->SetLineColor(kBlack);
  gr10_sys->Write();
  
  p2->Write();
  p4->Write();

}








