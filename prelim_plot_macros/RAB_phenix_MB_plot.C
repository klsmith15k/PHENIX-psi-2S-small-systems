
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


using namespace std;


void RAB_phenix_MB_plot()
{
    gROOT->ForceStyle();

  bool north_arm = false;
  int arm = 0;
  double FACTOR = 0.0005;

  static const int nbins = 2;  // north and south
  static const int narm = 2;

  if(north_arm)
    arm = 1;
    
  double eff_bbc = 0.55;
  double eff_jpsi = 0.79;

  double abs1 = 0;
  double abs2 = 0;
 
  
  const Double_t rap_array1[3] = {-3, 0, 3};
  const Double_t rap_array2[nbins] = {-1.5, 1.5};
  const Double_t rap_array3[nbins] = {-1.7, 1.7};
  const Double_t rap_array4[nbins] = {-1.9, 1.9};

  double delta_rap[nbins] = {1, 1};
  
  double bias = 0.86;
  double ncoll = 4.7;

  double pAu_MB[narm] = {2.09769926640*pow(10,11), 1.93066407424*pow(10,11)};
 
  double pp_MB_228[narm] = {1.092422549504*pow(10,12) ,1.073082815488*pow(10,12)};    // PPG228 good runlist

  double pp_MB[narm] = {9.16550312960000000*pow(10,11), 8.16541736960000000*pow(10,11)};   // PPG201 good runlist

 TFile *file_data1, *file_data2;

  TH1D *t1, *t2;

  TGraphAsymmErrors *f1;
  TGraphAsymmErrors *f2;

  TEfficiency *e1;
  TEfficiency *e2;

  //  for sng+dbl Jpsi
    //====================================================================================
//====================================================================================


  double sngdbl_jpsi[narm]= {11227, 11701};
  double sngdbl_jpsi_err[narm] = {144, 130};
    
    // correct
    double sngdbl_jpsi_trig_eff[narm] = {0};
   double sngdbl_jpsi_acc_eff[narm]  = {0};
      
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_jpsi_acc_eff[0] = e1->GetEfficiency(pt+1);

       sngdbl_jpsi_acc_eff[0] = 0.0269;
           
	 	              
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_jpsi_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_jpsi_acc_eff[1] = e2->GetEfficiency(pt+1);

       sngdbl_jpsi_acc_eff[1] = 0.0345;
     	 
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_S.root");
       file_data1->GetObject("h1",f1);
       
       sngdbl_jpsi_trig_eff[0] = 0.6617;
     	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       sngdbl_jpsi_trig_eff[1] = 0.7503;
    	            
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }


    double sngdbl_jpsi_pAu_inv_yield[narm] = {0};
    double sngdbl_jpsi_pAu_inv_yield_err[narm] = {0};

    // calculcate pAu inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	FACTOR = 0;
	sngdbl_jpsi_pAu_inv_yield[arm] = bias * sngdbl_jpsi[arm] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm] ) * ( FACTOR + sngdbl_jpsi_acc_eff[arm] ) * pAu_MB[arm]  );
	sngdbl_jpsi_pAu_inv_yield_err[arm] = bias * sngdbl_jpsi_err[arm] / ( ( FACTOR + sngdbl_jpsi_trig_eff[arm]) * ( FACTOR + sngdbl_jpsi_acc_eff[arm]) * pAu_MB[arm]  );
	cout << "sngdbl_jpsi_pAu inv yield: " << sngdbl_jpsi_pAu_inv_yield[arm] << "+\-" << sngdbl_jpsi_pAu_inv_yield_err[arm] << endl;
	
      }
    
    
    double pp_sngdbl_jpsi_inv_yield[narm] = {0}; 
    double pp_sngdbl_jpsi_inv_yield_err[narm] = {0}; 

    double pp_sngdbl_jpsi[narm] = {18712, 13251};
    double pp_sngdbl_jpsi_err[narm] = {350, 124};
    
    double pp_sngdbl_jpsi_acc_eff[narm] = {0};
    double pp_sngdbl_jpsi_trig_eff[narm] = {0};

     //RUN15pp
     {
          
       /////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_S.root");
       file_data1->GetObject("geff_y_int_S",t1);
       //pp_sngdbl_jpsi_acc_eff[0] = t1->GetBinContent(1);
       pp_sngdbl_jpsi_acc_eff[0] = 0.0343;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_acceff_rap_N.root");
       file_data2->GetObject("geff_y_int_N",t2);
       //pp_sngdbl_jpsi_acc_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_acc_eff[1] = 0.0286;
       
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
       pp_sngdbl_jpsi_trig_eff[0] = 0.7471;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //ypp_sngdbl_jpsi_trig_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[1] = 0.7219;
       
       file_data1->Close();
       file_data2->Close();
       /////////////////////////////////////////
     }



    // calculate pp inv. yield for jpsi fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	FACTOR = 0;  // no prob cut adjustment now
	    pp_sngdbl_jpsi_inv_yield[arm] = pp_sngdbl_jpsi[arm]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi  );
	    pp_sngdbl_jpsi_inv_yield_err[arm] = pp_sngdbl_jpsi_err[arm]*eff_bbc / ( ( FACTOR + pp_sngdbl_jpsi_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_jpsi_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
	    cout << "pp sngdbl_jpsi fvtx inv yield: " << pp_sngdbl_jpsi_inv_yield[arm] << "+\-" << pp_sngdbl_jpsi_inv_yield_err[arm] << endl;
	  
      }


    double RpAu_sngdbl_jpsi_cent[narm] = {0};
    double RpAu_sngdbl_jpsi_cent_err[narm] = {0};
    
    // jpsi sngdbl
    for(int arm = 0; arm < 2; arm++)
      {

	    abs1 = 0;
	    abs2 = 0;

	    RpAu_sngdbl_jpsi_cent[arm] =( sngdbl_jpsi_pAu_inv_yield[arm] ) / ( ncoll * pp_sngdbl_jpsi_inv_yield[arm] ) ;
	    abs1 = sngdbl_jpsi_pAu_inv_yield_err[arm] / sngdbl_jpsi_pAu_inv_yield[arm];
	    abs2 = pp_sngdbl_jpsi_inv_yield_err[arm] / pp_sngdbl_jpsi_inv_yield[arm];
	    RpAu_sngdbl_jpsi_cent_err[arm] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_jpsi_cent[arm];
	  cout << "RpAu jpsi[" << arm << "]: " <<  RpAu_sngdbl_jpsi_cent[arm] <<  endl;
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
       
       sngdbl_psi2s_trig_eff[0][0] = 0.7792;
       sngdbl_psi2s_trig_eff[0][1] = 0.6292;
	 sngdbl_psi2s_trig_eff[0][2] = 0;
       sngdbl_psi2s_trig_eff[0][3] = 0;
	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("h2",f2);
       
       sngdbl_psi2s_trig_eff[1][0] = 0.7604;
       sngdbl_psi2s_trig_eff[1][1] = 0.7708;
       sngdbl_psi2s_trig_eff[1][2] = 0.7690;
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
     double pp_sngdbl_psi2s[narm] = {571,451};
    double pp_sngdbl_psi2s_err[narm] = {36,37};

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
     
       pp_sngdbl_psi2s_trig_eff[0] = 0.7654;
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //pp_sngdbl_psi2s_trig_eff[1] = t2->GetBinContent(1);
     
       pp_sngdbl_psi2s_trig_eff[1] = 0.7390;
       
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

    // using J/psi and psi(2S) results from weighted average 


  //====================================================================================
  
    // dAu
    double RdAu_jpsi[3] = {10, 0.76, 10};  // PPG125, AN933 page 53
    double RdAu_jpsi_err[3] = {0,0.03, 0};

    double RdAu_mid[3] = {10, 0.54, 10};
    double RdAu_mid_err[3] = {0,0.11, 0};

    double RdAu_err_sys_up[3] = {0,0.19, 0};
    double RdAu_err_sys_down[3] = {0,0.16, 0};

    // PHENIX 3HE+Au
    double RHeAu_jpsi[2] = {0.89, 0.75};
    double RHeAu_jpsi_err[2] = {0.04, 0.04};

    double RHeAu_psi2s[2] = {0.56, 0.57};
    double RHeAu_psi2s_err[2] = {0.24, 0.16};

    double RHeAu_err_sys_up[2] = {0};
    double RHeAu_err_sys_down[2] = {0};

    // PHENIX p+Al
    //  double RpAl_jpsi[2] = {1.09, 0.98};
    //  double RpAl_jpsi_err[2] = {0.03, 0.04};
   
    //  double RpAl_jpsi[2] = {1.13, 1.13};
    //  double RpAl_jpsi_err[2] = {0.02, 0.03};
    double RpAl_jpsi[2] = {1.11, 1.13};
    double RpAl_jpsi_err[2] = {0.03, 0.03};

    // double RpAl_psi2s[2] = {0.57, 0.87};  // no mutr correction
    // double RpAl_psi2s[2] = {0.70, 1.03};
    // double RpAl_psi2s_err[2] = {0.19, 0.16};
    double RpAl_psi2s[2] = {0.65, 1.03};
    double RpAl_psi2s_err[2] = {0.21, 0.16};


    const Double_t err_x[2] = {0,0};
    const Double_t err_x4[3] = {0,0,0};
    const Double_t width = 0.145;
    const Double_t jpsi_ppg228_sys[narm] = {0.04, 0.06};

        
    double sys_err_228[2] = {0.06, 0.07};

    TGraphErrors *gr0 = new TGraphErrors(3,rap_array1,RdAu_jpsi,err_x4, RdAu_jpsi_err);
    TGraphErrors *gr1 = new TGraphErrors(3,rap_array1,RdAu_mid,err_x4, RdAu_mid_err);
    TGraphErrors *gr2 = new TGraphErrors(2,rap_array2,RpAl_jpsi, err_x, RpAl_jpsi_err);  
    TGraphErrors *gr3 = new TGraphErrors(2,rap_array2,RpAl_psi2s, err_x, RpAl_psi2s_err); 
    TGraphErrors *gr4 = new TGraphErrors(2,rap_array3,RpAu_sngdbl_jpsi_cent,err_x, RpAu_sngdbl_jpsi_cent_err);
    TGraphErrors *gr5 = new TGraphErrors(2,rap_array3,RpAu_sngdbl_psi2s_cent,err_x, RpAu_sngdbl_psi2s_cent_err);
    TGraphErrors *gr6 = new TGraphErrors(2,rap_array4,RHeAu_jpsi, err_x, RHeAu_jpsi_err);
    TGraphErrors *gr7 = new TGraphErrors(2,rap_array4,RHeAu_psi2s,err_x4, RHeAu_psi2s_err);  

    // TGraphErrors *gr5 = new TGraphErrors(18,rap_array_228,RpAu_jpsi_ppg228_cent,err_x, RpAu_jpsi_ppg228_cent_err );
   // TGraphErrors *gr7 = new TGraphErrors(4,rap_array4,RpAu_jpsi_sng_cent,err_x,RpAu_jpsi_sng_cent_err);
   // TGraphErrors *gr9 = new TGraphErrors(2,rap_array4,RpAu_sngdbl_psi2s_cent,err_x, RpAu_sngdbl_psi2s_cent_err);
   //  TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(2);
   //   TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(2);
  
   
   // TGraphErrors *gr2 = new TGraphErrors(2,rap_array1,RpAu_sngdbl_psi2s_cent,err_x, RpAu_sngdbl_psi2s_cent_err);
   // TGraphErrors *gr4 = new TGraphErrors(2,rap_array2,RpAu_sngdbl_jpsi_cent[0],err_x, RpAu_sngdbl_jpsi_cent_err[0]);
   // // TGraphErrors *gr6 = new TGraphErrors(18,rap_array_228,RpAu_jpsi_ppg228_cent[0],err_x, RpAu_jpsi_ppg228_cent_err[0] );;
   // //  TGraphErrors *gr8 = new TGraphErrors(4,rap_array4,RpAu_jpsi_sng_cent[0],err_x,RpAu_jpsi_sng_cent_err[0]);
   //  TGraphErrors *gr10 = new TGraphErrors(2,rap_array4,RpAu_sngdbl_psi2s_cent[0],err_x, RpAu_sngdbl_psi2s_cent_err[0]);
   // TGraphAsymmErrors *gr4_sys = new TGraphAsymmErrors(2);
   // TGraphAsymmErrors *gr6_sys = new TGraphAsymmErrors(2);
    

//    for(int i = 0; i < 2; i++)
//      {
//        // gr3_sys->SetPoint(i, rap_array2[i],RpAu_sngdbl_jpsi_cent[i]);   // need to also include array index here as well as arm index
//        // gr3_sys->SetPointError(i, width, width, jpsi_ppg228_sys[i], jpsi_ppg228_sys[i]);

//        gr5_sys->SetPoint(i, rap_array1[i],RpAu_jpsi_ppg228_cent[i]);   // need to also include array index here as well as arm index
//        gr5_sys->SetPointError(i, width, width, sys_err_228[i], sys_err_228[i]);

//        // gr4_sys->SetPoint(i, rap_array2[i],RpAu_sngdbl_jpsi_cent[0][i]);
//        // gr4_sys->SetPointError(i, width, width, jpsi_ppg228_sys[0][i], jpsi_ppg228_sys[0][i]);

//        gr6_sys->SetPoint(i, rap_array1[i],RpAu_jpsi_ppg228_cent[i]);
//        gr6_sys->SetPointError(i, width, width, sys_err_228[i], sys_err_228[i]);
	  
//      }

      
//     gr3_sys->SetMarkerColor(kMagneta+1);                     
//     gr3_sys->SetLineColor(kMagneta+1);                         
//     gr3_sys->SetFillStyle(0); 
//     gr3_sys->SetLineWidth(2);

//     gr4_sys->SetMarkerColor(kMagneta+1);                     
//     gr4_sys->SetLineColor(kMagneta+1);                         
//     gr4_sys->SetFillStyle(0); 
//     gr4_sys->SetLineWidth(2);

//     gr5_sys->SetMarkerColor(kBlack);                     
//     gr5_sys->SetLineColor(kBlack);                         
//     gr5_sys->SetFillStyle(0); 
//     gr5_sys->SetLineWidth(2);

//     gr6_sys->SetMarkerColor(kBlack);                     
//     gr6_sys->SetLineColor(kBlack);                         
//     gr6_sys->SetFillStyle(0); 
//     gr6_sys->SetLineWidth(2);


    TCanvas *c1 = new TCanvas("c1","c1",0,0,690,545);
  
    c1->SetTickx(1);
    c1->SetTicky(1);
   
    c1->cd();
    c1->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.04);
    gPad->SetBottomMargin(0.15);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetLineWidth(2);
   
 
    gr1->SetMarkerColor(kOrange-1);
    gr1->SetTitle("");
    gr1->SetLineWidth(3);
    gr1->SetLineColor(kOrange-1);
    gr1->SetMarkerStyle(23);
    gr1->SetMarkerSize(2.5);
  
  
    gr1->GetXaxis()->SetTitleOffset(0.89);
    gr1->GetXaxis()->SetTitleFont(102);
    gr1->GetXaxis()->SetLabelFont(102);
    gr1->GetXaxis()->SetTitle("y");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.06);
    
    gr1->GetXaxis()->SetRangeUser(-3, 3);
     gr1->GetXaxis()->SetNdivisions(6);
    // gr1->GetYaxis()->SetNdivisions(8);
   

    gr1->GetYaxis()->SetTitleOffset(0.65);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.075);
    gr1->GetYaxis()->SetLabelSize(0.06);
 
    gr1->GetYaxis()->SetTitle("R_{AB}");
 
    gr1->GetYaxis()->SetRangeUser(0,2);
    gr1->GetYaxis()->SetNdivisions(7);
 

    gr0->SetMarkerColor(kOrange-1);
    gr0->SetLineWidth(4);
    gr0->SetLineColor(kOrange-1);
    gr0->SetMarkerStyle(32);
    gr0->SetMarkerSize(2.5);                       


    gr2->SetMarkerColor(kRed-2);
    gr2->SetLineWidth(4);
    gr2->SetLineColor(kRed-2);
    gr2->SetMarkerStyle(24);
    gr2->SetMarkerSize(2);                       

    gr3->SetMarkerColor(kRed-2);
    gr3->SetLineWidth(3);
    gr3->SetLineColor(kRed-2);
    gr3->SetMarkerStyle(20); 
    gr3->SetMarkerSize(2);                       

    gr4->SetMarkerColor(kBlue-2);
    gr4->SetLineWidth(4);
    gr4->SetLineColor(kBlue-2);
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerSize(2);                       

    gr5->SetMarkerColor(kBlue-2);                     
    gr5->SetLineColor(kBlue-2);                         
    gr5->SetLineWidth(3);                           
    gr5->SetMarkerStyle(21);                       
    gr5->SetMarkerSize(2); 

    gr6->SetMarkerColor(kGray+3);                     
    gr6->SetLineColor(kGray+3);                         
    gr6->SetLineWidth(4);                           
    gr6->SetMarkerStyle(26);                       
    gr6->SetMarkerSize(2.5); 

    gr7->SetMarkerColor(kGray+3);                     
    gr7->SetLineColor(kGray+3);                         
    gr7->SetLineWidth(3);                           
    gr7->SetMarkerStyle(22);                       
    gr7->SetMarkerSize(2.5);                       
    
 
  

    gr1->Draw("AP");
    gr0->Draw("P");
    gr2->Draw("P");
    gr3->Draw("P");
    gr4->Draw("P");
    gr5->Draw("P");
    gr6->Draw("P");
    gr7->Draw("P");

    TLine *line = new TLine(-3,1,3,1);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");

   
     TLatex l2;
    l2.SetTextSize(0.0425);
    l2.SetTextAlign(12);
    l2.SetTextColor(kBlack);
  
    char text3[100];
    char text4[100];
  
    sprintf(text3,"Inclusive           |y|<0.35, 1.2<|y|<2.2");
    sprintf(text4,"");
    
    // no mid rapidity point
    // l2.SetTextFont(102);  
    // l2.DrawLatexNDC(0.17, 0.725, text3); // x0, y0
    // l2.DrawLatexNDC(0.17, 0.675, text4); // x0, y0

    // TLegend *leg1 = new TLegend(0.16, 0.77, 0.40, 0.94);  //(start x, start y, end x, end y)
    // //  leg1->SetFillColor(19); 
    // leg1->SetFillStyle(3003); 
    // leg1->SetLineWidth(1);
    // leg1->SetLineColor(0);
    // leg1->SetTextSize(0.0425); 
    // leg1->SetTextFont(102);

    // //  leg1->AddEntry(gr3,"R_{dAu} #psi(2S)","p");  
    // leg1->AddEntry(gr3,"R_{pAl} #psi(2S)","p");  
    // leg1->AddEntry(gr5,"R_{pAu} #psi(2S)","p");   
    // leg1->AddEntry(gr7,"R_{HeAu} #psi(2S)","p");  
    // leg1->Draw();

    // TLegend *leg2 = new TLegend(0.40, 0.77, 0.64, 0.94);  //(start x, start y, end x, end y)
    // leg2->SetFillColor(19); 
    // leg2->SetFillStyle(3003); 
    // leg2->SetLineWidth(1);
    // leg2->SetLineColor(0);
    // leg2->SetTextSize(0.0425); 
    // leg2->SetTextFont(102);

    // leg2->AddEntry(gr2,"R_{pAl} J/#psi, p+Al#sqrt{s_{NN}}=200GeV", "p");
    // leg2->AddEntry(gr4,"R_{pAu} J/#psi, p+Au#sqrt{s_{NN}}=200GeV", "p");
    // leg2->AddEntry(gr6,"R_{HeAu} J/#psi, ^{3}He+Au#sqrt{s_{NN}}=200GeV", "p");
    // leg2->Draw();

    // with mid rapidity point
 l2.SetTextFont(102);  
    l2.DrawLatexNDC(0.17, 0.69, text3); // x0, y0
    l2.DrawLatexNDC(0.17, 0.645, text4); // x0, y0

    TLegend *leg1 = new TLegend(0.16, 0.72, 0.40, 0.94);  //(start x, start y, end x, end y)
    //  leg1->SetFillColor(19); 
    leg1->SetFillStyle(3003); 
    leg1->SetLineWidth(1);
    leg1->SetLineColor(0);
    leg1->SetTextSize(0.0425); 
    leg1->SetTextFont(102);

    //  leg1->AddEntry(gr3,"R_{dAu} #psi(2S)","p");  
    leg1->AddEntry(gr3,"R_{pAl} #psi(2S)","p");  
    leg1->AddEntry(gr5,"R_{pAu} #psi(2S)","p");  
    leg1->AddEntry(gr1,"R_{dAu} #psi(2S)","p"); 
    leg1->AddEntry(gr7,"R_{HeAu} #psi(2S)","p");  
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.40, 0.72, 0.64, 0.94);  //(start x, start y, end x, end y)
    leg2->SetFillColor(19); 
    leg2->SetFillStyle(3003); 
    leg2->SetLineWidth(1);
    leg2->SetLineColor(0);
    leg2->SetTextSize(0.0425); 
    leg2->SetTextFont(102);

    leg2->AddEntry(gr2,"R_{pAl} J/#psi, p+Al#sqrt{s_{NN}}=200GeV", "p");
    leg2->AddEntry(gr4,"R_{pAu} J/#psi, p+Au#sqrt{s_{NN}}=200GeV", "p");
    leg2->AddEntry(gr0,"R_{dAu} J/#psi, d+Au#sqrt{s_{NN}}=200GeV", "p");
    leg2->AddEntry(gr6,"R_{HeAu} J/#psi, ^{3}He+Au#sqrt{s_{NN}}=200GeV", "p");
    leg2->Draw();

}




