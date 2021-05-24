
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


void RpAu_PHENIX_MB_plot()
{
    gROOT->ForceStyle();

  bool north_arm = true;
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
  double pi = 3.141593;

  const Double_t rap_array1[nbins] = {-1.7, 1.7};
  const Double_t rap_array2[nbins] = {-1.7, 1.7};
  const Double_t rap_array3[nbins] = {-3.5, 3.5};
  const Double_t rap_array4[nbins] = {-1.7, 1.7};

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
       
       sngdbl_jpsi_trig_eff[0] = 0.6799;  // bkgd hits, prob00
     	            
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("h1",f2);
       
       sngdbl_jpsi_trig_eff[1] = 0.7535;  // bkgd hits, prob00
    	            
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
       pp_sngdbl_jpsi_trig_eff[0] = 0.7450;  // bkgd hits, prob00
       
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/trigger_eff_sngdbl/sanghoon_simulations_final/Run15pp_sngdbl_jpsi_trigeff_N.root");
       file_data2->GetObject("geff_y_int",t2);
       //ypp_sngdbl_jpsi_trig_eff[1] = t2->GetBinContent(1);
       pp_sngdbl_jpsi_trig_eff[1] = 0.7185;  // bkgd hits, prob00
       
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

   
    double sngdbl_psi2s[narm] = {216, 321};
    double sngdbl_psi2s_err[narm] = {39, 34};
     
    double sngdbl_psi2s_acc_eff[narm] = {0};
    double sngdbl_psi2s_trig_eff[narm] = {0};
    
  {
    
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       // RUN15pAu
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_S.root");
       file_data1->GetObject("S_denom_clone",e1);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_psi2s_acc_eff[0] = e1->GetEfficiency(pt+1);

       sngdbl_psi2s_acc_eff[0] = 0.0332;
     
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/sng_dbl_tracks/sanghoon_simulations/Run15pAu_sngdbl_psi2s_acceff_cent_N.root");
       file_data2->GetObject("N_denom_clone",e2);
       for(int pt =0; pt < 4; pt++)
	 sngdbl_psi2s_acc_eff[1] = e2->GetEfficiency(pt+1);

       sngdbl_psi2s_acc_eff[1] = 0.0447;
    
       file_data1->Close();
       file_data2->Close();
    }
       
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RUN15pAu
     {
      
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       file_data1 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_S.root");
       file_data1->GetObject("h2",f1);
       
       sngdbl_psi2s_trig_eff[0] = 0.6938;  // bkgd hits, prob00
     
       file_data2 = TFile::Open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run15pAu/trigger_eff_sngdbl/sanghoon_simulations/Run15pAu_sngdbl_psi2s_trigeff_N.root");
       file_data2->GetObject("h2",f2);
       
       sngdbl_psi2s_trig_eff[1] = 0.7704;  // bkgd hits, prob00
   
	file_data1->Close();
	file_data2->Close();
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     }

    double sngdbl_psi2s_pAu_inv_yield[narm] = {0};
    double sngdbl_psi2s_pAu_inv_yield_err[narm] = {0};

    // calculcate pAu inv. yield for psi2s fvtx
    for(int arm = 0; arm < 2; arm++)
      {
	FACTOR = 0;

 	    sngdbl_psi2s_pAu_inv_yield[arm] = bias * sngdbl_psi2s[arm] / ( ( FACTOR + sngdbl_psi2s_trig_eff[arm] ) * ( FACTOR + sngdbl_psi2s_acc_eff[arm] ) * pAu_MB[arm] );
 	    sngdbl_psi2s_pAu_inv_yield_err[arm] = bias * sngdbl_psi2s_err[arm] / ( sngdbl_psi2s_trig_eff[arm]*sngdbl_psi2s_acc_eff[arm] * pAu_MB[arm]  );
 	    // cout << "sngdbl_psi2s_pAu inv yield: " << sngdbl_psi2s_pAu_inv_yield[arm] << "+\-" << sngdbl_psi2s_pAu_inv_yield_err[arm] << endl;
 	  
      }

  
    double pp_sngdbl_psi2s_inv_yield[narm] = {0}; 
    double pp_sngdbl_psi2s_inv_yield_err[narm] = {0}; 
    
    double pp_sngdbl_psi2s[narm] = {571, 451};
    double pp_sngdbl_psi2s_err[narm] = {37, 36};

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

	    pp_sngdbl_psi2s_inv_yield[arm] = pp_sngdbl_psi2s[arm] * eff_bbc / ( (FACTOR + pp_sngdbl_psi2s_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_psi2s_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
	    pp_sngdbl_psi2s_inv_yield_err[arm] = pp_sngdbl_psi2s_err[arm] * eff_bbc / ( ( FACTOR + pp_sngdbl_psi2s_trig_eff[arm] ) * ( FACTOR + pp_sngdbl_psi2s_acc_eff[arm] ) * pp_MB[arm] * eff_jpsi );
	    // cout << "pp sngdbl_psi2s fvtx inv yield: " << pp_sngdbl_psi2s_inv_yield[arm] << "+\-" << pp_sngdbl_psi2s_inv_yield_err[arm] << endl;
	  
      }
     double RpAu_sngdbl_psi2s_cent[narm] = {0};
     double RpAu_sngdbl_psi2s_cent_err[narm] = {0};
    
    // psi2s sngdbl
    for(int arm = 0; arm < 2; arm++)
      {
 
 	    abs1 = 0;
 	    abs2 = 0;

 	    RpAu_sngdbl_psi2s_cent[arm] =( sngdbl_psi2s_pAu_inv_yield[arm] ) / ( ncoll * pp_sngdbl_psi2s_inv_yield[arm] ) ;
 	    abs1 = sngdbl_psi2s_pAu_inv_yield_err[arm] / sngdbl_psi2s_pAu_inv_yield[arm];
 	    abs2 = pp_sngdbl_psi2s_inv_yield_err[arm] / pp_sngdbl_psi2s_inv_yield[arm];
 	    RpAu_sngdbl_psi2s_cent_err[arm] = sqrt( abs1 * abs1 + abs2 * abs2)*RpAu_sngdbl_psi2s_cent[arm];
	    cout << "RpAu psi2s[" << arm << "]: " <<  RpAu_sngdbl_psi2s_cent[arm] << ", ncoll: " << ncoll << ", raw psi2s yield: " <<  sngdbl_psi2s[arm] << ", psi2s pAu inv yield: " <<  sngdbl_psi2s_pAu_inv_yield[arm] <<  endl;
 	  
      }

  
  //====================================================================================
  
    // Shao MB predicitons

    TGraphAsymmErrors *p1; 
    TGraphAsymmErrors *p2; 
    TGraphAsymmErrors *p3; 
    TGraphAsymmErrors *p4; 
    TGraphAsymmErrors *p5; 
    TGraphAsymmErrors *p6; 
    TGraphAsymmErrors *p7; 
    TGraphAsymmErrors *p8;

    TFile *file_data3;
    TFile *file_data4;


 // Shao predciitons, jpsi  (filenames were mislabled)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
      file_data3 = TFile::Open("psi2S_root_filess/psi2s_RpAu_MB_predictions_rap.root");
     
      file_data3->GetObject("RpAu_epps_MB_predictions_bkwd",p1);  
      file_data3->GetObject("RpAu_epps_MB_predictions_fwd",p2);

      file_data3->GetObject("RpAu_ncteq_MB_predictions_bkwd",p3);
      file_data3->GetObject("RpAu_ncteq_MB_predictions_fwd",p4);
    
      file_data3->Close();
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    // Shao predciitons,psi2s  (filenames were mislabled)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
      file_data4 = TFile::Open("psi2S_root_files/jpsi_RpAu_MB_predictions_rap.root");
    
      file_data4->GetObject("RpAu_epps_MB_predictions_bkwd",p5);
      file_data4->GetObject("RpAu_epps_MB_predictions_fwd",p6);

      file_data4->GetObject("RpAu_ncteq_MB_predictions_bkwd",p7);
      file_data4->GetObject("RpAu_ncteq_MB_predictions_fwd",p8);
    
      file_data4->Close();
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

        
    // double RdAu_mid[3] = {10, 0.54, 10};
    double RdAu_mid[2] = {10, 10};
    double RdAu_mid_err[2] = {0,0.11};
    double RdAu_err_sys_up[3] = {0,0.19, 0};
    double RdAu_err_sys_down[3] = {0,0.16, 0};

    
    const Double_t err_x[2] = {0,0};
    const Double_t width = 0.14;
    const Double_t jpsi_ppg228_sys[narm] = {0.04, 0.06};    
    double sys_err_228[2] = {0.06, 0.07};


    double jpsi_sngdbl_percent[2] = {0.07, 0.07};    
    double psi2s_sngdbl_percent[2] = {0.15, 0.14};  

    double jpsi_sngdbl_sys[2] = {0};
    double psi2s_sngdbl_sys[2] = {0};
  

    for(int arm = 0; arm < 2; arm++)
     {
       jpsi_sngdbl_sys[arm] = jpsi_sngdbl_percent[arm]*RpAu_sngdbl_jpsi_cent[arm];
       psi2s_sngdbl_sys[arm] = psi2s_sngdbl_percent[arm]*RpAu_sngdbl_psi2s_cent[arm];
     }


    TGraphErrors *gr1 = new TGraphErrors(2,rap_array3,RdAu_mid,err_x, RdAu_mid_err);
    TGraphErrors *gr3 = new TGraphErrors(2,rap_array2,RpAu_sngdbl_jpsi_cent,err_x, RpAu_sngdbl_jpsi_cent_err);
    TGraphErrors *gr5 = new TGraphErrors(2,rap_array1,RpAu_sngdbl_psi2s_cent,err_x, RpAu_sngdbl_psi2s_cent_err);

    TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(2);
    TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(2);
  

    for(int i = 0; i < 2; i++)
      {
	gr3_sys->SetPoint(i, rap_array2[i],RpAu_sngdbl_jpsi_cent[i]); 
	gr3_sys->SetPointError(i, width, width, jpsi_sngdbl_sys[i], jpsi_sngdbl_sys[i]);

	gr5_sys->SetPoint(i, rap_array1[i],RpAu_sngdbl_psi2s_cent[i]);  
	gr5_sys->SetPointError(i, width, width, psi2s_sngdbl_sys[i], psi2s_sngdbl_sys[i]);
	  
      }

    gr3_sys->SetMarkerColor(kBlack);                     
    gr3_sys->SetLineColor(kBlack);                         
    gr3_sys->SetFillStyle(0); 
    gr3_sys->SetLineWidth(1);

    gr5_sys->SetMarkerColor(kBlue+1);                     
    gr5_sys->SetLineColor(kBlue+1);                         
    gr5_sys->SetFillStyle(0); 
    gr5_sys->SetLineWidth(1);


  
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
   
 
    gr1->SetMarkerColor(kBlack);
    gr1->SetTitle("");
    gr1->SetLineWidth(3);
    gr1->SetLineColor(kBlack);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(2);
  
  
    gr1->GetXaxis()->SetTitleOffset(1.05);
    gr1->GetXaxis()->SetTitleFont(102);
    gr1->GetXaxis()->SetLabelFont(102);
    gr1->GetXaxis()->SetTitle("y");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.05);

    gr1->GetXaxis()->SetRangeUser(-3, 3);
    gr1->GetXaxis()->SetNdivisions(6);

    gr1->GetYaxis()->SetTitleOffset(0.8);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.075);
    gr1->GetYaxis()->SetLabelSize(0.05);

    gr1->GetYaxis()->SetTitle("R_{pAu}");
 
    gr1->GetYaxis()->SetRangeUser(0,2.5);
    gr1->GetYaxis()->SetNdivisions(7);
 
  
     //////////////////////////////////////////////////     J/psi         
    gr3->SetMarkerColorAlpha(kBlack, 0.8);
    gr3->SetLineWidth(3);
    gr3->SetLineColorAlpha(kBlack, 0.8);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(2);        

    gr5->SetMarkerColorAlpha(kBlue+1, 0.8);
    gr5->SetLineWidth(3);
    gr5->SetLineColorAlpha(kBlue+1, 0.8);
    gr5->SetMarkerStyle(21);
    gr5->SetMarkerSize(2);                       
    //////////////////////////////////////////////////


   
    /////////////////// predicitons - J/psi 
    // epps south
    p1->SetFillStyle(3224);
    p1->SetLineColor(kBlack);
    p1->SetFillColorAlpha(kBlack, 0.5);
    p1->SetLineWidth(1);
    // epps north
    p2->SetFillStyle(3224);
    p2->SetLineColor(kBlack);
    p2->SetFillColorAlpha(kBlack, 0.5);
    p2->SetLineWidth(1);
    // ncteq south
    p3->SetFillStyle(3224);
    p3->SetLineColor(kBlack);
    p3->SetFillColorAlpha(kBlack, 0.5);
    p3->SetLineWidth(1);
    // ncteq north
    p4->SetFillStyle(3224);
    p4->SetLineColor(kBlack);
    p4->SetFillColorAlpha(kBlack, 0.5);
    p4->SetLineWidth(1);

    /////////////////// predicitons - psi(2S)
    // epps south
    p5->SetFillStyle(3001);
    p5->SetLineColor(kBlue+1);
    p5->SetFillColorAlpha(kBlue+1, 0.3);
    p5->SetLineWidth(1);
    // epps north
    p6->SetFillStyle(3001);
    p6->SetLineColor(kBlue+1);
    p6->SetFillColorAlpha(kBlue+1, 0.3);
    p6->SetLineWidth(1);
    // ncteq south
    p7->SetFillStyle(3001);
    p7->SetLineColor(kBlue+1);
    p7->SetFillColorAlpha(kBlue+1, 0.3);
    p7->SetLineWidth(1);
    // ncteq north
    p8->SetFillStyle(3001);
    p8->SetLineColor(kBlue+1);
    p8->SetFillColorAlpha(kBlue+1, 0.3);
    p8->SetLineWidth(1);
   
    

    gr1->Draw("AP");
    gr3->Draw("P");
    gr3_sys->Draw("e2same");
    gr5->Draw("P");
    gr5_sys->Draw("e2same");
    // p1->Draw("e3same");
    // p2->Draw("e3same");
    // p5->Draw("e3same");
    // p6->Draw("e3same");

    p3->Draw("e3same");
    p4->Draw("e3same");
    p7->Draw("e3same");
    p8->Draw("e3same");

      
    TLatex l2;
    l2.SetTextSize(0.042);
    l2.SetTextAlign(12);
    l2.SetTextColor(kBlack);
  
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
    l2.DrawLatexNDC(0.24, 0.91, text4); // x0, y0
    //  l2.DrawLatexNDC(0.76, 0.91, text4); // x0, y0
    l2.SetTextSize(0.04);
    l2.DrawLatexNDC(0.24, 0.805, text5); // x0, y0

      
    TLine *line = new TLine(-3.0, 1, 3.0, 1);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");

    double global_sys = 0.121;  // RpAu J/psi PPG228 in both Norht and South (checked)

    TBox *box = new TBox(2.8,1-global_sys,3.0,1+global_sys);
    box->SetFillStyle(1000);
    box->SetFillColorAlpha(1,0.5);
    box->Draw();
    
     TLegend *leg = new TLegend(0.18, 0.83, 0.44, 0.89);  //(start x, start y, end x, end y)
    leg->SetFillColor(19); 
    leg->SetFillStyle(3003); 
    leg->SetLineWidth(1);
    leg->SetLineColor(0);
    leg->SetTextSize(0.0425); 
    leg->SetTextFont(102);
    leg->AddEntry(gr3,"R_{pAu} J/#psi, p+Au #sqrt{s_{NN}}=200 GeV", "p");

    leg->Draw();
    
    TLegend *leg2 = new TLegend(0.18, 0.61, 0.44, 0.79);  //(start x, start y, end x, end y)
    leg2->SetFillColor(19); 
    leg2->SetFillStyle(3003); 
    leg2->SetLineWidth(1);
    leg2->SetLineColor(0);
    leg2->SetTextSize(0.0425); 
    leg2->SetTextFont(102);

    leg2->AddEntry(gr5,"R_{pAu} #psi(2S), p+Au #sqrt{s_{NN}}=200 GeV", "p");
    leg2->AddEntry(p3,"J/#psi nCTEQ15 (Shao et al)", "f");
    leg2->AddEntry(p7,"#psi(2S) nCTEQ15 (Shao et al)", "f");
    
    leg2->Draw();
    
    c1->SaveAs(Form("../psi2S_pdf/RpAu_jpsi_psi2s_MB_rap.pdf"));
   
}




