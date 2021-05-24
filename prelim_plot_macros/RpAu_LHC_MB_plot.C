
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


void RpAu_LHC_MB_plot()
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
 
  const Double_t rap_array1[nbins] = {-1.7, 1.7};
  const Double_t rap_array2[3] = {-4.5, 0, 4.5};
  const Double_t rap_array5[2] = {-3.25, 3.25};
  const Double_t rap_array3[nbins] = {-3.71, 2.78};  // same for psi(2S) and jpsi
  const Double_t rap_array4[nbins] = {-3.25, 3.25};

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
       
       sngdbl_jpsi_trig_eff[1] = 0.7535; // bkgd hits, prob00
    	            
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

  
    // using J/psi and psi(2S) results from weighted average 


  //====================================================================================
  
    // PHENIX mid psip   // PPG 151 PRL 111,202301 (2013)
    double RdAu_mid[3] = {10, 0.54, 10};
    double RdAu_mid_err[3] = {0,0.11, 0};

    double RdAu_err_sys_up[3] = {0,0.19, 0};
    double RdAu_err_sys_down[3] = {0,0.16, 0};

  // dAu
    double RdAu_jpsi[3] = {10, 0.76, 10};  // PPG109 text tables - PRL107, 142301
    double RdAu_jpsi_err[3] = {0,0.03, 0};
    double RdAu_jpsi_sys[3] = {10, 0.06, 0};  

    // PHENIX PPG228 Integrated result (confirmw ith Sanghoon about efficiencies)
    double RpAu_ppg228[2] = {0.908, 0.737};
    double RpAu_ppg228_err[2] = {0.014, 0.0087};

   // LHCb J/psi
    double LHCb_jpsi[2] = {0.93, 0.63};
    //double LHCb_jpsi_err[2] = {0.08, 0.07};
    double LHCb_jpsi_err[2] = {0.07, 0.04};
   
    // LHCb psi(2S)
    double LHCb_psi2s[2] = {0.51, 0.47};
    double LHCb_psi2s_err[2] = {0.16, 0.09}; // actaully total error bar (quadrature of all uncertainties)
    double psi2s_LHCb_sys[2] = {0.15, 0.07};  // actually only stat err

    //ALICE psi(2S)
    double ALICE_psi2s[2] = {0.56, 0.48};   //   figure 3 https://www.hepdata.net/record/ins1296307
    double ALICE_psi2s_err[2] = {0.1, 0.06};

    double ALICE_psi2s_sys[2] = {0.06, 0.05};

    double ALICE_jpsi[2] = {1.08, 0.70};   // https://www.hepdata.net/record/ins1251898    figure 1
    double ALICE_jpsi_err[2] = {0.01, 0.01};

    double ALICE_jpsi_sys[2] = {0.09, 0.05};  // boxes are uncorr sys only

    const Double_t err_x[2] = {0,0};
    const Double_t err_x4[3] = {0,0,0};
    const Double_t width = 0.15;
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

    TGraphErrors *gr6 = new TGraphErrors(2,rap_array1,RpAu_sngdbl_psi2s_cent,err_x, RpAu_sngdbl_psi2s_cent_err);
    TGraphErrors *gr2 = new TGraphErrors(2,rap_array1,RpAu_sngdbl_jpsi_cent,err_x, RpAu_sngdbl_jpsi_cent_err);
    TGraphErrors *gr22 = new TGraphErrors(2,rap_array1,RpAu_sngdbl_jpsi_cent,err_x, RpAu_sngdbl_jpsi_cent_err);
    TGraphErrors *gr7 = new TGraphErrors(3,rap_array2,RdAu_jpsi,err_x4, RdAu_jpsi_err);  // phenix mid
    TGraphErrors *gr77 = new TGraphErrors(3,rap_array2,RdAu_jpsi,err_x4, RdAu_jpsi_err);  // phenix mid
    TGraphErrors *gr1 = new TGraphErrors(3,rap_array2,RdAu_mid,err_x4, RdAu_mid_err);  // phenix mid
    TGraphErrors *gr4 = new TGraphErrors(2,rap_array4,LHCb_jpsi, err_x, LHCb_jpsi_err);  // LHCb jpsi
    TGraphErrors *gr44 = new TGraphErrors(2,rap_array4,LHCb_jpsi, err_x, LHCb_jpsi_err);  // LHCb jpsi
    TGraphErrors *gr5 = new TGraphErrors(2,rap_array5,LHCb_psi2s, err_x, LHCb_psi2s_err);  // LHCb psi(2S)
    TGraphErrors *gr3 = new TGraphErrors(2,rap_array3,ALICE_psi2s, err_x, ALICE_psi2s_err);  // LHCb psi(2S)
    TGraphErrors *gr8 = new TGraphErrors(2,rap_array3,ALICE_jpsi, err_x, ALICE_jpsi_err);  // LHCb psi(2S)
    TGraphErrors *gr88 = new TGraphErrors(2,rap_array3,ALICE_jpsi, err_x, ALICE_jpsi_err);  // LHCb psi(2S)

    TGraphAsymmErrors *gr7_sys = new TGraphAsymmErrors(3);  // dAu J/psi
    TGraphAsymmErrors *gr1_sys = new TGraphAsymmErrors(3);  // dAu psip
    TGraphAsymmErrors *gr2_sys = new TGraphAsymmErrors(2);
    TGraphAsymmErrors *gr3_sys = new TGraphAsymmErrors(2);
    TGraphAsymmErrors *gr5_sys = new TGraphAsymmErrors(2);
    TGraphAsymmErrors *gr6_sys = new TGraphAsymmErrors(2);
    TGraphAsymmErrors *gr8_sys = new TGraphAsymmErrors(2);
  

    for(int i = 0; i < 2; i++)
      {
	gr2_sys->SetPoint(i, rap_array1[i],RpAu_sngdbl_jpsi_cent[i]); 
	gr2_sys->SetPointError(i, width, width, jpsi_sngdbl_sys[i], jpsi_sngdbl_sys[i]);

	gr6_sys->SetPoint(i, rap_array1[i],RpAu_sngdbl_psi2s_cent[i]);  
	gr6_sys->SetPointError(i, width, width, psi2s_sngdbl_sys[i], psi2s_sngdbl_sys[i]);

	gr5_sys->SetPoint(i, rap_array4[i],LHCb_psi2s[i]);  
	gr5_sys->SetPointError(i, 0, 0, psi2s_LHCb_sys[i], psi2s_LHCb_sys[i]);  // LHCb does not use a box or gloabl error

	gr8_sys->SetPoint(i, rap_array3[i],ALICE_jpsi[i]);  
	gr8_sys->SetPointError(i, width, width, ALICE_jpsi_sys[i], ALICE_jpsi_sys[i]);

	gr3_sys->SetPoint(i, rap_array3[i],ALICE_psi2s[i]);  
	gr3_sys->SetPointError(i, width, width, ALICE_psi2s_sys[i], ALICE_psi2s_sys[i]);
	  
      }

    for(int i = 0; i < 2; i++)
      {
	gr1_sys->SetPoint(i, rap_array2[i],RdAu_mid[i]);  // psi2s
	gr1_sys->SetPointError(i, width, width,RdAu_err_sys_up[i], RdAu_err_sys_down[i]);

	gr7_sys->SetPoint(i, rap_array2[i],RdAu_jpsi[i]);  // psi2s
	gr7_sys->SetPointError(i, width, width,RdAu_jpsi_sys[i], RdAu_jpsi_sys[i]);

      }

    gr1_sys->SetMarkerColor(kOrange-1);                     
    gr1_sys->SetLineColor(kOrange-1);                         
    gr1_sys->SetFillStyle(0); 
    gr1_sys->SetLineWidth(1);

    gr7_sys->SetMarkerColor(kOrange-1);                     
    gr7_sys->SetLineColor(kOrange-1);                         
    gr7_sys->SetFillStyle(0); 
    gr7_sys->SetLineWidth(1);

    gr2_sys->SetMarkerColor(kBlue-2);                     
    gr2_sys->SetLineColor(kBlue-2);                         
    gr2_sys->SetFillStyle(0); 
    gr2_sys->SetLineWidth(1);

    gr6_sys->SetMarkerColor(kBlue-2);                     
    gr6_sys->SetLineColor(kBlue-2);                         
    gr6_sys->SetFillStyle(0); 
    gr6_sys->SetLineWidth(1);

    gr5_sys->SetMarkerColor(kRed-2);                     
    gr5_sys->SetLineColor(kRed-2);                         
    gr5_sys->SetFillStyle(0); 
    gr5_sys->SetLineWidth(1);

    gr8_sys->SetMarkerColor(kGray+3);                     
    gr8_sys->SetLineColor(kGray+3);                         
    gr8_sys->SetFillStyle(0); 
    gr8_sys->SetLineWidth(1);

    gr3_sys->SetMarkerColor(kGray+3);                     
    gr3_sys->SetLineColor(kGray+3);                         
    gr3_sys->SetFillStyle(0); 
    gr3_sys->SetLineWidth(1);


  
    
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
    gr1->GetXaxis()->SetTitle("y");
    gr1->GetXaxis()->SetTitleSize(0.065);
    gr1->GetXaxis()->SetLabelSize(0.05);

    gr1->GetXaxis()->SetRangeUser(-4.5,4.5);
    gr1->GetXaxis()->SetNdivisions(5);

    gr1->GetYaxis()->SetTitleOffset(0.8);
    gr1->GetYaxis()->SetTitleFont(102);
    gr1->GetYaxis()->SetLabelFont(102);
    gr1->GetYaxis()->SetTitleSize(0.075);
    gr1->GetYaxis()->SetLabelSize(0.05);

    gr1->GetYaxis()->SetTitle("R_{AB}");
 
    gr1->GetYaxis()->SetRangeUser(0,2.5);
    gr1->GetYaxis()->SetNdivisions(7);
 
   
    gr1->SetMarkerColor(kOrange-1);
    gr1->SetLineWidth(3);
    gr1->SetLineColor(kOrange-1);
    gr1->SetMarkerStyle(23); 
    gr1->SetMarkerSize(2);                       

   
    gr4->SetMarkerColor(kRed-2);  // J/psi LHCb
    gr4->SetLineWidth(4);
    gr4->SetLineColor(kRed-2);
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerSize(2); 

    gr44->SetMarkerColor(kRed-2);  // J/psi LHCb
    gr44->SetLineWidth(4);
    gr44->SetLineColor(kRed-2);
    gr44->SetMarkerStyle(25);
    gr44->SetMarkerSize(2);                       

    gr5->SetMarkerColor(kRed-2);
    gr5->SetLineWidth(3);
    gr5->SetLineColor(kRed-2);
    gr5->SetMarkerStyle(21); 
    gr5->SetMarkerSize(2);  

    gr2->SetMarkerColor(kBlue-2);  // J/psi pAu
    gr2->SetLineWidth(3);
    gr2->SetLineColor(kBlue-2);
    gr2->SetMarkerStyle(24);
    gr2->SetMarkerSize(1.75);                       
    
    gr22->SetMarkerColor(kBlue-2);  // J/psi pAu
    gr22->SetLineWidth(3);
    gr22->SetLineColor(kBlue-2);
    gr22->SetMarkerStyle(24);
    gr22->SetMarkerSize(1.75);                       

    gr6->SetMarkerColor(kBlue-2);                     
    gr6->SetLineColor(kBlue-2);                         
    gr6->SetLineWidth(3);                           
    gr6->SetMarkerStyle(20);                       
    gr6->SetMarkerSize(2);                        
   
    gr7->SetMarkerColor(kOrange-1);  // jpsi dAu
    gr7->SetLineWidth(3);
    gr7->SetLineColor(kWhite);
    gr7->SetMarkerStyle(32);
    gr7->SetMarkerSize(1.5); 

    gr77->SetMarkerColor(kOrange-1);  // jpsi dAu
    gr77->SetLineWidth(3);
    gr77->SetLineColor(kWhite);
    gr77->SetMarkerStyle(32);
    gr77->SetMarkerSize(1.5);        

    gr8->SetMarkerColor(kGray+3);   // J/spi ALICE                                      
    gr8->SetLineColor(kGray+3);                         
    gr8->SetLineWidth(3);                           
    gr8->SetMarkerStyle(26);                       
    gr8->SetMarkerSize(1.5); 

    gr88->SetMarkerColor(kGray+3);   // J/spi ALICE                                      
    gr88->SetLineColor(kGray+3);                         
    gr88->SetLineWidth(3);                           
    gr88->SetMarkerStyle(26);                       
    gr88->SetMarkerSize(1.5); 

    gr3->SetMarkerColor(kGray+3);  
    gr3->SetLineColor(kGray+3);                         
    gr3->SetLineWidth(3);                           
    gr3->SetMarkerStyle(22);                       
    gr3->SetMarkerSize(2.2);                       
    

    gr1->Draw("AP");
    gr2->Draw("P");
    gr22->Draw("P");
    gr3->Draw("P");
    gr4->Draw("P");
    // gr44->Draw("P");
    gr5->Draw("P");
    gr6->Draw("P");
    gr7->Draw("P");
    gr77->Draw("P");
    gr8->Draw("P");
    gr88->Draw("P");
    gr1_sys->Draw("e2same"); // dAu psip
    gr2_sys->Draw("e2same"); // pAu jpsi
    // gr5_sys->Draw("e1same"); //
    gr6_sys->Draw("e2same"); // pAu psip
    gr3_sys->Draw("e2same"); // pAu psip
    gr7_sys->Draw("e2same"); // pAu psip
    gr8_sys->Draw("e2same"); // pAu psip
  
     TLatex l2;
    l2.SetTextSize(0.0325);
    l2.SetTextAlign(12);
    l2.SetTextColor(kBlack);
  
    char text3[100];
    char text4[100];

  sprintf(text4,"Inclusive");
    l2.SetTextFont(102);  
    //  l2.DrawLatexNDC(0.24, 0.91, text3); // x0, y0
    //  l2.DrawLatexNDC(0.76, 0.55, text4); // x0, y0

    TLine *line = new TLine(-4.4,1,4.4,1);
    line->SetLineWidth(4);
    line->SetLineColor(kBlack);
    line->SetLineStyle(3);
    line->Draw("L");

   
    double global_sys_pAu = 0.121;  

    TBox *box = new TBox(4.35, 1-global_sys_pAu, 4.55, 1+global_sys_pAu);
 
    box->SetFillStyle(1000);
    box->SetFillColorAlpha(kBlue+2,0.5);
    box->Draw();

    double global_sys_pPb = 0.07;    // https://www.hepdata.net/record/ins1251898  figure 1

    TBox *box0 = new TBox(3.95, 1-global_sys_pPb, 4.15, 1+global_sys_pPb); // Alice 
    box0->SetFillStyle(1000);
    box0->SetFillColorAlpha(13,0.5);
    box0->Draw();


    double global_sys_dAu = 0.10;   // // 5.63% from Ncoll and then 10% form BBC.  No bias corrcetion in dAu centrality note

    TBox *box2 = new TBox(4.15, 1-global_sys_dAu, 4.35, 1+global_sys_dAu);
    box2->SetFillStyle(1000);
    box2->SetFillColorAlpha(kOrange-1,0.5);
    box2->Draw();

    // gStyle->SetEndErrorSize(1);
    

    TLegend *leg2 = new TLegend(0.16, 0.58, 0.44, 0.93);  //(start x, start y, end x, end y)
    leg2->SetFillColor(19); 
    leg2->SetFillStyle(3003); 
    leg2->SetLineWidth(1);
    leg2->SetLineColor(0);
    leg2->SetTextSize(0.0325); 
    leg2->SetTextFont(102);

  
    leg2->AddEntry(gr2,"R_{pAu} J/#psi PRC 102, 014902", "p");
    leg2->AddEntry(gr6,"R_{pAu} #psi(2S) PHENIX, p+Au #sqrt{s_{NN}}=200 GeV", "p"); // PHENIX p+Au #sqrt{s_{NN}}=200 GeV", "p");
    leg2->AddEntry(gr7,"R_{dAu} J/#psi PRL 111,202301 (2013)", "p");
    leg2->AddEntry(gr1,"R_{dAu} #psi(2S) PHENIX, d+Au #sqrt{s_{NN}}=200 GeV","p");   // PHENIX, d+Au #sqrt{s_{NN}}=200 GeV", "p");
    leg2->AddEntry(gr4,"R_{pPb} J/#psi JHEP 1603 (2016) 133", "p");
    leg2->AddEntry(gr5,"R_{pPb} #psi(2S) LHCb, p+Pb #sqrt{s_{NN}}=5 TeV", "p");  // LHCb, p+Pb #sqrt{s_{NN}}=5 TeV", "p");
    leg2->AddEntry(gr8,"R_{pPb} J/#psi JHEP 12 (2014) 073", "p");
    leg2->AddEntry(gr3,"R_{pPb} #psi(2S) ALICE, p+Pb #sqrt{s_{NN}}=5.02 TeV", "p"); // ALICE, p+Pb #sqrt{s_{NN}}=5 TeV", "p");
   
   
    leg2->Draw();



 if(north_arm)
      c1->SaveAs(Form("../psi2S_pdf/RpAu_LHC_MB_plot.pdf"));
    else
      c1->SaveAs(Form("../psi2S_pdf/RpAu_LHC_MB_plot.pdf"));



}




