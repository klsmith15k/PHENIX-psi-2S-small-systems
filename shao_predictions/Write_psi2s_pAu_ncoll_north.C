
#include <TF1.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <TVirtualFitter.h>
#include <TMatrixDSym.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLine.h>
#include <TRandom1.h>
#include <TPolyLine.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>

using namespace std;

void Write_psi2s_pAu_ncoll_north()
{

   bool north_arm = true;
  
  TFile *Run15pAu_prediction_files;

  const int bins_pAu = 3;  // cent bins
  const int columns = 5;

  char name100[500];

  const int narm = 2;

  double uepps_68[narm][columns][bins_pAu] = {0}; 
  double depps_68[narm][columns][bins_pAu] = {0};
  double cepps_68[narm][columns][bins_pAu] = {0};
	
  double epps16_env_68_min[narm][bins_pAu] = {0};
  double epps16_env_68_max[narm][bins_pAu] = {0};
  double epps16_env_68_RpAu_pred[narm][bins_pAu] = {0};

  double epps16_env_68_down[narm][bins_pAu] = {0};
  double epps16_env_68_up[narm][bins_pAu] = {0};
  
  double uncteq_68[narm][columns][bins_pAu] = {0}; 
  double dncteq_68[narm][columns][bins_pAu] = {0};
  double cncteq_68[narm][columns][bins_pAu] = {0};

  double ncteq15_env_68_min[narm][bins_pAu] = {0};
  double ncteq15_env_68_max[narm][bins_pAu] = {0};
  double ncteq15_env_68_RpAu_pred[narm][bins_pAu] = {0};

   double ncteq15_env_68_down[narm][bins_pAu] = {0};
  double ncteq15_env_68_up[narm][bins_pAu] = {0};

  double value_d  = 0;
  double value_c  = 0;
  double value_u = 0;

  double min_value = 0;
  double max_value = 0;

  double N_epps68_020_c = 0;   // using same c for central energy scale as in PPG228
 double N_epps68_2040_c = 0;
 double N_epps68_4085_c = 0;

 double N_ncteq68_020_c = 0;
 double N_ncteq68_2040_c = 0;
 double N_ncteq68_4085_c = 0;
 
 //  for(int arm = 0; arm < 2; arm++)
 {
      for(int i = 0; i < 18; i++)
	{
	  // if(arm == 0)
	    {
	     if(i == 0)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.Cen0To20.a1.dat"); // u_ncteq15_68_S;
	      if(i == 1)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.Cen0To20.a1.dat"); // 
	      if(i == 2)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.Cen0To20.a1.dat"); // u_epps_68_S;
	      if(i == 3)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.Cen0To20.a1.dat"); // u_epps_68_S;
	      if(i == 4)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.Cen0To20.a1.dat"); // d_epps_68_S;
	      if(i == 5)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.Cen0To20.a1.dat"); // c_epps_68_S;
	    }
	    {
	      if(i == 6)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.Cen20To40.a1.dat"); // u_ncteq15_68_S;
	      if(i == 7)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.Cen20To40.a1.dat"); // 
	      if(i == 8)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.Cen20To40.a1.dat"); // u_epps_68_S;
	      if(i == 9)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.Cen20To40.a1.dat"); // u_epps_68_S;
	      if(i == 10)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.Cen20To40.a1.dat"); // d_epps_68_S;
	      if(i == 11)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.Cen20To40.a1.dat"); // c_epps_68_S;
	    }
	    {
	      if(i == 12)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.Cen40To85.a1.dat"); // u_ncteq15_68_S;
	      if(i == 13)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.Cen40To85.a1.dat"); // 
	      if(i == 14)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.Cen40To85.a1.dat"); // u_epps_68_S;
	      if(i == 15)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.Cen40To85.a1.dat"); // u_epps_68_S;
	      if(i == 16)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.Cen40To85.a1.dat"); // d_epps_68_S;
	      if(i == 17)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.Cen40To85.a1.dat"); // c_epps_68_S;
	    }
	
      
	   double xmin = 0;
	  double xmax = 0;
	  double sigma_central = 0;
	  double sigma_min = 0;
	  double sigma_max = 0;
	  double sigma_pp = 0;
	  double RpAu_theory = 0;
	  double RpAu_theory_min = 0;
	  double RpAu_theory_max = 0;

	  int counter = 0;

	  ifstream Run15pAu_prediction_files(name100);
	  if(Run15pAu_prediction_files)
	    {
	      do 
		{
		  double f0, f1, f2, f3, f4;
		  Run15pAu_prediction_files >> f0 >> f1 >> f2 >> f3 >> f4;
		  //xmin = f0; xmax = f1; sigma_central = f2; sigma_min = f3; sigma_max = f4; sigma_pp = f5; RpAu_theory = f2/f5; RpAu_theory_min = f3/f5; RpAu_theory_max = f4/f5;
		  xmin = f0; xmax = f1; RpAu_theory = f2; RpAu_theory_min = f3; RpAu_theory_max = f4;
	     
		  cout << "RpAu Theory: " << RpAu_theory << endl;


		  if(counter == 1)
		    {
		      if(i == 0)
			{
			  uepps_68[1][0][0] = 0;
			  uepps_68[1][1][0] = RpAu_theory_min;
			  uepps_68[1][2][0] = RpAu_theory_max;
			  uepps_68[1][3][0] = 0;
			  uepps_68[1][4][0] = RpAu_theory;

			  //cout << uepps_68[1][0][0] << endl;
			}
		      if(i == 1)
			{
			  depps_68[1][0][0] = 0;
			  depps_68[1][1][0] = RpAu_theory_min;
			  depps_68[1][2][0] = RpAu_theory_max;
			  depps_68[1][3][0] = 0;
			  depps_68[1][4][0] = RpAu_theory;
		
			}
		      if(i == 2)
			{
			  cepps_68[1][0][0] = 0;
			  cepps_68[1][1][0] = RpAu_theory_min;
			  cepps_68[1][2][0] = RpAu_theory_max;
			  cepps_68[1][3][0] = 0;
			  cepps_68[1][4][0] = RpAu_theory;
			 
			  N_epps68_020_c = RpAu_theory;

			}
		      if(i == 3)
			{
			  uncteq_68[1][0][0] = 0;
			  uncteq_68[1][1][0] = RpAu_theory_min;
			  uncteq_68[1][2][0] = RpAu_theory_max;
			  uncteq_68[1][3][0] = 0;
			  uncteq_68[1][4][0] = RpAu_theory;
			}
		      if(i == 4)
			{
			  dncteq_68[1][0][0] = 0;
			  dncteq_68[1][1][0] = RpAu_theory_min;
			  dncteq_68[1][2][0] = RpAu_theory_max;
			  dncteq_68[1][3][0] = 0;
			  dncteq_68[1][4][0] = RpAu_theory;
			 
			}
		      if(i == 5)
			{
			  cncteq_68[1][0][0] = 0;
			  cncteq_68[1][1][0] = RpAu_theory_min;
			  cncteq_68[1][2][0] = RpAu_theory_max;
			  cncteq_68[1][3][0] = 0;
			  cncteq_68[1][4][0] = RpAu_theory;

			  N_ncteq68_020_c = RpAu_theory;

			}
		    } // counter == 0 is south arm for 0-20 cent

		  if(counter == 1)
		    {
		      if(i == 6)
			{
			  uepps_68[1][0][1] = 0;
			  uepps_68[1][1][1] = RpAu_theory_min;
			  uepps_68[1][2][1] = RpAu_theory_max;
			  uepps_68[1][3][1] = 0;
			  uepps_68[1][4][1] = RpAu_theory;
			}
		      if(i == 7)
			{
			  depps_68[1][0][1] = 0;
			  depps_68[1][1][1] = RpAu_theory_min;
			  depps_68[1][2][1] = RpAu_theory_max;
			  depps_68[1][3][1] = 0;
			  depps_68[1][4][1] = RpAu_theory;
			}
		      if(i == 8)
			{
			  cepps_68[1][0][1] = 0;
			  cepps_68[1][1][1] = RpAu_theory_min;
			  cepps_68[1][2][1] = RpAu_theory_max;
			  cepps_68[1][3][1] = 0;
			  cepps_68[1][4][1] = RpAu_theory;
			
			  N_epps68_2040_c = RpAu_theory;
			 
			}
		      if(i == 9)
			{
			  uncteq_68[1][0][1] = 0;
			  uncteq_68[1][1][1] = RpAu_theory_min;
			  uncteq_68[1][2][1] = RpAu_theory_max;
			  uncteq_68[1][3][1] = 0;
			  uncteq_68[1][4][1] = RpAu_theory;
			}
		      if(i == 10)
			{
			  dncteq_68[1][0][1] = 0;
			  dncteq_68[1][1][1] = RpAu_theory_min;
			  dncteq_68[1][2][1] = RpAu_theory_max;
			  dncteq_68[1][3][1] = 0;
			  dncteq_68[1][4][1] = RpAu_theory;
		
			}
		      if(i == 11)
			{
			  cncteq_68[1][0][1] = 0;
			  cncteq_68[1][1][1] = RpAu_theory_min;
			  cncteq_68[1][2][1] = RpAu_theory_max;
			  cncteq_68[1][3][1] = 0;
			  cncteq_68[1][4][1] = RpAu_theory;

			  N_ncteq68_2040_c = RpAu_theory;

			}
		    } // counter == 0 for cent 20-40

		  if(counter == 1)
		    {
		      if(i == 12)
			{
			  uepps_68[1][0][2] = 0;
			  uepps_68[1][1][2] = RpAu_theory_min;
			  uepps_68[1][2][2] = RpAu_theory_max;
			  uepps_68[1][3][2] = 0;
			  uepps_68[1][4][2] = RpAu_theory;
			}
		      if(i == 13)
			{
			  depps_68[1][0][2] = 0;
			  depps_68[1][1][2] = RpAu_theory_min;
			  depps_68[1][2][2] = RpAu_theory_max;
			  depps_68[1][3][2] = 0;
			  depps_68[1][4][2] = RpAu_theory;
			}
		      if(i == 14)
			{
			  cepps_68[1][0][2] = 0;
			  cepps_68[1][1][2] = RpAu_theory_min;
			  cepps_68[1][2][2] = RpAu_theory_max;
			  cepps_68[1][3][2] = 0;
			  cepps_68[1][4][2] = RpAu_theory;
			
			  N_epps68_4085_c = RpAu_theory;
			 
			}
		      if(i == 15)
			{
			  uncteq_68[1][0][2] = 0;
			  uncteq_68[1][1][2] = RpAu_theory_min;
			  uncteq_68[1][2][2] = RpAu_theory_max;
			  uncteq_68[1][3][2] = 0;
			  uncteq_68[1][4][2] = RpAu_theory;
			}
		      if(i == 16)
			{
			  dncteq_68[1][0][2] = 0;
			  dncteq_68[1][1][2] = RpAu_theory_min;
			  dncteq_68[1][2][2] = RpAu_theory_max;
			  dncteq_68[1][3][2] = 0;
			  dncteq_68[1][4][2] = RpAu_theory;
		
			}
		      if(i == 17)
			{
			  cncteq_68[1][0][2] = 0;
			  cncteq_68[1][1][2] = RpAu_theory_min;
			  cncteq_68[1][2][2] = RpAu_theory_max;
			  cncteq_68[1][3][2] = 0;
			  cncteq_68[1][4][2] = RpAu_theory;

			  N_ncteq68_4085_c = RpAu_theory;

			}
		    } // counter == 0 for cent 40-85

		  counter++;
		} while(counter <2);  //  two lines in each data file.  First line is south arm
	    }
	  else
	    cout << "file does not exist" << endl;     
	}
    }

  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 3; i++)  // 3 centrality bins
	{
	  value_c = cepps_68[arm][2][i];
	  value_d = depps_68[arm][2][i];
	  value_u = uepps_68[arm][2][i];
	      
	  // compare u and d to get max_value
	  if(value_d > value_u)
	    max_value = value_d;
	  if(value_u > value_d)
	    max_value = value_u;
	  if(value_d == value_u)
	    max_value = value_u;
	      
	  // compare max_value with c
	  if(value_c > max_value)
	    max_value = value_c;
	  if(max_value > value_c)
	    max_value = max_value;
	  if(value_c == max_value)
	    max_value = value_c;
	      
	  // fill max envope with max values of the three scales
	  epps16_env_68_max[arm][i] = max_value;
	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;


  // find min value [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 3; i++)  // 3 cent bins
	{
	  value_c = cepps_68[arm][1][i];
	  value_d = depps_68[arm][1][i];
	  value_u = uepps_68[arm][1][i];
	      
	  // compare u and d to get min_value
	  if(value_d < value_u)
	    min_value = value_d;
	  if(value_u < value_d)
	    min_value = value_u;
	  if(value_d == value_u)
	    min_value = value_u;
	      
	  // compare min_value with c
	  if(value_c < min_value)
	    min_value = value_c;
	  if(min_value < value_c)
	    min_value = min_value;
	  if(value_c == min_value)
	    min_value = value_c;
	      
	  // fill min envelope with min values of the three scales
	  epps16_env_68_min[arm][i] = min_value;
	   	      
	  epps16_env_68_up[arm][i] = epps16_env_68_max[arm][i] - cepps_68[arm][4][i];
	  epps16_env_68_down[arm][i] = cepps_68[arm][4][i] - epps16_env_68_min[arm][i];


	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;


  ///////////////////////////////////////////////
  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 3; i++) // 3 cent bins
	{
	  value_c = cncteq_68[arm][2][i];
	  value_d = dncteq_68[arm][2][i];
	  value_u = uncteq_68[arm][2][i];
	      
	  // compare u and d to get max_value
	  if(value_d > value_u)
	    max_value = value_d;
	  if(value_u > value_d)
	    max_value = value_u;
	  if(value_d == value_u)
	    max_value = value_u;
	      
	  // compare max_value with c
	  if(value_c > max_value)
	    max_value = value_c;
	  if(max_value > value_c)
	    max_value = max_value;
	  if(value_c == max_value)
	    max_value = value_c;
	      
	  // fill max envope with max values of the three scales
	  ncteq15_env_68_max[arm][i] = max_value;
	      
	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 3; i++)  // 3 cent bins
	{
	  value_c = cncteq_68[arm][1][i];
	  value_d = dncteq_68[arm][1][i];
	  value_u = uncteq_68[arm][1][i];
	      
	  // compare u and d to get min_value
	  if(value_d < value_u)
	    min_value = value_d;
	  if(value_u < value_d)
	    min_value = value_u;
	  if(value_d == value_u)
	    min_value = value_u;
	      
	  // compare min_value with c
	  if(value_c < min_value)
	    min_value = value_c;
	  if(min_value < value_c)
	    min_value = min_value;
	  if(value_c == min_value)
	    min_value = value_c;
	      
	  // fill min envope with min values of the three scales
	  ncteq15_env_68_min[arm][i] = min_value;
		      
	  ncteq15_env_68_up[arm][i] = ncteq15_env_68_max[arm][i] -  cncteq_68[arm][4][i];   // using central value energy scale, as in PPG228
	  ncteq15_env_68_down[arm][i] =  cncteq_68[arm][4][i] - ncteq15_env_68_min[arm][i];     

	}
    }

 
  bool draw_uscale = false;
  bool draw_dscale = false;
  bool draw_cscale = true;
  bool draw_all = true;
 
  char name[500];
  char name2[500];

  double x_errors[6] = {0};
  double y_errors[6] = {0};
 

  double ncoll_array_pAu[3] =  {8.2, 6.1, 3.4};

  TGraphAsymmErrors *gr_ncteq_N;
  TGraphAsymmErrors *gr_epps_N;
  TGraphAsymmErrors *gr_ncteq_S;
  TGraphAsymmErrors *gr_epps_S;


  //  gr_epps_S = new TGraphAsymmErrors(3,ncoll_array_pAu, cent_epps, x_errors, x_errors,cent_epps_down, cent_epps_up);
  // gr_ncteq_S = new TGraphAsymmErrors(3,ncoll_array_pAu, cent_ncteq, x_errors, x_errors, cent_ncteq_down, cent_ncteq_up);
  

 gr_epps_N = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cepps_68[1][4], x_errors, x_errors,epps16_env_68_down[1], epps16_env_68_up[1]);
 gr_ncteq_N = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cncteq_68[1][4], x_errors, x_errors,ncteq15_env_68_down[1], ncteq15_env_68_up[1]);
  
  //gr_epps_S = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cepps_68[0][4], x_errors, x_errors,epps16_env_68_down[0], epps16_env_68_up[0]);
  // gr_ncteq_S = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cncteq_68[0][4], x_errors, x_errors,ncteq15_env_68_down[0], ncteq15_env_68_up[0]);
  

  TCanvas *c1 = new TCanvas("c1","RpAu cent Predictions",200,10,700,500);
  c1->SetGrid();
  c1->Divide(2,2);
  c1->cd(1);  

  gr_epps_N->SetTitle("RpAu_epps_cent_predictions fwd");
  gr_epps_N->GetXaxis()->SetTitle("NColl");
  gr_epps_N->SetMarkerStyle(20);
  gr_epps_N->SetMarkerColor(kBlue+2);
  gr_epps_N->SetMarkerSize(1);
  gr_epps_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_N->GetXaxis()->SetLimits(0, 10);
  
  gr_epps_N->SetMarkerStyle(20);                                             
  gr_epps_N->Draw("AP");


  TLegend *leg = new TLegend(0.1, 0.7, 0.7, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.03); 
  leg->AddEntry(gr_epps_N, "EPPS16 Reweighted 68CL, combined error scale", "p");
  leg->Draw();

  c1->cd(3);  

  gr_ncteq_N->SetTitle("RpAu_ncteq_cent_predictions fwd");
  gr_ncteq_N->GetXaxis()->SetTitle("Ncoll");
  gr_ncteq_N->SetMarkerStyle(20);
  gr_ncteq_N->SetMarkerColor(kCyan+2);
  gr_ncteq_N->SetMarkerSize(1);
  gr_ncteq_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_N->GetXaxis()->SetLimits(0, 10);
  
  gr_ncteq_N->SetMarkerStyle(20);                                             
  gr_ncteq_N->Draw("AP");  

   
  TLegend *leg3 = new TLegend(0.1, 0.7, 0.7, 0.9);  //(start x, start y, end x, end y)
  leg3->SetFillColor(0); 
  leg3->SetTextSize(0.03); 

  leg3->AddEntry(gr_ncteq_N, "c-scale Reweighted 68CL, combined error scale", "p");
  leg3->Draw(); 

 
  //// write mulitple TH2D's to a root file
  TFile *h = new TFile("psi2s_RpAu_ncoll_shao_N.root", "RECREATE");
  h->cd();

  
  gr_ncteq_N->SetName("RpAu_ncteq_shao_ncoll_fwd");
  gr_epps_N->SetName("RpAu_epps_shao_ncoll_fwd");


  gr_ncteq_N->Write();
  gr_epps_N->Write();
  
}



  



