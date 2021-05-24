
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

void Write_jpsi_pAu_MB_predictions_rap()
{

   bool north_arm = false;
  
  TFile *Run15pAu_prediction_files;

  const int bins_pAu = 1;
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

 double N_epps68_c = 0;
 double S_epps68_c = 0;
 double N_ncteq68_c = 0;
 double S_ncteq68_c = 0;
 
 //  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < 6; i++)
	{
	  //	  if(arm == 0)
	    {
	     if(i == 0)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.MB.dat"); // u_ncteq15_68_S;
	      if(i == 1)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.MB.dat"); // 
	      if(i == 2)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.MB.dat"); // u_epps_68_S;
	      if(i == 3)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.MB.dat"); // u_epps_68_S;
	      if(i == 4)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.MB.dat"); // d_epps_68_S;
	      if(i == 5)
		sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.MB.dat"); // c_epps_68_S;
	    }
	  // if(arm == 1)
	  //   {
	  //     if(i == 0)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.MB.dat"); // u_ncteq15_68_S;
	  //     if(i == 1)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.MB.dat"); // 
	  //     if(i == 2)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.MB.dat"); // u_epps_68_S;
	  //     if(i == 3)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.MB.dat"); // u_epps_68_S;
	  //     if(i == 4)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.MB.dat"); // d_epps_68_S;
	  //     if(i == 5)
	  // 	sprintf(name100,"tables_RpAu/jpsi_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.MB.dat"); // c_epps_68_S;
	  //   }
      
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


		  if(counter == 0)
		    {
		      if(i == 0)
			{
			  uepps_68[0][0][counter] = 0;
			  uepps_68[0][1][counter] = RpAu_theory_min;
			  uepps_68[0][2][counter] = RpAu_theory_max;
			  uepps_68[0][3][counter] = 0;
			  uepps_68[0][4][counter] = RpAu_theory;

			  //cout << uepps_68[0][0][counter] << endl;
			}
		      if(i == 1)
			{
			  depps_68[0][0][counter] = 0;
			  depps_68[0][1][counter] = RpAu_theory_min;
			  depps_68[0][2][counter] = RpAu_theory_max;
			  depps_68[0][3][counter] = 0;
			  depps_68[0][4][counter] = RpAu_theory;
		
			}
		      if(i == 2)
			{
			  cepps_68[0][0][counter] = 0;
			  cepps_68[0][1][counter] = RpAu_theory_min;
			  cepps_68[0][2][counter] = RpAu_theory_max;
			  cepps_68[0][3][counter] = 0;
			  cepps_68[0][4][counter] = RpAu_theory;
			 
			  S_epps68_c = RpAu_theory;

			}
		      if(i == 3)
			{
			  uncteq_68[0][0][counter] = 0;
			  uncteq_68[0][1][counter] = RpAu_theory_min;
			  uncteq_68[0][2][counter] = RpAu_theory_max;
			  uncteq_68[0][3][counter] = 0;
			  uncteq_68[0][4][counter] = RpAu_theory;
			}
		      if(i == 4)
			{
			  dncteq_68[0][0][counter] = 0;
			  dncteq_68[0][1][counter] = RpAu_theory_min;
			  dncteq_68[0][2][counter] = RpAu_theory_max;
			  dncteq_68[0][3][counter] = 0;
			  dncteq_68[0][4][counter] = RpAu_theory;
			 
			}
		      if(i == 5)
			{
			  cncteq_68[0][0][counter] = 0;
			  cncteq_68[0][1][counter] = RpAu_theory_min;
			  cncteq_68[0][2][counter] = RpAu_theory_max;
			  cncteq_68[0][3][counter] = 0;
			  cncteq_68[0][4][counter] = RpAu_theory;

			  S_ncteq68_c = RpAu_theory;

			}
		    } // counter == 0 is south arm

		  if(counter == 1)
		    {
		      if(i == 0)
			{
			  uepps_68[1][0][0] = 0;
			  uepps_68[1][1][0] = RpAu_theory_min;
			  uepps_68[1][2][0] = RpAu_theory_max;
			  uepps_68[1][3][0] = 0;
			  uepps_68[1][4][0] = RpAu_theory;
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
			
			  N_epps68_c = RpAu_theory;
			 
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

			  N_ncteq68_c = RpAu_theory;

			}
		    } // counter == 1 is north arm

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
      for(int i =  0; i < 1; i++)
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

  double rap_epps[2][3] = {
    {S_epps68_c, S_epps68_c, S_epps68_c},
    {N_epps68_c, N_epps68_c, N_epps68_c}};
  
  double rap_epps_up_N = 0;
  double rap_epps_up_S = 0;
  double rap_epps_down_N = 0;
  double rap_epps_down_S = 0;

double rap_ncteq[2][3] = {
    {S_ncteq68_c, S_ncteq68_c, S_ncteq68_c},
    {N_ncteq68_c, N_ncteq68_c, N_ncteq68_c}};

  double rap_ncteq_up_N = 0;
  double rap_ncteq_up_S = 0;
  double rap_ncteq_down_N = 0;
  double rap_ncteq_down_S = 0;


  // find min value [1]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 1; i++)
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

	  rap_epps_up_S = epps16_env_68_up[0][0];
	  rap_epps_up_N =  epps16_env_68_up[1][0];
	  rap_epps_down_S =  epps16_env_68_down[0][0] ;
	  rap_epps_down_N =  epps16_env_68_down[1][0] ;

	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;


  ///////////////////////////////////////////////
  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 1; i++)
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
      for(int i =  0; i < 1; i++)
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
	
	  rap_ncteq_up_S = ncteq15_env_68_up[0][0];
	  rap_ncteq_up_N = ncteq15_env_68_up[1][0];
	  rap_ncteq_down_S = ncteq15_env_68_down[0][0] ;
	  rap_ncteq_down_N = ncteq15_env_68_down[1][0] ;

	}
    }

  double rap_ncteq_up[2][3] = {
    {rap_ncteq_up_S, rap_ncteq_up_S, rap_ncteq_up_S},
    {rap_ncteq_up_N, rap_ncteq_up_N, rap_ncteq_up_N}};

double rap_ncteq_down[2][3] = {
  {rap_ncteq_down_S, rap_ncteq_down_S, rap_ncteq_down_S},
  {rap_ncteq_down_N, rap_ncteq_down_N, rap_ncteq_down_N}};
 
double rap_epps_up[2][3] = {
  {rap_epps_up_S, rap_epps_up_S, rap_epps_up_S},
  {rap_epps_up_N, rap_epps_up_N, rap_epps_up_N}};

double rap_epps_down[2][3] = {
  {rap_epps_down_S, rap_epps_down_S, rap_epps_down_S},
  {rap_epps_down_N, rap_epps_down_N, rap_epps_down_N}};
  
  bool draw_uscale = false;
  bool draw_dscale = false;
  bool draw_cscale = true;
  bool draw_all = true;
 
  char name[500];
  char name2[500];

  double x_errors[6] = {0};
  double y_errors[6] = {0};
 
  
  double rap_array_pAu[2][3] =  {
    {-2.0, -1.7, -1.4},
    {1.4, 1.7, 2.0}};

  TGraphAsymmErrors *gr_ncteq_N;
  TGraphAsymmErrors *gr_epps_N;
  TGraphAsymmErrors *gr_ncteq_S;
  TGraphAsymmErrors *gr_epps_S;

  gr_epps_S = new TGraphAsymmErrors(3,rap_array_pAu[0], rap_epps[0], x_errors, x_errors,rap_epps_down[0], rap_epps_up[0]);
  gr_ncteq_S = new TGraphAsymmErrors(3,rap_array_pAu[0], rap_ncteq[0], x_errors, x_errors, rap_ncteq_down[0], rap_ncteq_up[0]);
  

  gr_epps_N = new TGraphAsymmErrors(3,rap_array_pAu[1], rap_epps[1], x_errors, x_errors,rap_epps_down[1], rap_epps_up[1]);
  gr_ncteq_N = new TGraphAsymmErrors(3,rap_array_pAu[1], rap_ncteq[1], x_errors, x_errors, rap_ncteq_down[1], rap_ncteq_up[1]);
  

  TCanvas *c1 = new TCanvas("c1","RpAu MB Predictions",200,10,700,500);
  c1->SetGrid();
  c1->Divide(2,2);
  c1->cd(1);  

  gr_epps_S->SetTitle("RpAu_epps_MB_predictions bkwd");
  gr_epps_S->GetXaxis()->SetTitle("y");
  gr_epps_S->SetMarkerStyle(20);
  gr_epps_S->SetMarkerColor(kBlue+2);
  gr_epps_S->SetMarkerSize(1);
  gr_epps_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_S->GetXaxis()->SetLimits(-3,3);
  
  gr_epps_S->SetMarkerStyle(20);                                             
  gr_epps_S->Draw("AP");

  gr_epps_N->SetTitle("RpAu_epps_MB_predictions fwd");
  gr_epps_N->GetXaxis()->SetTitle("y");
  gr_epps_N->SetMarkerStyle(20);
  gr_epps_N->SetMarkerColor(kBlue+2);
  gr_epps_N->SetMarkerSize(1);
  gr_epps_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_N->GetXaxis()->SetLimits(-3,3);
  
  gr_epps_N->Draw("P");
   
  TLegend *leg = new TLegend(0.1, 0.7, 0.7, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.03); 
  leg->AddEntry(gr_epps_N, "EPPS16 Reweighted 68CL, combined error scale", "p");
  leg->Draw();

  c1->cd(3);  

  gr_ncteq_S->SetTitle("RpAu_ncteq_MB_predictions bkwd");
  gr_ncteq_S->GetXaxis()->SetTitle("y");
  gr_ncteq_S->SetMarkerStyle(20);
  gr_ncteq_S->SetMarkerColor(kCyan+2);
  gr_ncteq_S->SetMarkerSize(1);
  gr_ncteq_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_S->GetXaxis()->SetLimits(-3,3);
  
  gr_ncteq_S->SetMarkerStyle(20);                                             
  gr_ncteq_S->Draw("AP");  

  gr_ncteq_N->SetTitle("RpAu_ncteq_MB_predictions fwd");
  gr_ncteq_N->GetXaxis()->SetTitle("y");
  gr_ncteq_N->SetMarkerStyle(20);
  gr_ncteq_N->SetMarkerColor(kCyan+2);
  gr_ncteq_N->SetMarkerSize(1);
  gr_ncteq_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_N->GetXaxis()->SetLimits(-3,3);                                    
  gr_ncteq_N->Draw("P");  


   
  TLegend *leg3 = new TLegend(0.1, 0.7, 0.7, 0.9);  //(start x, start y, end x, end y)
  leg3->SetFillColor(0); 
  leg3->SetTextSize(0.03); 

  leg3->AddEntry(gr_ncteq_N, "c-scale Reweighted 68CL, combined error scale", "p");
  leg3->Draw(); 

 
  //// write mulitple TH2D's to a root file
  TFile *h = new TFile("jpsi_RpAu_MB_predictions_rap.root", "RECREATE");
  h->cd();

  gr_ncteq_N->SetName("RpAu_ncteq_MB_predictions_fwd");
  gr_epps_N->SetName("RpAu_epps_MB_predictions_fwd");
  
  gr_ncteq_S->SetName("RpAu_ncteq_MB_predictions_bkwd");
  gr_epps_S->SetName("RpAu_epps_MB_predictions_bkwd");

  gr_ncteq_N->Write();
  gr_epps_N->Write();

  gr_ncteq_S->Write();
  gr_epps_S->Write();
  
}



  


