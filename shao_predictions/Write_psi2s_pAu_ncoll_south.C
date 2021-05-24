
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

void Write_psi2s_pAu_ncoll_south()
{

   bool north_arm = false;
  
  TFile *Run15pAu_prediction_files;

  const int bins_pAu = 2;  // cent bins
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

 double S_epps68_020_c = 0;
 double S_epps68_2085_c = 0;
 double S_ncteq68_020_c = 0;
 double S_ncteq68_2085_c = 0;
 
 //  for(int arm = 0; arm < 2; arm++)
 {
      for(int i = 0; i < 12; i++)
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
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFu.Rwgt.Cen20To85.a1.dat"); // u_ncteq15_68_S;
	      if(i == 7)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFd.Rwgt.Cen20To85.a1.dat"); // 
	      if(i == 8)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_EPPS16.68CL.MUFc.Rwgt.Cen20To85.a1.dat"); // u_epps_68_S;
	      if(i == 9)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFu.Rwgt.Cen20To85.a1.dat"); // u_epps_68_S;
	      if(i == 10)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFd.Rwgt.Cen20To85.a1.dat"); // d_epps_68_S;
	      if(i == 11)
		sprintf(name100,"tables_RpAu/psi2S_yCMS_PHENIX_nCTEQ15.68CL.MUFc.Rwgt.Cen20To85.a1.dat"); // c_epps_68_S;
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
			 
			  S_epps68_020_c = RpAu_theory;

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

			  S_ncteq68_020_c = RpAu_theory;

			}
		    } // counter == 0 is south arm for 0-20 cent

		  if(counter == 0)
		    {
		      if(i == 6)
			{
			  uepps_68[0][0][1] = 0;
			  uepps_68[0][1][1] = RpAu_theory_min;
			  uepps_68[0][2][1] = RpAu_theory_max;
			  uepps_68[0][3][1] = 0;
			  uepps_68[0][4][1] = RpAu_theory;
			}
		      if(i == 7)
			{
			  depps_68[0][0][1] = 0;
			  depps_68[0][1][1] = RpAu_theory_min;
			  depps_68[0][2][1] = RpAu_theory_max;
			  depps_68[0][3][1] = 0;
			  depps_68[0][4][1] = RpAu_theory;
			}
		      if(i == 8)
			{
			  cepps_68[0][0][1] = 0;
			  cepps_68[0][1][1] = RpAu_theory_min;
			  cepps_68[0][2][1] = RpAu_theory_max;
			  cepps_68[0][3][1] = 0;
			  cepps_68[0][4][1] = RpAu_theory;
			
			  S_epps68_2085_c = RpAu_theory;
			 
			}
		      if(i == 9)
			{
			  uncteq_68[0][0][1] = 0;
			  uncteq_68[0][1][1] = RpAu_theory_min;
			  uncteq_68[0][2][1] = RpAu_theory_max;
			  uncteq_68[0][3][1] = 0;
			  uncteq_68[0][4][1] = RpAu_theory;
			}
		      if(i == 10)
			{
			  dncteq_68[0][0][1] = 0;
			  dncteq_68[0][1][1] = RpAu_theory_min;
			  dncteq_68[0][2][1] = RpAu_theory_max;
			  dncteq_68[0][3][1] = 0;
			  dncteq_68[0][4][1] = RpAu_theory;
		
			}
		      if(i == 11)
			{
			  cncteq_68[0][0][1] = 0;
			  cncteq_68[0][1][1] = RpAu_theory_min;
			  cncteq_68[0][2][1] = RpAu_theory_max;
			  cncteq_68[0][3][1] = 0;
			  cncteq_68[0][4][1] = RpAu_theory;

			  S_ncteq68_2085_c = RpAu_theory;

			}
		    } // counter == 0 for cent 20-85

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
      for(int i =  0; i < 2; i++)  // 2 centrality bins
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

  double cent_epps[2] = {S_epps68_020_c, S_epps68_2085_c};
  
  double cent_epps_up_N = 0;
  double cent_epps_up_S = 0;
  double cent_epps_down_N = 0;
  double cent_epps_down_S = 0;

double cent_ncteq[2] = {S_ncteq68_020_c, S_ncteq68_2085_c};

  double cent_ncteq_up_N = 0;
  double cent_ncteq_up_S = 0;
  double cent_ncteq_down_N = 0;
  double cent_ncteq_down_S = 0;


  // find min value [1]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 2; i++)  // 2 cent bins
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

	  cent_epps_up_S = epps16_env_68_up[0][0];
	  cent_epps_up_N =  epps16_env_68_up[1][0];
	  cent_epps_down_S =  epps16_env_68_down[0][0] ;
	  cent_epps_down_N =  epps16_env_68_down[1][0] ;

	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;


  ///////////////////////////////////////////////
  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 2; i++) // 2 cent bins
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
      for(int i =  0; i < 2; i++)  // 2 cent bins
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
	
	  cent_ncteq_up_S = ncteq15_env_68_up[0][0];
	  cent_ncteq_up_N = ncteq15_env_68_up[1][0];
	  cent_ncteq_down_S = ncteq15_env_68_down[0][0] ;
	  cent_ncteq_down_N = ncteq15_env_68_down[1][0] ;

	}
    }

  double cent_ncteq_up[2] = {cent_ncteq_up_S, cent_ncteq_up_S};

double cent_ncteq_down[2] = {cent_ncteq_down_S, cent_ncteq_down_S};
 
double cent_epps_up[2] = {cent_epps_up_S, cent_epps_up_S};

double cent_epps_down[2] = {cent_epps_down_S, cent_epps_down_S};
  
  bool draw_uscale = false;
  bool draw_dscale = false;
  bool draw_cscale = true;
  bool draw_all = true;
 
  char name[500];
  char name2[500];

  double x_errors[6] = {0};
  double y_errors[6] = {0};
 

  double ncoll_array_pAu[2] =  {8.2, 4.3};

  TGraphAsymmErrors *gr_ncteq_N;
  TGraphAsymmErrors *gr_epps_N;
  TGraphAsymmErrors *gr_ncteq_S;
  TGraphAsymmErrors *gr_epps_S;


  // gr_epps_S = new TGraphAsymmErrors(2,ncoll_array_pAu, cent_epps, x_errors, x_errors,cent_epps_down, cent_epps_up);
  // gr_ncteq_S = new TGraphAsymmErrors(2,ncoll_array_pAu, cent_ncteq, x_errors, x_errors, cent_ncteq_down, cent_ncteq_up);
  

  gr_epps_S = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cepps_68[0][4], x_errors, x_errors,epps16_env_68_down[0], epps16_env_68_up[0]);
  gr_ncteq_S = new TGraphAsymmErrors(bins_pAu,ncoll_array_pAu, cncteq_68[0][4], x_errors, x_errors,ncteq15_env_68_down[0], ncteq15_env_68_up[0]);
  

  TCanvas *c1 = new TCanvas("c1","RpAu cent Predictions",200,10,700,500);
  c1->SetGrid();
  c1->Divide(2,2);
  c1->cd(1);  

  gr_epps_S->SetTitle("RpAu_epps_cent_predictions bkwd");
  gr_epps_S->GetXaxis()->SetTitle("NColl");
  gr_epps_S->SetMarkerStyle(20);
  gr_epps_S->SetMarkerColor(kBlue+2);
  gr_epps_S->SetMarkerSize(1);
  gr_epps_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_S->GetXaxis()->SetLimits(0, 10);
  
  gr_epps_S->SetMarkerStyle(20);                                             
  gr_epps_S->Draw("AP");


  TLegend *leg = new TLegend(0.1, 0.7, 0.7, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.03); 
  leg->AddEntry(gr_epps_S, "EPPS16 Reweighted 68CL, combined error scale", "p");
  leg->Draw();

  c1->cd(3);  

  gr_ncteq_S->SetTitle("RpAu_ncteq_cent_predictions bkwd");
  gr_ncteq_S->GetXaxis()->SetTitle("Ncoll");
  gr_ncteq_S->SetMarkerStyle(20);
  gr_ncteq_S->SetMarkerColor(kCyan+2);
  gr_ncteq_S->SetMarkerSize(1);
  gr_ncteq_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_S->GetXaxis()->SetLimits(0, 10);
  
  gr_ncteq_S->SetMarkerStyle(20);                                             
  gr_ncteq_S->Draw("AP");  

   
  TLegend *leg3 = new TLegend(0.1, 0.7, 0.7, 0.9);  //(start x, start y, end x, end y)
  leg3->SetFillColor(0); 
  leg3->SetTextSize(0.03); 

  leg3->AddEntry(gr_ncteq_S, "c-scale Reweighted 68CL, combined error scale", "p");
  leg3->Draw(); 

 
  //// write mulitple TH2D's to a root file
  TFile *h = new TFile("psi2s_RpAu_ncoll_shao_S.root", "RECREATE");
  h->cd();

  
  gr_ncteq_S->SetName("RpAu_ncteq_shao_ncoll_bkwd");
  gr_epps_S->SetName("RpAu_epps_shao_ncoll_bkwd");


  gr_ncteq_S->Write();
  gr_epps_S->Write();
  
}



  



