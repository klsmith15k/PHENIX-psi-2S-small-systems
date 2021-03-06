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

void Write_pAu_MB_predictions_pT()
{

   bool north_arm = false;
  
  TFile *Run15pAu_prediction_files;

  const int bins_pAu = 5;
  const int columns = 5;

  char name100[500];

  const int narm = 2;

  double uepps_90[narm][columns][bins_pAu] = {0}; 
  double depps_90[narm][columns][bins_pAu] = {0};
  double cepps_90[narm][columns][bins_pAu] = {0};
	
  double epps16_env_90_min[narm][bins_pAu] = {0};
  double epps16_env_90_max[narm][bins_pAu] = {0};
  double epps16_env_90_RpAu_pred[narm][bins_pAu] = {0};

  double epps16_env_90_down[narm][bins_pAu] = {0};
  double epps16_env_90_up[narm][bins_pAu] = {0};
  
  double uncteq_90[narm][columns][bins_pAu] = {0}; 
  double dncteq_90[narm][columns][bins_pAu] = {0};
  double cncteq_90[narm][columns][bins_pAu] = {0};

  double ncteq15_env_90_min[narm][bins_pAu] = {0};
  double ncteq15_env_90_max[narm][bins_pAu] = {0};
  double ncteq15_env_90_RpAu_pred[narm][bins_pAu] = {0};

   double ncteq15_env_90_down[narm][bins_pAu] = {0};
  double ncteq15_env_90_up[narm][bins_pAu] = {0};

  double value_d  = 0;
  double value_c  = 0;
  double value_u = 0;

  double min_value = 0;
  double max_value = 0;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < 6; i++)
	{
	  if(arm == 0)
	    {
	      if(i == 0)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_EPPS16.90CL.MUFu.Rwgt.MB.dat"); // u_epps_90_S;
	      if(i == 1)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_EPPS16.90CL.MUFd.Rwgt.MB.dat"); // d_epps_90_S;
	      if(i == 2)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_EPPS16.90CL.MUFc.Rwgt.MB.dat"); // c_epps_90_S;
	      if(i == 3)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_nCTEQ15.90CL.MUFu.Rwgt.MB.dat"); // u_epps_90_N;
	      if(i == 4)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_nCTEQ15.90CL.MUFd.Rwgt.MB.dat"); // d_epps_90_N;
	      if(i == 5)
		sprintf(name100,"tables_RpAu/psi2S_PT_bwd_PHENIX_nCTEQ15.90CL.MUFc.Rwgt.MB.dat"); // c_epps_90_N;
	    }
	  if(arm == 1)
	    {
	      if(i == 0)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_EPPS16.90CL.MUFu.Rwgt.MB.dat"); // u_epps_90_S;
	      if(i == 1)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_EPPS16.90CL.MUFd.Rwgt.MB.dat"); // d_epps_90_S;
	      if(i == 2)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_EPPS16.90CL.MUFc.Rwgt.MB.dat"); // c_epps_90_S;
	      if(i == 3)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_nCTEQ15.90CL.MUFu.Rwgt.MB.dat"); //  u_epps_90_N;
	      if(i == 4)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_nCTEQ15.90CL.MUFd.Rwgt.MB.dat"); // d_epps_90_N;
	      if(i == 5)
		sprintf(name100,"tables_RpAu/psi2S_PT_fwd_PHENIX_nCTEQ15.90CL.MUFc.Rwgt.MB.dat"); // c_epps_90_N;
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

		  if(i == 0)
		    {
		      uepps_90[arm][0][counter] = 0;
		      uepps_90[arm][1][counter] = RpAu_theory_min;
		      uepps_90[arm][2][counter] = RpAu_theory_max;
		      uepps_90[arm][3][counter] = 0;
		      uepps_90[arm][4][counter] = RpAu_theory;

		      //cout << uepps_90[arm][0][counter] << endl;
		    }
		  if(i == 1)
		    {
		      depps_90[arm][0][counter] = 0;
		      depps_90[arm][1][counter] = RpAu_theory_min;
		      depps_90[arm][2][counter] = RpAu_theory_max;
		      depps_90[arm][3][counter] = 0;
		      depps_90[arm][4][counter] = RpAu_theory;
		    }
		  if(i == 2)
		    {
		      cepps_90[arm][0][counter] = 0;
		      cepps_90[arm][1][counter] = RpAu_theory_min;
		      cepps_90[arm][2][counter] = RpAu_theory_max;
		      cepps_90[arm][3][counter] = 0;
		      cepps_90[arm][4][counter] = RpAu_theory;
		    }
		  if(i == 3)
		    {
		      uncteq_90[arm][0][counter] = 0;
		      uncteq_90[arm][1][counter] = RpAu_theory_min;
		      uncteq_90[arm][2][counter] = RpAu_theory_max;
		      uncteq_90[arm][3][counter] = 0;
		      uncteq_90[arm][4][counter] = RpAu_theory;
		    }
		  if(i == 4)
		    {
		      dncteq_90[arm][0][counter] = 0;
		      dncteq_90[arm][1][counter] = RpAu_theory_min;
		      dncteq_90[arm][2][counter] = RpAu_theory_max;
		      dncteq_90[arm][3][counter] = 0;
		      dncteq_90[arm][4][counter] = RpAu_theory;
		    }
		  if(i == 5)
		    {
		      cncteq_90[arm][0][counter] = 0;
		      cncteq_90[arm][1][counter] = RpAu_theory_min;
		      cncteq_90[arm][2][counter] = RpAu_theory_max;
		      cncteq_90[arm][3][counter] = 0;
		      cncteq_90[arm][4][counter] = RpAu_theory;
		    }

		  counter++;
		} while(counter <5);  // lines in each .dat file (omitting 4.5-7 data point)
	    }
	  else
	    cout << "file does not exist" << endl;     
	}
    }

  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 5; i++)
	{
	  value_c = cepps_90[arm][2][i];
	  value_d = depps_90[arm][2][i];
	  value_u = uepps_90[arm][2][i];
	      
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
	  epps16_env_90_max[arm][i] = max_value;
	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;

  // find min value [1]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 5; i++)
	{
	  value_c = cepps_90[arm][1][i];
	  value_d = depps_90[arm][1][i];
	  value_u = uepps_90[arm][1][i];
	      
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
	  epps16_env_90_min[arm][i] = min_value;
	   	      
	  epps16_env_90_up[arm][i] = epps16_env_90_max[arm][i] - cepps_90[arm][4][i];
	  epps16_env_90_down[arm][i] = cepps_90[arm][4][i] - epps16_env_90_min[arm][i];
	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;

  ///////////////////////////////////////////////
  // find max value  [2]
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 5; i++)
	{
	  value_c = cncteq_90[arm][2][i];
	  value_d = dncteq_90[arm][2][i];
	  value_u = uncteq_90[arm][2][i];
	      
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
	  ncteq15_env_90_max[arm][i] = max_value;
	      
	}
    }

  value_c = 0;
  value_d = 0;
  value_u = 0;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i =  0; i < 5; i++)
	{
	  value_c = cncteq_90[arm][1][i];
	  value_d = dncteq_90[arm][1][i];
	  value_u = uncteq_90[arm][1][i];
	      
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
	  ncteq15_env_90_min[arm][i] = min_value;
		      
	  ncteq15_env_90_up[arm][i] = ncteq15_env_90_max[arm][i] -  cncteq_90[arm][4][i];
	  ncteq15_env_90_down[arm][i] =  cncteq_90[arm][4][i] - ncteq15_env_90_min[arm][i];
	      
	}
    }
  bool draw_uscale = false;
  bool draw_dscale = true;
  bool draw_cscale = false;
  bool draw_all = true;
 
  char name[500];
  char name2[500];

  double x_errors[5] = {0};
  double y_errors[5] = {0};
 
  double pt_array_pAu[5] =  {0.25,1.0, 2.0, 3.5, 5.75};
  //double pt_array_pAu[5] =  {0.25,1.0, 2.0, 3.5};

  TGraphAsymmErrors *gr_ncteq_N;
  TGraphAsymmErrors *gr_epps_N;

  TGraphAsymmErrors *gr_ncteq_S;
  TGraphAsymmErrors *gr_epps_S;

  // gr_epps_N = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, cepps_90[1][4], x_errors, x_errors,epps16_env_90_down[1], epps16_env_90_up[1]);
  // gr_ncteq_N = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, cncteq_90[1][4], x_errors, x_errors,ncteq15_env_90_down[1], ncteq15_env_90_up[1]);
  
  // gr_epps_S = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, cepps_90[0][4], x_errors, x_errors,epps16_env_90_down[0], epps16_env_90_up[0]);
  // gr_ncteq_S = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, cncteq_90[0][4], x_errors, x_errors,ncteq15_env_90_down[0], ncteq15_env_90_up[0]);
  
  gr_epps_N = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, depps_90[1][4], x_errors, x_errors,epps16_env_90_down[1], epps16_env_90_up[1]);
  gr_ncteq_N = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, dncteq_90[1][4], x_errors, x_errors,ncteq15_env_90_down[1], ncteq15_env_90_up[1]);
  
  gr_epps_S = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, depps_90[0][4], x_errors, x_errors,epps16_env_90_down[0], epps16_env_90_up[0]);
  gr_ncteq_S = new TGraphAsymmErrors(bins_pAu,pt_array_pAu, dncteq_90[0][4], x_errors, x_errors,ncteq15_env_90_down[0], ncteq15_env_90_up[0]);



  TCanvas *c1 = new TCanvas("c1","RpAu MB Predictions",200,10,700,500);
  c1->SetGrid();
  c1->Divide(2,2);
  c1->cd(1);  

  gr_epps_S->SetTitle("RpAu_epps_MB_predictions_bkwd");
  gr_epps_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
  gr_epps_S->SetMarkerStyle(20);
  gr_epps_S->SetMarkerColor(kBlue+2);
  gr_epps_S->SetMarkerSize(1);
  gr_epps_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_S->GetXaxis()->SetLimits(0,7);
  
  gr_epps_S->SetMarkerStyle(20);                                             
  gr_epps_S->Draw("AP");
   
  TLegend *leg = new TLegend(0.1, 0.7, 0.7, 0.9);  
  leg->SetFillColor(0); 
  leg->SetTextSize(0.03); 
  leg->AddEntry(gr_epps_S, "EPPS16 Reweighted 90CL, combined error scale", "p");
  leg->Draw();
  
  c1->cd(2);  

  gr_epps_N->SetTitle("RpAu_epps_MB_predictions_fwd");
  gr_epps_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
  gr_epps_N->SetMarkerStyle(20);
  gr_epps_N->SetMarkerColor(kRed+2);
  gr_epps_N->SetMarkerSize(1);
  gr_epps_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_epps_N->GetXaxis()->SetLimits(0,7);
  
  gr_epps_N->SetMarkerStyle(20);                                             
  gr_epps_N->Draw("AP");
   
  TLegend *leg2 = new TLegend(0.1, 0.7, 0.7, 0.9);  
  leg2->SetFillColor(0); 
  leg2->SetTextSize(0.03); 
  leg2->AddEntry(gr_epps_N, "EPPS16 Reweighted 90CL, combined error scale", "p");
  leg2->Draw();

  c1->cd(3);  

  gr_ncteq_S->SetTitle("RpAu_ncteq_MB_predictions_bkwd");
  gr_ncteq_S->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
  gr_ncteq_S->SetMarkerStyle(20);
  gr_ncteq_S->SetMarkerColor(kCyan+2);
  gr_ncteq_S->SetMarkerSize(1);
  gr_ncteq_S->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_S->GetXaxis()->SetLimits(0,7);
  
  gr_ncteq_S->SetMarkerStyle(20);                                             
  gr_ncteq_S->Draw("AP");  
   
  TLegend *leg3 = new TLegend(0.1, 0.7, 0.7, 0.9);  //(start x, start y, end x, end y)
  leg3->SetFillColor(0); 
  leg3->SetTextSize(0.03); 

  leg3->AddEntry(gr_ncteq_S, "nCTEQ15 Reweighted 90CL, combined error scale", "p");
  leg3->Draw(); 

  c1->cd(4);  

  gr_ncteq_N->SetTitle("RpAu_ncteq_MB_predictions_fwd");
  gr_ncteq_N->GetXaxis()->SetTitle("Transverse Momentum (GeV/c)");
  gr_ncteq_N->SetMarkerStyle(20);
  gr_ncteq_N->SetMarkerColor(kMagenta+2);
  gr_ncteq_N->SetMarkerSize(1);
  gr_ncteq_N->GetYaxis()->SetRangeUser(0,2.5);
  gr_ncteq_N->GetXaxis()->SetLimits(0,7);
  
  gr_ncteq_N->SetMarkerStyle(20);                                             
  gr_ncteq_N->Draw("AP");
   
   
  TLegend *leg4 = new TLegend(0.1, 0.7, 0.7, 0.9);  //(start x, start y, end x, end y)
  leg4->SetFillColor(0); 
  leg4->SetTextSize(0.03); 

  leg4->AddEntry(gr_ncteq_N, "nCTEQ15 Reweighted 90CL, combined error scale", "p");
  leg4->Draw();  

  //// write mulitple TH2D's to a root file
 
  TFile *h = new TFile("psi2S_RpAu_MB_predictions_pT.root", "RECREATE");
  h->cd();

  gr_epps_N->SetName("RpAu_epps_MB_predictions_fwd");
  gr_epps_S->SetName("RpAu_epps_MB_predictions_bkwd");
  gr_ncteq_N->SetName("RpAu_ncteq_MB_predictions_fwd");
  gr_ncteq_S->SetName("RpAu_ncteq_MB_predictions_bkwd");

  gr_ncteq_N->Write();
  gr_ncteq_S->Write();
  gr_epps_N->Write();
  gr_epps_S->Write();
  
}



  
