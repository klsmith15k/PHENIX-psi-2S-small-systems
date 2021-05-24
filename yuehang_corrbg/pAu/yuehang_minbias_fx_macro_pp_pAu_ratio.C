

#include <TF1.h>
#include <TMinuit.h>
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
#include <TPaveText.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TNtuple.h>

using namespace std;

void yuehang_minbias_fx_macro_pp_pAu_ratio()
{

  double x_array[100] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,5.05,5.15,5.25,5.35,5.45,5.55,5.65,5.75,5.85,5.95,6.05,6.15,6.25,6.35,6.45,6.55,6.65,6.75,6.85,6.95,7.05,7.15,7.25,7.35,7.45,7.55,7.65,7.75,7.85,7.95,8.05,8.15,8.25,8.35,8.45,8.55,8.65,8.75,8.85,8.95,9.05,9.15,9.25,9.35,9.45,9.55,9.65,9.75,9.85,9.95};

//for each refit evaluated at each mass bin -> pAu refits
#include "pt_int/Run15pAu_reco_pt_int_N.C"
#include "pt_int/Run15pAu_reco_pt_int_S.C"

  // reconstructed pt integrated fit for each mass bin -> pp 
#include "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_pt_int/Run15pp_reco_pt_int_N.C"
#include "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/yuehang_pt_int/Run15pp_reco_pt_int_S.C"

  bool south_arm = false;
  //bool south_arm = true;
  
  int arm_low;

  if(south_arm)
    arm_low = 0;
  else
    arm_low = 1;

  double x_errors[2][100] = {0};  
  double y_errors[2][100] = {0};
  
  double fit_pAu[2][100] = {0};
  double fit_pp[2][100] = {0};   
  double pAu_pp_ratio[2][100] = {0};

  int bins = 100;
  
  double scale[2] = {0.00025,0.00075};
 
  for(int arm = 0; arm < 2; arm++)
   {
     for(int bin = 0; bin < 100; bin++)
       {
	 fit_pp[0][bin] = Run15pp_reco_pt_int_S[bin];
	 fit_pp[1][bin] = Run15pp_reco_pt_int_N[bin];

	 fit_pAu[0][bin] = Run15pAu_reco_pt_int_S[bin];
	 fit_pAu[1][bin] = Run15pAu_reco_pt_int_N[bin];

	 pAu_pp_ratio[arm][bin] = fit_pAu[arm][bin] / fit_pp[arm][bin];
	 y_errors[arm][bin] = sqrt(scale[arm]*pAu_pp_ratio[arm][bin]);
       }
   }

   cout << "South num pp: " << fit_pp[0][26] << endl;
   cout << "North num pp: " << fit_pp[1][26] << endl;
 
   cout << "South num pAu: " << fit_pAu[0][26] << endl;
   cout << "North num pAu: " << fit_pAu[1][26] << endl;
 
   double scale_peak[2] = {317.361/(152.251),358.544/(159.877)};
 
   for(int arm = 0; arm < 2; arm++)
     {
       for(int bin = 0; bin < 100; bin++)
	 {
	   pAu_pp_ratio[arm][bin] *= scale_peak[arm];
	 }
     }

  //PLOT RATIOS
  
  TCanvas *c20 = new TCanvas("c20","Ratio of pAu/pp Corr BG pT Integrated",200,10,700,500);
  gPad->SetGrid();
  
  TGraphErrors *gr1 = new TGraphErrors(bins,x_array,pAu_pp_ratio[arm_low],x_errors[arm_low],y_errors[arm_low]);
  
  if(south_arm)
    gr1->SetTitle("Ratio of pAu/pp Corr BG pT Integrated, South");  
  else
    gr1->SetTitle("Ratio of pAu/pp Corr BG pT Integrated, North");  
  gr1->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr1->GetXaxis()->SetLabelSize(0.04);
  gr1->GetXaxis()->SetTitleSize(0.04);
  gr1->GetXaxis()->SetTitleOffset(0.9);
  gr1->GetXaxis()->SetLimits(2,5);
  gr1->GetYaxis()->SetLabelSize(0.04);
  gr1->GetYaxis()->SetTitleSize(0.1); 
  gr1->GetYaxis()->SetTitleOffset(0.52);
  gr1->GetYaxis()->SetRangeUser(0.4,1.4);
  
  if(south_arm)
    gr1->SetMarkerColor(kCyan+1);  
  else
    gr1->SetMarkerColor(kCyan+1);  
  gr1->SetMarkerSize(0.75);
  gr1->SetMarkerStyle(20);
  gr1->Draw("AP");
  
  TF1 *fit;
  
  fit = new TF1("fit","[0]",2,5);
  fit->SetParameter(0,1);
  fit->SetLineStyle(10);
  
  gr1->Fit(fit,"R"); 
  
  Char_t message[200];
  
  sprintf(message,"y = %.3f +/- %.3f",fit->GetParameter(0),fit->GetParError(0));
  
  TPaveText *mytext = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
  mytext->SetTextSize(0.035);
  mytext->SetFillColor(0);
  mytext->SetTextAlign(12);
  mytext->AddText(message);
  mytext->Draw();
   

 
  
}// end void 


