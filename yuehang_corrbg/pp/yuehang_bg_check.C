// original binwidths were 1 GeV


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

void yuehang_bg_check()
{

  bool write = false;
  bool norm = true;


  ////////////////////////////////////////////////// 
 
  double scale[2] = {2.415*pow(10,11),0.045*pow(10,13)};
   
  int pt_binwidth;
  cout << "What pt binning are you running?  Enter '0' for Yue Hang binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth" << endl;
  cin >> pt_binwidth;

  int pt_slices[5] = {4,38,20,10,5};
  int num_bins_on_plots = 10;
  double pt_width[5][38] = {1.0,1.0,1.0,7.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
			    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    2.0,2.0,2.0,2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  double pt_center[5][38] = {0.5,1.5,2.5,6.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     1.0,3.0,5.0,7.0,9.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
 
  std::string filename[5][2] = {"yuehang_Run15pp_S_rebinned_1000_mev.root", 
				"yuehang_Run15pp_N_rebinned_1000_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_200_mev.root", 
				 "yuehang_Run15pp_N_rebinned_200_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_500_mev.root", 
				 "yuehang_Run15pp_N_rebinned_500_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_1_gev.root", 
				"yuehang_Run15pp_N_rebinned_1_gev.root",

 				 "yuehang_Run15pp_S_rebinned_2_gev.root", 
				 "yuehang_Run15pp_N_rebinned_2_gev.root"};
 
 
  //////////////////////////////////////////////////
  //everything below this line is created as a result of making an initial selection about which data you want to use

  TH1D *bg[6][2][38] ={0};
  TH1D *t_check[2]; 
  std::string obj_filename[6][2][38]; 
  
  char unique1[800];
  char unique2[800];
  char unique3[800];
  char unique4[800];
  char unique5[500];
  char unique6[500];

  TAxis *xaxis;
  TAxis *yaxis;	  
 
// corr. hadron binning is 150 along x-axis, but cc,bb,dy are 100 bins along x-axis
  int bins = 100;

  double pt_low;
  double pt_high;

  double bin_low;
  double bin_high;

  double delta_pt = 0.001;
  double a_array[38];
  double b_array[38];
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < pt_slices[pt_binwidth]; k++)
	{ 
	  pt_low = pt_center[pt_binwidth][k] - pt_width[pt_binwidth][k]/2 + delta_pt; 
	  pt_high = pt_center[pt_binwidth][k] + pt_width[pt_binwidth][k]/2 - delta_pt;
	  
	  double a = (pt_low - delta_pt)*1000; // for unique names in projected histograms, use meV
	  double b = (pt_high + delta_pt)*1000;

	  cout << "a: " << ", b: " << b << endl;
	   
	  sprintf(unique1,"cc_%.0f_%.0f",a,b); 
	  sprintf(unique2,"bb_UL_%.0f_%.0f",a,b); 
	  sprintf(unique3,"dy_%.0f_%.0f",a,b);  
	  sprintf(unique4,"bb_LS_%.0f_%.0f",a,b);
	  sprintf(unique5,"corr_had_LS_%.0f_%.0f",a,b);  
	  sprintf(unique6,"corr_had_UL_%.0f_%.0f",a,b); 
 
	  
	  a_array[k] = a/1000;
	  b_array[k] = b/1000;
	  cout << "a: " << a_array[k] << ", b: " << b_array[k] << endl;

	  obj_filename[0][arm][k] = unique1;
	  obj_filename[1][arm][k] = unique2;
	  obj_filename[2][arm][k] = unique3;
	  obj_filename[3][arm][k] = unique4;	
	  obj_filename[4][arm][k] = unique5;	
	  obj_filename[5][arm][k] = unique6;	
	}
    }
  
  for(int i_histo = 0; i_histo < 6; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      // cout << obj_filename[i_histo][arm][pt].c_str() << " HERE " << endl;
	    }
	}
    }
  
   TFile *file_data;
   TFile *file_check;
  
   std::string filename_ch[2] = {"for_krista_mass_arm0_12_5.root",
				 "for_krista_mass_arm1_12_5.root"};
   
   double x_array[100];  
   
  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Double_t numBinsX;
  
  for(int i_histo = 0; i_histo < 6; i_histo++)
    { 
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      file_data = TFile::Open(filename[pt_binwidth][arm].c_str()); 
	      file_data->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	      file_check = TFile::Open(filename_ch[arm].c_str()); 
	      file_check->GetObject("h_corr_mass",t_check[arm]);
	    }
	}	  
    }

  for(int i_histo = 0; i_histo < 6; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      
	      xaxis = bg[i_histo][arm][pt]->GetXaxis();  
	      binCenterX = xaxis->GetBinCenter(1);
	      binLowX = xaxis->GetBinLowEdge(1);
	      binWidthX = xaxis->GetBinWidth(1);
	      numBinsX = xaxis->GetNbins();
	      // cout << " ______________________________ " << endl;
	      // cout << obj_filename[i_histo][arm][pt].c_str() << endl;
	      // cout << "x axis starts at: " << binLowX <<  endl;
	      // cout << "x binwidth: " << binWidthX <<  endl;
	      // cout << "x axis bin 1 center: " << binCenterX << endl;
	      // cout << "x axis total bins: " << numBinsX << endl;
	      // cout << " ______________________________ " << endl;
	      
	    }
	}
    }
    
  xaxis = bg[0][1][1]->GetXaxis();   
  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
      //cout << "pt array: " << x_array[i] << " for bin " << i+1 << endl;
    }

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 6; i_histo++)
        {
  	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  	    {
  	      //cout << obj_filename[i_histo][arm][pt].c_str() << ", i_histo: " << i_histo << ", arm: " << arm << ", pt: " << pt << endl;
  	    }
  	} 
    }

  TH1D *h_sum[2][38];
 
  double cc[2][38][100];
  double bb_UL[2][38][100];
  double dy[2][38][100];
  double bb_LS[2][38][100];
  double corr_had_LS[2][38][100];
  double corr_had_UL[2][38][100];
  double corr_bg[2][38][100];
  double corr_total[2][100] = {0.0,0.0};
  double corr_mass[2][100];
  double errors[2][38][100];
  double x_errors[2][38][100];
  double counter = 0.0;

  // for bg check:
  double corr_bg_UL[2][38][100];
  double corr_bg_LS_v1[2][38][100];
  double corr_bg_LS_v2[2][38][100];
  double corr_bg_LS_v3[2][38][100];
  double corr_bg_UL_total[2][100];
  double corr_bg_LS_v1_total[2][100];
  double corr_bg_LS_v2_total[2][100];
  double corr_bg_LS_v3_total[2][100];

	  
  for(int i = 0; i < 6; i++)
    {
      for(int j = 0; j < 2; j++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      bg[i][j][pt]->Scale(scale[j]);
	      t_check[j]->Scale(scale[j]);
	    }
	}
    }
  
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  	{
  	  for(int i = 0; i < bins; i++)
  	    {
	      cc[arm][pt][i] =  (bg[0][arm][pt]->GetBinContent(i+1));
	      bb_UL[arm][pt][i] = (bg[1][arm][pt]->GetBinContent(i+1));
	      dy[arm][pt][i] = (bg[2][arm][pt]->GetBinContent(i+1));
	      bb_LS[arm][pt][i] = (bg[3][arm][pt]->GetBinContent(i+1));
	      corr_had_LS[arm][pt][i] = (bg[4][arm][pt]->GetBinContent(i+1));
	      corr_had_UL[arm][pt][i] = (bg[5][arm][pt]->GetBinContent(i+1));

	      corr_bg[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i] + corr_had_UL[arm][pt][i] - corr_had_LS[arm][pt][i];
	      
	      corr_mass[arm][i] = (t_check[arm]->GetBinContent(i+1));
	      corr_total[arm][i]+= corr_bg[arm][pt][i];

	      // for bg check:
	      corr_bg_UL[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] + corr_had_UL[arm][pt][i]; 
	      corr_bg_LS_v1[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i];
	      corr_bg_LS_v2[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] + corr_had_UL[arm][pt][i]; // mixed events
	      corr_bg_LS_v3[arm][pt][i] =  cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] + corr_had_UL[arm][pt][i] - corr_had_LS[arm][pt][i];
	     
	      corr_bg_UL_total[arm][i] += corr_bg_UL[arm][pt][i];
	      corr_bg_LS_v1_total[arm][i] += corr_bg_LS_v1[arm][pt][i];
	      corr_bg_LS_v2_total[arm][i] += corr_bg_LS_v2[arm][pt][i];
	      corr_bg_LS_v3_total[arm][i] += corr_bg_LS_v3[arm][pt][i];
	
	    }
	 	   
  	}
    }
	  
  char uniqueh[38][500];
  int i = 0;
  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      i++;
      sprintf(uniqueh[pt],"h_sum_%d",i); 
      //cout << uniqueh[pt] << endl;
    }

 TH1D *h1 = new TH1D("h1", "Distribution", 100, 0, 10);

  for(int arm = 0; arm < 2; arm++)  // cannot use add method when histogramsdo not have exact same x-axis properties
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {
	      h1->SetBinContent(i+1,(bg[0][arm][pt]->GetBinContent(i+1) + bg[1][arm][pt]->GetBinContent(i+1) + bg[2][arm][pt]->GetBinContent(i+1) - bg[3][arm][pt]->GetBinContent(i+1))); 

	      if(h1->GetBinContent(i+1) > 0)
	      	h1->SetBinError(i+1,h1->GetBinContent(i+1)*0.15);

	      h_sum[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_sum[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));
	      h_sum[arm][pt]->SetBinError(i+1,h1->GetBinError(i+1));



	      // h_sum[arm][pt] = (TH1D *) bg[0][arm][pt]->Clone();
	      // h_sum[arm][pt]->SetBinContent(i+1,(bg[0][arm][pt]->GetBinContent(i+1) + bg[1][arm][pt]->GetBinContent(i+1) + bg[2][arm][pt]->GetBinContent(i+1) - bg[3][arm][pt]->GetBinContent(i+1))); 
	      // x_errors[arm][pt][i] = 0.0;
	      // //h_sum[arm][pt]->SetBinError(i+1,(bg[0][arm][pt]->GetBinError(i+1) + bg[1][arm][pt]->GetBinError(i+1) + bg[2][arm][pt]->GetBinError(i+1) + h_sum[arm][pt]->GetBinContent(i+1)/4.0));

	      // // yue hang email says conervatively cc ~10%, bb ~ 15% and DY ~ 10%
	      // h_sum[arm][pt]->SetBinError(i+1,(h_sum[arm][pt]->GetBinContent(i+1)*0.1 + h_sum[arm][pt]->GetBinContent(i+1)*0.15 + h_sum[arm][pt]->GetBinContent(i+1)*0.1 + h_sum[arm][pt]->GetBinContent(i+1)*0.15));
	      //h_sum[arm][pt]->SetBinError(i+1,(bg[0][arm][pt]->GetBinContent(i+1)*0.1 + bg[1][arm][pt]->GetBinContent(i+1)*0.15 + bg[2][arm][pt]->GetBinContent(i+1)*0.1 + bg[3][arm][pt]->GetBinContent(i+1)*0.15));
	      	      
	      //	cout << "bin content : " << h_sum[arm][pt]->GetBinContent(i+1) << ", and error : " << h_sum[arm][pt]->GetBinError(i+1) << endl;
	    }
	}
    }

  double area_h[2][38];
  double norm_h[2][38];
    
  if(norm)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      area_h[arm][pt] = h_sum[arm][pt]->Integral( h_sum[0][0]->FindBin(2.001),  h_sum[0][0]->FindBin(4.999) );
	     
	      if(area_h[arm][pt] != 0)
		{
		  norm_h[arm][pt] = 1/(area_h[arm][pt]);
		  h_sum[arm][pt]->Scale(norm_h[arm][pt]);
		}
	      //norm_area[arm][pt] = h_sum[arm][pt]->Integral( h_sum[0][0]->FindBin(2.001),  h_sum[0][0]->FindBin(4.999) );
	    }
	}
    }

 
  TGraph *gr_corrbg_UL = new TGraph(bins,x_array,corr_bg_UL_total[1]); 
  TGraph *gr_corrbg_v1  = new TGraph(bins,x_array,corr_bg_LS_v1_total[1]); 
  TGraph *gr_corrbg_v2  = new TGraph(bins,x_array,corr_bg_LS_v2_total[1]);
  TGraph *gr_corrbg_v3  = new TGraph(bins,x_array,corr_bg_LS_v3_total[1]);
  TGraph *corrbg_total = new TGraph(bins,x_array,corr_total[1]);
  
  TCanvas *c3 = new TCanvas("c3","Corr BG summed comparisons",200,10,700,500);
  gPad->SetGrid();
  gPad->SetLogy();
    
  gr_corrbg_UL->SetTitle("Corr BG Comparison");  
  gr_corrbg_UL->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr_corrbg_UL->GetXaxis()->SetLabelSize(0.04);
  gr_corrbg_UL->GetXaxis()->SetTitleSize(0.04);
  gr_corrbg_UL->GetXaxis()->SetTitleOffset(0.9);
  gr_corrbg_UL->GetXaxis()->SetLimits(2,5);
  gr_corrbg_UL->GetYaxis()->SetLabelSize(0.04);
  gr_corrbg_UL->GetYaxis()->SetTitleSize(0.1); 
  gr_corrbg_UL->GetYaxis()->SetTitleOffset(0.52);
  gr_corrbg_UL->SetMaximum(pow(10,5));
  gr_corrbg_UL->SetMinimum(pow(10,2));

  /// assign different colors and draw each bg
  gr_corrbg_UL->SetLineColor(kAzure);  
  gr_corrbg_UL->SetLineWidth(5.0);
  gr_corrbg_UL->SetLineStyle(2);
  gr_corrbg_UL->Draw("AL");
  
  gr_corrbg_v1->SetLineColor(kSpring);  
  gr_corrbg_v1->SetLineWidth(5.0);
  gr_corrbg_v1->SetLineStyle(2);
  gr_corrbg_v1->Draw("L");  

  gr_corrbg_v2->SetLineColor(kCyan);  
  gr_corrbg_v2->SetLineWidth(5.0);
  gr_corrbg_v2->SetLineStyle(2);
  //gr_corrbg_v2->Draw("L");

  gr_corrbg_v3->SetLineColor(kMagenta);  
  gr_corrbg_v3->SetLineWidth(5.0);
  gr_corrbg_v3->SetLineStyle(2);
  gr_corrbg_v3->Draw("L");

  corrbg_total->SetLineColor(kBlack);  
  corrbg_total->SetLineWidth(5.0);
  corrbg_total->SetLineStyle(2);
  corrbg_total->Draw("L");

 TLegend *leg;
 leg = new TLegend(0.51, 0.6, 0.8, 0.9);  
 leg->SetFillColor(0); 
 leg->SetTextSize(0.035);
 leg->AddEntry(gr_corrbg_UL,"corr bg mixed events", "l");
 leg->AddEntry(gr_corrbg_v1,"corr bg like-sign", "l");
 // leg->AddEntry(gr_corrbg_v2,"corr bg like-sign test 1", "l");
 leg->AddEntry(gr_corrbg_v3,"corr bg like-sign test 2", "l");
 leg->AddEntry(corrbg_total,"corr bg like-sign with hadrons", "l");
 
 leg->Draw();


/*


  double a = 3;

  //TGraph *gr_corrbg[2][2];
  TGraph *gr[2][38];

  for(int arm = 0; arm < 2; arm++) 
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	   gr[arm][pt] = new TGraph(bins,x_array,corr_bg[arm][pt]); //100 bins is 0.0 - 10.0 GeV/c^2
	   gr[arm][pt]->SetLineColor(pt+1);
	   if(pt+1 > 9)
	     gr[arm][pt]->SetLineColor(pt+2);
	   gr[arm][pt]->SetLineWidth(3.0);
	   // gr_corrbg[arm][0] = new TGraph(bins,x_array,corr_total[arm]); 
	   // gr_corrbg[arm][1] = new TGraph(bins,x_array,corr_mass[arm]); 
      
	  if(pt == 0)
	    {	     
	      gr[arm][pt]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
	      gr[arm][pt]->GetXaxis()->SetLabelSize(0.04);
	      gr[arm][pt]->GetXaxis()->SetTitleSize(0.04);
	      gr[arm][pt]->GetXaxis()->SetTitleOffset(0.9);
	      gr[arm][pt]->GetXaxis()->SetLimits(2,5);
	      gr[arm][pt]->GetYaxis()->SetLabelSize(0.04);
	      gr[arm][pt]->GetYaxis()->SetTitleSize(0.1); 
	      gr[arm][pt]->GetYaxis()->SetTitleOffset(0.52);
	      gr[arm][pt]->SetMaximum(pow(10,4));
	      gr[arm][pt]->SetMinimum(pow(10,0));
	    }
	
	}
    }

  gr[1][0]->SetTitle("Rebinned Yue Hang Results, North");

 
  //for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  // for(int pt = 0; pt < num_bins_on_plots; pt++)
  //   {
  //     if(pt == 0)
  // 	gr_corrbg[1][pt]->Draw("AL");
  //     //  else
  //     //	gr[1][pt]->Draw("L");
  //   }
  
  gr_corrbg[1][0]->SetLineColor(kAzure);  
  gr_corrbg[1][0]->SetLineWidth(5.0);
  gr_corrbg[1][0]->Draw("AL");
 
  gr_corrbg[1][1]->SetLineColor(kCyan);  
  gr_corrbg[1][1]->SetLineWidth(3.0);
  gr_corrbg[1][1]->SetLineStyle(10);
  gr_corrbg[1][1]->Draw("L");
   
  
  char unique[800];
  TLegend *leg[4];

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_graph = 0; i_graph < 4; i_graph++)
	{
	  leg[i_graph] = new TLegend(0.51, 0.6, 0.8, 0.9);  
	  leg[i_graph]->SetFillColor(0); 
	  leg[i_graph]->SetTextSize(0.035);
	  // leg[i_graph]->AddEntry(gr_corrbg[arm][0],"corr bg (sum of p_T bins)", "l");
	  //leg[i_graph]->AddEntry(gr_corrbg[arm][1],"corr bg (p_T integrated)", "l");
	 
	  //for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	  for(int pt = 0; pt < num_bins_on_plots; pt++)
	    {
	      cout << "a: " << a_array[pt] << ", b: " << b_array[pt] << endl;
	      sprintf(unique,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	      leg[i_graph]->AddEntry(gr[arm][pt],unique, "l");
	    }
	}
    }

  leg[0]->Draw();

  // Plot South TGraphs
  TCanvas *c4 = new TCanvas("c4","TGraph 2",200,10,700,500);
  gPad->SetLogy();
  gPad->SetGrid();
  gr[0][0]->SetTitle("Rebinned Yue Hang Results, South");
  
  // for(int pt = 0; pt < num_bins_on_plots; pt++)
  
  //   {
  //     if(pt == 0)
  // 	gr[0][pt]->Draw("AL");
  //     else
  // 	gr[0][pt]->Draw("L");
  //   }

  gr_corrbg[0][0]->SetLineColor(kAzure);  
  gr_corrbg[0][0]->SetLineWidth(5.0);
  gr_corrbg[0][0]->Draw("AL");
  
  gr_corrbg[0][1]->SetLineColor(kCyan);  
  gr_corrbg[0][1]->SetLineWidth(3.0);
  gr_corrbg[0][1]->SetLineStyle(10);
  gr_corrbg[0][1]->Draw("L");
 
  leg[1]->Draw();

  TCanvas *c1 = new TCanvas("c1","Histogram 1",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetGrid();
  
  for(int arm = 0; arm < 2; arm++)
    {
      h_sum[arm][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm][0]->GetXaxis()->SetLabelSize(0.04);
      h_sum[arm][0]->GetXaxis()->SetTitleSize(0.04);
      h_sum[arm][0]->GetXaxis()->SetTitleOffset(0.9);
      h_sum[arm][0]->SetAxisRange(2,5); 
      h_sum[arm][0]->SetMaximum(pow(10,4));
      h_sum[arm][0]->SetMinimum(pow(10,0));
      if(norm)
	{
	  h_sum[arm][0]->SetMaximum(pow(10,0));
	  h_sum[arm][0]->SetMinimum(pow(10,-3));
	}
    }

  h_sum[0][0]->SetTitle("Yue Hang Corr BG (cc+bb_UL+dy-bb_LS), South");
  
  //for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  for(int pt = 0; pt < num_bins_on_plots; pt++)
    {
      h_sum[0][pt]->SetLineColor(pt+1);  
      h_sum[0][pt]->SetLineWidth(3.0);
      if(pt == 0)
	h_sum[0][pt]->Draw();
      else
	h_sum[0][pt]->Draw("SAME");
    }

  // gr_corrbg[0][0]->Draw("L");
  // gr_corrbg[0][1]->Draw("L");
  
  leg[2]->Draw();


  // North Histogram Results
  TCanvas *c2 = new TCanvas("c2","Histogram 2",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetGrid();
  h_sum[1][0]->SetTitle("Yue Hang Corr BG (cc+bb_UL+dy-bb_LS), North");
   
  //for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  for(int pt = 0; pt < num_bins_on_plots; pt++)
    {
      h_sum[1][pt]->SetLineColor(pt+1);  
      h_sum[1][pt]->SetLineWidth(3.0);
      if(pt == 0)
	h_sum[1][pt]->Draw();
      else
	h_sum[1][pt]->Draw("SAME");
    }

  // gr_corrbg[1][0]->Draw("L");
  // gr_corrbg[1][1]->Draw("L");
  
  leg[3]->Draw();

  //////////////////////////// Fit histograms //////////////////////////
  
  char name1[800];
  char name2[800];

  // write out bestfit parameters for the fits that CONVERGE
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;
  int bins_fail[bins];

  TF1 *mix_ul_fit;
  
  mix_ul_fit = new TF1("mix_ul_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);
 
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < num_bins_on_plots - 1; pt++)
	{
	  cout << "Arm " << arm << ", Bin " << pt+1 << endl;
	  mix_ul_fit->SetParameter(0, 0.295358);  
	  // mix_ul_fit->SetParLimits(0,1,100);
	  mix_ul_fit->SetParameter(1, -0.040983); 
	  mix_ul_fit->SetParameter(2, 5.0); 
	  mix_ul_fit->SetParameter(3, 12.9981); 
	  mix_ul_fit->SetParameter(4, 18.7372); 
	  
	  mix_ul_fit->SetLineColor(kRed);
	  mix_ul_fit->SetLineStyle(2);
	  
	   
	  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 
	  
	 
	  
	  TFitResultPtr r = h_sum[arm][pt]->Fit(mix_ul_fit,"S R"); // fit already defined for range 2,5 in defintion

	  h_sum[arm][pt]->SetAxisRange(2,5); 
	  h_sum[arm][0]->SetMaximum(pow(10,4));
	  h_sum[arm][0]->SetMinimum(pow(10,0));
	  if(norm)
	    {
	      h_sum[arm][0]->SetMaximum(pow(10,0));
	      h_sum[arm][0]->SetMinimum(pow(10,-3));
	    }
	  
	  cout << "Chi Square/NDF: " << mix_ul_fit->GetChisquare() << "/" << mix_ul_fit->GetNDF() << "" << endl;
	  p0 = r->Parameter(0);
	  p1 = r->Parameter(1);
	  p2 = r->Parameter(2);
	  p3 = r->Parameter(3);
	  p4 = r->Parameter(4);

	  e0 = r->ParError(0);
	  e1 = r->ParError(1);
	  e2 = r->ParError(2);
	  e3 = r->ParError(3);
	  e4 = r->ParError(4);

	  // r->Print("V");
	  Int_t fitStatus = r;
	  // cout << "status: " << fitStatus << " for pt bin " << pt+1 << endl;
	  if(fitStatus == 0) // 0 indicates it converged, and 4 indicates it failed 
	    // https://root.cern.ch/doc/master/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
	    {
	      if(write)	     
		{
		  if(arm == 1)
		    {
		      sprintf(name1,"yuehang_bestfit_1gev/bestfit_parameters_N_%i.dat",pt+1); 
		      std::fstream bestfit_parameters_1(name1,std::ofstream::out); 
		      bestfit_parameters_1 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_1.close();
		    }
		  else
		    {
		      sprintf(name2,"yuehang_bestfit_1gev/bestfit_parameters_S_%i.dat",pt+1); 
		      std::fstream bestfit_parameters_2(name2,std::ofstream::out); 
		      bestfit_parameters_2 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_2.close();
		    }
		} // write
	      
	    }// fit converged
	  else
	    {
	      counter+= 1;
	      cout << "number of fits that failed : " << counter << endl;
	    }
	} // for loop pt

    } // for loop arm

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i = 0; i < bins;  i++)
	{
	  // cout << "corr mass : " << corr_mass[arm][i] << ", corr total : " << corr_total[arm][i] << " for bin : " << i+1 << endl;
	}
    } 
   
	
  // Draw pT Integrated plot only
  TCanvas *c6 = new TCanvas("c6","Corr BG pT Int",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetGrid();

  gr_corrbg[0][0]->SetTitle("pT Int. Corr BG (cc+bb_UL+dy-bb_LS), South");
  gr_corrbg[0][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr_corrbg[0][0]->GetXaxis()->SetLabelSize(0.04);
  gr_corrbg[0][0]->GetXaxis()->SetTitleSize(0.04);
  gr_corrbg[0][0]->GetXaxis()->SetTitleOffset(0.9);
  gr_corrbg[0][0]->GetXaxis()->SetLimits(2,5);
  gr_corrbg[0][0]->GetYaxis()->SetLabelSize(0.04);
  gr_corrbg[0][0]->GetYaxis()->SetTitleSize(0.1); 
  gr_corrbg[0][0]->GetYaxis()->SetTitleOffset(0.52);
  
  gr_corrbg[0][0]->SetLineColor(kCyan);  
  gr_corrbg[0][0]->SetLineWidth(3.0);
  gr_corrbg[0][0]->SetLineStyle(2);

  gr_corrbg[0][0]->Draw("AL");

  gr_corrbg[0][0]->SetMaximum(pow(10,6));
  gr_corrbg[0][0]->SetMinimum(pow(10,0));

  */
   

}// void end macro
  
  
  
  
