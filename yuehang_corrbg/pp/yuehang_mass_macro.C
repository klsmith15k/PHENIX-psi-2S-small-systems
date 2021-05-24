
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
#include <sstream>
#include <cstring>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TNtuple.h>

using namespace std;

void yuehang_mass_macro()
{
  TFile *file_data;

  TH1D *bg[2][5]; 
  TH1F *bg_fx[2];

  //double pt_array_fx[28] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.25,6.75,0.0,0.0};

  std::string filename[2] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/for_krista_mass_arm0_12_5.root",
			     "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/for_krista_mass_arm1_12_5.root"};

  std::string obj_filename[2][5] = {"h_corr_mass","h_bb_mass","h_cc_mass","h_dy_mass","h_pp_mass_MX12_tot",
				    "h_corr_mass","h_bb_mass","h_cc_mass","h_dy_mass","h_pp_mass_MX12_tot"};
  
  std::string obj_filename_fx[2] = {"h_sim_mass_FG12_jet_tot",	
				    "h_sim_mass_FG12_jet_tot"};
  double x_array[150];  
  
  TAxis *xaxis;
  TAxis *yaxis;	  
  
  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Double_t numBinsX;

  for(int arm = 0; arm < 2; arm++)
    { 
      for(int i_histo = 0; i_histo < 5; i_histo++)
	{
	  file_data = TFile::Open(filename[arm].c_str()); 
	  file_data->GetObject(obj_filename[arm][i_histo].c_str(),bg[arm][i_histo]);
	  //file_data->cd();
	  file_data->GetObject(obj_filename_fx[arm].c_str(),bg_fx[arm]);
	}
    }	  
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int j = 0; j < 5; j++)
	{
	  xaxis = bg[arm][j]->GetXaxis();  
	  
	  binCenterX = xaxis->GetBinCenter(1);
	  binLowX = xaxis->GetBinLowEdge(1);
	  binWidthX = xaxis->GetBinWidth(1);
	  numBinsX = xaxis->GetNbins();
	  
	  cout << " ______________________________ " << endl;
	  cout << "x axis starts at: " << binLowX <<  endl;
	  cout << "x binwidth: " << binWidthX <<  endl;
	  cout << "x axis bin 1 center: " << binCenterX << endl;
	  cout << "x axis total bins: " << numBinsX << endl;
	  cout << " ______________________________ " << endl;
	  
	}
    }
  
  // make x_array using center of x axis bins
 for(int arm = 0; arm < 2; arm++)
    {
      for(int j = 0; j < 5; j++)
	{
	  xaxis = bg[arm][j]->GetXaxis();  
	  for(int i = 0; i < 150; i++)
	    {
	      x_array[i]  = xaxis->GetBinCenter(i+1);
	      //cout << "pt array: " << x_array[i] << " for bin " << i+1 << endl;
	    }
	}
    }

  double pt;
  int bins = numBinsX;
  double counter = 0.0;

  double corr_mass[2][150] = {0}; 
  double bb_mass[2][150] = {0};   
  double cc_mass[2][150] = {0}; 
  double drell_mass[2][150] = {0};  
  double pp_mass[2][150] = {0}; 
  double sim_mass[2][150] = {0}; 

  double corr_sum[2][150] = {0};
  double difference[2][150];
    
  // for(int arm = 0; arm < 2; arm++)
  //   {
  //     for(int k = 0; k < 5; k++)
  //       {
  // 	  cout << obj_filename[arm][k].c_str() << endl;
  // 	}
  //   } // prints out the desired system, not all systems (don't have pAl centrality or HeAu centrality files yet)

  TH1D *h_sum[2];

  for(int arm = 0; arm < 2; arm++)
    { 
      h_sum[arm] = (TH1D *) bg[arm][1]->Clone();
      h_sum[arm]->Add(bg[arm][2]);
      h_sum[arm]->Add(bg[arm][3]);
      h_sum[arm]->Add(bg_fx[arm]);

      for(int i = 0; i < bins; i++)
	{
	  corr_mass[arm][i] = bg[arm][0]->GetBinContent(i+1);
	  bb_mass[arm][i] = bg[arm][1]->GetBinContent(i+1);
	  cc_mass[arm][i] = bg[arm][2]->GetBinContent(i+1);
	  drell_mass[arm][i] = bg[arm][3]->GetBinContent(i+1);
	  pp_mass[arm][i] = bg[arm][4]->GetBinContent(i+1);
	  sim_mass[arm][i] = bg_fx[arm]->GetBinContent(i+1);
	  corr_sum[arm][i] = bb_mass[arm][i] + cc_mass[arm][i] + drell_mass[arm][i] + sim_mass[arm][i];
	  
	  difference[arm][i] = corr_sum[arm][i] - corr_mass[arm][i];
      	  //cout << "difference: " << difference[arm][i] << endl;

	  if(difference[arm][i] == 0)
	    counter++;

	  if(i+1 == 102)
	   {
	     cout << "__________________________________ " << endl;     
	     cout << "Bin:" << i+1 << endl;
	     cout << "bb mass " << bb_mass[arm][101] << " , arm " << arm << ", cc mass " << cc_mass[arm][102-1] << endl;
	     cout << "drell mass " << drell_mass[arm][101] << " , arm " << arm << ", corr mass " << corr_mass[arm][101] << endl;
	     cout << "pp mass " << pp_mass[arm][101] << " , arm " << arm << ", sim mass " << sim_mass[arm][101] << endl;
	   }
	 }
    }
  


  TF1 *mix_ul_fit;
  
  mix_ul_fit = new TF1("mix_ul_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);
  mix_ul_fit->SetParameter(0, 0.1);  // fix
  mix_ul_fit->SetParameter(1, -0.017);
  mix_ul_fit->SetParameter(2, 2*pow(10,-8)); // for North Arm
  mix_ul_fit->SetParameter(3, 39);
  mix_ul_fit->SetParameter(4, 21);

  mix_ul_fit->SetLineColor(kRed);
  mix_ul_fit->SetLineStyle(2);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 


  // cout << "" << counter << "/150 have zero difference" << endl;
  /////////// PLOT ///////////////
  
  TCanvas *c1 = new TCanvas("c1","title",200,10,700,500);
  c1->SetGrid();
  c1->SetLogy();
  
  TGraph *gr[10];
  
  gr[0] = new TGraph(150,x_array,corr_mass[0]); 
  gr[1] = new TGraph(150,x_array,corr_mass[1]);
  gr[2] = new TGraph(150,x_array,bb_mass[0]); 
  gr[3] = new TGraph(150,x_array,bb_mass[1]); 
  gr[4] = new TGraph(150,x_array,cc_mass[0]); 
  gr[5] = new TGraph(150,x_array,cc_mass[1]); 
  gr[6] = new TGraph(150,x_array,drell_mass[0]); 
  gr[7] = new TGraph(150,x_array,drell_mass[1]);
  gr[8] = new TGraph(150,x_array,sim_mass[0]); 
  gr[9] = new TGraph(150,x_array,sim_mass[1]);
  gr[10]= new TGraph(150,x_array,corr_mass[0]);
  gr[11]= new TGraph(150,x_array,corr_mass[1]);
  
  // Plot for North pAu,pAl,HeAu
  gr[1]->SetTitle("Yue Hang Results, North");
  
  gPad->SetLeftMargin(0.11); // 0.3
  gPad->SetBottomMargin(0.15);
  gr[1]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr[1]->GetXaxis()->SetLabelSize(0.06);
  gr[1]->GetXaxis()->SetTitleSize(0.07);
  gr[1]->GetXaxis()->SetTitleOffset(0.9);
  // gr[1]->GetYaxis()->SetTitle("Ncounts");
  gr[1]->GetYaxis()->SetLabelSize(0.06);
  gr[1]->GetYaxis()->SetTitleSize(0.1); 
  gr[1]->GetYaxis()->SetTitleOffset(0.52);
  gr[1]->SetLineColor(kBlue); 
  gr[1]->SetLineWidth(3.0);
  //gr[1]->SetMarkerStyle(20); 
  gr[1]->Draw("AL"); 
  // gr[1]->Fit(mix_ul_fit,"LL","",2,5); 

  gr[1]->GetXaxis()->SetLimits(2,5);
  gr[1]->GetYaxis()->SetRangeUser(0.0,25*pow(10,-9));
  
  gr[3]->SetMarkerColor(kRed);  
  gr[3]->SetMarkerSize(1.0);
  gr[3]->SetMarkerStyle(20); 
  // gr[3]->Draw("P");
 
  gr[5]->SetMarkerColor(kGreen);  
  gr[5]->SetMarkerSize(1.0);
  gr[5]->SetMarkerStyle(20); 
  // gr[5]->Draw("P");
 
  gr[7]->SetMarkerColor(kBlack);  
  gr[7]->SetMarkerSize(1.0);
  gr[7]->SetMarkerStyle(20); 
  // gr[7]->Draw("P");

  gr[9]->SetMarkerColor(kCyan);  
  gr[9]->SetMarkerSize(1.0);
  gr[9]->SetMarkerStyle(20); 
  //gr[9]->Draw("P");

  gr[11]->SetMarkerColor(kCyan);  
  gr[11]->SetMarkerSize(1.0);
  gr[11]->SetMarkerStyle(20); 
  gr[11]->Draw("P");

  
  h_sum[1]->SetMarkerColor(kYellow);  
  h_sum[1]->SetMarkerSize(1.0);
  h_sum[1]->SetMarkerStyle(20); 
  h_sum[1]->Draw("SAME");
  
 
  TLegend *leg = new TLegend(0.51, 0.6, 0.8, 0.9);   //0.11
  leg->SetFillColor(0); 
  leg->SetTextSize(0.035);
  leg->AddEntry(gr[1],"corr bg, N", "l");
  // leg->AddEntry(gr[3],"bb bg, N", "p");
  // leg->AddEntry(gr[5],"cc bg, N", "p");
  // leg->AddEntry(gr[7],"drell yan bg, N", "p"); 
  // leg->AddEntry(gr[9],"had bg, N", "p"); 
  leg->AddEntry(gr[11],"sum bb+cc+drell+hadrons, N", "p"); 
  leg->AddEntry(h_sum[1],"histogram sum, N", "p"); 
  leg->Draw();

   // Plot for South pAu,pAl,HeAu

  TCanvas *c2 = new TCanvas("c2","title",200,10,700,500);
  c2->SetGrid();
  c2->SetLogy();

  gr[0]->SetTitle("Yue Hang Results, South");
  gPad->SetLeftMargin(0.11); // 0.3
  gPad->SetBottomMargin(0.15);
  gr[0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr[0]->GetXaxis()->SetLabelSize(0.06);
  gr[0]->GetXaxis()->SetTitleSize(0.07);
  gr[0]->GetXaxis()->SetTitleOffset(0.9);
  //gr[0]->GetYaxis()->SetTitle("Ncounts");
  gr[0]->GetYaxis()->SetLabelSize(0.06);
  gr[0]->GetYaxis()->SetTitleSize(0.1);
  gr[0]->GetYaxis()->SetTitleOffset(0.52);
  gr[0]->SetLineColor(kBlue); 
  gr[0]->SetLineWidth(3.0);
  //gr[0]->SetMarkerStyle(20); 
  gr[0]->Draw("AL"); 
  // mix_ul_fit->SetParameter(2, 1.5*pow(10,-8)); // for South arm
  // gr[0]->Fit(mix_ul_fit,"LL","",2,5); 
  cout << "Chi Square/NDF: " << mix_ul_fit->GetChisquare() << "/" << mix_ul_fit->GetNDF() << "" << endl;
 
  gr[0]->GetXaxis()->SetLimits(2,5);
  gr[0]->GetYaxis()->SetRangeUser(0.0,25*pow(10,-9));
  

  gr[2]->SetMarkerColor(kRed);  
  gr[2]->SetMarkerSize(1.0);
  gr[2]->SetMarkerStyle(20); 
  // gr[2]->Draw("P");
 
  gr[4]->SetMarkerColor(kViolet);  
  gr[4]->SetMarkerSize(1.0);
  gr[4]->SetMarkerStyle(20); 
  // gr[4]->Draw("P");


  gr[6]->SetMarkerColor(kBlack);  
  gr[6]->SetMarkerSize(1.0);
  gr[6]->SetMarkerStyle(20); 
  // gr[6]->Draw("P");

  gr[8]->SetMarkerColor(kCyan);  
  gr[8]->SetMarkerSize(1.0);
  gr[8]->SetMarkerStyle(20); 
  // gr[8]->Draw("P");

  gr[10]->SetMarkerColor(kCyan);  
  gr[10]->SetMarkerSize(1.0);
  gr[10]->SetMarkerStyle(20); 
  gr[10]->Draw("P");

  h_sum[0]->SetMarkerColor(kYellow);  
  h_sum[0]->SetMarkerSize(1.0);
  h_sum[0]->SetMarkerStyle(20); 
  h_sum[0]->Draw("SAME");
   
  TLegend *leg2 = new TLegend(0.51, 0.6, 0.8, 0.9);   //0.11
  leg2->SetFillColor(0); 
  leg2->SetTextSize(0.035);
  leg2->AddEntry(gr[0], "corr bg, S", "l"); 
  // leg2->AddEntry(gr[2],"bb bg, S", "p"); 
  // leg2->AddEntry(gr[4],"cc bg, S", "p"); 
  //  leg2->AddEntry(gr[6],"drell yan bg, S", "p");
  // leg2->AddEntry(gr[8],"had bg, S", "p");
  leg2->AddEntry(gr[10], "sum bb+cc+drell+had, S", "p"); 
  leg2->AddEntry(h_sum[0],"histogram sum, S", "p");
  leg2->Draw();

  // plot corr_mass with fit
  TCanvas *c3 = new TCanvas("c3","title",200,10,700,500);
  c3->Divide(2,1);
  c3->cd(2);
  gPad->SetLogy();
  c3->SetGrid();

  gr[0]->SetTitle("Yue Hang Corr Mass, South");
  gPad->SetLeftMargin(0.11); // 0.3
  gPad->SetBottomMargin(0.15);
  gr[0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr[0]->GetXaxis()->SetLabelSize(0.06);
  gr[0]->GetXaxis()->SetTitleSize(0.05);
  gr[0]->GetXaxis()->SetTitleOffset(0.9);
  //gr[0]->GetYaxis()->SetTitle("Ncounts");
  gr[0]->GetYaxis()->SetLabelSize(0.06);
  gr[0]->GetYaxis()->SetTitleSize(0.06);
  gr[0]->GetYaxis()->SetTitleOffset(0.99);
  gr[0]->SetLineColor(kViolet); 
  gr[0]->SetLineWidth(3.0);
  //gr[0]->SetMarkerStyle(20); 
  gr[0]->Draw("AL"); 
  gr[0]->GetXaxis()->SetLimits(2,5);
  gr[0]->GetYaxis()->SetRangeUser(pow(10,-10),25*pow(10,-9));
  gr[0]->Draw();
 
  c3->cd(1);
  gPad->SetLogy();
  c3->SetGrid();

  gr[1]->SetTitle("Yue Hang Corr Mass, North");
  gPad->SetLeftMargin(0.11); // 0.3
  gPad->SetBottomMargin(0.15);
  gr[1]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr[1]->GetXaxis()->SetLabelSize(0.06);
  gr[1]->GetXaxis()->SetTitleSize(0.05);
  gr[1]->GetXaxis()->SetTitleOffset(0.9);
  //gr[1]->GetYaxis()->SetTitle("Ncounts");
  gr[1]->GetYaxis()->SetLabelSize(0.06);
  gr[1]->GetYaxis()->SetTitleSize(0.06);
  gr[1]->GetYaxis()->SetTitleOffset(0.99);
  gr[1]->SetLineColor(kBlue); 
  gr[1]->SetLineWidth(3.0);
  gr[1]->Draw("AL"); 
  gr[1]->GetXaxis()->SetLimits(2,5);
  gr[1]->GetYaxis()->SetRangeUser(pow(10,-10),25*pow(10,-9));
  gr[1]->SetLineColor(kBlue); 
  gr[1]->SetLineWidth(3.0);
  gr[1]->Draw("AL");

} // void end macro



