// this macro opens and plots yue hang's TH1D histograms that he projected himself of the four bg contributions, in the file "for_krista_mass_arm0.root" and "for_krista_mass_arm1.root".  The other macro of similar name, "yuehang_pt_macro.C" opens and plots the TH2D histograms that I projected myself using the rebinning macro "yuehang_rebinning_macro.C"


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

void yuehang_pt_macro_sqbr()
{
  
  bool tgraph = true;

  TFile *file_data;
  TFile *file_check;


  TH1D *bg[4][2][4];
  TH1D *t_check[2]; 
  
  std::string filename[2] = {"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/for_krista_mass_arm0_12_5.root",
			     "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/tony_bestfit_parameters/for_krista_mass_arm1_12_5.root"};

  std::string obj_filename[4][2][4] = {
    "h_cc_mass_ptslice[0]","h_cc_mass_ptslice[1]","h_cc_mass_ptslice[2]","h_cc_mass_ptslice[3]",
    "h_cc_mass_ptslice[0]","h_cc_mass_ptslice[1]","h_cc_mass_ptslice[2]","h_cc_mass_ptslice[3]",
    
    "h_bb_mass_ptslice[0]","h_bb_mass_ptslice[1]","h_bb_mass_ptslice[2]","h_bb_mass_ptslice[3]",
    "h_bb_mass_ptslice[0]","h_bb_mass_ptslice[1]","h_bb_mass_ptslice[2]","h_bb_mass_ptslice[3]",
    
    "h_dy_mass_ptslice[0]","h_dy_mass_ptslice[1]","h_dy_mass_ptslice[2]","h_dy_mass_ptslice[3]",
    "h_dy_mass_ptslice[0]","h_dy_mass_ptslice[1]","h_dy_mass_ptslice[2]","h_dy_mass_ptslice[3]",
    
    "h_sim_jet_mass_ptslice_FG12_tot[0]","h_sim_jet_mass_ptslice_FG12_tot[1]","h_sim_jet_mass_ptslice_FG12_tot[2]","h_sim_jet_mass_ptslice_FG12_tot[3]",
    "h_sim_jet_mass_ptslice_FG12_tot[0]","h_sim_jet_mass_ptslice_FG12_tot[1]","h_sim_jet_mass_ptslice_FG12_tot[2]","h_sim_jet_mass_ptslice_FG12_tot[3]"};
  

 double x_array[150];  
  
  TAxis *xaxis;
  TAxis *yaxis;	  
  
  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Double_t numBinsX;

for(int i_histo = 0; i_histo < 4; i_histo++)
  {
    for(int arm = 0; arm < 2; arm++)
      { 
	for(int pt = 0; pt < 4; pt++)
	  {
	    file_data = TFile::Open(filename[arm].c_str()); 
	    file_data->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	    file_check = TFile::Open(filename[arm].c_str()); 
	    file_check->GetObject("h_corr_mass",t_check[arm]);
	  }
      }	  
  }

 for(int i_histo = 0; i_histo < 4; i_histo++)
   {
     for(int arm = 0; arm < 2; arm++)
       {
	 for(int pt = 0; pt < 4; pt++)
	   {
	     xaxis = bg[i_histo][arm][pt]->GetXaxis();  
	     
	     binCenterX = xaxis->GetBinCenter(1);
	     binLowX = xaxis->GetBinLowEdge(1);
	     binWidthX = xaxis->GetBinWidth(1);
	     numBinsX = xaxis->GetNbins();
	     
	     cout << " ______________________________ " << endl;
	     cout << obj_filename[i_histo][arm][pt].c_str() << endl;
	     cout << "x axis starts at: " << binLowX <<  endl;
	     cout << "x binwidth: " << binWidthX <<  endl;
	     cout << "x axis bin 1 center: " << binCenterX << endl;
	     cout << "x axis total bins: " << numBinsX << endl;
	     cout << " ______________________________ " << endl;
	     
	   }
       }
   }

  // make x_array using center of x axis bins
 
for(int i_histo = 0; i_histo < 4; i_histo++)
  {
    for(int arm = 0; arm < 2; arm++)
      {
	for(int pt = 0; pt < 4; pt++)
	  {
	    xaxis = bg[i_histo][arm][pt]->GetXaxis();  
	    for(int i = 0; i < 150; i++)
	      {
		x_array[i]  = xaxis->GetBinCenter(i+1);
		if(arm == 0)
		  {
		    cout << obj_filename[i_histo][arm][pt].c_str() << endl;
		    cout << "pt array: " << x_array[i] << " for bin " << i+1 << " in histo " << i_histo <<  endl;
		 }
	      }
	  }
      }
  }

  double pt;
  static const int bins = 150;
  double counter = 0.0;
  double scale[2] = {0.868*pow(10,15),1.78*pow(10,15)}; // for LS background 

  TH1D *h_sum[2][4]; // [arm][pt] (adding the 4 histograms togther of the same arm and pt)

  for(int arm = 0; arm < 2; arm++)
    { 
      for(int pt = 0; pt < 4; pt++)
      	{
	  //cout << " " << bg[0][arm][pt]->GetNbinsX() << endl;
	  h_sum[arm][pt] = (TH1D *) bg[0][arm][pt]->Clone(); // different histogram for each arm and pt
	  // h_sum[arm][pt]->Sumw2();
	  h_sum[arm][pt]->Add(bg[1][arm][pt]);
	  h_sum[arm][pt]->Add(bg[2][arm][pt]);
	  h_sum[arm][pt]->Add(bg[3][arm][pt]);
	  
	  bg[0][arm][pt]->Scale(scale[arm]);
	  bg[1][arm][pt]->Scale(scale[arm]);
	  bg[2][arm][pt]->Scale(scale[arm]);
	  bg[3][arm][pt]->Scale(scale[arm]);
	}
    } 

  //********************************************
  double cc[2][4][bins];
  double bb[2][4][bins];
  double dy[2][4][bins];
  double had[2][4][bins];
  double corr_sum[2][4][bins];
  double corr_total[2][bins];
  double corr_mass[2][bins];
  double errors[2][4][bins];
  double x_errors[2][bins];
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < 4; pt++)
	{
	  for(int i = 0; i < bins; i++)
	    {
	      cc[arm][pt][i] = bg[0][arm][pt]->GetBinContent(i+1); 
	      bb[arm][pt][i] = bg[1][arm][pt]->GetBinContent(i+1);
	      dy[arm][pt][i] = bg[2][arm][pt]->GetBinContent(i+1);
	      had[arm][pt][i] = bg[3][arm][pt]->GetBinContent(i+1);
	      corr_mass[arm][i] = t_check[arm]->GetBinContent(i+1);
	      // cout << " " << arm << " " << pt << " " << i << " " <<   "cc arm " << cc[arm][pt][i] << endl;
	      corr_sum[arm][pt][i] = cc[arm][pt][i] + bb[arm][pt][i] + dy[arm][pt][i] + had[arm][pt][i]; // pt dependent
	      //cout << " " << arm << " " << pt << " " << i << " " <<   "sum arm " << corr_sum[arm][pt][i] << endl;
	      corr_total[arm][i] = corr_sum[arm][0][i] + corr_sum[arm][1][i] + corr_sum[arm][2][i] + corr_sum[arm][3][i]; // pt integrated
	      //  cout << "Bin Content and errors:" << endl;
	      //  cout << "cc: " <<  cc[arm][pt][i] << ", and error: " <<  bg[0][arm][pt]->GetBinError(i+1) << ", bb: " <<  bb[arm][pt][i] << ", and error: " <<   bg[1][arm][pt]->GetBinError(i+1) <<  ", dy: " << dy[arm][pt][i] << ", and error: " <<   bg[2][arm][pt]->GetBinError(i+1) << ", had: "  << had[arm][pt][i] << ", and error: " <<   bg[3][arm][pt]->GetBinError(i+1) << endl;
	      errors[arm][0][i] = bg[0][arm][0]->GetBinError(i+1); // cc errors for initial slicing
	      errors[arm][1][i] = bg[1][arm][0]->GetBinError(i+1); // bb errors
	      errors[arm][2][i] = bg[2][arm][0]->GetBinError(i+1); // dy errors
	      errors[arm][3][i] = bg[3][arm][0]->GetBinContent(i+1)*0.10; // had errors
	      x_errors[arm][i] = 0.0;
	    }
	}
    }
  //********************************************	

      TGraphErrors *gr[12];
      // plot of individual components for south arm
      gr[0] = new TGraphErrors(150,x_array,cc[0][0],x_errors[0],errors[0][0]);  //[arm][pt]
      gr[1] = new TGraphErrors(150,x_array,bb[0][0],x_errors[0],errors[0][1]); 
      gr[2] = new TGraphErrors(150,x_array,dy[0][0],x_errors[0],errors[0][2]); 
      gr[3] = new TGraphErrors(150,x_array,had[0][0],x_errors[0],errors[0][3]); 

      // plot of individual components for north arm
      // gr[0] = new TGraphErrors(150,x_array,cc[1][0],x_errors[1],errors[1][0]);   // x_errors is just an empty array for both arms
      // gr[1] = new TGraphErrors(150,x_array,bb[1][0],x_errors[1],errors[1][1]); 
      // gr[2] = new TGraphErrors(150,x_array,dy[1][0],x_errors[1],errors[1][2]); 
      // gr[3] = new TGraphErrors(150,x_array,had[1][0],x_errors[1],errors[1][3]); 




      gr[4] = new TGraphErrors(150,x_array,corr_total[0]); 
      gr[5] = new TGraphErrors(150,x_array,corr_mass[0]); 
      
      gr[6] = new TGraphErrors(150,x_array,corr_sum[1][0]); 
      gr[7] = new TGraphErrors(150,x_array,corr_sum[1][1]); 
      gr[8] = new TGraphErrors(150,x_array,corr_sum[1][2]); 
      gr[9] = new TGraphErrors(150,x_array,corr_sum[1][3]); 
      gr[10] = new TGraphErrors(150,x_array,corr_total[1]); 
      gr[11] = new TGraphErrors(150,x_array,corr_mass[1]); 
      // gr[0] = new TGraphErrors(150,x_array,cc[0][0]); 
      
      
      TCanvas *c1 = new TCanvas("c1","0-1 GeV",200,10,700,500);
      gPad->SetLogy();
      gPad->SetGrid();
      
      gr[0]->SetTitle("Yue Hang Corr BG, South");
      gPad->SetLeftMargin(0.15); // 0.3
      gPad->SetBottomMargin(0.15);
      gr[0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      gr[0]->GetXaxis()->SetLabelSize(0.06);
      gr[0]->GetXaxis()->SetTitleSize(0.05);
      gr[0]->GetXaxis()->SetTitleOffset(1.00);
      gr[0]->GetYaxis()->SetTitle("");
      gr[0]->GetYaxis()->SetLabelSize(0.06);
      gr[0]->GetYaxis()->SetTitleSize(0.06);
      gr[0]->GetYaxis()->SetTitleOffset(0.99);
      gr[0]->GetYaxis()->SetRangeUser(pow(10,-1),pow(10,1));
      
      gr[0]->SetMarkerColor(kTeal);  
      gr[0]->SetMarkerStyle(20);
      gr[0]->SetMarkerSize(0.75);
      gr[1]->SetMarkerColor(kAzure);  
      gr[1]->SetMarkerStyle(20);
      gr[1]->SetMarkerSize(0.75);
      gr[2]->SetMarkerColor(kViolet);  
      gr[2]->SetMarkerStyle(20);
      gr[2]->SetMarkerSize(0.75);
      gr[3]->SetMarkerColor(kMagenta);  
      gr[3]->SetMarkerStyle(20);
      gr[3]->SetMarkerSize(0.75);

      // gr[0]->GetXaxis()->SetLimits(2,5);
      // gr[0]->GetYaxis()->SetRangeUser(0.0,pow(10,-8));
      gr[0]->Draw("AP"); 
      gr[1]->Draw("P");
      gr[2]->Draw("P");
      gr[3]->Draw("P");
   
  TLegend *leg = new TLegend(0.51, 0.6, 0.8, 0.9);   //0.11
  leg->SetFillColor(0); 
  leg->SetTextSize(0.025);
  leg->AddEntry(gr[0],"cc 0-1 GeV/c","p");
  leg->AddEntry(gr[1],"bb 0-1 GeV/c","p");
  leg->AddEntry(gr[2],"dy 0-1 GeV/c","p");
  leg->AddEntry(gr[3],"had 0-1 GeV/c","p");
  //leg->Draw();

  /////////// PLOT ///////////////
  
  
  gStyle->SetOptStat(0);
  
  //TCanvas *c4 = new TCanvas("c4","title",200,10,700,500);
  //gPad->SetLogy();
  gPad->SetGrid();

  h_sum[0][0]->SetTitle("Yue Hang Corr BG (cc+bb+dy+had), S");
  gPad->SetLeftMargin(0.15); // 0.3
  gPad->SetBottomMargin(0.15);
  h_sum[0][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  h_sum[0][0]->GetXaxis()->SetLabelSize(0.06);
  h_sum[0][0]->GetXaxis()->SetTitleSize(0.05);
  h_sum[0][0]->GetXaxis()->SetTitleOffset(1.00);
  //h_sum[0][0]->GetXaxis()->SetLimits(2,5);
  h_sum[0][0]->GetYaxis()->SetTitle("");
  h_sum[0][0]->GetYaxis()->SetLabelSize(0.06);
  h_sum[0][0]->GetYaxis()->SetTitleSize(0.06);
  h_sum[0][0]->GetYaxis()->SetTitleOffset(0.99);

  //  h_sum[0][0]->SetLineColor(kTeal-6);  
  // h_sum[0][0]->SetLineWidth(3.0);
  h_sum[0][0]->SetLineColor(kTeal);  
  h_sum[0][0]->SetLineWidth(3.0);
  h_sum[0][0]->SetMinimum(pow(10,-20));
  h_sum[0][1]->SetLineColor(kAzure);  
  h_sum[0][1]->SetLineWidth(3.0);
  h_sum[0][2]->SetLineColor(kViolet);  
  h_sum[0][2]->SetLineWidth(3.0);
  h_sum[0][3]->SetLineColor(kBlack);  
  h_sum[0][3]->SetLineWidth(3.0);
    
  // h_sum[0][0]->Draw(); 
  // // h_sum[0][0]->GetXaxis()->SetLimits(2,5);
  // h_sum[0][1]->Draw("SAME"); 
  // h_sum[0][2]->Draw("SAME");
  // h_sum[0][3]->Draw("SAME"); 
  
  //h_sum[0][0]->GetYaxis()->SetRangeUser(0,pow(10,3));


  TLegend *leg2 = new TLegend(0.51, 0.6, 0.8, 0.9);   //0.11
  leg2->SetFillColor(0); 
  leg2->SetTextSize(0.025);
  leg2->AddEntry(h_sum[0][0],"Corr bg 0-1 GeV", "l"); 
  leg2->AddEntry(h_sum[0][1],"Corr bg 1-2 GeV", "l"); 
  leg2->AddEntry(h_sum[0][2],"Corr bg 2-3 GeV", "l"); 
  leg2->AddEntry(h_sum[0][3],"Corr bg 3-7 GeV", "l");  
  // leg2->Draw();

} // void end macro
