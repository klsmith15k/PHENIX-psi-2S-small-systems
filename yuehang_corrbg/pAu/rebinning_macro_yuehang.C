
// need to change name of output file to reflect binwidth of the first bin (change at line 50)
// also need to select the corresponding pt array (binwidth, number of bins and center of bin)

#include <TF1.h>
#include <TMath.h>
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

using namespace std;

void rebinning_macro_yuehang()
{
  
  bool write = false;
 
  TFile *fout; 
  TFile *file_data;
 
  TAxis *xaxis;
  TAxis *yaxis;	  
  

  std::string anataxi_filename[2] = {"for_krista_mass_arm0_12_5.root",   // these have been rescaled and rebinned
				     "for_krista_mass_arm1_12_5.root"};
  
  std::string rebinned_filename[2] = {"yuehang_Run15pp_S_rebinned_2_gev.root","yuehang_Run15pp_N_rebinned_2_gev.root"};

  for(int arm = 0; arm < 2; arm++)
    {
      cout << rebinned_filename[arm] << endl;
      
    }

  double sum1[2] = {0,0};
  double sum2[2] = {0,0};
  double sum3[2] = {0,0};
  double sum4[2] = {0,0};
  
  char unique1[500];
  char unique2[500];
  char unique3[500];
  char unique4[500];
  char unique5[500];
  char unique6[500];

  TH2F *t[2][6];

  file_data = TFile::Open(anataxi_filename[1].c_str()); 	 
  file_data->GetObject("h_cc_mass_pt",t[1][0]);
  file_data->GetObject("h_bb_mass_pt",t[1][1]);
  file_data->GetObject("h_dy_mass_pt",t[1][2]); 
  //file_data->GetObject("h_sim_mass_pt_FG12_jet_tot",t[1][3]); // correlated hadrons are for UL sign only
  file_data->GetObject("h_bb_FGLS_mass_pt",t[1][3]);
  file_data->GetObject("h_sim_mass_pt_FGLS_jet_tot",t[1][4]); // LS
  file_data->GetObject("h_sim_mass_pt_FG12_jet_tot",t[1][5]); // UL

  file_data = TFile::Open(anataxi_filename[0].c_str()); 	 
  file_data->GetObject("h_cc_mass_pt",t[0][0]);
  file_data->GetObject("h_bb_mass_pt",t[0][1]);
  file_data->GetObject("h_dy_mass_pt",t[0][2]);
  //file_data->GetObject("h_sim_mass_pt_FG12_jet_tot",t[0][3]); // correlated hadrons are for UL sign only
  file_data->GetObject("h_bb_FGLS_mass_pt",t[0][3]);
  file_data->GetObject("h_sim_mass_pt_FGLS_jet_tot",t[0][4]);
  file_data->GetObject("h_sim_mass_pt_FG12_jet_tot",t[0][5]);

  //double scale[2] = {2.415*pow(10,11),2.3*pow(10,11)}; // try not using scale factors
  double scale[2] = {1,1};
  
  for(int i = 0; i < 2;i++)
    {
      for(int j = 0; j < 6;j++)
	{
	  t[i][j]->Scale(scale[i]);

	}
    }

    // Original pT binwidths (1 GeV)

  // int pt_slices = 4;
  // double pt_width[4] = {1.0,1.0,1.0,7.0}; // in GeV/c
  // double pt_center[4] = {0.5,1.5,2.5,6.5}; // in GeV/c



  //////////////////////////////////////// Rebinned pT binwidths (0.1 GeV)

  // int pt_slices = 10;  // the pT binwidth in Yue Hang's original TH2D is 0.1 GeV/c
  // double pt_width[10] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  // double pt_center[10] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95};
 
  //////////////////////////////////////// Rebinned pT binwidths (0.2 GeV)

  // int pt_slices = 38;  // the pT binwidth in Yue Hang's original TH2D is 0.1 GeV/c
  // double pt_width[38] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};  // 0.0 - 10.0 GeV/c have binwidths of 0.2 geV (Matt's has binwidths 0.25 GeV/c) to 6.0 GeV and then 0.5 GeV/c to 10.0 Gev
  // double pt_center[38] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
 
  //////////////////////////////////////// Rebinned pT binwidths (0.3 GeV)

  // int pt_slices = 21;  // the pT binwidth in Yue Hang's original TH2D is 0.1 GeV/c
  // double pt_width[21] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.1};
  // double pt_center[21] = {0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.25,5.55,5.85,6.5};

  ///////////////////////////////////////////////  Rebinned pT binwidths (0.5 GeV)

  // int pt_slices = 20;  
  // double pt_width[20] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  // double pt_center[20] = {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
 
  ///////////////////////////////////////////// Rebinned pT binwidths (1 GeV)
 
  // int pt_slices = 10;  
  // double pt_width[10] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  // double pt_center[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};

 ///////////////////////////////////////////// Rebinned pT binwidths (2 GeV)
 
  int pt_slices = 5;  
  double pt_width[5] = {1,2,2,2,2};
  double pt_center[5] = {0.5,2.0,4.0,6.0,8.0};

  // double mass[200];
  
  double pt_low;
  double pt_high;

  double bin_low;
  double bin_high;

  double delta_pt = 0.001;

  Double_t binCenterX;
  Double_t binCenterY;
  Double_t binLowX;
  Double_t binLowY;
  Double_t binWidthX;
  Double_t binWidthY;
  Double_t numBinsX;
  Double_t numBinsY;

  // // Fill x axis mass array
  // for(int i = 0; i < 200; i ++)
  //   {
  //     mass[i] = 0.05*(i+1);
  //     cout << "mass array: " << mass[i] << endl;
  //   }

  // Extract information from TH2F
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 6; i_histo++)
	{
	  xaxis = t[arm][i_histo]->GetXaxis();  
	  yaxis = t[arm][i_histo]->GetYaxis();  
	
	  binCenterX = xaxis->GetBinCenter(1);
	  binCenterY = yaxis->GetBinCenter(1);
	  binLowX = xaxis->GetBinLowEdge(1);
	  binLowY = yaxis->GetBinLowEdge(1);
	  binWidthX = xaxis->GetBinWidth(1);
	  binWidthY = yaxis->GetBinWidth(1);
	  numBinsX = xaxis->GetNbins();
	  numBinsY = yaxis->GetNbins();
	  
	  cout << " ______________________________ " << endl;
	  cout << "t[" << arm << "][" << i_histo << "]" << endl;
	  cout << "x axis starts at: " << binLowX << ", y axis starts at: " << binLowY << endl;
	  cout << "x binwidth: " << binWidthX << ", y binwidth: " << binWidthY << endl;
	  cout << "x axis bin 1 center: " << binCenterX << ", y axis bin 1 center: " << binCenterY << endl;
	  cout << "x axis total bins: " << numBinsX << ", y axis total bins: " << numBinsY << endl;
	  cout << " ______________________________ " << endl;
	
	}
    }
 
  int bins = numBinsY;

  TH1D *cc[2][pt_slices];  
  TH1D *bb_UL[2][pt_slices];
  TH1D *dy[2][pt_slices];
  TH1D *bb_LS[2][pt_slices];
  TH1D *corr_had_LS[2][pt_slices];
  TH1D *corr_had_UL[2][pt_slices];

  double cc_entries; 
  double bb_UL_entries;
  double dy_entries;
  double bb_LS_entries;
  double corr_had_LS_entries;
  double corr_had_US_entries;
  
  // Project onto X axis
  for(int arm = 0; arm < 2; arm++)
    {
      fout = new TFile(rebinned_filename[arm].c_str(), "RECREATE"); //resulting root file with name specified in array above 
 
      for(int k = 0; k < pt_slices; k++)
	{ 
	  pt_low = pt_center[k] - pt_width[k]/2 + delta_pt; 
	  pt_high = pt_center[k] + pt_width[k]/2 - delta_pt;
	 
	  bin_low = yaxis->FindBin(pt_low);
	  bin_high = yaxis->FindBin(pt_high);
	  
	  cout << "bin low: " << bin_low << ", bin high: " << bin_high << endl;


	  double a = (pt_low - delta_pt)*1000; // convert to MeV/c for unique labeling 
	  double b = (pt_high + delta_pt)*1000;

	  sprintf(unique1,"cc_%.0f_%.0f",a,b); 
	  sprintf(unique2,"bb_UL_%.0f_%.0f",a,b); 
	  sprintf(unique3,"dy_%.0f_%.0f",a,b);  
	  sprintf(unique4,"bb_LS_%.0f_%.0f",a,b); 
	  sprintf(unique5,"corr_had_LS_%.0f_%.0f",a,b);  
	  sprintf(unique6,"corr_had_UL_%.0f_%.0f",a,b); 
 

	  // cc[arm][k] = (TH1F *) t[arm][0]->ProjectionX(unique1,bin_low,bin_high); 
	  // bb_UL[arm][k] = (TH1F *) t[arm][1]->ProjectionX(unique2,bin_low,bin_high); 
	  // dy[arm][k] = (TH1F *) t[arm][2]->ProjectionX(unique3,bin_low,bin_high); 
	  // bb_LS[arm][k] = (TH1F *) t[arm][3]->ProjectionX(unique4,bin_low,bin_high); 
	  
	  cc[arm][k] =  (TH1D *) t[arm][0]->ProjectionX(unique1,bin_low,bin_high); 
	  bb_UL[arm][k] =  (TH1D *) t[arm][1]->ProjectionX(unique2,bin_low,bin_high); 
	  dy[arm][k] =  (TH1D *) t[arm][2]->ProjectionX(unique3,bin_low,bin_high); 
	  bb_LS[arm][k] =  (TH1D *) t[arm][3]->ProjectionX(unique4,bin_low,bin_high); 
	  corr_had_LS[arm][k] =  (TH1D *) t[arm][4]->ProjectionX(unique5,bin_low,bin_high); 
	  corr_had_UL[arm][k] =  (TH1D *) t[arm][5]->ProjectionX(unique6,bin_low,bin_high); 
	  
	  fout->cd();  // indicates we do not want to use the most recent file
	  
	  if(write == true)
	    {
	      cc[arm][k]->Write(); 
	      bb_UL[arm][k]->Write(); 
	      dy[arm][k]->Write(); 
	      bb_LS[arm][k]->Write();
	      corr_had_LS[arm][k]->Write();
	      corr_had_UL[arm][k]->Write();
	    }

	  sum1[arm]+= cc[arm][k]->GetBinContent(k+1);
	  cout << "cc count: " << cc[arm][k]->GetBinContent(k+1) << ", and bin error: " << cc[arm][k]->GetBinError(k+1) << endl;

	  sum2[arm]+= bb_UL[arm][k]->GetBinContent(k+1);
	  cout << "bb count: " << bb_UL[arm][k]->GetBinContent(k+1) << ", and bin error: " << bb_UL[arm][k]->GetBinError(k+1) << endl;

	  sum3[arm]+= dy[arm][k]->GetBinContent(k+1);
	  cout << "dy count: " << dy[arm][k]->GetBinContent(k+1) << ", and bin error: " << dy[arm][k]->GetBinError(k+1) << endl;

	  sum4[arm]+= bb_LS[arm][k]->GetBinContent(k+1);
	  cout << "had count: " << bb_LS[arm][k]->GetBinContent(k+1) << ", and bin error: " << bb_LS[arm][k]->GetBinError(k+1) << endl;
	  
	}// end for loop k
      cout << "_________________________ " << endl;
      cout << "     SUMMARY " << endl;
      cout << "Arm " << arm << ", total cc entries: " << sum1[arm] << endl;
      cout << "Arm " << arm << ", total bb_UL entries: " << sum2[arm] << endl;
      cout << "Arm " << arm << ", total dy entries: " << sum3[arm] << endl;
      cout << "Arm " << arm << ", total bb_LS entries: " << sum4[arm] << endl;
      cout << "_________________________ " << endl;

    } // end for loop arm
  
  fout->Close();
  
} // end void macro

//*/












