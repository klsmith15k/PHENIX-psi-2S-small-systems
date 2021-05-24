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

void cc_rebinning_macro()
{


  bool write = true;

  TFile *file_cc_pp;
  TFile *file_cc_powheg;
  TFile *file_cc_suppression;

  TH2F *cc_pp[2] = {0,0}; // north and south arms
  TH2D *cc_0 = 0;
  TH2D *cc_dAu[2] = {0,0};  // north and south arms
  
  std::string filename_pp[2] = {"for_krista_mass_arm0_12_5.root",
				"for_krista_mass_arm1_12_5.root"};
  
  std::string filename_powheg = "powheg_cc_vf_leptons_14664_ccrpa0_vkrista.root";

  std::string filename_suppression[2] = {"powheg_cc_vf_leptons_14664_ccrpa2_vkrista.root",
					 "powheg_cc_vf_leptons_14664_ccrpa1_vkrista.root"};
  
  for(int arm = 0; arm < 2; arm++)
    {
      file_cc_pp = TFile::Open(filename_pp[arm].c_str()); 
      file_cc_pp->GetObject("h_cc_mass_pt",cc_pp[arm]);
      
      file_cc_powheg = TFile::Open(filename_powheg.c_str()); 
      file_cc_powheg->GetObject("h_11_mass_pt_mm_cc_FG12",cc_0); // unlike cc histograms
     
      file_cc_suppression = TFile::Open(filename_suppression[arm].c_str()); 
      file_cc_suppression->GetObject("h_11_mass_pt_mm_cc_FG12",cc_dAu[arm]);
      
    }
	  
  TAxis *xaxis;
  TAxis *yaxis;

  double bin_two;
  double bin_five;

  double x_array[100];
  double x_low[100];
  double x_high[100];

  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Int_t numBinsX;
 
  Double_t binCenterY;
  Double_t binLowY;
  Double_t binWidthY;
  Int_t numBinsY;

  cc_0->RebinX(10); // x axis has 1500 bins ---> 150 bins
  cc_0->RebinY(10); // y axis has 800 bins ---> 80 bins

  for(int arm = 0; arm < 2; arm++)
    {
      cc_dAu[arm]->RebinX(10);
      cc_dAu[arm]->RebinY(10);

      xaxis = cc_0->GetXaxis(); 
      binCenterX = xaxis->GetBinCenter(1);
      binLowX = xaxis->GetBinLowEdge(1);
      binWidthX = xaxis->GetBinWidth(1);
      numBinsX = xaxis->GetNbins();
	     
      bin_two = xaxis->FindBin(2.001);
      bin_five = xaxis->FindBin(4.999);
      cout << " ______________________________ " << endl;
      cout << "powheg mass bin 2.0: " << bin_two << ", and bin 5.0: " << bin_five << endl;
      cout << "x axis starts at: " << binLowX <<  endl;
      cout << "x binwidth: " << binWidthX <<  endl;
      cout << "x axis bin 1 center: " << binCenterX << endl;
      cout << "x axis total bins: " << numBinsX << endl;
      cout << " ______________________________ " << endl;

      yaxis = cc_dAu[arm]->GetYaxis();
      binCenterY = yaxis->GetBinCenter(1);
      binLowY = yaxis->GetBinLowEdge(1);
      binWidthY = yaxis->GetBinWidth(1);
      numBinsY = yaxis->GetNbins();

      cout << " ______________________________ " << endl;
      cout << "RdAu suppression mass bin 2.0: " << bin_two << ", and bin 5.0: " << bin_five << endl;
      cout << "y axis starts at: " << binLowY <<  endl;
      cout << "y binwidth: " << binWidthY <<  endl;
      cout << "y axis bin 1 center: " << binCenterY << endl;
      cout << "y axis total bins: " << numBinsY << endl;
      cout << " ______________________________ " << endl;
	      
    }
        
  xaxis = cc_0->GetXaxis();  
 
  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
      cout << "x array: " << x_array[i] << endl;
    }


  /// now the 1500 x 800 TH2D histograms have been correctly rebinned to have 150 x 80 bins
  // start projection now in 500 mev bins


  TFile *fout; 

  std::string rebinned_filename[2] = {"yuehang_Run15pAu_S_rebinned_300mev_cc.root","yuehang_Run15pAu_N_rebinned_300mev_cc.root"};
 
  char unique1[500];
  char unique2[500];
  char unique3[500];


  //////////////////////////////////////// Rebinned pT binwidths (200 MeV)

  // int pt_slices = 38;  
  // double pt_width[38] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};  
  // double pt_center[38] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};

 //////////////////////////////////////// Rebinned pT binwidths (0.3 GeV)

  int pt_slices = 21;  // the pT binwidth in Yue Hang's original TH2D is 0.1 GeV/c
  double pt_width[21] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,1.0};
  double pt_center[21] = {0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.25,5.55,5.85,6.5};


 ///////////////////////////////////////////// Rebinned pT binwidths (500 MeV)
 
  // int pt_slices = 20;  
  // double pt_width[20] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  // double pt_center[20] = {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};

  ///////////////////////////////////////////// Rebinned pT binwidths (1 GeV)
 
  // int pt_slices = 10;  
  // double pt_width[10] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  // double pt_center[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};

 ///////////////////////////////////////////// Rebinned pT binwidths (2 GeV)
 
  // int pt_slices = 5;  
  // double pt_width[5] = {1,2,2,2,2};
  // double pt_center[5] = {0.5,2.0,4.0,6.0,8.0};

///////////////////////////////////////////// Rebinned pT binwidths (3 GeV)
 
  // int pt_slices = 3;  
  // double pt_width[3] = {1,3,3};
  // double pt_center[3] = {0.5,2.5,5.5};




  for(int i = 0; i < 1; i++)
    {
      for(int j = 0; j < 21; j++)
	{
	  cout << " pt center: " << pt_center[j] << ", for j = " << j << endl;
	}
    }



 
  double pt_low;
  double pt_high;

  double bin_low;
  double bin_high;

  double delta_pt = 0.001;

  TH1D *cc_1[2][pt_slices];  
  TH1D *cc_2[pt_slices];
  TH1D *cc_3[2][pt_slices];

  double cc_entries_pp; 
  double cc_entries_powheg;
  double cc_entries_suppression;
 
  double sum1[2] = {0,0};
  double sum2 = 0;
  double sum3[2] = {0,0};
  
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
	  
	  double a = (pt_low - delta_pt)*1000; // convert to MeV/c for unique labeling 
	  double b = (pt_high + delta_pt)*1000;

	  sprintf(unique1,"cc_pp%.0f_%.0f",a,b); 
	  sprintf(unique2,"cc_powheg_%.0f_%.0f",a,b); 
	  sprintf(unique3,"cc_suppression_%.0f_%.0f",a,b);  
	  
	  cc_1[arm][k] =  (TH1D *) cc_pp[arm]->ProjectionX(unique1,bin_low,bin_high); 
	  cc_2[k] =  (TH1D *) cc_0->ProjectionX(unique2,bin_low,bin_high); 
	  cc_3[arm][k] =  (TH1D *) cc_dAu[arm]->ProjectionX(unique3,bin_low,bin_high); 
		  
	  fout->cd();  // indicates we do not want to use the most recent file
	  
	  if(write == true)
	    {
	      cc_1[arm][k]->Write(); 
	      cc_2[k]->Write(); 
	      cc_3[arm][k]->Write(); 
	 
	    }

	  sum1[arm]+= cc_1[arm][k]->GetBinContent(k+1);
	  cout << "cc_pp count: " << cc_1[arm][k]->GetBinContent(k+1) << ", and bin error: " << cc_1[arm][k]->GetBinError(k+1) << endl;

	  sum2+= cc_2[k]->GetBinContent(k+1);
	  cout << "cc_powheg count: " << cc_2[k]->GetBinContent(k+1) << ", and bin error: " << cc_2[k]->GetBinError(k+1) << endl;

	  sum3[arm]+= cc_3[arm][k]->GetBinContent(k+1);
	  cout << "cc suppression count: " << cc_3[arm][k]->GetBinContent(k+1) << ", and bin error: " << cc_3[arm][k]->GetBinError(k+1) << endl;

	}// end for loop k

      cout << "_________________________ " << endl;
      cout << "     SUMMARY " << endl;
      cout << "Arm " << arm << ", total cc from pp entries: " << sum1[arm] << endl;
      cout << "Arm " << arm << ", total cc from powheg entries: " << sum2 << endl;
      cout << "Arm " << arm << ", total cc from suppression entries: " << sum3[arm] << endl;
      cout << "_________________________ " << endl;

    } // end for loop arm
  
  fout->Close();
  

} // void end macro
