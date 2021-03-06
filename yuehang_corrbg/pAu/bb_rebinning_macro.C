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

void bb_rebinning_macro()
{

  bool write = true;

  TFile *file_bb_pp;
  TFile *file_bb_powheg;
  TFile *file_bb_suppression;

  // Run15 pp 
  TH2F *bb_UL_run15pp[2] = {0,0}; // north and south arms
  TH2F *bb_LS_run15pp[2] = {0,0}; // north and south arms

  // powheg unmodified
  TH2D *bb_UL_0 = 0; // powheg unmodified bb  UL 
  TH2D *bb_LS_pp = 0; // powheg unmodified bb LS pp
  TH2D *bb_LS_mm = 0; // powheg unmodified bb LS mm

  // powheg modified
  TH2D *bb_dAu_UL = 0;  // powheg modified UL for north and south arms
  TH2D *bb_dAu_LS_mm = 0;  // powheg modified LS pp for north and south arms
  TH2D *bb_dAu_LS_pp = 0;  // powheg modified LS mm for north and south arms
   
  std::string filename_pp[2] = {"for_krista_mass_arm0_12_5.root",
				"for_krista_mass_arm1_12_5.root"};
 
  std::string filename_powheg = "powheg_bbparent_vrpa_leptons_5000_ccrpa0_vkrista.root";

  std::string filename_suppression = "powheg_bbparent_vrpa_leptons_5000_ccrpa30_vf.root";
					 
				    					
  for(int arm = 0; arm < 2; arm++)
    {
      file_bb_pp = TFile::Open(filename_pp[arm].c_str()); 
      file_bb_pp->GetObject("h_bb_mass_pt",bb_UL_run15pp[arm]); // unlike bb from run15pp
      file_bb_pp->GetObject("h_bb_FGLS_mass_pt",bb_LS_run15pp[arm]); // like bb from run15pp
      
      // the same file is used for north and south arms
      file_bb_powheg = TFile::Open(filename_powheg.c_str()); 
      file_bb_powheg->GetObject("h_11_mass_pt_mm_bb_FG12",bb_UL_0); // unlike sign bb from powheg_0
      file_bb_powheg->GetObject("h_11_mass_pt_mm_bb_FG22",bb_LS_pp); // like sign bb from powheg_0
      file_bb_powheg->GetObject("h_11_mass_pt_mm_bb_FG11",bb_LS_mm); // like sign bb from powheg_0
      
      file_bb_suppression = TFile::Open(filename_suppression.c_str()); 
      file_bb_suppression->GetObject("h_11_mass_pt_mm_bb_FG12",bb_dAu_UL); // UL sign modifications from powheg
      file_bb_suppression->GetObject("h_11_mass_pt_mm_bb_FG22",bb_dAu_LS_pp); // UL sign modifications from powheg
      file_bb_suppression->GetObject("h_11_mass_pt_mm_bb_FG11",bb_dAu_LS_mm); // UL sign modifications from powheg
      
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

  // all powheg unmodified files
  bb_UL_0->RebinX(10); // x axis has 1500 bins ---> 150 bins
  bb_UL_0->RebinY(10); // y axis has 800 bins ---> 80 bins

  bb_LS_mm->RebinX(10); 
  bb_LS_mm->RebinY(10);

  bb_LS_pp->RebinX(10); 
  bb_LS_pp->RebinY(10); 

  // all powheg modified files
  bb_dAu_UL->RebinX(10);
  bb_dAu_UL->RebinY(10);

  bb_dAu_LS_pp->RebinX(10);
  bb_dAu_LS_pp->RebinY(10);

  bb_dAu_LS_mm->RebinX(10);
  bb_dAu_LS_mm->RebinY(10);
  
  for(int arm = 0; arm < 2; arm++)
    {
  
      xaxis = bb_UL_0->GetXaxis(); 
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

      yaxis = bb_dAu_UL->GetYaxis();
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
        
  xaxis = bb_UL_0->GetXaxis();  
 
  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
      cout << "x array: " << x_array[i] << endl;
    }


  /// now the 1500 x 800 TH2D histograms have been correctly rebinned to have 150 x 80 bins
  // start projection now in 500 mev bins


  TFile *fout; 

  std::string rebinned_filename[2] = {"yuehang_Run15pAu_S_rebinned_300mev_bb.root","yuehang_Run15pAu_N_rebinned_300mev_bb.root"};
 
  char unique1[500];
  char unique2[500];
  char unique3[500];
  char unique4[500];
  char unique5[500];
  char unique6[500];
  char unique7[500];
  char unique8[500];
  
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
 
  double pt_low;
  double pt_high;

  double bin_low;
  double bin_high;

  double delta_pt = 0.001;

  TH1D *bb_1[2][pt_slices];  // bb UL from Run15pp
  TH1D *bb_2[pt_slices]; 
  TH1D *bb_3[pt_slices]; 
  TH1D *bb_4[pt_slices];
  TH1D *bb_5[pt_slices];
  TH1D *bb_6[2][pt_slices]; // bb LS from Run15pp
  TH1D *bb_7[pt_slices]; 
  TH1D *bb_8[pt_slices];
  
  double bb_entries_pp; 
  double bb_entries_powheg;
  double bb_entries_suppression;
 
  double sum1[2] = {0,0};
  double sum2 = 0;
  double sum3 = 0;
  
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

	  sprintf(unique1,"bb_UL_run15pp_%.0f_%.0f",a,b); 
	  sprintf(unique2,"bb_powheg_%.0f_%.0f",a,b); 
	  sprintf(unique3,"bb_UL_suppression_%.0f_%.0f",a,b);  
	  sprintf(unique4,"bb_LS_pp_%.0f_%.0f",a,b);  
	  sprintf(unique5,"bb_LS_mm_%.0f_%.0f",a,b);  
	  sprintf(unique6,"bb_LS_run15pp_%.0f_%.0f",a,b);  
	  sprintf(unique7,"bb_LS_pp_suppression_%.0f_%.0f",a,b);  
	  sprintf(unique8,"bb_LS_mm_suppression_%.0f_%.0f",a,b);  
	  
	  bb_1[arm][k] =  (TH1D *) bb_UL_run15pp[arm]->ProjectionX(unique1,bin_low,bin_high); 
	  bb_2[k] =  (TH1D *) bb_UL_0->ProjectionX(unique2,bin_low,bin_high); 
	  bb_3[k] =  (TH1D *) bb_dAu_UL->ProjectionX(unique3,bin_low,bin_high); 
	  bb_4[k] =  (TH1D *) bb_LS_pp->ProjectionX(unique4,bin_low,bin_high); 
	  bb_5[k] =  (TH1D *) bb_LS_mm->ProjectionX(unique5,bin_low,bin_high); 
	  bb_6[arm][k] =  (TH1D *) bb_LS_run15pp[arm]->ProjectionX(unique6,bin_low,bin_high); 
	  bb_7[k] =  (TH1D *) bb_dAu_LS_pp->ProjectionX(unique7,bin_low,bin_high); 
	  bb_8[k] =  (TH1D *) bb_dAu_LS_mm->ProjectionX(unique8,bin_low,bin_high); 
	  
	  fout->cd();  // indicates we do not want to use the most recent file
	  
	  if(write == true)
	    {
	      bb_1[arm][k]->Write(); 
	      bb_2[k]->Write(); 
	      bb_3[k]->Write(); 
	      bb_4[k]->Write(); 
	      bb_5[k]->Write(); 
	      bb_6[arm][k]->Write(); 
	      bb_7[k]->Write(); 
	      bb_8[k]->Write(); 
	    }

	  sum1[arm]+= bb_1[arm][k]->GetBinContent(k+1);
	  cout << "bb_pp count: " << bb_1[arm][k]->GetBinContent(k+1) << ", and bin error: " << bb_1[arm][k]->GetBinError(k+1) << endl;

	  sum2+= bb_2[k]->GetBinContent(k+1);
	  cout << "bb_powheg count: " << bb_2[k]->GetBinContent(k+1) << ", and bin error: " << bb_2[k]->GetBinError(k+1) << endl;

	  sum3+= bb_3[k]->GetBinContent(k+1);
	  cout << "bb suppression count: " << bb_3[k]->GetBinContent(k+1) << ", and bin error: " << bb_3[k]->GetBinError(k+1) << endl;

	}// end for loop k

      cout << "_________________________ " << endl;
      cout << "     SUMMARY " << endl;
      cout << "Arm " << arm << ", total bb from pp entries: " << sum1[arm] << endl;
      cout << "Arm " << arm << ", total bb from powheg entries: " << sum2 << endl;
      cout << "Arm " << arm << ", total bb from suppression entries: " << sum3 << endl;
      cout << "_________________________ " << endl;

    } // end for loop arm
  
  fout->Close();
  

} // void end macro
