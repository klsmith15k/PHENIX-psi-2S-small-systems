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

void yuehang_bg_macro_300mev_refit()
{

#include "yuehang_fit_coeff_300/fit_coeff_S_array_nineteen_point.C"   
#include "yuehang_fit_coeff_300/fit_coeff_N_array_nineteen_point.C"   

  bool all_fits = true; // set to true for ratio plots
  bool write = false;

  //bool south_arm = true;
  bool south_arm = false;
   bool pt_int_fit = false;

  bool initial_fit = false;  // set to false for ratio plots
  bool refit = true;  // set to true for ratio plots
  bool normalization = false;
  
  double weight[2] = {1.03,0.985};
			 
  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; 

  int pt_binwidth = 3;

  // cout << "What pt binning are you running?  Enter '0' for 100 MeV binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 300 MeV binwidth, '4' for 1 GeV binwidth" << endl;
  // cin >> pt_binwidth;

  int pt_slices[5] = {10,6,10,21,10}; 
 
  std::string filename_comp[2][2] = {"yuehang_Run15pAu_S_rebinned_300mev_cc.root",
  				     "yuehang_Run15pAu_N_rebinned_300mev_cc.root",

				     "yuehang_Run15pAu_S_rebinned_300mev_bb.root",
  				     "yuehang_Run15pAu_N_rebinned_300mev_bb.root"};

std::string filename_corrhad[2] = {"yuehang_Run15pAu_S_rebinned_corr_had.root",
  				     "yuehang_Run15pAu_N_rebinned_corr_had.root"};

std::string filename_Run15pp[2] = {"yuehang_Run15pp_S_rebinned_300mev.root",
				  "yuehang_Run15pp_N_rebinned_300mev.root"};
   
 double pt_width[5][38] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
			   0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
 
 double pt_center[5][38] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.25,5.55,5.85,6.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

 double pt_par_plot[21] = {0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.25,5.55,5.85,6.5};

  TFile *file_comp[2][2];
  TFile *file_corrhad[2];
  TFile *file_Run15pp[2];
  
  TH1D *bg[11][2][38] = {0};
  TH1D *corr_had[2][2][38] = {0}; 
  TH1D *pp[3][2][38] = {0}; 



  for(int i = 0; i < 5; i++)
    {
      for(int j = 0; j < 38; j++)
	{
	  cout << " pt center: " << pt_center[i][j] << ", for i = " <<  i << ", j = " << j << endl;
	}
    }


 
  std::string obj_filename[11][2][38];  //[number of files][arm][pt slices]
  std::string obj_filename_corrhad[2][2][4];  
  std::string obj_filename_pp[3][2][38]; 

  char unique0[800];
  char unique1[800];
  char unique2[800];
  char unique3[800];
  char unique4[800];
  char unique5[800];
  char unique6[800];
  char unique7[800];
  char unique8[800];
  char unique9[800];
  char unique10[800];
  char unique11[800]; 
  char unique12[800];
  char unique13[800];
  char unique14[800];
  char unique15[800];
  
  int bins = 100;
  
  double pt_low;
  double pt_high;
  
  double bin_low;
  double bin_high;
  
  double delta_pt = 0.001;

  double a_array[38];
  double b_array[38];

  int corr_a;
  int corr_b;
  
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{ 
	  pt_low = pt_center[pt_binwidth][pt] - pt_width[pt_binwidth][pt]/2 + delta_pt; 
	  pt_high = pt_center[pt_binwidth][pt] + pt_width[pt_binwidth][pt]/2 - delta_pt;
  
	  double a = (pt_low - delta_pt)*1000; 
	  double b = (pt_high + delta_pt)*1000;

	  sprintf(unique0,"cc_pp%.0f_%.0f",a,b);                  // cc unmodified(run15pp) is [0], which is included in the rebinned pAu root file
	  sprintf(unique1,"cc_powheg_%.0f_%.0f",a,b);        // pp unmodified is [1]
	  sprintf(unique2,"cc_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) is [2]  
		  
	  sprintf(unique3,"bb_UL_run15pp_%.0f_%.0f",a,b);                  // bb unmodified (run15pp) 
	  sprintf(unique4,"bb_powheg_%.0f_%.0f",a,b);        // pp unmodified 
	  sprintf(unique5,"bb_UL_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) 
	  sprintf(unique6,"bb_LS_pp_%.0f_%.0f",a,b);      // pp like sign unmodified for bb plus plus
	  sprintf(unique7,"bb_LS_mm_%.0f_%.0f",a,b);     // pp like sign unmodified for bb minus minus
	  sprintf(unique8,"bb_LS_run15pp_%.0f_%.0f",a,b);     // pp like sign unmodified for bb minus minus
	  sprintf(unique9,"bb_LS_pp_suppression_%.0f_%.0f",a,b);     // bb(LS) modified 
	  sprintf(unique10,"bb_LS_mm_suppression_%.0f_%.0f",a,b);     // bb(LS) modified
	      
	  sprintf(unique11,"dy_%.0f_%.0f",a,b);     // run15pp components
	  sprintf(unique12,"corr_had_UL_%.0f_%.0f",a,b);      // run15pp components
	  sprintf(unique13,"corr_had_LS_%.0f_%.0f",a,b);       // run15pp components
	  
	  a_array[pt] = a/1000;
	  b_array[pt] = b/1000;

	  obj_filename[0][arm][pt] = unique0;
	  obj_filename[1][arm][pt] = unique1;
	  obj_filename[2][arm][pt] = unique2;
	  obj_filename[3][arm][pt] = unique3;
	  obj_filename[4][arm][pt] = unique4;
	  obj_filename[5][arm][pt] = unique5;
	  obj_filename[6][arm][pt] = unique6;
	  obj_filename[7][arm][pt] = unique7;
	  obj_filename[8][arm][pt] = unique8;
	  obj_filename[9][arm][pt] = unique9;
	  obj_filename[10][arm][pt] = unique10;

	  obj_filename_pp[0][arm][pt] = unique11;
	  obj_filename_pp[1][arm][pt] = unique12;
	  obj_filename_pp[2][arm][pt] = unique13;
	  
	}
    }
  
  for(int i_histo = 0; i_histo < 11; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      //cout << "obj 1: " << obj_filename[i_histo][arm][pt].c_str() << " " << endl;
	      //cout << pt << endl;
	    }
	}
    }

  for(int i_histo = 0; i_histo < 2; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {

	      if(pt < 3)
		{
		  corr_a = pt;
		  corr_b = pt+1;
		  
		  sprintf(unique14,"corr_had_UL_%i_%i",corr_a,corr_b); // modified
		  sprintf(unique15,"corr_had_LS_%i_%i",corr_a,corr_b); // modified
		}
	      if(pt == 3)
		{
		  corr_a = pt;
		  corr_b = pt+5;
		  
		  sprintf(unique14,"corr_had_UL_%i_%i",corr_a,corr_b);  // modifed
		  sprintf(unique15,"corr_had_LS_%i_%i",corr_a,corr_b); // modified
		}
	      
	      obj_filename_corrhad[0][arm][pt] = unique14; 
	      obj_filename_corrhad[1][arm][pt] = unique15;
	      //cout << "obj 2: " << obj_filename_corrhad[i_histo][arm][pt].c_str() << " " << endl;
	    }
	}
    }
  
  TH1 *read_this = 0;
   
  // Fill bg TH1D array
  for(int arm = 0; arm < 2; arm++)      
    {
      file_comp[0][arm] = TFile::Open(filename_comp[0][arm].c_str()); 
      file_comp[1][arm] = TFile::Open(filename_comp[1][arm].c_str()); 

      for(int i_histo = 0; i_histo < 11; i_histo++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      if(i_histo < 3)
		file_comp[0][arm]->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]); // cc UL 
	     
	      if(i_histo > 2) 
		file_comp[1][arm]->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]); // bb series

	      if(!bg[i_histo][arm][pt])
		{
		  cout << "Get object failed for i_hist " << i_histo <<  " arm " << arm << " pt " << pt << " quit!" << endl;
		  exit(1);
		}
	      
	    }
	}
      // file_comp->Close();
    }
 
  ///////////////////////////////////////////////////////// check files were retrieved ////////////////////////////////////////////////////////

  // i_histo == 0 is cc 
 file_comp[0][0]->GetObject("cc_powheg_0_300", read_this);
 cout << "read_this : " << read_this << endl;
 read_this = 0;

 // i_histo == 1 is bb
 file_comp[1][0]->GetObject("bb_powheg_0_300", read_this);
 cout << "read_this : " << read_this << endl;

 ///////////////////////////////////////////////////////// check files were retrieved ////////////////////////////////////////////////////////


  // Fill bg TH1D array
  
  for(int arm = 0; arm < 2; arm++)
    {
      file_Run15pp[arm] = TFile::Open(filename_Run15pp[arm].c_str()); 

      for(int i_histo = 0; i_histo < 3; i_histo++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      file_Run15pp[arm]->GetObject(obj_filename_pp[i_histo][arm][pt].c_str(),pp[i_histo][arm][pt]); // pp dy(UL),hd(UL),had(LS)
	    }
	}
      // file_Run15pp->Close();
    }
  
  // Fill corr had TH1D array
   for(int arm = 0; arm < 2; arm++)
     {
       file_corrhad[arm] = TFile::Open(filename_corrhad[arm].c_str()); 

       for(int i_histo = 0; i_histo < 2; i_histo++)
	 {
	   for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      if(pt < 3)
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][0].c_str(),corr_had[i_histo][arm][pt]); // first 3 elements are 0-1 GeV/c [0]
	      if( (pt > 2) && (pt < 6))
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][1].c_str(),corr_had[i_histo][arm][pt]); // next 3 elements are 1-2 GeV/c [1]
	      if( (pt > 5) && (pt < 10))
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][2].c_str(),corr_had[i_histo][arm][pt]); // next 4 elements are 2-3 GeV/c [2]
	      if( pt > 9)
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][3].c_str(),corr_had[i_histo][arm][pt]); // remaining elements are 3-8 GeV/c [3]
	    }
	}
       // file_corrhad->Close();
    }	  

  TAxis *xaxis;
  TAxis *yaxis;
  
  double bin_two;
  double bin_five;
  
  double x_array[100];
  double x_array_mass_ratio[100];
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
 
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      xaxis = bg[6][0][0]->GetXaxis();          
	      binCenterX = xaxis->GetBinCenter(1);
	      binLowX = xaxis->GetBinLowEdge(1);
	      binWidthX = xaxis->GetBinWidth(1);
	      numBinsX = xaxis->GetNbins();
	      
	      // find bin number for each histogram that corresponds to 2.0 Gev/c^2
	      bin_two = xaxis->FindBin(2.001);
	      bin_five = xaxis->FindBin(4.999);
	    
	      //  cout << bin_two << bin_five << endl;

	      // cout << binWidthX << endl;
	    }
	}
    }

  xaxis = corr_had[1][0][0]->GetXaxis();          
  // cout <<"bg: " << bg[1][0][0]->GetNbinsX() << endl;

  for(int i = 0; i < 100; i++)
    {

      x_array[i]  = xaxis->GetBinCenter(i+1);
      //cout << "x array: " << x_array[i] << endl;
    }
  

  int arm_low;
  int arm_high;
  
  if(south_arm == true)
    {
      arm_low = 0;
      arm_high = 1;
    }
  else
    {
      arm_low = 1;
      arm_high = 2;
    }
  
  // scale cc components from pp histograms (these were invariant yields, not the counts)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  
 	  bg[0][arm][pt]->Scale(scale[arm]); // cc component of Run15 pp
	  bg[3][arm][pt]->Scale(scale[arm]); // bb UL component of Run15 pp
	  bg[8][arm][pt]->Scale(scale[arm]); // bb LS component of Run15 pp
	  pp[0][arm][pt]->Scale(scale[arm]); // dy component of Run15pp
	  pp[1][arm][pt]->Scale(scale[arm]); // corr had UL  component of Run15pp
	  pp[2][arm][pt]->Scale(scale[arm]); // corr had LS component of Run15pp
	}
    }

  char uniqueh[38][500];
 
  int i = 0;
  
    for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      i++;
      sprintf(uniqueh[pt],"h_sum_%d",i); 
    }
  
  TH1D *h1 = new TH1D("h1", "cc distribution", 100, 0, 10);
  
  TH1D *h_sum[2][38];
    
  double cc_pAu_UL[2][38][100] ={0};
  double bb_pAu_UL[2][38][100] ={0};
  double dy_pAu_UL[2][38][100] ={0};
  double bb_pAu_LS_pp[2][38][100] ={0};
  double bb_pAu_LS_mm[2][38][100] ={0};
  double bb_pAu_LS_ave[2][38][100] ={0};
  double bb_pAu_LS[2][38][100] ={0};

 
  double corrhad_pAu_LS[2][38][100] ={0};
  double corrhad_pAu_UL[2][38][100] ={0};
  double pAu_corrbg[2][38][100] ={0};
  double corr_total[2][100] ={0};
  double corr_total_err[2][100] ={0};

  double cc_UL_err[2][38][100];
  double bb_UL_err[2][38][100];
  double bb_LS_err[2][38][100];
  double dy_UL_err[2][38][100];
  double corr_had_UL_err[2][38][100];
  double corr_had_LS_err[2][38][100];
    
  double quad[2][38][100];
  double quad_2[2][38][100];
  double quad_3[2][38][100];
  double quad_norm[2][38][100];
  double temp[2][100];
  double x_errors[2][100];
 


  for(int i_histo = 0; i_histo < 11; i_histo++)
    {
      for(int arm = 0; arm < 1; arm++)  
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      for(int i = 0; i < bins; i++)
		{
		  
		  
		  //cout << "bg[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << bg[i_histo][arm][pt]->GetBinContent(i+1) << endl;
		  // cout << "corr_had[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << corr_had[i_histo][arm][pt]->GetBinContent(i+1) << endl;
		  //cout << "pp[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << pp[i_histo][arm][pt]->GetBinContent(i+1) << endl;

		}
	    }
	}
    }

  // Z
  int pt_l = 4;
  int pt_h = 5;

  // For 0.25 - 5 GeV/c^2fit range
  double scale_factor[2][21] = {0.01,0.01,0.01,0.01,0.542,0.27,0.195,0.245,0.39,0.49,0.69,0.72,0.99,1.5,1.5,1.89,1.37,1.8,1.7,1.55,1.98,
				0.01,0.01,0.01,0.01,0.68,0.53,0.474,0.31,0.41,0.65,0.81,0.61,0.87,1.29,1.4,1.1,1.49,1.57,1.8,1.07,1.55};
  
   for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	{
	  for(int i = 0; i < bins; i++)
	    {
   
	      if (bg[1][arm][pt]->GetBinContent(i+1) != 0)
		cc_pAu_UL[arm][pt][i] = bg[0][arm][pt]->GetBinContent(i+1)*( bg[2][arm][pt]->GetBinContent(i+1)/bg[1][arm][pt]->GetBinContent(i+1) );
	      
	      if( bg[4][arm][pt]->GetBinContent(i+1) != 0)
		bb_pAu_UL[arm][pt][i] = bg[3][arm][pt]->GetBinContent(i+1)*( bg[5][arm][pt]->GetBinContent(i+1)/bg[4][arm][pt]->GetBinContent(i+1) ) *weight[arm];
	      
	      if( (bg[6][arm][pt]->GetBinContent(i+1) + bg[7][arm][pt]->GetBinContent(i+1) ) != 0)
	      	bb_pAu_LS_ave[arm][pt][i] =  ( bg[9][arm][pt]->GetBinContent(i+1) + bg[10][arm][pt]->GetBinContent(i+1) ) / (bg[6][arm][pt]->GetBinContent(i+1) + bg[7][arm][pt]->GetBinContent(i+1) );
	 
	      bb_pAu_LS[arm][pt][i] = bg[8][arm][pt]->GetBinContent(i+1)*( bb_pAu_LS_ave[arm][pt][i] );
	      
	      ///////////////////NEW
	      
	      dy_pAu_UL[arm][pt][i] = pp[0][arm][pt]->GetBinContent(i+1);
	      corrhad_pAu_UL[arm][pt][i] = pp[1][arm][pt]->GetBinContent(i+1)*corr_had[0][arm][pt]->GetBinContent(i+1);
	      corrhad_pAu_LS[arm][pt][i] = pp[2][arm][pt]->GetBinContent(i+1)*corr_had[1][arm][pt]->GetBinContent(i+1);
	      
	      pAu_corrbg[arm][pt][i] = cc_pAu_UL[arm][pt][i] + bb_pAu_UL[arm][pt][i]  + dy_pAu_UL[arm][pt][i] +  corrhad_pAu_UL[arm][pt][i] - corrhad_pAu_LS[arm][pt][i] - bb_pAu_LS[arm][pt][i];
	      
	      h1->SetBinContent(i+1, scale_factor[arm][pt]*pAu_corrbg[arm][pt][i] );
	      
	      corr_total[arm][i] += pAu_corrbg[arm][pt][i];
	      corr_total_err[arm][i] = corr_total[arm][i]*0.05;
	      
	      if(h1->GetBinContent(i+1) > 0)
		{
		  h1->SetBinError(i+1, sqrt( h1->GetBinContent(i+1)) );
		}
	      else
		{
		    h1->SetBinError(i+1,0);
		}
	      
	      x_errors[arm][i] = 0.0;


	      h_sum[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_sum[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));
	      
	    }
	}
    }
  
  
  
  ////////////////////////////////////////////////
  // Histogram comp
  TCanvas *c6 = new TCanvas("c6","Corr BG comp ",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
   gPad->SetGrid();
   // gPad->SetLogY();

   
  if(all_fits)
    {
      h_sum[0][0]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) + corr_had(UL) - bb(LS) - corr_had(LS) ), South");
      h_sum[1][0]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) + corr_had(UL) - bb(LS) - corr_had(LS) ), North");
      h_sum[arm_low][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm_low][0]->GetXaxis()->SetLabelSize(0.04);
      h_sum[arm_low][0]->GetXaxis()->SetTitleSize(0.04);
      h_sum[arm_low][0]->GetXaxis()->SetTitleOffset(0.9);
      h_sum[arm_low][0]->GetYaxis()->SetLabelSize(0.04);
      h_sum[arm_low][0]->GetYaxis()->SetTitleSize(0.1); 
      h_sum[arm_low][0]->GetYaxis()->SetTitleOffset(0.52);
      h_sum[arm_low][0]->SetMinimum(0);
      // if(south_arm == false)
	//	h_sum[arm_low][0]->SetMaximum(1500);
      // else
	h_sum[arm_low][0]->SetMaximum(2500);
	h_sum[arm_low][0]->SetAxisRange(0,5);
    }
  else
    {
      h_sum[0][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) + corr_had(UL) - corr_had(LS) ), South");
      h_sum[1][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) + corr_had(UL) - corr_had(LS) ), North");
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm_low][pt_l]->GetXaxis()->SetLabelSize(0.04);
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitleSize(0.04);
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitleOffset(0.9);
      h_sum[arm_low][pt_l]->GetYaxis()->SetLabelSize(0.04);
      h_sum[arm_low][pt_l]->GetYaxis()->SetTitleSize(0.1); 
      h_sum[arm_low][pt_l]->GetYaxis()->SetTitleOffset(0.52);
      h_sum[arm_low][pt_l]->SetMinimum(0);
      // h_sum[arm_low][pt_l]->SetMaximum(1200);
      h_sum[arm_low][pt_l]->SetAxisRange(0,5);
    }

  for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
    {
      h_sum[arm_low][pt]->SetMarkerColor(pt+1);  
      if(pt+1==10)
	h_sum[arm_low][pt]->SetMarkerColor(30); 
      if(pt+1==2)
	h_sum[arm_low][pt]->SetMarkerColor(28); 
      if(pt+1 == 5)
	h_sum[arm_low][pt]->SetMarkerColor(40); 
      h_sum[arm_low][pt]->SetMarkerSize(0.7);
      h_sum[arm_low][pt]->SetMarkerStyle(20);

      if(all_fits)
	{
	  if(pt == 0)
	    h_sum[arm_low][pt]->Draw();
	  else
	    h_sum[arm_low][pt]->Draw("SAME");
	}
      else
	{
	  if(pt == pt_l)
	    h_sum[arm_low][pt_l]->Draw();
	}
    }

 /////////////////////////////////////////////////////////////////////////////////
 /// Make legend for h_sum histogram

  char unique_l[800];
  TLegend *leg_sum[2];
 
  for(int arm = 0; arm < 2; arm++)
    {
      leg_sum[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
      leg_sum[arm]->SetFillColor(0); 
      leg_sum[arm]->SetTextSize(0.035);
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  sprintf(unique_l,"corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_sum[arm]->AddEntry(h_sum[arm][pt],unique_l, "p");
	}
    }
  
  leg_sum[arm_low]->Draw();



  TLatex l3;
  l3.SetTextSize(0.06);
  l3.SetTextAlign(13);
  l3.SetTextColor(4);

  char text4[100];
  sprintf(text4,"p + Au");
  l3.SetTextAlign(12);
  l3.DrawLatexNDC(0.75, 0.92, text4); //4.4,150
 

 

  ////////////////////////////////////////////////////////////// END HERE FOR BG PLOTS ONLY ///////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////
  /// Make legend for histogram

  char unique_l2[800];
  TLegend *leg_h[2];
 
  int for_begin; 
  int for_end; 

  if(all_fits)
    {
      for_begin = 0;
      for_end = pt_slices[pt_binwidth];
    }
  else
    {
      for_begin = pt_l;
      for_end = pt_h;
    }
 
  for(int arm = 0; arm < 2; arm++)
    {
      leg_h[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
      leg_h[arm]->SetFillColor(0); 
      leg_h[arm]->SetTextSize(0.035);
	  	 
      for(int pt = for_begin; pt < for_end; pt++)
	{
	  sprintf(unique_l2,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_h[arm]->AddEntry(h_sum[arm][pt],unique_l2, "p");
	}
    }
  
   //////////////////////////// Fit histograms //////////////////////////
  
  char name1[800];
  char name2[800];
  
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;
  int bins_fail[bins];
  
  TF1 *corr_bg_fit;
  TF1 *corr_total_fit;
  double f5,f6,f10,f11,f12;
  TFitResultPtr r;

  // BG Shape #2 For 500 MeV bins
  double par_0[14] = {1,1,1,1,1,1,1,1,0.1,1,1,1,1,1};  // original
  double par_1[14] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2[14] = {10,10,10,100,100,10,1,10,10,1,1,1,1,1};
  double par_3[14] = {10,10,10,10,10,10,100,10,10,10,10,10,10,10}; 
  double par_4[14] = {10,10,10,10,10,10,0.01,10,10,10,10,10,10,10};  // original
  
  double par_0_S[14] = {1,1,1,1,1,1,10,1,1,1,1,0.1,1,1};  
  double par_2_S[14] = {10,10,100,10,100,10,1,10,10,10,1,1,1,1};
  double par_3_S[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10}; 
  double par_4_S[14] = {10,10,10,10,10,10,0.01,10,10,10,10,10,10,10}; 

  double limit_a[14] = {1.5,1.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.25,0.25,0.25,0.25};   //1.15
  double limit_a_S[14] = {1.5,1.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};
 
  double limit_b[14] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5};
 
  double par_2_refit[14] = {10,10,100,100,100,10,1,10,10,10,1,1,1,1};
 
  //double limit_a_refit[14] = {1.5,1.5,2,2,2,2,2,2,2,2,2,2,2,2};
  double limit_a_refit[14] = {1.5,1.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};   //1.15
 
  // // For 300 MeV bins
  double par_0_300[21] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,1,1,1,1,1};
  double par_1_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2_300[21] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  double par_3_300[21] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,10,10,10,10,10};
  double par_4_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,10,10,10,10,10};
  
  double limit_a_300[21] = {0,0,0,0.25,2,2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.85,0.85,0.85}; 
  double limit_b_300[21] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};

  double limit_a_300_North[21] = {0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4,0.85}; 
  double limit_a_300_North_refit[21] = {0,0,0,0,0.25,2,2,2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4,0.85}; 
  //double limit_a_300_North_refit[21] = {0,0,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4,0.85}; 







  
  // For 1 GeV bins
  double par_0_1000[10] = {0,0,0,0,1,1,1,0,0,0};
  double par_1_1000[10] = {0,0,0,0,0.01,0.01,0.01,0,0,0};
  double par_2_1000[10] = {0,0,0,0,1,1,1,0,0,0};
  double par_3_1000[10] = {0,0,0,0,10,10,10,0,0,0};
  double par_4_1000[10] = {0,0,0,0,10,10,10,0,0,0};
  
  double limit_a_1000[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 
  double limit_b_1000[10] = {0,0,0,0,5,5,5,0,0,0};

  double limit_a_1000_North[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 

 double refit_par_values[2][5][38];

 // for pT Integrated

 double par_0_total = 1;
 double par_1_total = 0.01;
 double par_2_total = 100;
 double par_3_total = 10;
 double par_4_total = 10;
 double counter = 0;
 double x = 0;
 

 //////////////////////////////////////////////////////////////////////////  All code above this point is for both North and South 
 //cout << "for begin: " << for_begin << "for end" << for_end << endl;  


 for(int arm = arm_low; arm < arm_high; arm++)  
   {
     for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) // ALL FITS TO DISPLAY ON SAME PLOT
       {

	  x = pt_center[pt_binwidth][pt];

	  cout << "x = pT center: " << x << ", and pt: " << pt << endl;
	
	  if(pt_binwidth == 3)
	    {
	      par_0[pt] = par_0_300[pt];
	      par_1[pt] = par_1_300[pt];
	      par_2[pt] = par_2_300[pt];
	      par_3[pt] = par_3_300[pt];
	      par_4[pt] = par_4_300[pt];
	      
	      limit_a[pt] = limit_a_300[pt];
	      limit_b[pt] = limit_b_300[pt];
	      
	      if(south_arm == false)
		{
		  limit_a[pt] = limit_a_300_North[pt];
		  if(refit)
		    limit_a[pt] = limit_a_300_North_refit[pt];
		}
	    }
      

	  corr_bg_fit = new TF1("corr_bg_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",limit_a[pt],limit_b[pt]);

	  	  
	  cout <<" TF1 fit is from " << limit_a[pt] << "  to " << limit_b[pt] << endl;

	  // For fitting bg as function of pT and then plotting parameter values as function of pT
	  //**************************************************************************************
	  if(initial_fit == true)
	    {
	      corr_bg_fit->SetParameter(0, par_0[pt]);
	      corr_bg_fit->SetParameter(1, par_1[pt]);
	      corr_bg_fit->SetParameter(2, par_2[pt]);
	      corr_bg_fit->SetParameter(3, par_3[pt]);
	      corr_bg_fit->SetParameter(4, par_4[pt]);
       
	    }	 

	  // for calculating the parameter from a function plugging in the pt value
	  //************************************************************************

	  double par_a[10];
	  double par_b[10];
	  double par_d[10];
	  double par_e[10];
	  
	  double par_a_fx;
	  double par_b_fx;
	  double par_d_fx;
	  double par_e_fx;

	  if(refit == true)
	    {
	     	    	     
	      f10 = par_2[pt];
	      
	      /// cout << "f10: " << f10 << endl;
	      
	      for(int par = 0; par < 10; par++)
	  	{	 
	  	  if(south_arm)
	  	    {
	  	      par_a[par] = coeff_par_S[0][par];
	  	      par_b[par] = coeff_par_S[1][par];
	  	      par_d[par] = coeff_par_S[2][par];
	  	      par_e[par] = coeff_par_S[3][par];
	  	    }
	  	  else
	  	    {
	  	      par_a[par] = coeff_par_N[0][par];
	  	      par_b[par] = coeff_par_N[1][par];
	  	      par_d[par] = coeff_par_N[2][par];
	  	      par_e[par] = coeff_par_N[3][par];
		      
	  	    }
	  	}

	      
	      par_a_fx = par_a[0] + par_a[1]*x + par_a[2]*x*x + par_a[3]*x*x*x + par_a[4]*x*x*x*x + par_a[5]*x*x*x*x*x + par_a[6]*x*x*x*x*x*x + par_a[7]*x*x*x*x*x*x*x + par_a[8]*x*x*x*x*x*x*x*x + par_a[9]*x*x*x*x*x*x*x*x*x;
	      par_b_fx = par_b[0] + par_b[1]*x + par_b[2]*x*x + par_b[3]*x*x*x + par_b[4]*x*x*x*x + par_b[5]*x*x*x*x*x + par_b[6]*x*x*x*x*x*x + par_b[7]*x*x*x*x*x*x*x + par_b[8]*x*x*x*x*x*x*x*x + par_b[9]*x*x*x*x*x*x*x*x*x;
	      par_d_fx = par_d[0] + par_d[1]*x + par_d[2]*x*x + par_d[3]*x*x*x + par_d[4]*x*x*x*x + par_d[5]*x*x*x*x*x + par_d[6]*x*x*x*x*x*x + par_d[7]*x*x*x*x*x*x*x + par_d[8]*x*x*x*x*x*x*x*x + par_d[9]*x*x*x*x*x*x*x*x*x;
	      par_e_fx = par_e[0] + par_e[1]*x + par_e[2]*x*x + par_e[3]*x*x*x + par_e[4]*x*x*x*x + par_e[5]*x*x*x*x*x + par_e[6]*x*x*x*x*x*x + par_e[7]*x*x*x*x*x*x*x + par_e[8]*x*x*x*x*x*x*x*x + par_e[9]*x*x*x*x*x*x*x*x*x;
	      
	      //  cout << "a0: " << par_a[0] << ", a1: " << par_a[1] << ", a2: " << par_a[2] << ", a3: " << par_a[3] << ", a4: " << par_a[4] << ", a5: " << par_a[5] << ", a6: " << par_a[6] << ", a7: " << par_a[7] << ", a8: " << par_a[8] << ", a9: " << par_a[9] << endl;
	      
	      // cout << "b0: " << par_b[0] << ", b1: " << par_b[1] << ", b2: " << par_b[2] << ", b3: " << par_b[3] << ", b4: " << par_b[4] << ", b5: " << par_b[5] << ", b6: " << par_b[6] << ", b7: " << par_b[7] << ", b8: " << par_b[8] << ", b9: " << par_b[9] << endl;
	      
	      // cout << "d0: " << par_d[0] << ", d1: " << par_d[1] << ", d2: " << par_d[2] << ", d3: " << par_d[3] << ", d4: " << par_d[4] << ", d5: " << par_d[5] << ", d6: " << par_d[6] << ", d7: " << par_d[7] << ", d8: " << par_d[8] << ", d9: " << par_d[9] << endl;
	      
	      // cout << "e0: " << par_e[0] << ", e1: " << par_e[1] << ", e2: " << par_e[2] << ", e3: " << par_e[3] << ", e4: " << par_e[4] << ", e5: " << par_e[5] << ", e6: " << par_e[6] << ", e7: " << par_e[7] << ", e8: " << par_e[8] << ", e9: " << par_e[9] << endl;
	      
	      corr_bg_fit->FixParameter(0, par_a_fx);
	      corr_bg_fit->FixParameter(1, par_b_fx);
	      corr_bg_fit->SetParameter(2, f10);
	      corr_bg_fit->FixParameter(3, par_d_fx);
	      corr_bg_fit->FixParameter(4, par_e_fx);
  
	    }  
	  
	  //****************************************************************************
	  
	  corr_bg_fit->SetLineColor(kRed);
	  corr_bg_fit->SetLineStyle(5);
	  corr_bg_fit->SetLineWidth(2);
	  
		  
	  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
	  
	  //r = gr_corr_pt[pt]->Fit(corr_bg_fit,"S R");	  // for TGraphErrors
	  
	  
	  if(pt_int_fit == false)
	    {
	      r = h_sum[arm_low][pt]->Fit(corr_bg_fit,"S R");
	      cout << "Chi Square/NDF: " << corr_bg_fit->GetChisquare() << "/" << corr_bg_fit->GetNDF() << "" << endl;
	    }
		  
	  p0 = r->Parameter(0);
	  p1 = r->Parameter(1);
	  p2 = r->Parameter(2);
	  p3 = r->Parameter(3);
	  p4 = r->Parameter(4);
	  
	  
	  // fill array with refit parameter values calculated by functions
	  refit_par_values[arm][0][pt] = p0;
	  refit_par_values[arm][1][pt] = p1;
	  refit_par_values[arm][2][pt] = p2;
	  refit_par_values[arm][3][pt] = p3;
	  refit_par_values[arm][4][pt] = p4;
	  
	  e0 = r->ParError(0);
	  e1 = r->ParError(1);
	  e2 = r->ParError(2);
	  e3 = r->ParError(3);
	  e4 = r->ParError(4);
	  
	  Int_t fitStatus = r;  



	  // cout << "status: " << fitStatus << " for pt bin " << pt+1 << endl;

	  //  if(fitStatus == 0) 
	  //   {
	      if(write)	     
		{
		  if(arm == 1)
		    {
		      if(pt_binwidth == 2)
			sprintf(name1,"yuehang_bestfit_500mev/bestfit_parameters_N_%i.dat",pt+1); 
		      if(pt_binwidth == 1)
			sprintf(name1,"yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_%i.dat",pt+1); 
		      if(pt_binwidth == 3)
			sprintf(name1,"yuehang_bestfit_300mev/bestfit_parameters_N_%i.dat",pt+1); 
		      if(pt_binwidth == 0)
			sprintf(name1,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_N_%i.dat",pt+1);
		      if(pt_binwidth == 4)
			sprintf(name1,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_N_%i.dat",pt+1);
		      
		      std::fstream bestfit_parameters_1(name1,std::ofstream::out); 
		      bestfit_parameters_1 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_1.close();
		    }
		  else
		    {
		      if(pt_binwidth == 2)
			sprintf(name2,"yuehang_bestfit_500mev/bestfit_parameters_S_%i.dat",pt+1); 
		      if(pt_binwidth == 1)
			sprintf(name2,"yuehang_bestfit_200mev/bestfit_parameters_low_pt_S_%i.dat",pt+1); 	
		      if(pt_binwidth == 3)
			sprintf(name2,"yuehang_bestfit_300mev/bestfit_parameters_S_%i.dat",pt+1);
		      if(pt_binwidth == 0)
			sprintf(name2,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_%i.dat",pt+1);
		      if(pt_binwidth == 4)
			sprintf(name2,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);
		      
		      
		      std::fstream bestfit_parameters_2(name2,std::ofstream::out); 
		      bestfit_parameters_2 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_2.close();
		    }
		} // write
	      
	      // }// fit converged
	  // else
	  //   {
	  //     counter+= 1;
	  //     cout << "number of fits that failed : " << counter << endl;
	  //   }
	   

	      if(x == 6.5)
	      	break;
	      
       } // for loop pt
     
   } // for loop arm
 
 
 cout << "HERE" << endl;
 

 Char_t message[200];
 //  Char_t message2[200];
 
 if(pt_int_fit == false)
   sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_bg_fit->GetChisquare(),corr_bg_fit->GetNDF());
 else
   sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_total_fit->GetChisquare(),corr_total_fit->GetNDF());
 // sprintf(message2,"J/Psi Counts =  %.1f +/- %.3f",NJpsi,err_NJpsi);
 
 TPaveText *mytext = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
 mytext->SetTextSize(0.035);
 mytext->SetFillColor(0);
 mytext->SetTextAlign(12);
 mytext->AddText(message);
 // mytext->AddText(message2);
 mytext->Draw();
 
 
 
 
 
   
////////////////////////////////////// for ratio plots ///////////////////////////////////////////////

int a = 2;
int b = 3;
  double original_par_values[2][5][21] = {0}; // [arm][parameter][pt]   need [4][2][5][14]
  double ratio_par_values[2][5][21] = {0};
  double pt_error[21];
  double par_array[5];
  double par_error[2][5][38] = {0};
 double pt_slices2[4] = {6,14,21,7};

 pt_binwidth = 2; // overwrite previous decleration of pt_binwidth = 3
  
  char basename4[8][500] = {"yuehang_bestfit_100mev/bestfit_parameters_S_","yuehang_bestfit_200mev/bestfit_parameters_S_","yuehang_bestfit_300mev/bestfit_parameters_S_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_",

			    "yuehang_bestfit_100mev/bestfit_parameters_N_","yuehang_bestfit_200mev/bestfit_parameters_N_","yuehang_bestfit_300mev/bestfit_parameters_N_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_N_"};
  

  // fill 300 mev
  vector < vector < vector <std::string> > > bestfit_filename4;
    
  for(int j_arm = 0; j_arm < 2; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm4;
  
    
      for(int j_series = 0; j_series < 4; j_series++)
	{
	  vector < std::string >  filename_series4;
	  int j = j_arm*4 + j_series;
	  // cout << "j: " << j << endl;
	  
	  for(int i = 0; i < pt_slices2[j_series]; i++)
	    {
	      char name4[500];
	      sprintf(name4,"%s%i%s",basename4[j],i+1,".dat");	     
	      std::string filename_tmp = name4;
	      filename_series4.push_back(filename_tmp);
	    } 

	  filename_arm4.push_back(filename_series4);
	}

      bestfit_filename4.push_back(filename_arm4);

    } 

  // read out 300 mev
  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = 0; series < 4; series++)
	{
	  for(int pt = 0; pt < pt_slices2[series]; pt++)
	    {
	      //cout << bestfit_filename4[arm][series][pt].c_str() << endl; 
	    }
	}
    }


  // get 300 parameters
  for(int arm = arm_low; arm < arm_high; arm++)
    {
      for(int series = a; series < b; series++) 
	{
	  for(int pt = 0; pt < 19; pt++)  // for 300 meV 
	    {
	      for(int i = 0; i < 5; i++)
		{
		  
		  ifstream bestfit_parameters4(bestfit_filename4[arm][series][pt].c_str()); 
		 
		  if(bestfit_parameters4)
		    {
		      bestfit_parameters4 >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      par_array[4] = p4;
		      
		      cout << "par e: " << p4 << ", for bin: " << pt+1 << endl;
		      
		      par_error[arm][0][pt] = e0;
		      par_error[arm][1][pt] = e1;
		      par_error[arm][2][pt] = e2;
		      par_error[arm][3][pt] = e3;
		      par_error[arm][4][pt] = e4;
		      
		      bestfit_parameters4.close();
		    }
		  // else
		  //   {
		  //     par_array[i] = 0;
		  //     par_error[arm][i][pt] = 0.0;
		  //   }
		  
		  original_par_values[arm][i][pt] = par_array[i];
		 
		}
	    }
	} // for loop arm
      
    } // for series loop
 

  char name5[800];
  

  /*
  // get par values  5-6 GeV/c and 6-7 GeV/c
  for(int arm = arm_low; arm < arm_high; arm++)
    {
      for(int series = a+1; series < b+1; series++) 
	{
	  for(int pt = 5; pt < 7; pt++) 
	    {
	      for(int i = 0; i < 5; i++)
		{
		  
		  ifstream bestfit_parameters5(bestfit_filename4[arm][series][pt].c_str()); 

		  if(bestfit_parameters5)
		    {
		      bestfit_parameters5 >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      par_array[4] = p4;

		      cout << "par e: " << p4 << ", for bin: " << pt+1 << endl;

		      if(pt == 5)
			{
			  par_error[arm][0][17] = e0;
			  par_error[arm][1][17] = e1;
			  par_error[arm][2][17] = e2;
			  par_error[arm][3][17] = e3;
			  par_error[arm][4][17] = e4;
			  original_par_values[arm][i][17] = par_array[i];
			}
		      if(pt == 6)
			{
			  par_error[arm][0][18] = e0;
			  par_error[arm][1][18] = e1;
			  par_error[arm][2][18] = e2;
			  par_error[arm][3][18] = e3;
			  par_error[arm][4][18] = e4;
			  original_par_values[arm][i][18] = par_array[i];
			}
		   
		      bestfit_parameters5.close();
		    }
		  else
		    {
		      par_array[i] = 0;
		      par_error[arm][i][pt] = 0.0;
		    }
       
		  pt_error[pt] = 0;

		}
	    }
	} // for loop arm
      
    } // for series loop
  
  */
  
  for(int arm = arm_low; arm < arm_high; arm++)
    {
      for(int i = 0; i < 5; i++)
	{
	  for(int pt = 0; pt < pt_slices2[pt_binwidth]; pt++) // pt = 0
	    {
	      ratio_par_values[arm][i][pt] = refit_par_values[arm][i][pt]/original_par_values[arm][i][pt];
	      if( (ratio_par_values[arm][i][pt] < -100) || (ratio_par_values[arm][i][pt] > 100) )
		ratio_par_values[arm][i][pt] = 0;
	     
	      //cout << ratio_par_values[arm][i][pt] << ", for arm " << arm << ", and i: " << i << ", and pt: " << pt << endl;
	    }
	}
    }
 
  int par_value = 4;
 
  TGraph *ratio[2][5]; //1-5
  double points;
 
  for(int arm = arm_low; arm < arm_high; arm++)  
    {
      for(int i = 0; i < 5; i++) 
	{
	  points = 17;
	  ratio[arm][i] = new TGraph(points,pt_par_plot,ratio_par_values[arm][i]);//,pt_error,par_error[arm][i]); // should be quadratic par error
	  ratio[arm][i]->SetMarkerStyle(20);	      
	  ratio[arm][i]->SetMarkerSize(2);	     
	  ratio[arm][i]->SetMarkerColor(i+1);
	  ratio[arm][i]->GetXaxis()->SetTitle("p_{T} GeV/c");
	  ratio[arm][i]->GetXaxis()->SetLabelSize(0.06);
	  ratio[arm][i]->GetXaxis()->SetTitleSize(0.05);
	  ratio[arm][i]->GetXaxis()->SetTitleOffset(0.9);
	  ratio[arm][i]->GetYaxis()->SetLabelSize(0.06);
	  ratio[arm][i]->GetYaxis()->SetTitleSize(0.1); 
	  ratio[arm][i]->GetYaxis()->SetTitleOffset(0.52);
	  ratio[arm][i]->GetXaxis()->SetLimits(1.1,5);
	  ratio[arm][i]->SetMinimum(-10);
	  ratio[arm][i]->SetMaximum(10);
	 
	}
    }
 
  TCanvas *c12 = new TCanvas("c12","Ratio Corr BG Parameter Refit/Initial Fit",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetGrid();
 
  if(south_arm == true)
    {
      if(par_value == 0)
	ratio[0][par_value]->SetTitle("Ratio Parameter a Refit/Initial Fit, Corr BG South");
      if(par_value == 1)
	ratio[0][par_value]->SetTitle("Ratio Parameter b Refit/Initial Fit, Corr BG South");
      if(par_value == 3)
	ratio[0][par_value]->SetTitle("Ratio Parameter d Refit/Initial Fit, Corr BG South");
      if(par_value == 4)
	ratio[0][par_value]->SetTitle("Ratio Parameter e Refit/Initial Fit, Corr BG South");
    }
  else
    {
      if(par_value == 0)
	ratio[1][par_value]->SetTitle("Parameter a Refit/Initial Ratio, Corr BG North");
      if(par_value == 1)
	ratio[1][par_value]->SetTitle("Parameter b Refit/Initial Ratio, Corr BG North");
      if(par_value == 3)
	ratio[1][par_value]->SetTitle("Parameter d Refit/Initial Ratio, Corr BG North");
      if(par_value == 4)
	ratio[1][par_value]->SetTitle("Parameter e Refit/Initial Ratio, Corr BG North");
    }

  ratio[arm_low][par_value]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ratio[arm_low][par_value]->GetXaxis()->SetLabelSize(0.04);
  ratio[arm_low][par_value]->GetXaxis()->SetTitleSize(0.04);
  ratio[arm_low][par_value]->GetXaxis()->SetTitleOffset(0.9);
  ratio[arm_low][par_value]->GetYaxis()->SetLabelSize(0.04);
  ratio[arm_low][par_value]->GetYaxis()->SetTitleSize(0.1); 
  ratio[arm_low][par_value]->GetYaxis()->SetTitleOffset(0.52);
  ratio[arm_low][par_value]->Draw("AP");
 
  TF1 *linear_fit;
  linear_fit = new TF1("linear_fit","[0]",1.1,5);
 
  //linear_fit->FixParameter(0, 1);
  linear_fit->SetParameter(0, 1);
  linear_fit->SetLineColor(kRed);
  linear_fit->SetLineStyle(5);
  linear_fit->SetLineWidth(2);
 
  ratio[arm_low][par_value]->Fit(linear_fit,"R");
 
  Char_t message2[200];
 
  sprintf(message2,"y = %.3f  +/- %.3f",linear_fit->GetParameter(0),linear_fit->GetParError(0));
  //sprintf(message2,"y = %.1f ",linear_fit->GetParameter(0));
 
  TPaveText *mytext2 = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
  mytext2->SetTextSize(0.035);
  mytext2->SetFillColor(0);
  mytext2->SetTextAlign(12);
  mytext2->AddText(message2);
  mytext2->Draw();
 
 
  
}// void end macro



      
  



      
  
