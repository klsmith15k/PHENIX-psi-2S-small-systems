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


using namespace std;


void yuehang_bg_macro_1gev()
{

#include "yuehang_fit_coeff_300/fit_coeff_S_array_nineteen_point.C"   
#include "yuehang_fit_coeff_300/fit_coeff_N_array_nineteen_point.C"   

  /// to run fit_array setting ->  all_fits/write_minbias/refit = true, everything else false   ---> Chekc line 50 matches 1370
  // to run area under pt curve setting  -> all_fits/write_area/initial_fit = true, everything else false
  
  bool all_fits = true; // set to true for ratio plots
  bool write = false;  // writes bestfit parameters
  bool write_minbias = false; // writes f(x|pt)
  bool write_area = false;
  //bool south_arm = true;
  bool south_arm = false;
  bool initial_fit = true;  // set to false for ratio plots
  bool refit = false;  // set to true for ratio plots
  bool write_xuan = true;  

  bool pt_int_fit = false;
  bool normalization = false;
  
  double weight[2] = {1.03,0.985};
  
  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; // for LS background  ~ Dec 2018
  //double scale[2] = {0.868*pow(10,15),1.78*pow(10,15)}; // for LS background ?  Feb 2019
  //double scale[2] = {1,1};
  
  int pt_binwidth = 3;
  // cout << "What pt binning are you running?  Enter '0' for 100 MeV binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth" << endl;
  // cin >> pt_binwidth;
  
  // int pt_slices[5] = {10,6,14,10,5};  // the first 6 points are 0.0 - 1.2 GeV/c
  int pt_slices[5] = {10,6,10,10,5};  // the first 6 points are 0.0 - 1.2 GeV/c
 
  std::string filename_comp[2][2] = {"yuehang_Run15pAu_S_rebinned_1_gev_cc.root",
  				     "yuehang_Run15pAu_N_rebinned_1_gev_cc.root",

				     "yuehang_Run15pAu_S_rebinned_1_gev_bb.root",
  				     "yuehang_Run15pAu_N_rebinned_1_gev_bb.root"};

std::string filename_corrhad[2] = {"yuehang_Run15pAu_S_rebinned_corr_had.root",
  				     "yuehang_Run15pAu_N_rebinned_corr_had.root"};

std::string filename_Run15pp[2] = {"yuehang_Run15pp_S_rebinned_1_gev.root",
				  "yuehang_Run15pp_N_rebinned_1_gev.root"};
   
 double pt_width[5][38] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
			    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  double pt_center[5][38] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,2.0,4.0,6.0,8.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};





for(int i = 0; i < 5; i++)
  {
    for(int j = 0; j < 38; j++)
      {
	cout << " pt center: " << pt_center[i][j] << ", for i = " <<  i << ", j = " << j << endl;
      }
  }



TFile *file_comp[2][2];
TFile *file_corrhad[2];
TFile *file_Run15pp[2];

  TH1D *bg[11][2][38] = {0};
  TH1D *corr_had[2][2][38] = {0}; 
  TH1D *pp[3][2][38] = {0}; 

 
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
	      cout << "obj 1: " << obj_filename[i_histo][arm][pt].c_str() << " " << endl;
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
	      cout << "obj 2: " << obj_filename_corrhad[i_histo][arm][pt].c_str() << " " << endl;
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
 file_comp[0][0]->GetObject("cc_powheg_0_1000", read_this);
 cout << "read_this : " << read_this << endl;
 read_this = 0;

 // i_histo == 1 is bb
 file_comp[1][0]->GetObject("bb_powheg_0_1000", read_this);
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
	      if(pt < 1)
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][0].c_str(),corr_had[i_histo][arm][pt]); // first 2 elements are 0-1 GeV/c [0]
	      if( (pt > 0) && (pt < 2))
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][1].c_str(),corr_had[i_histo][arm][pt]); // next 2 elements are 1-2 GeV/c [1]
	      if( (pt > 1) && (pt < 3))
		file_corrhad[arm]->GetObject(obj_filename_corrhad[i_histo][arm][2].c_str(),corr_had[i_histo][arm][pt]); // next 2 elements are 2-3 GeV/c [2]
	      if( pt > 2)
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

  xaxis = bg[1][0][0]->GetXaxis();          
  // cout <<"bg: " << bg[1][0][0]->GetNbinsX() << endl;

  for(int i = 0; i < 100; i++)
    {

      x_array[i]  = xaxis->GetBinCenter(i+1);
      cout <<  x_array[i] << endl;
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
  double bb_pAu_LS[2][38][100] ={0};
  double bb_pAu_LS_ave[2][38][100] ={0};

 
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
  
  
  
  // for(int i_histo = 0; i_histo < 11; i_histo++) // for bg
  for(int i_histo = 0; i_histo < 2; i_histo++) // for corr had
    //for(int i_histo = 0; i_histo < 3; i_histo++) // for pp
      {
      for(int arm = 0; arm < 2; arm++)  
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      for(int i = 0; i < bins; i++)
		{
		  
		  
		  //cout << "bg[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << bg[i_histo][arm][pt]->GetBinContent(i+1) << endl;
		  // cout << "corr_had[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << corr_had[i_histo][arm][pt]->GetBinContent(i+1) << endl;
		  // cout << "pp[" << i_histo << "][" << arm << "][" << pt << "][" << i << "]: " << pp[i_histo][arm][pt]->GetBinContent(i+1) << endl;
		  
		}
	    }
	}
    }

  // Z
  int pt_l = 3;
  int pt_h = 4;

  double scale_factor[2][10] = {0.39,0.22,0.15,0.5,0.49,0.865,1.77,2.19,2.19,2.19,
				0.1,0.1,0.1,0.1,0.1,1.45,1.65,0.1,0.1,0.1};
  

  //////////////////////////     for pt integrated:
  TH1D *h_pt_int[2];
  h_pt_int[0] = new TH1D("h_pt_int", "pT Integrated South", 100, 0, 10);
  h_pt_int[1] = new TH1D("h_pt_int_2", "pT Integrated North", 100, 0, 10);

  double h_sum_pt_int[2][100];
  double area[2][38];
  ////////////////////////////
  
 
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
	      
	      ///////////////////NEW

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

	      if((pt > 4) && (pt < 7)) // only two pt bins are used for background
		 h_sum_pt_int[arm][i] += h_sum[arm][pt]->GetBinContent(i+1);
	    }

	  area[arm][pt] = h_sum[arm][pt]->Integral( h_sum[0][0]->FindBin(2.001),  h_sum[0][0]->FindBin(4.999) ); // for pt integrated
	}
    }
  
  int a; 
  int b; 

  if(pt_binwidth == 2)
    {
      a = 2;
      b = 3;
    }   
 if(pt_binwidth == 1)
    {
      a = 1;
      b = 2;
    }  
 if(pt_binwidth == 3)
    {
      a = 3;
      b = 4;
    }  

  char name100[800];
  char name200[800];

  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	{
	  if(write_area)	     
	    {
	      if(arm == 1)
		{
		  if(pt_binwidth == 0)
		    sprintf(name100,"pt_int/weighted_fits_100mev_N_%i.dat",pt+1); 
		  if(pt_binwidth == 1)
		    sprintf(name100,"pt_int/weighted_fits_200mev_N_%i.dat",pt+1); 
		  if(pt_binwidth == 2)
		    sprintf(name100,"pt_int/weighted_fits_500mev_N_%i.dat",pt+1); 
		  if(pt_binwidth == 3)
		    sprintf(name100,"pt_int/weighted_fits_1gev_N_%i.dat",pt+1); 
		  if(pt_binwidth == 4)
		    sprintf(name100,"pt_int/weighted_fits_2gev_N_%i.dat",pt+1); 
		  
		  std::fstream bestfit_parameters_100(name100,std::ofstream::out); 
		  bestfit_parameters_100 <<  area[arm][pt]  <<  " "  << endl;
		  bestfit_parameters_100.close();
		}
	      if(arm == 0)
		{
		  if(pt_binwidth == 0)
		    sprintf(name200,"pt_int/weighted_fits_100mev_S_%i.dat",pt+1); 
		  if(pt_binwidth == 1)
		    sprintf(name200,"pt_int/weighted_fits_200mev_S_%i.dat",pt+1); 
		  if(pt_binwidth == 2)
		    sprintf(name200,"pt_int/weighted_fits_500mev_S_%i.dat",pt+1); 
		  if(pt_binwidth == 3)
		    sprintf(name200,"pt_int/weighted_fits_1gev_S_%i.dat",pt+1); 
		  if(pt_binwidth == 4)
		    sprintf(name200,"pt_int/weighted_fits_2gev_S_%i.dat",pt+1); 
		  
		  std::fstream bestfit_parameters_200(name200,std::ofstream::out); 
		  bestfit_parameters_200 <<  area[arm][pt]  <<  " "  << endl;
		  bestfit_parameters_200.close();
		} // arm == 0
	    }// write
	} // for pt
    } // for arm
  
  ///for pt f(x) integrated


 
  if(write_xuan && !south_arm)
    {

      std::fstream fit_array_file0("pt_int/yuehang_direct_pt_int_1gev_N.C", std::ofstream::out);
       fit_array_file0 << "double direct_pt_int_1gev_N[100] = { " ;
        
       for(int bin = 0; bin < 100; bin++)
	 {
	   if(bin < 99)
	     fit_array_file0 << h_sum_pt_int[1][bin]  <<  ", ";
	   if((bin == 99))
	     fit_array_file0 << h_sum_pt_int[1][bin] << "};" << endl;
	 }
       
       fit_array_file0 << endl;

    } // write_xuan north
  
  if(write_xuan && south_arm)
    {

      std::fstream fit_array_file11("pt_int/yuehang_direct_pt_int_1gev_S.C", std::ofstream::out);
       fit_array_file11 << "double direct_pt_int_1gev_S[100] = { " ;
        
       for(int bin = 0; bin < 100; bin++)
	 {
	   if(bin < 99)
	     fit_array_file11 << h_sum_pt_int[0][bin]  <<  ", ";
	   if((bin == 99))
	     fit_array_file11 << h_sum_pt_int[0][bin] << "};" << endl;
	 }
       
       fit_array_file11 << endl;

    } // write_xuan south
  
  
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

  

  if(normalization)
    h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,0));    
  else
    {
      if(all_fits)
	{
	  //h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-1),1600);
	  h_sum[arm_low][0]->SetAxisRange(0,5);
	}
      else
	{
	  if(pt_binwidth == 2)
	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),1500);
	  //h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),180);
	  // h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),800);
	  if(pt_binwidth == 1)
	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),4*pow(10,2));
	  h_sum[arm_low][pt_l]->SetAxisRange(0,5);
	}
    }
  /////////////////////////////////////////////////////////////////////////////////
  /// Make legend for histogram

  char unique_l2[800];
  TLegend *leg_h[2];
 
  double for_begin; 
  double for_end; 

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
  
  //leg_h[arm_low]->Draw();
  
  
  
  
  
  // pT Integrated corrletaed background plot
  //  TGraphErrors *gr_corrbg_total = new TGraphErrors(bins,x_array,corr_total[arm_low],x_errors[arm_low],corr_total_err[arm_low]);
  //  TGraphErrors *gr_corrbg_total_N = new TGraphErrors(bins,x_array,corr_total[1],x_errors[1],corr_total_err[1]);
  
  
  //   TCanvas *c0 = new TCanvas("c0","Yue Hang's pAu Corr BG, pT Integrated",200,10,700,500);
//   gPad->SetGrid();
//   // gPad->SetLogy();
//   //  if(arm_low == 0)
//   gr_corrbg_total->SetTitle("Yue Hang's pAu Corr BG, pT Integrated");  
//   //  if(arm_low == 1)
//   //   gr_corrbg_total->SetTitle("Corr BG Test pT Integrated, North");  

//   gr_corrbg_total->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
//   gr_corrbg_total->GetXaxis()->SetLabelSize(0.04);
//   gr_corrbg_total->GetXaxis()->SetTitleSize(0.04);
//   gr_corrbg_total->GetXaxis()->SetTitleOffset(0.9);
//   gr_corrbg_total->GetXaxis()->SetLimits(0,5);
//   gr_corrbg_total->GetYaxis()->SetLabelSize(0.04);
//   gr_corrbg_total->GetYaxis()->SetTitleSize(0.1); 
//   gr_corrbg_total->GetYaxis()->SetTitleOffset(0.52);
//   gr_corrbg_total->SetMaximum(5*pow(10,4));
//   gr_corrbg_total->SetMinimum(pow(10,2));

//   /// assign different colors and draw each bg
//   // if(arm_low == 0)
//     gr_corrbg_total->SetMarkerColor(kMagenta);  
//     // else
//     gr_corrbg_total_N->SetMarkerColor(kAzure);  
//   gr_corrbg_total->SetMarkerSize(0.75);
//   gr_corrbg_total_N->SetMarkerSize(0.75);
//   gr_corrbg_total->SetMarkerStyle(20);
//   gr_corrbg_total_N->SetMarkerStyle(20);
//   gr_corrbg_total->Draw("AP");
//   gr_corrbg_total_N->Draw("P");
  
//  TLegend *leg0;
//  leg0 = new TLegend(0.51, 0.6, 0.8, 0.9);  
//  leg0->SetFillColor(0); 
//  leg0->SetTextSize(0.035);
//  leg0->AddEntry(gr_corrbg_total,"corr bg, South", "p");
//  leg0->AddEntry(gr_corrbg_total_N,"corr bg, North", "p");
 
//  leg0->Draw();
// //////// must comment our everything below this line to generate the total corr bg plots



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


  double limit_a[14] = {1.5,1.5,0.25,0.25,0.25,1,1,0.25,0.25,1.15,0.5,0.5,0.5,0.5}; 
  double limit_a_S[14] = {1.5,1.5,0.25,0.25,1,1,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5};
  double limit_b[14] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5};
 
  double par_2_refit[14] = {10,10,100,100,100,100,100,10,10,100,1,1,1,1}; 
  double limit_a_refit[14] = {1,1,1,1,1,1,1,1,1,1};


      // For 200 MeV bins (only interested in first 6 bins: 0.0 - 1.2 GeV/c, which is low pT corr BG shape)
  double par_0_200[6] = {10,10,10,10,10,10};
  double par_1_200[6] = {0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2_200[6] = {1,1,1,1,1,1};
  double par_3_200[6] = {100,100,100,100,100,100};
  double par_4_200[6] = {0.01,0.01,0.01,0.01,0.01,0.01};
 
  double limit_a_200[6] = {1.4,1.4,1.4,1.4,1.4,1.4};  // for south arm


  //double limit_a_200_North[6] = {1.5,1.4,1.5,1.4,1.5,1.5}; // for noth arm
  // double limit_a_200_North[6] = {1.45,1.45,1.45,1.45,1.45,1.45}; // for noth ar
  // double limit_a_200_North[6] = {1.4,1.75,1.25,1.25,1.25,1.25};  // all converge at these limits
  double limit_a_200_North[6] = {1.4,1.4,1.4,1.4,1.4,1.4};


  double limit_b_200[6] = {5,5,5,5,5,5};

  // // For 300 MeV bins
  double par_0_300[21] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,1,1,1,1,1};
  double par_1_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2_300[21] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  double par_3_300[21] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,10,10,10,10,10};
  double par_4_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,10,10,10,10,10};
  
  double limit_a_300[21] = {1.4,1.4,1.4,1.4,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.85,1.4,0.85}; 
  //double limit_a_300[21] = {1.7,1.7,1.7,1.7,1.7,1.7,1.7,0.25,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0}; 
  double limit_b_300[21] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};

  // For 100 MeV bins (only interested in first 5 bins: 0.0 - 1.0 GeV/c, which is low pT corr BG shape)
  double par_0_100[10] = {10,10,10,10,10,10,10,10,10,10}; // only try 400 - 800
  double par_1_100[10] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2_100[10] = {1,1,1,1,1,1,1,1,10,10};
  double par_3_100[10] = {100,100,100,100,100,100,100,100,100,100};
  double par_4_100[10] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
 
  double limit_a_100[10] = {1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.6};
  double limit_b_100[10] = {5,5,5,5,5,5,5,5,5,5};
  
  // For 1 GeV bins
  double par_0_1000[10] = {1,1,1,1,1,1,1,1,1,1};
  double par_1_1000[10] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_2_1000[10] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};


  double par_3_1000[10] = {10,10,10,10,10,10,10,10,10,10};
  double par_4_1000[10] = {10,10,10,10,10,10,10,10,10,10};
  
  double limit_a_1000[10] = {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};
  double limit_b_1000[10] = {5,5,5,5,5,5,5,5,5,5};

  double limit_a_1000_North[10] = {0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};

  // For 2 GeV bins
  double par_0_2000[5] = {1,1,1,1,1};
  double par_1_2000[5] = {0.01,0.01,0.01,0.01,0.01};
  double par_2_2000[5] = {1,1,1,1,1};
  double par_3_2000[5] = {10,10,10,10,10};
  double par_4_2000[5] = {10,10,10,10,10};
  
  double limit_a_2000[5] = {0.75,0.25,0.75,0.75,0.75};
  double limit_b_2000[5] = {5,5,5,5,5};





 double refit_par_values[2][5][38];

 // for pT Integrated

 double par_0_total = 1;
 double par_1_total = 0.01;
 double par_2_total = 100;
 double par_3_total = 10;
 double par_4_total = 10;
 double counter = 0;


 double fit_array[14][2][100] = {0}; 
 

 //////////////////////////////////////////////////////////////////////////  All code above this point is for both North and South 
  
	

 for(int arm = arm_low; arm < arm_high; arm++)  
   {
     for(int pt = for_begin; pt < for_end; pt++) // ALL FITS TO DISPLAY ON SAME PLOT
       {
	   for(int bin = 0; bin < 100; bin++)
	    {
	
	      if(south_arm && (pt_binwidth == 2))
		{
		  limit_a[pt] = limit_a_S[pt];
		  par_0[pt] = par_0_S[pt];
		  par_2[pt] = par_2_S[pt];
		  par_3[pt] = par_3_S[pt];
		  par_4[pt] = par_4_S[pt];
		}
	 
	 	 
	      if(pt_binwidth == 3)
		{
		  par_0[pt] = par_0_1000[pt];
		  par_1[pt] = par_1_1000[pt];
		  par_2[pt] = par_2_1000[pt];
		  par_3[pt] = par_3_1000[pt];
		  par_4[pt] = par_4_1000[pt];
	     
		  limit_a[pt] = limit_a_1000[pt];
		  limit_b[pt] = limit_b_1000[pt];
	     
       		}
		 
	      corr_bg_fit = new TF1("corr_bg_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",limit_a[pt],limit_b[pt]);

	      corr_total_fit = new TF1("corr_total_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,5);
	  	  
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

		  corr_total_fit->SetParameter(0, par_0_total);
		  corr_total_fit->SetParameter(1, par_1_total);
		  corr_total_fit->SetParameter(2, par_2_total);
		  corr_total_fit->SetParameter(3, par_3_total);
		  corr_total_fit->SetParameter(4, par_4_total);
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
	  
	      double x;
	  

	      x = pt_center[3][pt];
	      cout << "x = pT center: " << x << ", and pt: " << pt << endl;


	      if(refit == true)
		{
	      
		  if(pt < pt_l)
		    continue;
	     
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
	  
	      corr_total_fit->SetLineColor(kRed);
	      corr_total_fit->SetLineStyle(5);
	      corr_total_fit->SetLineWidth(2);
	  
	      ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
	  
	      //r = gr_corr_pt[pt]->Fit(corr_bg_fit,"S R");	  // for TGraphErrors
	  
	  
	      if(pt_int_fit == false)
		{
		  r = h_sum[arm_low][pt]->Fit(corr_bg_fit,"S R");
		  fit_array[pt][arm][bin] = corr_bg_fit->Eval(x_array[bin]);
		  cout << "Chi Square/NDF: " << corr_bg_fit->GetChisquare() << "/" << corr_bg_fit->GetNDF() << "" << endl;
		}
	      else
		{
		  //  r = gr_corrbg_total->Fit(corr_total_fit,"S R");
		   cout << "Chi Square/NDF: " << corr_total_fit->GetChisquare() << "/" << corr_total_fit->GetNDF() << "" << endl;
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

	      if(write)	     
		{
		  if(arm == 1)
		    {
		      if(pt_binwidth == 2)
			sprintf(name1,"yuehang_bestfit_500mev/bestfit_parameters_N_%i.dat",pt+1); 
		      if(pt_binwidth == 1)
			sprintf(name1,"yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_%i.dat",pt+1); 
		      if(pt_binwidth == 3)
			sprintf(name1,"yuehang_bestfit_1gev/bestfit_parameters_N_%i.dat",pt+1); 
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
			sprintf(name2,"yuehang_bestfit_1gev/bestfit_parameters_S_%i.dat",pt+1);
		      if(pt_binwidth == 0)
			sprintf(name2,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_%i.dat",pt+1);
		      if(pt_binwidth == 4)
			sprintf(name2,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);
		      
		      
		      std::fstream bestfit_parameters_2(name2,std::ofstream::out); 
		      bestfit_parameters_2 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_2.close();
		    }

		} // write

	    }// for bin 

       } // for loop pt
     
   } // for loop arm
    
 
  
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



//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // calulate ratio of each refit parameter to the original parameter
  
//   double original_par_values[2][5][14] = {0}; // [arm][parameter][pt]   need [4][2][5][14]
//   double ratio_par_values[2][5][14] = {0};
//   double pt_error[14];
//   double par_array[5];
//   double par_error[2][5][38] = {0};
//   // int pt_slices2[4] = {6,14,10,5};
  
//   char basename2[8][500] = {"yuehang_bestfit_100mev/bestfit_parameters_S_","yuehang_bestfit_200mev/bestfit_parameters_S_","yuehang_bestfit_500mev/bestfit_parameters_S_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_",

// 			    "yuehang_bestfit_100mev/bestfit_parameters_N_","yuehang_bestfit_200mev/bestfit_parameters_N_","yuehang_bestfit_500mev/bestfit_parameters_N_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_N_"};
  
//   vector < vector < vector <std::string> > > bestfit_filename;
    
//   for(int j_arm = 0; j_arm < 2; j_arm++)
//     {
//       vector  < vector <std::string> >  filename_arm;
  
    
//       for(int j_series = 0; j_series < 4; j_series++)
// 	{
// 	  vector < std::string >  filename_series;
// 	  int j = j_arm*4 + j_series;
// 	  // cout << "j: " << j << endl;
	  
// 	  for(int i = 0; i < pt_slices[j_series]; i++)
// 	    {
// 	      char name[500];
// 	      sprintf(name,"%s%i%s",basename2[j],i+1,".dat");	     
// 	      std::string filename_tmp = name;
// 	      filename_series.push_back(filename_tmp);
// 	    } 

// 	  filename_arm.push_back(filename_series);
// 	}

//       bestfit_filename.push_back(filename_arm);

//     } 

//   for(int arm = 0; arm < 2; arm++)
//     {
//       for(int series = 0; series < 4; series++)
// 	{
// 	  for(int pt = 0; pt < pt_slices[series]; pt++)
// 	    {
// 	      cout << bestfit_filename[arm][series][pt].c_str() << endl; 
// 	    }
// 	}
//     }

//   int a; 
//   int b; 

//   if(pt_binwidth == 2)
//     {
//       a = 2;
//       b = 3;
//     }

//   for(int arm = arm_low; arm < arm_high; arm++)
//     {
//       for(int series = a; series < b; series++) 
// 	{
// 	  for(int pt = 0; pt < 10; pt++)  // for 500 meV 4 - 5
// 	    {
// 	      for(int i = 0; i < 5; i++)
// 		{
		 
// 		  ifstream bestfit_parameters(bestfit_filename[arm][series][pt].c_str()); 
		 
// 		  if(bestfit_parameters)
// 		    {
// 		      bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
// 		      par_array[0] = p0;
// 		      par_array[1] = p1;
// 		      par_array[2] = p2;
// 		      par_array[3] = p3;
// 		      par_array[4] = p4;
		      
// 		      cout << "par e: " << p4 << ", for bin: " << pt+1 << endl;
		      
// 		      par_error[arm][0][pt] = e0;
// 		      par_error[arm][1][pt] = e1;
// 		      par_error[arm][2][pt] = e2;
// 		      par_error[arm][3][pt] = e3;
// 		      par_error[arm][4][pt] = e4;
		      
// 		      bestfit_parameters.close();
// 		    }
// 		  else
// 		    {
// 		      par_array[i] = 0;
// 		      par_error[arm][i][pt] = 0.0;
// 		    }
		  
// 		  original_par_values[arm][i][pt] = par_array[i];
		 
// 		}
// 	    }
// 	} // for loop arm
      
//     } // for series loop
 

//   char name5[800];
  
//   //For reading in values from 5-6 GeV/c and 6-7 GeV/c
//   for(int arm = arm_low; arm < arm_high; arm++)
//     {
//       for(int series = a+1; series < b+1; series++) 
// 	{
// 	  for(int pt = 5; pt < 7; pt++) 
// 	    {
// 	      for(int i = 0; i < 5; i++)
// 		{
		  
// 		  ifstream bestfit_parameters2(bestfit_filename[arm][series][pt].c_str()); 

// 		  if(bestfit_parameters2)
// 		    {
// 		      bestfit_parameters2 >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
// 		      par_array[0] = p0;
// 		      par_array[1] = p1;
// 		      par_array[2] = p2;
// 		      par_array[3] = p3;
// 		      par_array[4] = p4;

// 		      cout << "par e: " << p4 << ", for bin: " << pt+1 << endl;

// 		      if(pt == 5)
// 			{
// 			  par_error[arm][0][10] = e0;
// 			  par_error[arm][1][10] = e1;
// 			  par_error[arm][2][10] = e2;
// 			  par_error[arm][3][10] = e3;
// 			  par_error[arm][4][10] = e4;
// 			  original_par_values[arm][i][10] = par_array[i];
// 			}
// 		      if(pt == 6)
// 			{
// 			  par_error[arm][0][11] = e0;
// 			  par_error[arm][1][11] = e1;
// 			  par_error[arm][2][11] = e2;
// 			  par_error[arm][3][11] = e3;
// 			  par_error[arm][4][11] = e4;
// 			  original_par_values[arm][i][11] = par_array[i];
// 			}
		   
// 		      bestfit_parameters2.close();
// 		    }
// 		  else
// 		    {
// 		      par_array[i] = 0;
// 		      par_error[arm][i][pt] = 0.0;
// 		    }
       
// 		  pt_error[pt] = 0;

// 		}
// 	    }
// 	} // for loop arm
      
//     } // for series loop



 
 // CREATE MINBIAS 100 ENTRY ARRAY
   
  if(write_minbias && !south_arm)
    {

      std::fstream fit_array_file("pt_int/fit_array_1gev_N.C", std::ofstream::out);
      fit_array_file << "double fit_array_1gev_N[17][100] = { " ;
      
      for(int pt = 0; pt < 7; pt++)
	{
	  for(int bin = 0; bin < 100; bin++)
	    {
	      
	      if(pt < 6) // prevents comma 
		fit_array_file <<  fit_array[pt][1][bin]  <<  ",";
	      if(pt == 6 && bin < 99)
		fit_array_file <<  fit_array[pt][1][bin]  <<  ",";
	      if((bin == 99) && (pt == 6))
		fit_array_file <<  fit_array[pt][1][bin] << "};" << endl;    

	    }
	      
	}

      fit_array_file << endl;

    } // write_minbias
  
  if(write_minbias && south_arm)
    {

      std::fstream fit_array_file3("pt_int/fit_array_1gev_S.C", std::ofstream::out);
       fit_array_file3 << "double fit_array_1gev_S[17][100] = { " ;
      
      for(int pt = 0; pt < 7; pt++)
	{
	  for(int bin = 0; bin < 100; bin++)
	    {
	      
	      if(pt < 6) // prevents comma 
		fit_array_file3 <<  fit_array[pt][0][bin]  <<  ",";
	      if(pt == 6 && bin < 99)
		fit_array_file3 <<  fit_array[pt][0][bin]  <<  ",";
	      if((bin == 99) && (pt == 6))
		fit_array_file3 <<  fit_array[pt][0][bin] << "};" << endl;   	 
	      
	    }
	      
	}

      fit_array_file3 << endl;

    } // write_minbias
  

}// void end macro



      
  
