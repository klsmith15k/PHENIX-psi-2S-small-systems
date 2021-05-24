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

void yuehang_bg_macro_hard_coded()
{

  bool all_fits = true;
  bool write = false;
  bool south_arm = true;
  //bool south_arm = false;

  int pt_l = 2;
  int pt_h = 3;

  double weight[2][2] = {1,1.03,
			 1,0.985};
			 
  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; // for LS background 

  int pt_binwidth;
  cout << "What pt binning are you running?  Enter '0' for 100 MeV binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth" << endl;
  cin >> pt_binwidth;
  
  // int pt_slices[5] = {10,6,14,10,5};  // the first 6 points are 0.0 - 1.2 GeV/c
  int pt_slices[5] = {10,6,14,10,5};  // the first 6 points are 0.0 - 1.2 GeV/c
 
  std::string filename_comp[3][2] = {"yuehang_Run15pAu_S_rebinned_500mev_cc.root",
  				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root",

				     "yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  				     "yuehang_Run15pAu_N_rebinned_500mev_bb.root",

				     "yuehang_Run15pp_S_rebinned_500_mev.root",
  				     "yuehang_Run15pp_N_rebinned_500_mev.root"};

   
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


  TFile *file_comp;
  TH1D *bg[7][2][38] = {0};
  TH1D *RpAu_had_UL[2][38] = {0};
  TH1D *temp_had_UL[2][38] = {0};
  TH1D *RpAu_had_LS[2][38] = {0};
  TH1D *temp_had_LS[2][38] = {0};
  TH1D *corr_had_UL[2][38] = {0}; 
  TH1D *corr_had_LS[2][38] = {0}; 

  std::string obj_filename[9][2][38];  //[number of files][arm][pt slices]

  char unique0[800];
  char unique1[800];
  char unique2[800];
  char unique3[800];
  char unique4[800];
  char unique5[800];
  char unique6[800];
  char unique7[800];
  char unique8[800];
 
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
  
	  double a = (pt_low - delta_pt)*1000; 
	  double b = (pt_high + delta_pt)*1000;

	  sprintf(unique0,"cc_pp%.0f_%.0f",a,b);                  // cc unmodified is [0]
	  sprintf(unique1,"cc_powheg_%.0f_%.0f",a,b);        // pp unmodified is [1]
	  sprintf(unique2,"cc_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) is [2]  
		  
	  sprintf(unique3,"bb_pp%.0f_%.0f",a,b);                  // bb unmodified is [0]
	  sprintf(unique4,"bb_powheg_%.0f_%.0f",a,b);        // pp unmodified is [1]
	  sprintf(unique5,"bb_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) is [2]  
	
	  sprintf(unique6,"dy_%.0f_%.0f",a,b); 
	  sprintf(unique7,"corr_had_UL_%.0f_%.0f",a,b);
	  sprintf(unique8,"corr_had_LS_%.0f_%.0f",a,b);
               
	  a_array[k] = a/1000;
	  b_array[k] = b/1000;

	  obj_filename[0][arm][k] = unique0;
	  obj_filename[1][arm][k] = unique1;
	  obj_filename[2][arm][k] = unique2;
	  obj_filename[3][arm][k] = unique3;
	  obj_filename[4][arm][k] = unique4;
	  obj_filename[5][arm][k] = unique5;
	  obj_filename[6][arm][k] = unique6;
	  obj_filename[7][arm][k] = unique7;
	  obj_filename[8][arm][k] = unique8;
	}
    }
  
  for(int i_histo = 0; i_histo < 9; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      cout << obj_filename[i_histo][arm][pt].c_str() << " " << endl;
	    }
	}
    }
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 9; i_histo++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      file_comp = TFile::Open(filename_comp[0][arm].c_str()); 
	      if(i_histo < 3)
		file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	      file_comp = TFile::Open(filename_comp[1][arm].c_str()); 
	      if((i_histo > 2) && (i_histo < 6))
		file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	      file_comp = TFile::Open(filename_comp[2][arm].c_str()); 
	      if(i_histo == 6)
		file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	      if(i_histo == 7)
		file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),corr_had_UL[arm][pt]);
	      if(i_histo == 8)
		file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),corr_had_LS[arm][pt]);
	    }
	}	  
    }
    
  // fill the arrays for corr hadrons
  // TFile *file_corrhad;
 
  std::string obj_filename_2[2][2][4] = {
"r_rpA2_ccrpA2_sim_jet_mass_ptslice_FG12[0][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FG12[1][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FG12[2][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FG12[3][0]",
    "r_rpA2_ccrpA2_sim_jet_mass_ptslice_FGLS[0][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FGLS[1][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FGLS[2][0]","r_rpA2_ccrpA2_sim_jet_mass_ptslice_FGLS[3][0]",
    
    "r_rpA1_ccrpA1_sim_jet_mass_ptslice_FG12[0][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FG12[1][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FG12[2][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FG12[3][0]",
    "r_rpA1_ccrpA1_sim_jet_mass_ptslice_FGLS[0][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FGLS[1][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FGLS[2][0]","r_rpA1_ccrpA1_sim_jet_mass_ptslice_FGLS[3][0]"};
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      TFile *file_corrhad = TFile::Open("hbgrpAu_16_ccrpAu_1_anavf.root");
	      file_corrhad->GetObject(obj_filename_2[arm][0][pt].c_str(),temp_had_UL[arm][pt]);
	      file_corrhad->GetObject(obj_filename_2[arm][1][pt].c_str(),temp_had_LS[arm][pt]);
	    }
	}	  
    }

 for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      cout << "obj_filename_2[" << arm << "][" << i_histo << "][" << pt << "]: " << obj_filename_2[arm][i_histo][pt].c_str() << endl;
	    }
	}	  
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
  double numbins_array[2][2][4];
  double bin_2gev_array[2][2][4];
  double bin_5gev_array[2][2][4];

  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      if(i_histo == 0)
		xaxis = temp_had_UL[arm][pt]->GetXaxis();  // 
	      if(i_histo == 1)
		xaxis = temp_had_LS[arm][pt]->GetXaxis();  // 

	      binCenterX = xaxis->GetBinCenter(1);
	      binLowX = xaxis->GetBinLowEdge(1);
	      binWidthX = xaxis->GetBinWidth(1);
	      numBinsX = xaxis->GetNbins();
	      
	      numbins_array[arm][i_histo][pt] = numBinsX;
	      // cout << "num bins array: " << numbins_array[arm][i_histo][pt] << ", for arm: " << arm << ", i_histo: " << i_histo << ", and pt: " << pt << endl;	      

	      // find bin number for each histogram that corresponds to 2.0 Gev/c^2
	      bin_two = xaxis->FindBin(2.001);
	      bin_2gev_array[arm][i_histo][pt] = bin_two;

	      bin_five = xaxis->FindBin(4.999);
	      bin_5gev_array[arm][i_histo][pt] = bin_five;

	      cout << "bin 2gev array: " << bin_2gev_array[arm][i_histo][pt] << ", for arm: " << arm << ", i_histo: " << i_histo << ", and pt: " << pt << endl;	      
	      cout << "bin 5gev array: " << bin_5gev_array[arm][i_histo][pt] << ", for arm: " << arm << ", i_histo: " << i_histo << ", and pt: " << pt << endl;	      

	    }
	}
    }

  xaxis = bg[1][0][0]->GetXaxis();          

  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
      //  cout << "x array: " << x_array[i] << endl;
    }
  
  /// NEW FOR CORR HADRONS pAu

  double initial_binning[2][2][4][100] = {0};  //arm/i_histo/pt/mass bin


  TAxis *xaxis_UL[2][4];
  TAxis *xaxis_LS[2][4];

  double x_array_UL[2][4][100] = {0};
  double x_array_width_UL[2][4][100] = {0};
  double x_array_LS[2][4][100] = {0};
  double x_array_width_LS[2][4][100] = {0};
  int count = 0;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < 4; pt++)
	{  
	  for(int i = 0; i < 100; i++)
	    {
	      xaxis_UL[arm][pt] = temp_had_UL[arm][pt]->GetXaxis();
	      xaxis_LS[arm][pt] = temp_had_LS[arm][pt]->GetXaxis();
	      
	      x_array_UL[arm][pt][i] = xaxis_UL[arm][pt]->GetBinCenter(i+1);
	      x_array_LS[arm][pt][i] = xaxis_LS[arm][pt]->GetBinCenter(i+1);

	      x_array_width_UL[arm][pt][i] = xaxis_UL[arm][pt]->GetBinWidth(i+1);
	      x_array_width_LS[arm][pt][i] = xaxis_LS[arm][pt]->GetBinWidth(i+1);
	      
	      //cout << "x array UL: " << x_array_UL[arm][pt][i] << ", for [" << arm << "][" << pt << "][" << i << "]" << endl;
	      // cout << "x array LS: " << x_array_LS[arm][pt][i] << ", for [" << arm << "][" << pt << "][" << i << "]" << endl;
	      //cout << "x array width UL: " << x_array_width_UL[arm][pt][i] << ", for [" << arm << "][" << pt << "][" << i << "]" << endl;
	      // cout << "x array width LS: " << x_array_width_LS[arm][pt][i] << ", for [" << arm << "][" << pt << "][" << i << "]" << endl;
	    }
	}
    }

  // populate UL 0-1
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 0; pt < 1; pt++)
	    {  
	      for(int i = 10; i < 27; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_UL[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
  
  // populate UL 1-2
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 1; pt < 2; pt++)
	    {  
	      for(int i = 20; i < 37; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_UL[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
    
  // populate UL 2-3
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 2; pt < 3; pt++)
	    {  
	      for(int i = 10; i < 22; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_UL[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
    
  // populate UL 3-7
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 3; pt < 4; pt++)
	    {  
	      for(int i = 10; i < 22; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_UL[arm][pt]->GetBinContent(i+1); 
		  if(i == 15)
		    initial_binning[arm][i_histo][pt][i] =  ( temp_had_UL[arm][pt]->GetBinContent(i+3) + temp_had_UL[arm][pt]->GetBinContent(i+2) + initial_binning[arm][i_histo][pt][i-1] + initial_binning[arm][i_histo][pt][i-2] ) / 4;
		}
	    }
	}
    }
    
  // populate LS 0-1
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 1; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 1; pt++)
	    {  
	      for(int i = 10; i < 22; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_LS[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
    
  // populate LS 1-2
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 1; i_histo < 2; i_histo++)
	{
	  for(int pt = 1; pt < 2; pt++)
	    {  
	      for(int i = 10; i < 22; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_LS[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
    
  // populate LS 2-3
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 1; i_histo < 2; i_histo++)
	{
	  for(int pt = 2; pt < 3; pt++)
	    {  
	      for(int i = 9; i < 16; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_LS[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
   // populate LS 3-7
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 1; i_histo < 2; i_histo++)
	{
	  for(int pt = 3; pt < 4; pt++)
	    {  
	      for(int i = 9; i < 16; i++)
		{
		  initial_binning[arm][i_histo][pt][i] = temp_had_LS[arm][pt]->GetBinContent(i+1); 
		}
	    }
	}
    }
  
  // check the initial binning matrix

 for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      for(int i = 0; i < 100; i++)
		{
		  cout << "binning matrix: " << initial_binning[arm][i_histo][pt][i] << ", for [" << arm << "][" << i_histo << "][" << pt << "][" << i << "]" << endl;
		}
	    }
	}
    }
 //////////////////////////////////////////////////////////////////

 double final_binning[2][2][4][100] = {0};
 double temp1[2][2][4][100];
 double temp[2][100] = {0};
 double ave[2][100];

 // Now need to interpolate 

  // interpolate UL 0-1 Step 1 (3.2,3.4,3.6,3.8)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 0; pt < 1; pt++)
	    {  
	      for(int i = 10; i < 27; i++)
		{
		  for(int k = 20; k < 50; k++)
		    {
		      if( (k < 30))
			final_binning[arm][i_histo][pt][k] = initial_binning[arm][i_histo][pt][k-10]; // already in 100 mev bins

		      if(( i > 19) && (i < 24) )
			{
			  ave[arm][i] = ( initial_binning[arm][i_histo][pt][i] + initial_binning[arm][i_histo][pt][i + 1] )/2;
			  cout << ave[arm][i] << endl;
			}
		    }
		}
	    }
	}
    }
  
  // hard code this  UL 0-1
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 0; pt < 1; pt++)
	    {  
	      for(int i = 10; i < 27; i++)
		{
		  for(int k = 30; k < 39; k++)
		    {
		      
		      temp[arm][31] = ( initial_binning[arm][i_histo][pt][20] + ave[arm][20] )/2; // 3.15
		      temp[arm][32] = ( initial_binning[arm][i_histo][pt][21] + ave[arm][20] )/2; // 3.25
		      temp[arm][33] = ( initial_binning[arm][i_histo][pt][21] + ave[arm][21] )/2; // 3.35
		      temp[arm][34] = ( initial_binning[arm][i_histo][pt][22] + ave[arm][21] )/2; //3.45
		      temp[arm][35] = ( initial_binning[arm][i_histo][pt][22] + ave[arm][22] )/2; //3.55
		      temp[arm][36] = ( initial_binning[arm][i_histo][pt][23] + ave[arm][22] )/2; //3.65
		      temp[arm][37] = ( initial_binning[arm][i_histo][pt][23] + ave[arm][23] )/2;// 3.75
		      temp[arm][38] = ( initial_binning[arm][i_histo][pt][24] + ave[arm][23] )/2;// 3.85
		      temp[arm][30] = ( initial_binning[arm][i_histo][pt][19] + temp[arm][31] )/2; // 3.15

		      final_binning[arm][i_histo][pt][k] = temp[arm][k];
		      // cout << "final_binning[" << arm << "][" << i_histo << "][" << pt << "][" << k << "]: " <<  final_binning[arm][i_histo][pt][k] << endl;
		    }
		}
	    }
	}
    }

  ave[2][100] = {0};

  // interpolate UL 1-2 Step 1 (3.2,3.4,3.6,3.8)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 1; pt < 2; pt++)
	    {  
	      for(int i = 20; i < 37; i++)
		{
		  for(int k = 20; k < 50; k++)
		    {
		      if( (k < 30))
			final_binning[arm][i_histo][pt][k] = initial_binning[arm][i_histo][pt][k]; // already in 100 mev bins

		      if(( i > 29) && (i < 34) )
			{
			  ave[arm][i] = ( initial_binning[arm][i_histo][pt][i] + initial_binning[arm][i_histo][pt][i + 1] )/2;
			  cout << ave[arm][i] << endl;
			}
		    }
		}
	    }
	}
    }
  
  temp[2][100] = {0};

  // hard code this  UL 1-2
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 1; pt < 2; pt++)
	    {  
	      for(int i = 20; i < 37; i++)
		{
		  for(int k = 30; k < 39; k++)
		    {
		      
		      temp[arm][31] = ( initial_binning[arm][i_histo][pt][30] + ave[arm][30] )/2; // 3.15
		      temp[arm][32] = ( initial_binning[arm][i_histo][pt][31] + ave[arm][30] )/2; // 3.25
		      temp[arm][33] = ( initial_binning[arm][i_histo][pt][31] + ave[arm][31] )/2; // 3.35
		      temp[arm][34] = ( initial_binning[arm][i_histo][pt][32] + ave[arm][31] )/2; //3.45
		      temp[arm][35] = ( initial_binning[arm][i_histo][pt][32] + ave[arm][32] )/2; //3.55
		      temp[arm][36] = ( initial_binning[arm][i_histo][pt][33] + ave[arm][32] )/2; //3.65
		      temp[arm][37] = ( initial_binning[arm][i_histo][pt][33] + ave[arm][33] )/2;// 3.75
		      temp[arm][38] = ( initial_binning[arm][i_histo][pt][34] + ave[arm][33] )/2;// 3.85
		      temp[arm][30] = ( initial_binning[arm][i_histo][pt][29] + temp[arm][31] )/2; // 3.15
  
		      final_binning[arm][i_histo][pt][k] = temp[arm][k];
		      // cout << "final_binning[" << arm << "][" << i_histo << "][" << pt << "][" << k << "]: " <<  final_binning[arm][i_histo][pt][k] << endl;
		    }
		}
	    }
	}
    }


  ave[2][100] = {0};

  // UL 2-3 (2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 2; pt < 3; pt++)
	    {  
	      for(int i = 10; i < 20; i++)
		{
		  ave[arm][i] = ( initial_binning[arm][i_histo][pt][i] + initial_binning[arm][i_histo][pt][i + 1] )/2;
		  cout << ave[arm][i] << endl;
		}
	    }
	}
    }
  
  temp[2][100] = {0};

  // hard code this  UL 2-3
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 2; pt < 3; pt++)
	    {  
	      for(int i = 10; i < 21; i++)
		{
		  for(int k = 20; k < 39; k++)
		    {
		      
		      temp[arm][21] = ( initial_binning[arm][i_histo][pt][10] + ave[arm][10] )/2;// 2.15
		      temp[arm][22] = ( initial_binning[arm][i_histo][pt][11] + ave[arm][10] )/2;// 2.25
		      temp[arm][23] = ( initial_binning[arm][i_histo][pt][11] + ave[arm][11] )/2;// 2.35
		      temp[arm][24] = ( initial_binning[arm][i_histo][pt][12] + ave[arm][11] )/2;//2.45
		      temp[arm][25] = ( initial_binning[arm][i_histo][pt][12] + ave[arm][12] )/2;//2.55
		      temp[arm][26] = ( initial_binning[arm][i_histo][pt][13] + ave[arm][12] )/2;//2.65
		      temp[arm][27] = ( initial_binning[arm][i_histo][pt][13] + ave[arm][13] )/2;// 2.75
		      temp[arm][28] = ( initial_binning[arm][i_histo][pt][14] + ave[arm][13] )/2;// 2.85
		      temp[arm][29] = ( initial_binning[arm][i_histo][pt][14] + ave[arm][14] )/2;// 2.95
		      temp[arm][30] = ( initial_binning[arm][i_histo][pt][15] + ave[arm][14] )/2;// 3.05
		      temp[arm][31] = ( initial_binning[arm][i_histo][pt][15] + ave[arm][15] )/2;// 3.15
		      temp[arm][32] = ( initial_binning[arm][i_histo][pt][16] + ave[arm][15] )/2;// 3.25
		      temp[arm][33] = ( initial_binning[arm][i_histo][pt][16] + ave[arm][16] )/2;// 3.35
		      temp[arm][34] = ( initial_binning[arm][i_histo][pt][17] + ave[arm][16] )/2;// 3.45
		      temp[arm][35] = ( initial_binning[arm][i_histo][pt][17] + ave[arm][17] )/2;// 3.55
		      temp[arm][36] = ( initial_binning[arm][i_histo][pt][18] + ave[arm][17] )/2;// 3.65
		      temp[arm][37] = ( initial_binning[arm][i_histo][pt][18] + ave[arm][18] )/2;// 3.75
		      temp[arm][38] = ( initial_binning[arm][i_histo][pt][19] + ave[arm][18] )/2;// 3.85
		      
		      //  temp[arm][20] = ( initial_binning[arm][i_histo][pt][9] + temp[arm][21] )/2; // 2.05
  
		      final_binning[arm][i_histo][pt][k] = temp[arm][k];
		      // cout << "final_binning[" << arm << "][" << i_histo << "][" << pt << "][" << k << "]: " <<  final_binning[arm][i_histo][pt][k] << endl;
		    }
		}
	    }
	}
    }



  ave[2][100] = {0};

  // UL 3-7 (2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 3; pt < 4; pt++)
	    {  
	      for(int i = 10; i < 20; i++)
		{
		  ave[arm][i] = ( initial_binning[arm][i_histo][pt][i] + initial_binning[arm][i_histo][pt][i + 1] )/2;
		  cout << ave[arm][i] << endl;
		}
	    }
	}
    }
  
  temp[2][100] = {0};

  // hard code this  UL 3-7
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 3; pt < 4; pt++)
	    {  
	      for(int i = 10; i < 21; i++)
		{
		  for(int k = 20; k < 39; k++)
		    {
		      
		      temp[arm][21] = ( initial_binning[arm][i_histo][pt][10] + ave[arm][10] )/2;// 2.15
		      temp[arm][22] = ( initial_binning[arm][i_histo][pt][11] + ave[arm][10] )/2;// 2.25
		      temp[arm][23] = ( initial_binning[arm][i_histo][pt][11] + ave[arm][11] )/2;// 2.35
		      temp[arm][24] = ( initial_binning[arm][i_histo][pt][12] + ave[arm][11] )/2;//2.45
		      temp[arm][25] = ( initial_binning[arm][i_histo][pt][12] + ave[arm][12] )/2;//2.55
		      temp[arm][26] = ( initial_binning[arm][i_histo][pt][13] + ave[arm][12] )/2;//2.65
		      temp[arm][27] = ( initial_binning[arm][i_histo][pt][13] + ave[arm][13] )/2;// 2.75
		      temp[arm][28] = ( initial_binning[arm][i_histo][pt][14] + ave[arm][13] )/2;// 2.85
		      temp[arm][29] = ( initial_binning[arm][i_histo][pt][14] + ave[arm][14] )/2;// 2.95
		      temp[arm][30] = ( initial_binning[arm][i_histo][pt][15] + ave[arm][14] )/2;// 3.05
		      temp[arm][31] = ( initial_binning[arm][i_histo][pt][15] + ave[arm][15] )/2;// 3.15
		      temp[arm][32] = ( initial_binning[arm][i_histo][pt][16] + ave[arm][15] )/2;// 3.25
		      temp[arm][33] = ( initial_binning[arm][i_histo][pt][16] + ave[arm][16] )/2;// 3.35
		      temp[arm][34] = ( initial_binning[arm][i_histo][pt][17] + ave[arm][16] )/2;// 3.45
		      temp[arm][35] = ( initial_binning[arm][i_histo][pt][17] + ave[arm][17] )/2;// 3.55
		      temp[arm][36] = ( initial_binning[arm][i_histo][pt][18] + ave[arm][17] )/2;// 3.65
		      temp[arm][37] = ( initial_binning[arm][i_histo][pt][18] + ave[arm][18] )/2;// 3.75
		      temp[arm][38] = ( initial_binning[arm][i_histo][pt][19] + ave[arm][18] )/2;// 3.85
		      
		      //  temp[arm][20] = ( initial_binning[arm][i_histo][pt][9] + temp[arm][21] )/2; // 2.05
  
		      final_binning[arm][i_histo][pt][k] = temp[arm][k];
		      // cout << "final_binning[" << arm << "][" << i_histo << "][" << pt << "][" << k << "]: " <<  final_binning[arm][i_histo][pt][k] << endl;
		    }
		}
	    }
	}
    }



  ///////////////////////// print out final binning matrix modification coefficients ///////////////////////////////
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 1; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {  
	      for(int k = 0; k < 100; k++)
		{
		  cout << "final_binning[" << arm << "][" << i_histo << "][" << pt << "][" << k << "]: " <<  final_binning[arm][i_histo][pt][k] << endl;
		}
	    }
	}
    }
  
  



















} // end void macro

  /*





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
	  bg[3][arm][pt]->Scale(scale[arm]); // bb component of Run15 pp
	  bg[6][arm][pt]->Scale(scale[arm]); // dy compoennt of Run15pp
	  corr_had_UL[arm][pt]->Scale(scale[arm]);
	  corr_had_LS[arm][pt]->Scale(scale[arm]);
	}
    }

  // char uniquer[38][500];
 
  // int i = 0;
  
  //   for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  //   {
  //     i++;
  //     sprintf(uniquer[pt],"h_RpAu_%d",i); 
  //   }
  
  // set all errors in RpAu hadrons histograms to zero because they are too large
 
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < 4; pt++)
 	{
	  for(int i = 0; i < bins/2; i++)
	    {
	      //cout << "temp LS: " << temp_had_LS[arm][pt]->GetBinContent(i+1) << ", for arm: " << arm << ", pt: " << pt << ", and i: " << i << endl;
	      //RpAu_had_UL[arm][pt]->SetBinContent(i+1,temp_had_UL[arm][pt]->GetBinContent(i+1));
	      //RpAu_had_UL[arm][pt]->SetBinError(i+1, temp_had_UL[arm][pt]->GetBinContent(i+1)*0.01);
	    }
	}
    }

 cout << "HERE" << endl; 

  // RpAu histograms are binned in 200 Mev binwidths (50 bins: 0 -10 GeV/c).  Need 100 bins -> Can take average of two consective bins to get a new midpoint.  An array is fine here

  double RpAu_UL_array_final[2][38][100];
  double ave_UL_array[2][38][100] = {0};
  double RpAu_UL_array_initial[2][38][100];

  // For UL sign corr hadrons initial points
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < 2; i++)
	    {
	      RpAu_UL_array_initial[arm][pt][i] = RpAu_had_LS[arm][pt]->GetBinContent(i+1);
	      // ave_UL_array[arm][pt][i] = ( RpAu_UL_array_initial[arm][pt][i] + RpAu_UL_array_initial[arm][pt][i+1] ) / 2;
	      // RpAu_UL_array_final[arm][pt][2*i] = ave_UL_array[arm][pt][i]; // all even bins are 2k
	      // RpAu_UL_array_final[arm][pt][2*i+1] = RpAu_UL_array_initial[arm][pt][i];  // all odd bins are 2k+1
	    }
	}
    }
 // For LS sign corr hadrons average points
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {
	      cout << "RpAu final array for UL: " << RpAu_UL_array_final[arm][pt][i] << ", for element: " << i <<  endl;
	    }
	}
    }


  char uniqueh[38][500];
 
  int i = 0;
  
    for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      i++;
      sprintf(uniqueh[pt],"h_comp_%d",i); 
    }
  
  TH1D *h1 = new TH1D("h1", "cc distribution", 100, 0, 10);
  
  TH1D *h_comp[2][38];
    
  double cc_pAu[2][38][100];
  double bb_pAu[2][38][100];
  double dy_pAu[2][38][100];
  double corrhad_UL_pAu[2][38][100];
  double corrhad_LS_pAu[2][38][100];
   
  double quad[2][38][100];
  double quad_2[2][38][100];
  double quad_3[2][38][100];
  double quad_norm[2][38][100];
  double temp[2][100];
  double x_errors[2][100];
 
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {

	      // cc_pAu[arm][pt][bin] = bg[0][arm][pt]->GetBinContent(i+1)*( bg[2][arm][pt]->GetBinContent(i+1)/bg[1][arm][pt]->GetBinContent(i+1) ) *weight[arm][0] ;
	      // bb_pAu[arm][pt][bin] = bg[3][arm][pt]->GetBinContent(i+1)*( bg[5][arm][pt]->GetBinContent(i+1)/bg[4][arm][pt]->GetBinContent(i+1) ) *weight[arm][1];
	      // dy[arm][pt][bin] = bg[6][arm][pt]->GetBinContent(i+1);
	      // corrhad_UL_pAu[arm][pt][bin] = RpAu_had_UL[arm][pt]->GetBinContent(i+1)*corr_had_UL[arm][pt]->GetBinContent(i+1);
	      // corrhad_LS_pAu[arm][pt][bin]  = RpAu_had_LS[arm][pt]->GetBinContent(i+1)*corr_had_LS[arm][pt]->GetBinContent(i+1);

	      if( (bg[1][arm][pt]->GetBinContent(i+1) != 0) &&  (bg[4][arm][pt]->GetBinContent(i+1) != 0) )
	      	h1->SetBinContent(i+1, bg[0][arm][pt]->GetBinContent(i+1)*( bg[2][arm][pt]->GetBinContent(i+1)/bg[1][arm][pt]->GetBinContent(i+1) ) *weight[arm][0]  + bg[6][arm][pt]->GetBinContent(i+1) + bg[3][arm][pt]->GetBinContent(i+1)*( bg[5][arm][pt]->GetBinContent(i+1)/bg[4][arm][pt]->GetBinContent(i+1) ) *weight[arm][1] );
	      else
	      	h1->SetBinContent(i+1, 0);
	
	      // if( (bg[1][arm][pt]->GetBinContent(i+1) != 0) &&  (bg[4][arm][pt]->GetBinContent(i+1) != 0) )
	      // 	h1->SetBinContent(i+1, cc_pAu[arm][pt][i] + bb_pAu[arm][pt][i] + dy[arm][pt][i] + corrhad_UL_pAu[arm][pt][i] - corrhad_LS_pAu[arm][pt][i] );
	      // else
	      // 	h1->SetBinContent(i+1, 0);
	

	   
	      // comp_0[arm][pt][i] = bg[0][arm][pt]->GetBinContent(i+1);
	      // comp_1[arm][pt][i] = bg[1][arm][pt]->GetBinContent(i+1);
	      // comp_2[arm][pt][i] = bg[2][arm][pt]->GetBinContent(i+1);
	      // comp_3[arm][pt][i] = bg[3][arm][pt]->GetBinContent(i+1);
	      // comp_4[arm][pt][i] = bg[4][arm][pt]->GetBinContent(i+1);
	      // comp_5[arm][pt][i] = bg[5][arm][pt]->GetBinContent(i+1);
	      // comp_6[arm][pt][i] = bg[6][arm][pt]->GetBinContent(i+1);
	    
	      // comp_err_0[arm][pt][i] = (bg[0][arm][pt]->GetBinError(i+1)); 
	      // comp_err_1[arm][pt][i] = (bg[1][arm][pt]->GetBinError(i+1)); 
	      // comp_err_2[arm][pt][i] = (bg[2][arm][pt]->GetBinError(i+1));  
	      // comp_err_3[arm][pt][i] = (bg[3][arm][pt]->GetBinError(i+1)); 
	      // comp_err_4[arm][pt][i] = (bg[4][arm][pt]->GetBinError(i+1));  
	      // comp_err_5[arm][pt][i] = (bg[5][arm][pt]->GetBinError(i+1)); 
	     
     	      // quad[arm][pt][i] = ( pow(comp_err_0[arm][pt][i]/comp_0[arm][pt][i],2) +  pow(comp_err_1[arm][pt][i]/comp_1[arm][pt][i],2) +  pow(comp_err_2[arm][pt][i]/comp_2[arm][pt][i],2) )  * pow( (comp_0[arm][pt][i]*comp_2[arm][pt][i])/comp_1[arm][pt][i] ,2) ;

	   
	      if(h1->GetBinContent(i+1) > 0)
	     	{
		  h1->SetBinError(i+1,h1->GetBinContent(i+1)*0.05);
	     	}
	     
	      x_errors[arm][i] = 0.0;
	      
	      h_comp[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_comp[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));
	   	    	      
	    }
	}
    }

cout << "HERE" << endl;
 
  ////////////////////////////////////////////////
  // Histogram comp
  TCanvas *c6 = new TCanvas("c6","Corr BG comp ",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
   gPad->SetGrid();



  if(all_fits)
    {
      h_comp[0][0]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) ), South");
      h_comp[1][0]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) ), North");
      h_comp[arm_low][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_comp[arm_low][0]->GetXaxis()->SetLabelSize(0.04);
      h_comp[arm_low][0]->GetXaxis()->SetTitleSize(0.04);
      h_comp[arm_low][0]->GetXaxis()->SetTitleOffset(0.9);
      h_comp[arm_low][0]->GetYaxis()->SetLabelSize(0.04);
      h_comp[arm_low][0]->GetYaxis()->SetTitleSize(0.1); 
      h_comp[arm_low][0]->GetYaxis()->SetTitleOffset(0.52);
      h_comp[arm_low][0]->SetMinimum(0);
      if(south_arm == false)
	h_comp[arm_low][0]->SetMaximum(1000);
      else
	h_comp[arm_low][0]->SetMaximum(1500);
      h_comp[arm_low][0]->SetAxisRange(0,5);
    }
  else
    {
      h_comp[0][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) ), South");
      h_comp[1][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) + bb(UL) + dy(UL) ), North");
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_comp[arm_low][pt_l]->GetXaxis()->SetLabelSize(0.04);
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitleSize(0.04);
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitleOffset(0.9);
      h_comp[arm_low][pt_l]->GetYaxis()->SetLabelSize(0.04);
      h_comp[arm_low][pt_l]->GetYaxis()->SetTitleSize(0.1); 
      h_comp[arm_low][pt_l]->GetYaxis()->SetTitleOffset(0.52);
      h_comp[arm_low][pt_l]->SetMinimum(0);
      h_comp[arm_low][pt_l]->SetMaximum(1200);
      h_comp[arm_low][pt_l]->SetAxisRange(0,5);
    }

  for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
    {
      h_comp[arm_low][pt]->SetMarkerColor(pt+1);  
      if(pt+1==10)
	h_comp[arm_low][pt]->SetMarkerColor(30); 
      if(pt+1==2)
	h_comp[arm_low][pt]->SetMarkerColor(28); 
      if(pt+1 == 5)
	h_comp[arm_low][pt]->SetMarkerColor(40); 
      h_comp[arm_low][pt]->SetMarkerSize(0.7);
      h_comp[arm_low][pt]->SetMarkerStyle(20);

      if(all_fits)
	{
	  if(pt == 0)
	    h_comp[arm_low][pt]->Draw();
	  else
	    h_comp[arm_low][pt]->Draw("SAME");
	}
      else
	{
	  if(pt == pt_l)
	    h_comp[arm_low][pt_l]->Draw();
	}
    }

 /////////////////////////////////////////////////////////////////////////////////
 /// Make legend for h_comp histogram

  char unique_l[800];
  TLegend *leg_sum[2];
 
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
      leg_sum[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
      leg_sum[arm]->SetFillColor(0); 
      leg_sum[arm]->SetTextSize(0.035);
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  sprintf(unique_l,"cc + bb + dy corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_sum[arm]->AddEntry(h_comp[arm][pt],unique_l, "p");
	}
    }
  
  leg_sum[arm_low]->Draw();
  

}// void end macro




//   if(normalization)
//     h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,0));    
//   else
//     {
//       if(all_fits)
// 	{
// 	  h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-1),pow(10,4));
// 	  h_sum[arm_low][0]->SetAxisRange(0,5);
// 	}
//       else
// 	{
// 	  if(pt_binwidth == 2)
// 	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),1*pow(10,3));
// 	  //h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),180);
// 	  // h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),800);
// 	  if(pt_binwidth == 1)
// 	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),4*pow(10,2));
// 	  h_sum[arm_low][pt_l]->SetAxisRange(0,5);
// 	}
//     }
//   /////////////////////////////////////////////////////////////////////////////////
//   /// Make legend for histogram

//   char unique_l[800];
//   TLegend *leg_h[2];
 
//   double for_begin; 
//   double for_end; 

//   if(all_fits)
//     {
//       for_begin = 0;
//       for_end = pt_slices[pt_binwidth];
//     }
//   else
//     {
//       for_begin = pt_l;
//       for_end = pt_h;
//     }
 
//   for(int arm = 0; arm < 2; arm++)
//     {
//       leg_h[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
//       leg_h[arm]->SetFillColor(0); 
//       leg_h[arm]->SetTextSize(0.035);
	  	 
//       for(int pt = for_begin; pt < for_end; pt++)
// 	{
// 	  sprintf(unique_l,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
// 	  leg_h[arm]->AddEntry(h_sum[arm][pt],unique_l, "p");
// 	}
//     }

//   leg_h[arm_low]->Draw();





//   // pT Integrated corrletaed background plot
//   TGraphErrors *gr_corrbg_total = new TGraphErrors(bins,x_array,corr_total[arm_low],x_errors[arm_low],corr_total_err[arm_low]);
  
//   /*
//   TCanvas *c0 = new TCanvas("c0","Corr BG test pT Integrated",200,10,700,500);
//   gPad->SetGrid();
//   //gPad->SetLogy();
    
//   gr_corrbg_total->SetTitle("Corr BG Test pT Integrated");  
//   gr_corrbg_total->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
//   gr_corrbg_total->GetXaxis()->SetLabelSize(0.04);
//   gr_corrbg_total->GetXaxis()->SetTitleSize(0.04);
//   gr_corrbg_total->GetXaxis()->SetTitleOffset(0.9);
//   gr_corrbg_total->GetXaxis()->SetLimits(0,5);
//   gr_corrbg_total->GetYaxis()->SetLabelSize(0.04);
//   gr_corrbg_total->GetYaxis()->SetTitleSize(0.1); 
//   gr_corrbg_total->GetYaxis()->SetTitleOffset(0.52);
//   gr_corrbg_total->SetMaximum(pow(10,5));
//   gr_corrbg_total->SetMinimum(pow(10,2));

//   /// assign different colors and draw each bg
//   gr_corrbg_total->SetMarkerColor(kCyan);  
//   gr_corrbg_total->SetMarkerSize(0.75);
//   gr_corrbg_total->SetMarkerStyle(20);
//   gr_corrbg_total->Draw("AP");
  
//  TLegend *leg0;
//  leg0 = new TLegend(0.51, 0.6, 0.8, 0.9);  
//  leg0->SetFillColor(0); 
//  leg0->SetTextSize(0.035);
//  leg0->AddEntry(gr_corrbg_total,"corr bg like-sign test", "p");
 
//  leg0->Draw();
// ////////
// */

//   //////////////////////////// Fit histograms //////////////////////////
  
//   char name1[800];
//   char name2[800];

//   double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;
//   int bins_fail[bins];

//   TF1 *corr_bg_fit;
//   TF1 *corr_total_fit;
//   double f5,f6,f10,f11,f12;
//   TFitResultPtr r;

//   // For 500 MeV bins
//   double par_0[14] = {1,1,1,1,1,1,1,1,1,1,1,0.1,1,1};  // original
//   double par_1[14] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
//   double par_2[14] = {10,10,100,100,100,100,100,10,10,10,1,1,1,1};
//   double par_3[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10}; 
//   double par_4[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10};  // original

//   double par_2_refit[14] = {10,10,100,100,100,100,100,10,10,100,1,1,1,1};   //10
 
//   double limit_a_refit[14] = {1.4,1.4,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4};   //(10,0.01,10,100,0.01)
//   //double limit_a_refit[14] = {1.4,1.4,2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4};   //(10,0.01,10,100,0.01)

//   double limit_a[14] = {1.75,1.75,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.75,0.75,0.75,0.75};
//   double limit_b[14] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5};

//       // For 200 MeV bins (only interested in first 6 bins: 0.0 - 1.2 GeV/c, which is low pT corr BG shape)
//   double par_0_200[6] = {10,10,10,10,10,10};
//   double par_1_200[6] = {0.01,0.01,0.01,0.01,0.01,0.01};
//   double par_2_200[6] = {1,1,1,1,1,1};
//   double par_3_200[6] = {100,100,100,100,100,100};
//   double par_4_200[6] = {0.01,0.01,0.01,0.01,0.01,0.01};
 
//   double limit_a_200[6] = {1.4,1.4,1.4,1.4,1.4,1.4};  // for south arm


//   //double limit_a_200_North[6] = {1.5,1.4,1.5,1.4,1.5,1.5}; // for noth arm
//   // double limit_a_200_North[6] = {1.45,1.45,1.45,1.45,1.45,1.45}; // for noth ar
//   // double limit_a_200_North[6] = {1.4,1.75,1.25,1.25,1.25,1.25};  // all converge at these limits
//   double limit_a_200_North[6] = {1.4,1.4,1.4,1.4,1.4,1.4};


//   double limit_b_200[6] = {5,5,5,5,5,5};

//   // // For 300 MeV bins
//   double par_0_300[21] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,1,1,1,1,1};
//   double par_1_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
//   double par_2_300[21] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
//   double par_3_300[21] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,10,10,10,10,10};
//   double par_4_300[21] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,10,10,10,10,10};
  
//   double limit_a_300[21] = {1.4,1.4,1.4,1.4,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.85,1.4,0.85}; 
//   //double limit_a_300[21] = {1.7,1.7,1.7,1.7,1.7,1.7,1.7,0.25,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0}; 
//   double limit_b_300[21] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};

//   // For 100 MeV bins (only interested in first 5 bins: 0.0 - 1.0 GeV/c, which is low pT corr BG shape)
//   double par_0_100[10] = {10,10,10,10,10,10,10,10,10,10}; // only try 400 - 800
//   double par_1_100[10] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
//   double par_2_100[10] = {1,1,1,1,1,1,1,1,10,10};
//   double par_3_100[10] = {100,100,100,100,100,100,100,100,100,100};
//   double par_4_100[10] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
 
//   double limit_a_100[10] = {1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.6};
//   double limit_b_100[10] = {5,5,5,5,5,5,5,5,5,5};
  
//   // For 1 GeV bins
//   double par_0_1000[10] = {0,0,0,0,1,1,1,0,0,0};
//   double par_1_1000[10] = {0,0,0,0,0.01,0.01,0.01,0,0,0};
//   double par_2_1000[10] = {0,0,0,0,1,1,1,0,0,0};
//   double par_3_1000[10] = {0,0,0,0,10,10,10,0,0,0};
//   double par_4_1000[10] = {0,0,0,0,10,10,10,0,0,0};
  
//   double limit_a_1000[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 
//   double limit_b_1000[10] = {0,0,0,0,5,5,5,0,0,0};

//   double limit_a_1000_North[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 

//   // For 2 GeV bins
//   double par_0_2000[5] = {1,1,1,1,1};
//   double par_1_2000[5] = {0.01,0.01,0.01,0.01,0.01};
//   double par_2_2000[5] = {1,1,1,1,1};
//   double par_3_2000[5] = {10,10,10,10,10};
//   double par_4_2000[5] = {10,10,10,10,10};
  
//   double limit_a_2000[5] = {0.75,0.25,0.75,0.75,0.75};
//   double limit_b_2000[5] = {5,5,5,5,5};





//  double refit_par_values[2][5][38];

//  // for pT Integrated

//  double par_0_total = 1;
//  double par_1_total = 0.01;
//  double par_2_total = 100;
//  double par_3_total = 10;
//  double par_4_total = 10;

 

//  //////////////////////////////////////////////////////////////////////////  All code above this point is for both North and South 
  

//  for(int arm = arm_low; arm < arm_high; arm++)  
//    // for(int arm = 0; arm < 2; arm++)  
//     {
//       for(int pt = for_begin; pt < for_end; pt++) // ALL FITS TO DISPLAY ON SAME PLOT
// 	{

// 	  if(refit && pt_binwidth == 2)
// 	    {
// 	      par_2[pt] = par_2_refit[pt];
// 	      limit_a[pt] = limit_a_refit[pt];
// 	    }
	      
// 	  if(pt_binwidth == 1)
// 	    {
// 	      par_0[pt] = par_0_200[pt];
// 	      par_1[pt] = par_1_200[pt];
// 	      par_2[pt] = par_2_200[pt];
// 	      par_3[pt] = par_3_200[pt];
// 	      par_4[pt] = par_4_200[pt];
// 	      if(south_arm)
// 		limit_a[pt] = limit_a_200[pt];
// 	      else
// 		limit_a[pt] = limit_a_200_North[pt];
// 	      limit_b[pt] = limit_b_200[pt];
// 	    }
	  
// 	  if(pt_binwidth == 3)
//  	    {
// 	      par_0[pt] = par_0_1000[pt];
// 	      par_1[pt] = par_1_1000[pt];
// 	      par_2[pt] = par_2_1000[pt];
// 	      par_3[pt] = par_3_1000[pt];
// 	      par_4[pt] = par_4_1000[pt];
	      
// 	      limit_a[pt] = limit_a_1000[pt];
// 	      limit_b[pt] = limit_b_1000[pt];

// 	      if(south_arm == false)
// 		limit_a[pt] = limit_a_1000_North[pt];

// 	    }
// 	  if(pt_binwidth == 0)
//  	    {
// 	      par_0[pt] = par_0_100[pt];
// 	      par_1[pt] = par_1_100[pt];
// 	      par_2[pt] = par_2_100[pt];
// 	      par_3[pt] = par_3_100[pt];
// 	      par_4[pt] = par_4_100[pt];
	      
// 	      limit_a[pt] = limit_a_100[pt];
// 	      limit_b[pt] = limit_b_100[pt];
// 	    }
// 	  if(pt_binwidth == 0) //4
//  	    {
// 	      par_0[pt] = par_0_300[pt];
// 	      par_1[pt] = par_1_300[pt];
// 	      par_2[pt] = par_2_300[pt];
// 	      par_3[pt] = par_3_300[pt];
// 	      par_4[pt] = par_4_300[pt];
	      
// 	      limit_a[pt] = limit_a_300[pt];
// 	      limit_b[pt] = limit_b_300[pt];
// 	    }
// 	  if(pt_binwidth == 4)
// 	    {
// 	      par_0[pt] = par_0_2000[pt];
// 	      par_1[pt] = par_1_2000[pt];
// 	      par_2[pt] = par_2_2000[pt];
// 	      par_3[pt] = par_3_2000[pt];
// 	      par_4[pt] = par_4_2000[pt];
	      
// 	      limit_a[pt] = limit_a_2000[pt];
// 	      limit_b[pt] = limit_b_2000[pt];
// 	    }
	  	  
// 	  corr_bg_fit = new TF1("corr_bg_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",limit_a[pt],limit_b[pt]);

// 	  corr_total_fit = new TF1("corr_total_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,5);
	  	  
// 	  cout <<" TF1 fit is from " << limit_a[pt] << "  to " << limit_b[pt] << endl;

// 	  // For fitting bg as function of pT and then plotting parameter values as function of pT
// 	  //**************************************************************************************
// 	  if(initial_fit == true)
// 	    {
// 	      corr_bg_fit->SetParameter(0, par_0[pt]);
// 	      corr_bg_fit->SetParameter(1, par_1[pt]);
// 	      corr_bg_fit->SetParameter(2, par_2[pt]);
// 	      corr_bg_fit->SetParameter(3, par_3[pt]);
// 	      corr_bg_fit->SetParameter(4, par_4[pt]);

// 	      corr_total_fit->SetParameter(0, par_0_total);
// 	      corr_total_fit->SetParameter(1, par_1_total);
// 	      corr_total_fit->SetParameter(2, par_2_total);
// 	      corr_total_fit->SetParameter(3, par_3_total);
// 	      corr_total_fit->SetParameter(4, par_4_total);
// 	    }	 

// 	  // for calculating the parameter from a function plugging in the pt value
// 	  //************************************************************************

// 	  double par_a[10];
// 	  double par_b[10];
// 	  double par_d[10];
// 	  double par_e[10];
	  
// 	  double par_a_fx;
// 	  double par_b_fx;
// 	  double par_d_fx;
// 	  double par_e_fx;
	  
// 	  double x;
	  

// 	  x = pt_center[pt_binwidth][pt];
// 	  cout << "x = pT center: " << x << endl;
	  
// 	  if(refit == true)
// 	    {
  
// 	      f10 = par_2[pt];

// 	      cout << "f10: " << f10 << endl;
	      
	      
// 	      for(int par = 0; par < 10; par++)
// 		{	 
// 		  if(south_fx)
// 		    {
// 		      par_a[par] = coeff_par_S[0][par];
// 		      par_b[par] = coeff_par_S[1][par];
// 		      par_d[par] = coeff_par_S[2][par];
// 		      par_e[par] = coeff_par_S[3][par];
// 		    }
// 		  else
// 		    {
// 		      par_a[par] = coeff_par_N[0][par];
// 		      par_b[par] = coeff_par_N[1][par];
// 		      par_d[par] = coeff_par_N[2][par];
// 		      par_e[par] = coeff_par_N[3][par];
		  
// 		    }
// 		}
	      
// 	      par_a_fx = par_a[0] + par_a[1]*x + par_a[2]*x*x + par_a[3]*x*x*x + par_a[4]*x*x*x*x + par_a[5]*x*x*x*x*x + par_a[6]*x*x*x*x*x*x + par_a[7]*x*x*x*x*x*x*x + par_a[8]*x*x*x*x*x*x*x*x + par_a[9]*x*x*x*x*x*x*x*x*x;
// 	      par_b_fx = par_b[0] + par_b[1]*x + par_b[2]*x*x + par_b[3]*x*x*x + par_b[4]*x*x*x*x + par_b[5]*x*x*x*x*x + par_b[6]*x*x*x*x*x*x + par_b[7]*x*x*x*x*x*x*x + par_b[8]*x*x*x*x*x*x*x*x + par_b[9]*x*x*x*x*x*x*x*x*x;
// 	      par_d_fx = par_d[0] + par_d[1]*x + par_d[2]*x*x + par_d[3]*x*x*x + par_d[4]*x*x*x*x + par_d[5]*x*x*x*x*x + par_d[6]*x*x*x*x*x*x + par_d[7]*x*x*x*x*x*x*x + par_d[8]*x*x*x*x*x*x*x*x + par_d[9]*x*x*x*x*x*x*x*x*x;
// 	      par_e_fx = par_e[0] + par_e[1]*x + par_e[2]*x*x + par_e[3]*x*x*x + par_e[4]*x*x*x*x + par_e[5]*x*x*x*x*x + par_e[6]*x*x*x*x*x*x + par_e[7]*x*x*x*x*x*x*x + par_e[8]*x*x*x*x*x*x*x*x + par_e[9]*x*x*x*x*x*x*x*x*x;
	      
// 	      cout << "a0: " << par_a[0] << ", a1: " << par_a[1] << ", a2: " << par_a[2] << ", a3: " << par_a[3] << ", a4: " << par_a[4] << ", a5: " << par_a[5] << ", a6: " << par_a[6] << ", a7: " << par_a[7] << ", a8: " << par_a[8] << ", a9: " << par_a[9] << endl;
	      
// 	      cout << "b0: " << par_b[0] << ", b1: " << par_b[1] << ", b2: " << par_b[2] << ", b3: " << par_b[3] << ", b4: " << par_b[4] << ", b5: " << par_b[5] << ", b6: " << par_b[6] << ", b7: " << par_b[7] << ", b8: " << par_b[8] << ", b9: " << par_b[9] << endl;
	      
// 	      cout << "d0: " << par_d[0] << ", d1: " << par_d[1] << ", d2: " << par_d[2] << ", d3: " << par_d[3] << ", d4: " << par_d[4] << ", d5: " << par_d[5] << ", d6: " << par_d[6] << ", d7: " << par_d[7] << ", d8: " << par_d[8] << ", d9: " << par_d[9] << endl;
	      
// 	      cout << "e0: " << par_e[0] << ", e1: " << par_e[1] << ", e2: " << par_e[2] << ", e3: " << par_e[3] << ", e4: " << par_e[4] << ", e5: " << par_e[5] << ", e6: " << par_e[6] << ", e7: " << par_e[7] << ", e8: " << par_e[8] << ", e9: " << par_e[9] << endl;
	      
// 	      corr_bg_fit->FixParameter(0, par_a_fx);
// 	      corr_bg_fit->FixParameter(1, par_b_fx);
// 	      corr_bg_fit->SetParameter(2, f10);
// 	      corr_bg_fit->FixParameter(3, par_d_fx);
// 	      corr_bg_fit->FixParameter(4, par_e_fx);
  
// 	    }
	      
// 	  //****************************************************************************
	      
// 	  corr_bg_fit->SetLineColor(kRed);
// 	  corr_bg_fit->SetLineStyle(5);
// 	  corr_bg_fit->SetLineWidth(2);

// 	  corr_total_fit->SetLineColor(kRed);
// 	  corr_total_fit->SetLineStyle(5);
// 	  corr_total_fit->SetLineWidth(2);
	      
// 	  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
	     
// 	  //r = gr_corr_pt[pt]->Fit(corr_bg_fit,"S R");	  // for TGraphErrors


// 	  if(pt_int_fit == false)
// 	    {
// 	      r = h_sum[arm_low][pt]->Fit(corr_bg_fit,"S R");
// 	      cout << "Chi Square/NDF: " << corr_bg_fit->GetChisquare() << "/" << corr_bg_fit->GetNDF() << "" << endl;
// 	    }
// 	  else
// 	    {
// 	      r = gr_corrbg_total->Fit(corr_total_fit,"S R");
// 	      cout << "Chi Square/NDF: " << corr_total_fit->GetChisquare() << "/" << corr_total_fit->GetNDF() << "" << endl;
// 	    }
	  
// 	  p0 = r->Parameter(0);
// 	  p1 = r->Parameter(1);
// 	  p2 = r->Parameter(2);
// 	  p3 = r->Parameter(3);
// 	  p4 = r->Parameter(4);


// 	  // fill array with refit parameter values calculated by functions
// 	  refit_par_values[arm][0][pt] = p0;
// 	  refit_par_values[arm][1][pt] = p1;
// 	  refit_par_values[arm][2][pt] = p2;
// 	  refit_par_values[arm][3][pt] = p3;
// 	  refit_par_values[arm][4][pt] = p4;
	  
// 	  e0 = r->ParError(0);
// 	  e1 = r->ParError(1);
// 	  e2 = r->ParError(2);
// 	  e3 = r->ParError(3);
// 	  e4 = r->ParError(4);
	    
// 	  Int_t fitStatus = r;  


// 	  cout << "status: " << fitStatus << " for pt bin " << pt+1 << endl;

// 	  if(fitStatus == 0) 
// 	    {
// 	      if(write)	     
// 		{
// 		  if(arm == 1)
// 		    {
// 		      if(pt_binwidth == 2)
// 			sprintf(name1,"yuehang_bestfit_500mev/bestfit_parameters_N_%i.dat",pt+1); 
// 		      if(pt_binwidth == 1)
// 			sprintf(name1,"yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_%i.dat",pt+1); 
// 		      if(pt_binwidth == 3)
// 			sprintf(name1,"yuehang_bestfit_1gev/bestfit_parameters_high_pt_N_%i.dat",pt+1); 
// 		      if(pt_binwidth == 0)
// 			sprintf(name1,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_N_%i.dat",pt+1);
// 		      if(pt_binwidth == 4)
// 			sprintf(name1,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_N_%i.dat",pt+1);
   
// 		      std::fstream bestfit_parameters_1(name1,std::ofstream::out); 
// 		      bestfit_parameters_1 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
// 		      bestfit_parameters_1.close();
// 		    }
// 		  else
// 		    {
// 		      if(pt_binwidth == 2)
// 			sprintf(name2,"yuehang_bestfit_500mev/bestfit_parameters_S_%i.dat",pt+1); 
// 		      if(pt_binwidth == 1)
// 			sprintf(name2,"yuehang_bestfit_200mev/bestfit_parameters_low_pt_S_%i.dat",pt+1); 	
// 		      if(pt_binwidth == 3)
// 			sprintf(name2,"yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);
// 		      if(pt_binwidth == 0)
// 			sprintf(name2,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_%i.dat",pt+1);
// 		      if(pt_binwidth == 4)
// 			sprintf(name2,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);

 	      
// 		      std::fstream bestfit_parameters_2(name2,std::ofstream::out); 
// 		      bestfit_parameters_2 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
// 		      bestfit_parameters_2.close();
// 		    }
// 		} // write
	     
// 	    }// fit converged
// 	  else
// 	    {
// 	      counter+= 1;
// 	      cout << "number of fits that failed : " << counter << endl;
// 	    }
	  
// 	} // for loop pt
      
//     } // for loop arm
  
   

//  cout << "HERE" << endl;
  
//   Char_t message[200];
//   //  Char_t message2[200];
  
//   if(pt_int_fit == false)
//     sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_bg_fit->GetChisquare(),corr_bg_fit->GetNDF());
//   else
//     sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_total_fit->GetChisquare(),corr_total_fit->GetNDF());
//   // sprintf(message2,"J/Psi Counts =  %.1f +/- %.3f",NJpsi,err_NJpsi);
  
//   TPaveText *mytext = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
//   mytext->SetTextSize(0.035);
//   mytext->SetFillColor(0);
//   mytext->SetTextAlign(12);
//   mytext->AddText(message);
//   // mytext->AddText(message2);
//   mytext->Draw();












  




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
  

 
//}// void end macro



      
  
