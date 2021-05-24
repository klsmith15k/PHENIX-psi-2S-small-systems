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

void yuehang_pt_integrated()
{
  // for f(x) parameter functions
#include "yuehang_fit_coeff_500/fit_coeff_S_array_twelve_point.C"    // 4-4.5,4.5-5,5-6,6-7 ---> 12 points
#include "yuehang_fit_coeff_500/fit_coeff_N_array_twelve_point.C"  // 4-4.5,4.5-5,5-6,6-7 ---> 12 points

//for each fit evaluated at each mass bin
#include "yuehang_pt_int/fit_array_200mev_N.C" 
#include "yuehang_pt_int/fit_array_200mev_S.C" 
  #include "yuehang_pt_int/fit_array_500mev_N.C" 
#include "yuehang_pt_int/fit_array_500mev_S.C" 
  #include "yuehang_pt_int/fit_array_1gev_N.C" 
#include "yuehang_pt_int/fit_array_1gev_S.C" 

  bool south_arm = false;
  //bool south_arm = true;
  bool south_fx = false; // need to set to true ONLY if refitting south arm
  bool pt_int_fit = false;

  bool initial_fit = false;  // set to false for ratio plots
  bool refit = true;  // set to true for ratio plots
  bool normalization = false;
  bool all_fits = true;  // must be set to true to generate parameter ratio plots at end of macro
  bool write = false;

 
 
  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; 
  
  int pt_binwidth;
  cout << "What pt binning are you running?  Enter '0' for 100 MeV binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth" << endl;
  cin >> pt_binwidth;
  
  int pt_slices[5] = {10,6,14,10,5}; 
  
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

  TAxis *xaxis;
  TAxis *yaxis;	  
  
  int bins = 100;
  
  double pt_low;
  double pt_high;
  
  double bin_low;
  double bin_high;
  
  double delta_pt = 0.001;
  double a_array[38];
  double b_array[38];
   
  char basename2[8][500] = {"yuehang_pt_int/weighted_fits_100mev_S_","yuehang_pt_int/weighted_fits_200mev_S_","yuehang_pt_int/weighted_fits_500mev_S_","yuehang_pt_int/weighted_fits_1gev_S_",
			    "yuehang_pt_int/weighted_fits_100mev_N_","yuehang_pt_int/weighted_fits_200mev_N_","yuehang_pt_int/weighted_fits_500mev_N_","yuehang_pt_int/weighted_fits_1gev_N_"};
		  
  vector < vector < vector <std::string> > > bestfit_filename;
    
  for(int j_arm = 0; j_arm < 2; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm;
  
      for(int j_series = 0; j_series < 4; j_series++)
	{
	  vector < std::string >  filename_series;
	  int j = j_arm*4 + j_series;
	  // cout << "j: " << j << endl;
	  
	  for(int i = 0; i < pt_slices[j_series]; i++)
	    {
	      char name[500];
	      sprintf(name,"%s%i%s",basename2[j],i+1,".dat");	     
	      std::string filename_tmp = name;
	      filename_series.push_back(filename_tmp);
	    } 

	  filename_arm.push_back(filename_series);
	}

      bestfit_filename.push_back(filename_arm);

    } 

  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = 0; series < 4; series++)
	{
	  for(int pt = 0; pt < pt_slices[series]; pt++)
	    {
	      cout << bestfit_filename[arm][series][pt].c_str() << endl; 
	    }
	}
    }


  double area_corrbg[2][14];
  int counter1 = 0;
  int counter2 = 0;
  
  for(int arm = 0; arm < 2; arm++)
    {
      counter1 = 4; 
      counter2 = 12;
     
      for(int series = 0; series < 4; series++) 
	{
	  for(int pt = 0; pt < 14; pt++)  
	    {
	     		 
	      ifstream bestfit_parameters(bestfit_filename[arm][series][pt].c_str()); 
	      
	      if(bestfit_parameters)
		{
		  double a0;
		  bestfit_parameters >> a0;
		  if(series == 1 && pt < 5) 
		    {
		      if(pt < 2)
			area_corrbg[arm][pt] = a0;
		      if( (pt == 3) || (pt == 4) )
			area_corrbg[arm][pt-1] = a0; 
		    }
		  if( (series == 2) && (pt < 10) )
		    {
		      if( (pt > 1) && (pt < 10) ) 
			{
			  area_corrbg[arm][counter1] = a0;
			  counter1++;
			}
		    }
		  if((series == 3) && ( (pt == 5) ||  (pt == 6) ) )
		    {
		      area_corrbg[arm][counter2] = a0;
		      counter2++;
		    } 

		  bestfit_parameters.close();
		}
	      else
		continue;
	     	     
	    }// for pt

	} // for loop series
      
      counter1 = 0;
      counter2 = 0;
    } // for loop arm
  
  // counter1 = 0;
  // counter2 = 0;

  // for(int arm = 1; arm < 2; arm++)
  //   {
  //     counter1 = 4; 
  //     counter2 = 12;
     
  //     for(int series = 0; series < 4; series++) 
  // 	{
  // 	  for(int pt = 0; pt < 14; pt++)  
  // 	    {
	     		 
  // 	      ifstream bestfit_parameters(bestfit_filename[arm][series][pt].c_str()); 
	      
  // 	      if(bestfit_parameters)
  // 		{
  // 		  double a0;
  // 		  bestfit_parameters >> a0;
  // 		  if(series == 1 && pt < 5) 
  // 		    {
  // 		      if(pt < 2)
  // 			area_corrbg[arm][pt] = a0;
  // 		      if( (pt == 3) || (pt == 4) )
  // 			area_corrbg[arm][pt-1] = a0; 
  // 		    }
  // 		  if( (series == 2) && (pt < 10) )
  // 		    {
  // 		      if( (pt > 1) && (pt < 10) ) 
  // 			{
  // 			  area_corrbg[arm][counter1] = a0;
  // 			  counter1++;
  // 			}
  // 		    }
  // 		  if((series == 3) && ( (pt == 5) ||  (pt == 6) ) )
  // 		    {
  // 		      area_corrbg[arm][counter2] = a0;
  // 		      counter2++;
  // 		    } 

  // 		  bestfit_parameters.close();
  // 		}
  // 	      else
  // 		continue;
	     	     
  // 	    }
  // 	} // for loop arm
      
  //   } // for series loop
  

  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < 14; k++)
	{
	  cout << "arm " << arm << ", bin: " << k << endl;
	  cout << area_corrbg[arm][k] << endl;
	}
    }



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  std::string filename[5][2] = {"yuehang_Run15pp_S_rebinned_100_mev.root", 
				"yuehang_Run15pp_N_rebinned_100_mev.root",
				
				"yuehang_Run15pp_S_rebinned_200_mev.root", 
				"yuehang_Run15pp_N_rebinned_200_mev.root",
				
				"yuehang_Run15pp_S_rebinned_500_mev.root", 
				"yuehang_Run15pp_N_rebinned_500_mev.root",
				
				"yuehang_Run15pp_S_rebinned_1000_mev.root", 
				"yuehang_Run15pp_N_rebinned_1000_mev.root",
				
				"yuehang_Run15pp_S_rebinned_2_gev.root", 
				"yuehang_Run15pp_N_rebinned_2_gev.root"};
 

  //////////////////////////////////////////////////
  //everything below this line is created as a result of making an initial selection about which data you want to use
  
  TH1D *bg[6][2][38];
  TH1D *minbias[5][2][38];

  TH1D *t_check[2]; 
  TH2F *yuehang_2d = 0;
  
  std::string obj_filename[6][2][38]; 
  
  char unique1[800];
  char unique2[800];
  char unique3[800];
  char unique4[800];
  char unique5[800];
  char unique6[800];
  
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < pt_slices[pt_binwidth]; k++)
	{ 
	  pt_low = pt_center[pt_binwidth][k] - pt_width[pt_binwidth][k]/2 + delta_pt; 
	  pt_high = pt_center[pt_binwidth][k] + pt_width[pt_binwidth][k]/2 - delta_pt;
  
	  double a = (pt_low - delta_pt)*1000; 
	  double b = (pt_high + delta_pt)*1000;

	  sprintf(unique1,"cc_%.0f_%.0f",a,b); 
	  sprintf(unique2,"bb_UL_%.0f_%.0f",a,b); 
	  sprintf(unique3,"dy_%.0f_%.0f",a,b);  
	  sprintf(unique4,"bb_LS_%.0f_%.0f",a,b);
	  sprintf(unique5,"corr_had_LS_%.0f_%.0f",a,b);
	  sprintf(unique6,"corr_had_UL_%.0f_%.0f",a,b);
	  
	  a_array[k] = a/1000;
	  b_array[k] = b/1000;

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
	      //cout << obj_filename[i_histo][arm][pt].c_str() << " HERE " << endl;
	      //cout << obj_filename2[i_histo][arm][pt].c_str() << " HERE " << endl;
	    }
	}
    }
  
  TFile *file_data;
  TFile *file_weights;
  TFile *file_check;
  
  std::string filename_ch[2] = {"for_krista_mass_arm0_12_5.root",
				"for_krista_mass_arm1_12_5.root"};
   
  double x_array[100];
  double x_array_mass_ratio[100];
  double x_low[100];
  double x_high[100];
  int bin_two;
  int bin_five;  
   
  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Int_t numBinsX;
 
  Double_t binCenterY;
  Double_t binLowY;
  Double_t binWidthY;
  Int_t numBinsY;

  double bin_35gev;
  double bin_4gev;
  double bin_5gev;
  

  // for rebinned data files
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
	      file_check->GetObject("h_cc_mass_pt",yuehang_2d);
	      
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
	     
	      // bin_two = xaxis->FindBin(2.001);
	      // bin_five = xaxis->FindBin(4.999);
	      // cout << " ______________________________ " << endl;
	      // cout << obj_filename[i_histo][arm][pt].c_str() << endl;
	      // cout << "x axis starts at: " << binLowX <<  endl;
	      // cout << "x binwidth: " << binWidthX <<  endl;
	      // cout << "x axis bin 1 center: " << binCenterX << endl;
	      // cout << "x axis total bins: " << numBinsX << endl;
	      // cout << " ______________________________ " << endl;

	      // yaxis = yuehang_2d->GetYaxis();
	      // binCenterY = yaxis->GetBinCenter(1);
	      // binLowY = yaxis->GetBinLowEdge(1);
	      // binWidthY = yaxis->GetBinWidth(1);
	      // numBinsY = yaxis->GetNbins();

	      // cout << " ______________________________ " << endl;
	      // cout << "y axis starts at: " << binLowY <<  endl;
	      // cout << "y binwidth: " << binWidthY <<  endl;
	      // cout << "y axis bin 1 center: " << binCenterY << endl;
	      // cout << "y axis total bins: " << numBinsY << endl;
	      // cout << " ______________________________ " << endl;
	      
	    }
	}
    }
    
  xaxis = bg[0][1][1]->GetXaxis();  
  int bin_counter = 0;

  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
  
      // cout << "bin : " << i+1 << ", and bin 3.5: " << bin_35gev << endl;
       // cout << "mass ratio array: " << x_array_mass_ratio[i] << " for bin " << i+1 << endl;
    }

  // for(int arm = 0; arm < 2; arm++)
  //   {
  //     for(int i_histo = 0; i_histo < 4; i_histo++)
  //       {
  // 	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  // 	    {
  // 	      cout << obj_filename[i_histo][arm][pt].c_str() << ", i_histo: " << i_histo << ", arm: " << arm << ", pt: " << pt << endl;
  // 	    }
  // 	} 
  //   }

  TH1D *h_sum[2][38];
  TH1D *h_int[2];
 
  double cc[2][38][100];
  double bb_UL[2][38][100];
  double dy[2][38][100];
  double bb_LS[2][38][100];
  double corr_had_LS[2][38][100];
  double corr_had_UL[2][38][100];
  double corr_bg[2][38][100];
  double corr_total[2][100] = {0.0,0.0};
  double corr_total_err[2][100] = {0.0,0.0};
  double corr_total_norm[2][100] = {0.0,0.0};
  double corr_bg_err[2][38][100] = {0.0,0.0};
  double pt_int[2][38];
  double corr_mass[2][100];
  double errors[2][38][100];
  double x_errors[2][100];
  double y_errors[2][100];
  double counter = 0.0;
  double mass_ratio[2][38][100];
	   
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
    
  int pt_l = 6;
  int pt_h = 7;
  //double lim_a = 0;
  //double lim_b = 0;
  

  double per_err[14] = {0.11,0.0537,0.0989,0.10,0.11,0.105,0.11,0.14,0.19,0.245,0.36,0.419,0.395,0.40};
  //double per_err[14] = {0.11,0.0537,0.0989,0.1119,0.11,0.105,0.11,0.14,0.19,0.245,0.36,0.419,0.395,0.40};
 double per_err_North[14] = {0.15,0.0537,0.086,0.121,0.114,0.127,0.138,0.201,0.189,0.216,0.395,0.495,0.49,0.50};

  double per_err_200[6] = {0.238,0.155,0.10,0.086,0.088,0.064};
  double per_err_200_North[6] = {0.286,0.127,0.13,0.132,0.103,0.074};

  double per_err_1000[10] = {0,0,0,0,0.141,0.458,0.455,0,0,0};
  double per_err_1000_North[10] = {10,10,10,10,0.145,0.273,0.435,10,10,10};

  double per_err_100[10] = {0,0,0,0.155,0.148,0.145,0.125,0.123,0.29,0};
  double per_err_300[21] = {0.163,0.089,0.073,0.061,0.11,0.123,0.106,0.11,0.111,0.132,0.12,0.157,0.18,0.25,0.274,0.287,0.432,0.479,0.522,0.58,0.47};
  
  double per_err_2000[5] = {0.1,0.1,0.1,0.322,0.1};
  double per_err_2000_North[5] = {0.1,0.1,0.1,0.246,0.1};

  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  if((pt_binwidth == 1))
	    per_err[pt] = per_err_200[pt];
	  if(pt_binwidth == 3)
	    per_err[pt] = per_err_1000[pt];
	  if(pt_binwidth == 3 && arm_low == 1)
	    per_err[pt] = per_err_1000_North[pt];
	  if((pt_binwidth == 0))
	    per_err[pt] = per_err_100[pt];
	  if((pt_binwidth == 4))
	    per_err[pt] = per_err_2000[pt];
	  if(pt_binwidth == 4 && arm_low == 1)
	    per_err[pt] = per_err_2000_North[pt];
	  if(( arm_low == 1) && (pt_binwidth == 2) )
	    per_err[pt] = per_err_North[pt];
	  if(( arm_low == 1) && (pt_binwidth == 3) )
	    per_err[pt] = per_err_1000_North[pt];
	  if((arm_low == 1) && (pt_binwidth == 1))
	    per_err[pt] = per_err_200_North[pt];


	  for(int i = 0; i < bins; i++)
  	    {
	      cc[arm][pt][i] =  (bg[0][arm][pt]->GetBinContent(i+1));
	      bb_UL[arm][pt][i] = (bg[1][arm][pt]->GetBinContent(i+1));
	      dy[arm][pt][i] = (bg[2][arm][pt]->GetBinContent(i+1));
	      bb_LS[arm][pt][i] = (bg[3][arm][pt]->GetBinContent(i+1));
	      corr_had_LS[arm][pt][i] = (bg[4][arm][pt]->GetBinContent(i+1));
	      corr_had_UL[arm][pt][i] = (bg[5][arm][pt]->GetBinContent(i+1));
	      corr_mass[arm][i] = (t_check[arm]->GetBinContent(i+1));

	      //YUE HANG SHAPE
	      // corr_bg[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i];  // sanghoon said to change WWND

	      // corr bg SANGHOON:
	      corr_bg[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i] +corr_had_UL[arm][pt][i] - corr_had_LS[arm][pt][i];
	      corr_bg_err[arm][pt][i] = per_err[pt]*corr_bg[arm][pt][i]; 
	      corr_total[arm][i]+= corr_bg[arm][pt][i];
	      corr_total_err[arm][i] += corr_bg_err[arm][pt][i]; 
	      y_errors[arm][i] = 0.00000000035*corr_total[arm][i];
	    }
	}
    }

  char uniqueh[38][500];
  int i = 0;
  
  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      i++;
      sprintf(uniqueh[pt],"h_sum_%d",i); 
    }

  TH1D *h1 = new TH1D("h1", "Distribution", 100, 0, 10);
  TH1D *h_pt_int[2];
  h_pt_int[0] = new TH1D("h_pt_int", "pT Integrated South", 100, 0, 10);
  h_pt_int[1] = new TH1D("h_pt_int_2", "pT Integrated North", 100, 0, 10);

  double cc_err[2][38][100];
  double bb_UL_err[2][38][100];
  double dy_err[2][38][100];
  double bb_LS_err[2][38][100];
  double corr_had_LS_err[2][38][100];
  double corr_had_UL_err[2][38][100];

  double quad[2][38][100];
  double quad_norm[2][38][100];
  double temp[2][100];
 
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {
	      // YUE HANG CORR BG SHAPE
	      //  h1->SetBinContent(i+1,(bg[0][arm][pt]->GetBinContent(i+1) + bg[1][arm][pt]->GetBinContent(i+1) + bg[2][arm][pt]->GetBinContent(i+1) - bg[3][arm][pt]->GetBinContent(i+1)));   // sanghoon changed WWND

	      // SANGHOON CORR BG SHAPE
	      h1->SetBinContent(i+1,(bg[0][arm][pt]->GetBinContent(i+1) + bg[1][arm][pt]->GetBinContent(i+1) + bg[2][arm][pt]->GetBinContent(i+1) - bg[3][arm][pt]->GetBinContent(i+1) + bg[5][arm][pt]->GetBinContent(i+1) - bg[4][arm][pt]->GetBinContent(i+1)));

	      cc_err[arm][pt][i] = (bg[0][arm][pt]->GetBinContent(i+1))*0.1;  //  ~10% 
	      bb_UL_err[arm][pt][i] = (bg[1][arm][pt]->GetBinContent(i+1))*0.1;  //  ~10% 
	      dy_err[arm][pt][i] = (bg[2][arm][pt]->GetBinContent(i+1))*0.1;  //  ~10% 
	      bb_LS_err[arm][pt][i] = (bg[3][arm][pt]->GetBinContent(i+1))*0.15;  //  ~15% 
	      corr_had_LS_err[arm][pt][i] = (bg[4][arm][pt]->GetBinContent(i+1))*0.1;  //  ~10%  ??
	      corr_had_UL_err[arm][pt][i] = (bg[5][arm][pt]->GetBinContent(i+1))*0.1;  //  ~10%  ??

	      quad[arm][pt][i] = sqrt( pow(cc_err[arm][pt][i],2) +  pow(bb_UL_err[arm][pt][i],2) +  pow(dy_err[arm][pt][i],2) +  pow(bb_LS_err[arm][pt][i],2) +  pow(corr_had_LS_err[arm][pt][i],2) +  pow(corr_had_UL_err[arm][pt][i],2));
	       
	      if(h1->GetBinContent(i+1) > 0)
		{
		  //h1->SetBinError(i+1,quad[arm][pt][i]);
		  h1->SetBinError(i+1,(h1->GetBinContent(i+1))*per_err[pt]); 
		}
	     
	      //  temp[arm][i] = 0.02*corr_total[arm][i]; 
	      x_errors[arm][i] = 0.0;

	      h_sum[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_sum[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));
	      
	    }
	}
    }
   
  ///////////// NORMALIZE HISTOGRAMS ////////////////// 

  double area[2][38];
  Double_t norm[2][38];
  double norm_area[2][38];
  int num_pt_bins = pt_slices[pt_binwidth];
  
  if(normalization)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      area[arm][pt] = h_sum[arm][pt]->Integral( h_sum[0][0]->FindBin(2.001),  h_sum[0][0]->FindBin(4.999) );
	     	      
	      if(area[arm][pt] != 0)
		{
		  norm[arm][pt] = 1/(area[arm][pt]);
		  h_sum[arm][pt]->Scale(norm[arm][pt]);
		}
	    		
	      norm_area[arm][pt] = h_sum[arm][pt]->Integral( h_sum[0][0]->FindBin(2.001),  h_sum[0][0]->FindBin(4.999) );
	    }
	 
	}
    }

  //Normalize mass ratio components 
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  	{
  	  for(int i = 0; i < bins; i++)
  	    {
	      h_pt_int[arm]->Fill(corr_total[arm][i]);
	      corr_bg[arm][pt][i] *= norm[arm][pt];
	      corr_total_norm[arm][i] += corr_bg[arm][pt][i] ;
	    }
	}
    }

  // Calculate mass ratio after corr_total is filled and normalized
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
  	{
	  for(int i = 0; i < bins; i++)
  	    {
	      
	      if(corr_total_norm[arm][i] != 0)
		mass_ratio[arm][pt][i] = corr_bg[arm][pt][i]/(corr_total_norm[arm][i]/num_pt_bins);
	      else
		mass_ratio[arm][pt][i] = 0;

	    }
	}
    }

  ///////////// CUSTOMIZED PLOTS ///////////////////////
  // initialize

  int pt_lo = 0;
  int pt_hi = 12;

  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      //cout << "a: " << a_array[pt] << ", b: " << b_array[pt] << endl;
    }

  ////////////////////////////////////////////////
  // Histogram
  TCanvas *c6 = new TCanvas("c6","Corr BG Fits",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  // gPad->SetLogy();
  gPad->SetGrid();

  if(all_fits)
    {
      h_sum[arm_low][0]->SetTitle("pT Int. Corr BG (cc(UL)+dy(UL)+bb(UL)+corr_had(UL)-bb(LS)-corr_had(LS)), South");
      h_sum[arm_low][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm_low][0]->GetXaxis()->SetLabelSize(0.04);
      h_sum[arm_low][0]->GetXaxis()->SetTitleSize(0.04);
      h_sum[arm_low][0]->GetXaxis()->SetTitleOffset(0.9);
      h_sum[arm_low][0]->GetYaxis()->SetLabelSize(0.04);
      h_sum[arm_low][0]->GetYaxis()->SetTitleSize(0.1); 
      h_sum[arm_low][0]->GetYaxis()->SetTitleOffset(0.52);
    }
  else
    {
      h_sum[0][pt_l]->SetTitle("pT Int. Corr BG (cc(UL)+dy(UL)+bb(UL)+corr_had(UL)-bb(LS)-corr_had(LS)), South");
      h_sum[1][pt_l]->SetTitle("pT Int. Corr BG (cc(UL)+dy(UL)+bb(UL)+corr_had(UL)-bb(LS)-corr_had(LS)), North");
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm_low][pt_l]->GetXaxis()->SetLabelSize(0.04);
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitleSize(0.04);
      h_sum[arm_low][pt_l]->GetXaxis()->SetTitleOffset(0.9);
      h_sum[arm_low][pt_l]->GetYaxis()->SetLabelSize(0.04);
      h_sum[arm_low][pt_l]->GetYaxis()->SetTitleSize(0.1); 
      h_sum[arm_low][pt_l]->GetYaxis()->SetTitleOffset(0.52);
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
  
  if(normalization)
    h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-3),pow(10,0));    
  else
    {
      if(all_fits)
	{
	  h_sum[arm_low][0]->GetYaxis()->SetRangeUser(pow(10,-1),pow(10,4));
	  h_sum[arm_low][0]->SetAxisRange(0,5);
	}
      else
	{
	  if(pt_binwidth == 2)
	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),1*pow(10,3));
	  //h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),180);
	  // h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,0),800);
	  if(pt_binwidth == 1)
	    h_sum[arm_low][pt_l]->GetYaxis()->SetRangeUser(pow(10,-1),4*pow(10,2));
	  h_sum[arm_low][pt_l]->SetAxisRange(0,5);
	}
    }
  /////////////////////////////////////////////////////////////////////////////////
  /// Make legend for histogram

  char unique_l[800];
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
	  sprintf(unique_l,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_h[arm]->AddEntry(h_sum[arm][pt],unique_l, "p");
	}
    }

  leg_h[arm_low]->Draw();

  // pT Integrated corrletaed background plot
  TGraphErrors *gr_corrbg_total = new TGraphErrors(bins,x_array,corr_total[arm_low],x_errors[arm_low],corr_total_err[arm_low]);
  
  //////////////////////////// Fit histograms //////////////////////////
  
  char name1[800];
  char name2[800];

  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;
  int bins_fail[bins];

  TF1 *corr_bg_fit;
  TF1 *corr_total_fit;
  TF1 *int_fit_pt;
  double f5,f6,f10,f11,f12;
 

  // For 500 MeV bins
  double par_0[14] = {1,1,1,1,1,1,1,1,1,1,1,0.1,1,1};  // original
  double par_1[14] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  // double par_2[14] = {10,10,100,100,100,100,100,10,10,10,1,1,1,1};
  double par_2[14] = {10,10,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double par_3[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10}; 
  double par_4[14] = {10,10,10,10,10,10,10,10,10,10,10,10,10,10};  // original

  double par_2_refit[14] = {10,10,100,100,100,100,100,10,10,100,1,1,1,1};   //10
 
  double limit_a_refit[14] = {1.4,1.4,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4};   //(10,0.01,10,100,0.01)
  //double limit_a_refit[14] = {1.4,1.4,2,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,1.4,1.4};   //(10,0.01,10,100,0.01)

  double limit_a[14] = {1.75,1.75,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.75,0.75,0.75,0.75};
  double limit_b[14] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5};

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
  double par_0_1000[10] = {0,0,0,0,1,1,1,0,0,0};
  double par_1_1000[10] = {0,0,0,0,0.01,0.01,0.01,0,0,0};
  double par_2_1000[10] = {0.01,0.01,0.01,0.01,1,1,1,0,0,0};
  double par_3_1000[10] = {0,0,0,0,10,10,10,0,0,0};
  double par_4_1000[10] = {0,0,0,0,10,10,10,0,0,0};
  
  double limit_a_1000[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 
  double limit_b_1000[10] = {0,0,0,0,5,5,5,0,0,0};

  double limit_a_1000_North[10] = {0,0,0,0,0.5,0.5,0.5,0,0,0}; 

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

 double fit_array[14][2][100] = {0};
 TFitResultPtr r;

 //////////////////////////////////////////////////////////////////////////  All code above this point is for both North and South 
  

 for(int arm = arm_low; arm < arm_high; arm++)  
   // for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = for_begin; pt < for_end; pt++) // ALL FITS TO DISPLAY ON SAME PLOT
	//	for(int pt = 1; pt < for_end; pt++)
	{
	  for(int bin = 0; bin < 100; bin++)
	    {
	      if(refit && pt_binwidth == 2)
		{
		  par_2[pt] = par_2_refit[pt];
		  limit_a[pt] = limit_a_refit[pt];
		}
	      
	      if(pt_binwidth == 1)
		{
		  par_0[pt] = par_0_200[pt];
		  par_1[pt] = par_1_200[pt];
		  par_2[pt] = par_2_200[pt];
		  par_3[pt] = par_3_200[pt];
		  par_4[pt] = par_4_200[pt];
		  if(south_arm)
		    limit_a[pt] = limit_a_200[pt];
		  else
		    limit_a[pt] = limit_a_200_North[pt];
		  limit_b[pt] = limit_b_200[pt];
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
		  
		  if(south_arm == false)
		    limit_a[pt] = limit_a_1000_North[pt];
		  
		}
	      if(pt_binwidth == 0)
		{
		  par_0[pt] = par_0_100[pt];
		  par_1[pt] = par_1_100[pt];
		  par_2[pt] = par_2_100[pt];
		  par_3[pt] = par_3_100[pt];
		  par_4[pt] = par_4_100[pt];
		  
		  limit_a[pt] = limit_a_100[pt];
		  limit_b[pt] = limit_b_100[pt];
		}
	      if(pt_binwidth == 0) //4
		{
		  par_0[pt] = par_0_300[pt];
		  par_1[pt] = par_1_300[pt];
		  par_2[pt] = par_2_300[pt];
		  par_3[pt] = par_3_300[pt];
		  par_4[pt] = par_4_300[pt];
		  
		  limit_a[pt] = limit_a_300[pt];
		  limit_b[pt] = limit_b_300[pt];
		}
	      if(pt_binwidth == 4)
		{
		  par_0[pt] = par_0_2000[pt];
		  par_1[pt] = par_1_2000[pt];
		  par_2[pt] = par_2_2000[pt];
		  par_3[pt] = par_3_2000[pt];
		  par_4[pt] = par_4_2000[pt];
		  
		  limit_a[pt] = limit_a_2000[pt];
		  limit_b[pt] = limit_b_2000[pt];
		}
	      
	      corr_bg_fit = new TF1("corr_bg_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",limit_a[pt],limit_b[pt]);
	      
	      corr_total_fit = new TF1("corr_total_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",0,5);
	      
	      //cout <<" TF1 fit is from " << limit_a[pt] << "  to " << limit_b[pt] << endl;
	      
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
	      
	      
	      x = pt_center[pt_binwidth][pt];
	      //cout << "x = pT center: " << x << endl;
	      
	      if(refit == true)
		{
		  
		  f10 = par_2[pt];
		  
		  // cout << "f10: " << f10 << endl;
		  
		  
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
		  
		  // cout << "a0: " << par_a[0] << ", a1: " << par_a[1] << ", a2: " << par_a[2] << ", a3: " << par_a[3] << ", a4: " << par_a[4] << ", a5: " << par_a[5] << ", a6: " << par_a[6] << ", a7: " << par_a[7] << ", a8: " << par_a[8] << ", a9: " << par_a[9] << endl;
		  
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

	      if(pt_int_fit == false)
		{
		  //  if(pt < 2 && south_arm)
		  //  continue;
		  r = h_sum[arm_low][pt]->Fit(corr_bg_fit,"S R");
		  fit_array[pt][arm][bin] = corr_bg_fit->Eval(x_array[bin]);
		  // cout << "Chi Square/NDF: " << corr_bg_fit->GetChisquare() << "/" << corr_bg_fit->GetNDF() << "" << endl;
		}
	      else
		{
		  r = gr_corrbg_total->Fit(corr_total_fit,"S R");
		  // cout << "Chi Square/NDF: " << corr_total_fit->GetChisquare() << "/" << corr_total_fit->GetNDF() << "" << endl;
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

	      //  cout << "status: " << fitStatus << " for pt bin " << pt+1 << endl;
	      /*
		if(fitStatus == 0) 
		{
		if(write)	     
		{
		if(arm == 1)
		{
		if(pt_binwidth == 2)
		sprintf(name1,"yuehang_pt_int/f(x)_fit_value_N_%i.dat",pt+1); 
		if(pt_binwidth == 1)
		sprintf(name1,"yuehang_pt_int/bestfit_parameters_low_pt_N_%i.dat",pt+1); 
		if(pt_binwidth == 3)
		sprintf(name1,"yuehang_pt_int/bestfit_parameters_high_pt_N_%i.dat",pt+1); 
		if(pt_binwidth == 0)
		sprintf(name1,"yuehang_pt_int/bestfit_parameters_extra_low_pt_N_%i.dat",pt+1);
	   
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
		sprintf(name2,"yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);
		if(pt_binwidth == 0)
		sprintf(name2,"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_%i.dat",pt+1);
		if(pt_binwidth == 4)
		sprintf(name2,"yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_%i.dat",pt+1);

 	      
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
	  
	      */
	    } // for loop bin
	} // for loop pt
      
    } // for loop arm
 
 Char_t message[200];
 
 if(pt_int_fit == false)
   sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_bg_fit->GetChisquare(),corr_bg_fit->GetNDF());
 else
   sprintf(message,"#chi^{2}/NDF = %.1f / %.d",corr_total_fit->GetChisquare(),corr_total_fit->GetNDF());
 
 TPaveText *mytext = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
 mytext->SetTextSize(0.035);
 mytext->SetFillColor(0);
 mytext->SetTextAlign(12);
 mytext->AddText(message);
 mytext->Draw();
 
 
 
 // CREATE MINBIAS SUM
 //pt_high = 0; 
 // pt_high = 5;
 //pt_high = 10; 
 //pt_high = 7;

 // write fit value at each pt bin
 
 // std::fstream fit_array_file("yuehang_pt_int/fit_array_200mev_S.C", std::ofstream::out);
  //  fit_array_file << "double fit_array_200mev_S[4][100] = { " ;
 
  
  //std::fstream fit_array_file("yuehang_pt_int/fit_array_200mev_N.C", std::ofstream::out);
  //fit_array_file << "double fit_array_200mev_N[4][100] = { " ;
 
 //std::fstream fit_array_file("yuehang_pt_int/fit_array_500mev_S.C", std::ofstream::out);
 // fit_array_file << "double fit_array_500mev_S[10][100] = { " ;

 //std::fstream fit_array_file("yuehang_pt_int/fit_array_500mev_N.C", std::ofstream::out);
 // fit_array_file << "double fit_array_500mev_N[10][100] = { " ;

 //std::fstream fit_array_file("yuehang_pt_int/fit_array_1gev_S.C", std::ofstream::out);
 // fit_array_file << "double fit_array_1gev_S[7][100] = { " ;
 
 //std::fstream fit_array_file("yuehang_pt_int/fit_array_1gev_N.C", std::ofstream::out);
 // fit_array_file << "double fit_array_1gev_N[7][100] = { " ;
 
 /*
 if(write)
   {
     for(int pt = 0; pt < pt_high; pt++)
       {
	 for(int bin = 0; bin < 100; bin++)
	   {
	     if(pt_binwidth == 1)
	       {
		 if((bin < 99) && (pt != 2))
		   fit_array_file <<  fit_array[pt][arm_low][bin]  <<  ", ";
		 if((bin == 99) && (pt == 4))
		   fit_array_file <<  fit_array[pt][arm_low][bin] << "};" << endl;

		 // try this one to not change numbering of array when skipping 400 - 600
		 // if(bin < 99)
		 //   {
		 //     if(pt < 2)
		 //       fit_array_file <<  fit_array[pt][arm_low][bin]  <<  ", ";
		 //     if(pt == 3)
		 //       fit_array_file <<  fit_array[pt-1][arm_low][bin]  <<  ", ";
		 //   }
		 //     if((bin == 99) && (pt == 4))
		 //       fit_array_file <<  fit_array[pt-1][arm_low][bin] << "};" << endl;
		   
		 
	       }
	     if(pt_binwidth == 2)
	       {
		 if(bin < 99)
		   fit_array_file <<  fit_array[pt][arm_low][bin]  <<  ", ";
		 if((bin == 99) && (pt == 9))
		   fit_array_file <<  fit_array[pt][arm_low][bin] << "};" << endl;
	       }
	     if(pt_binwidth == 3)
	       {
		 if(bin < 99)
		   fit_array_file <<  fit_array[pt][arm_low][bin]  <<  ", ";
		 if((bin == 99) && (pt == 6))
		   fit_array_file <<  fit_array[pt][arm_low][bin] << "};" << endl;
	       }
	   }
	 fit_array_file << endl;
       }
   } // write
 
 */


 for( int pt = 0; pt < 4; pt++)
   {
     for(int arm = 1; arm < 2; arm++)
       {
	 for(int i = 0; i < 100; i ++)
	   {
	     // cout << "pt = " << pt << endl;
	     //cout << fit_array_200mev_N[pt][i] << ", for bin: " << i << endl;
	   }
       }
   }
	    
 	  
	  
	    
		    
 
 // NOW ACCESS FIT ARRAY FILES AFTER ALL HAVE BEEN WRITTEN
 //double fit_array_200mev_N[4][100]
 //double fit_array_500mev_N[10][100]
 //double fit_array_1gev_N[7][100] 
 
 //area under each curve for all 14 points is in array area_corrbg[arm][14];
 double num[2][100];
 fit_array[14][2][100] = {0};
 double sum_area[2];
  double Run15pp_minbias[2][100];

 for(int arm = 0; arm < 2; arm++)
   {
     for(int k = 0; k < 14; k++)
       {
	 sum_area[arm] += area_corrbg[arm][k];
       }
   }
 
 for(int k = 0; k < 14; k++)
   {
     for(int arm = 0; arm < 2; arm++)
       {
	 for(int bin = 0; bin < 100; bin++)
	   {
	     if(k < 4)
	       {
		 fit_array[k][0][bin] = fit_array_200mev_S[k][bin];
		 fit_array[k][1][bin] = fit_array_200mev_N[k][bin];
	       }
	     if((k < 12) && (k > 3))
	       {
		 fit_array[k][0][bin] = fit_array_500mev_S[k][bin];
		 fit_array[k][1][bin] = fit_array_500mev_N[k][bin];
	       }
	     if(k > 11)
	       {
		 fit_array[k][0][bin] = fit_array_1gev_S[k][bin];
		 fit_array[k][1][bin] = fit_array_1gev_N[k][bin];
	       }
	     
	     cout << fit_array[k][arm][bin] << endl;
	     
	     //  cout << "Sum of area: " << sum_area[arm] << endl;
	     // cout << sum_area[arm] << endl;
	     num[arm][bin] += fit_array[k][arm][bin] * area_corrbg[arm][k];
	     //  cout << "numerator: " << num[arm][bin] << endl;
	     Run15pp_minbias[arm][bin] =   num[arm][bin] /sum_area[arm];
	     // cout << "Result: " << Run15pp_minbias[arm][bin] << endl;
	   }
       }
   }
 
 
 //PLOT MINBIAS AND FIT
 
   TCanvas *c20 = new TCanvas("c20","Corr BG pT Integrated",200,10,700,500);
   gPad->SetGrid();
   //gPad->SetLogy();
   
   // pT Integrated corrletaed background plot
   //TGraphErrors *corrbg_int = new TGraphErrors(bins,x_array,Run15pp_minbias[arm_low],x_errors[arm_low],corr_total_err[arm_low]);
   TGraphErrors *corrbg_int = new TGraphErrors(bins,x_array,Run15pp_minbias[arm_low],x_errors[arm_low],y_errors[arm_low]);
   
   if(south_arm)
     corrbg_int->SetTitle("Corr BG pT Integrated, South");  
   else
     corrbg_int->SetTitle("Corr BG pT Integrated, North");  
   corrbg_int->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
   corrbg_int->GetXaxis()->SetLabelSize(0.04);
   corrbg_int->GetXaxis()->SetTitleSize(0.04);
   corrbg_int->GetXaxis()->SetTitleOffset(0.9);
   corrbg_int->GetXaxis()->SetLimits(0,5);
   corrbg_int->GetYaxis()->SetLabelSize(0.04);
   corrbg_int->GetYaxis()->SetTitleSize(0.1); 
   corrbg_int->GetYaxis()->SetTitleOffset(0.52);
   // corrbg_int->SetMaximum(pow(10,5));
   // corrbg_int->SetMinimum(pow(10,2));
   
   /// assign different colors and draw each bg
   corrbg_int->SetMarkerColor(kViolet);  
   corrbg_int->SetMarkerSize(0.75);
   corrbg_int->SetMarkerStyle(20);
   corrbg_int->Draw("AP");
   

   
   ////////
   double f0,f1,f2,f3,f4;

   if(south_arm == false)
     {
       // f0 = 10;
       // f1 = 0.01;
       // f2 = 0.01;
       // f3 = 100;
       // f4 = 0.01;
       // used Feb 2019
       // f0 = 1;
       f0 = 10;
       f1 = 0.1;
       f2 = 0.01;
       f3 = 1;
       f4 = 1;
   
     }
   else
     {
       f0 = 10;
       f1 = 0.01;
       f2 = 0.01;
       f3 = 100;
       f4 = 0.01;
     }

   int_fit_pt = new TF1("int_fit_pt"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",1,5);

   int_fit_pt->SetParameter(0, f0);
   int_fit_pt->SetParameter(1, f1);
   int_fit_pt->SetParameter(2, f2);
   int_fit_pt->SetParameter(3, f3);
   int_fit_pt->SetParameter(4, f4);
   
   int_fit_pt->SetLineColor(kRed);
   int_fit_pt->SetLineStyle(5);
   int_fit_pt->SetLineWidth(2);
 
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
 
   r = corrbg_int->Fit(int_fit_pt,"S R");
   
   // r->Parameter(0);
   // r->Parameter(1);
   // r->Parameter(2);
   // p3 = r->Parameter(3);
   // p4 = r->Parameter(4);
   
   //  mytext->Draw();
   
 
    
} // void end macro
