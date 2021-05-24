// original binwidths were 1 GeV
// change pt slices plotted on graph at line 388

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

void yuehang_pt_macro()
{

  #include "yuehang_fit_coeff_500/fit_coeff_S_array.C"

  bool normalization = false;
  bool write = false;
   
  //double scale[5] = {25*pow(10,11),25*pow(10,10),25*pow(10,9),25*pow(10,9),25*pow(10,10)};  // rebinned hisotgrams already scaled.  This is just for pT integrated


  //double scale[2] = {2.415*pow(10,11),2.3*pow(10,11)}; // for mixed event BG
 

  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; // for LS background 
  //double scale[2] = {1,1};
  
  
  int pt_binwidth;
  cout << "What pt binning are you running?  Enter '0' for Yue Hang binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 300 MeV binwidth" << endl;
  cin >> pt_binwidth;

  int pt_slices[5] = {4,38,20,10,21};

  int num_bins_on_plots[5] = {4,10,6,5,21};

  double pt_width[5][38] = {1.0,1.0,1.0,7.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
			    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  double pt_center[5][38] = {0.5,1.5,2.5,6.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.25,5.55,5.85,6.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
 
  std::string filename[5][2] = {"yuehang_Run15pp_S_rebinned_1000_mev.root", 
				"yuehang_Run15pp_N_rebinned_1000_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_200_mev.root", 
				 "yuehang_Run15pp_N_rebinned_200_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_500_mev.root", 
				 "yuehang_Run15pp_N_rebinned_500_mev.root",
				 
				 "yuehang_Run15pp_S_rebinned_1_gev.root", 
				"yuehang_Run15pp_N_rebinned_1_gev.root",

 				 "yuehang_Run15pp_S_rebinned_300_mev.root", 
				 "yuehang_Run15pp_N_rebinned_300_mev.root"};
 
 
  //////////////////////////////////////////////////
  //everything below this line is created as a result of making an initial selection about which data you want to use

  TH1D *bg[6][2][38];
  TH1D *t_check[2]; 
  TH2F *yuehang_2d = 0;

  std::string obj_filename[6][2][38]; 
  
  char unique1[800];
  char unique2[800];
  char unique3[800];
  char unique4[800];
  char unique5[800];
  char unique6[800];

  TAxis *xaxis;
  TAxis *yaxis;	  
 
// corr. hadron binning is 150 along x-axis, but cc,bb,dy are 100 bins along x-axis
  int bins = 100;

  double pt_low;
  double pt_high;

  double bin_low;
  double bin_high;

  double delta_pt = 0.001;
  double a_array[pt_slices[pt_binwidth]];
  double b_array[pt_slices[pt_binwidth]];
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < pt_slices[pt_binwidth]; k++)
	{ 
	  pt_low = pt_center[pt_binwidth][k] - pt_width[pt_binwidth][k]/2 + delta_pt; 
	  pt_high = pt_center[pt_binwidth][k] + pt_width[pt_binwidth][k]/2 - delta_pt;
	  
	  double a = (pt_low - delta_pt)*1000; // for unique names in projected histograms, use meV
	  double b = (pt_high + delta_pt)*1000;

	  // cout << "a: " << a << ", b: " << b << endl;
	   
	  sprintf(unique1,"cc_%.0f_%.0f",a,b); 
	  sprintf(unique2,"bb_UL_%.0f_%.0f",a,b); 
	  sprintf(unique3,"dy_%.0f_%.0f",a,b);  
	  sprintf(unique4,"bb_LS_%.0f_%.0f",a,b);
	  sprintf(unique5,"corr_had_LS_%.0f_%.0f",a,b);
	  sprintf(unique6,"corr_had_UL_%.0f_%.0f",a,b);

	  
	  a_array[k] = a/1000;
	  b_array[k] = b/1000;
	  // cout << "a: " << a_array[k] << ", b: " << b_array[k] << endl;

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
	    }
	}
    }
  
   TFile *file_data;
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
     
      bin_35gev = xaxis->FindBin(3.501);
      bin_4gev = xaxis->FindBin(3.999);
      bin_5gev = xaxis->FindBin(4.999);

      cout << "bin : " << i+1 << ", and bin 3.5: " << bin_35gev << endl;

      if(i+1 < bin_35gev)
	{
	  x_array_mass_ratio[bin_counter] = xaxis->GetBinCenter(i+1);
	  bin_counter++;
	}
      if(i+1 == bin_35gev)
	{
	  bin_counter++;
	  x_array_mass_ratio[bin_counter] = 3.75;
	}
      if((i+1 > bin_35gev) && (i+1 < bin_4gev))
	continue;
      if(i+1 == bin_5gev)
	{
	  bin_counter++;
	  x_array_mass_ratio[bin_counter] = 4.5;
	}
      if((i+1 > bin_4gev) && (i+1 < bin_5gev))
	continue;
      if(i+1 > bin_5gev)
	{
	  bin_counter++;
	  x_array_mass_ratio[bin_counter] = xaxis->GetBinCenter(i+1);
	}

      //cout << "mass array: " << x_array[i] << " for bin " << i+1 << endl;
      cout << "mass ratio array: " << x_array_mass_ratio[i] << " for bin " << i+1 << endl;
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
 
  double cc[2][pt_slices[pt_binwidth]][bins]; 
  double bb_UL[2][pt_slices[pt_binwidth]][bins];
  double dy[2][pt_slices[pt_binwidth]][bins];
  double bb_LS[2][pt_slices[pt_binwidth]][bins]; 
  double corr_had_LS[2][pt_slices[pt_binwidth]][bins];  // NEW
  double corr_had_UL[2][pt_slices[pt_binwidth]][bins];  // NEW
  double corr_bg[2][pt_slices[pt_binwidth]][bins];
  double corr_total[2][100] = {0.0,0.0};
  double corr_total_err[2][100] = {0.0,0.0};
  double corr_total_norm[2][100] = {0.0,0.0};
  double corr_bg_err[2][38][100] = {0.0,0.0};
  double pt_int[2][38];
  double corr_mass[2][bins];
  double errors[2][pt_slices[pt_binwidth]][bins];
  double x_errors[2][bins];
  double counter = 0.0;
  double mass_ratio[2][pt_slices[pt_binwidth]][bins];
	  
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
    
  int pt_l = 0;
  int pt_h = 1;
  //double lim_a = 2;
  //double lim_b = 5;
  double per_err[21] = {0.05,0.065,0.065,0.13,0.115,0.13,0.16,0.175,0.2,0.4,0.3,0.17,0.17,0.17,0,0,0,0,0,0,0};

  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  if((pt  < pt_l) || (pt > pt_h))
	    continue;
  	  for(int i = 0; i < bins; i++)
  	    {
	      cc[arm][pt][i] =  (bg[0][arm][pt]->GetBinContent(i+1));
	      //cout << "cc:" << cc[arm][pt][i] << ", bg[0] : " << bg[0][arm][pt]->GetBinContent(i+1) <<  endl;
	      bb_UL[arm][pt][i] = (bg[1][arm][pt]->GetBinContent(i+1));
	      dy[arm][pt][i] = (bg[2][arm][pt]->GetBinContent(i+1));
	      bb_LS[arm][pt][i] = (bg[3][arm][pt]->GetBinContent(i+1));
	      corr_had_LS[arm][pt][i] = (bg[4][arm][pt]->GetBinContent(i+1));
	      corr_had_UL[arm][pt][i] = (bg[5][arm][pt]->GetBinContent(i+1));
	      corr_mass[arm][i] = (t_check[arm]->GetBinContent(i+1));
	      //corr_bg[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i];  // sanghoon said to change WWND
	      corr_bg[arm][pt][i] = cc[arm][pt][i] + bb_UL[arm][pt][i] + dy[arm][pt][i] - bb_LS[arm][pt][i] +corr_had_UL[arm][pt][i] - corr_had_LS[arm][pt][i];
	      corr_bg_err[arm][pt][i] = per_err[pt_l]*corr_bg[arm][pt][i]; 
	      corr_total[arm][i]+= corr_bg[arm][pt][i];
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
 TH1D *h_pt_int[2];
 h_pt_int[0] = new TH1D("h_pt_int", "pT Integrated South", 100, 0, 10);
 h_pt_int[1] = new TH1D("h_pt_int_2", "pT Integrated North", 100, 0, 10);

 double cc_err[2][38][bins];
 double bb_UL_err[2][38][bins];
 double dy_err[2][38][bins];
 double bb_LS_err[2][38][bins];
 double corr_had_LS_err[2][38][bins];
 double corr_had_UL_err[2][38][bins];

 double quad[2][38][bins];
 double quad_norm[2][38][bins];
 double temp[2][bins];
 
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {
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
		  h1->SetBinError(i+1,quad[arm][pt][i]);
		}
	     
	      temp[arm][i] = 0.02*corr_total[arm][i]; // Tony said to overwrite errors to 2%
	      corr_total_err[arm][i] += temp[arm][i]; 
	      x_errors[arm][i] = 0.0;

	      //cout << "bin: " << i+1 << ", corr total err: " << corr_total_err[arm][i] << endl;

	      h_sum[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_sum[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));
	      h_sum[arm][pt]->SetBinError(i+1,h1->GetBinError(i+1));
	      //h_sum[arm][pt]->SetBinError(i+1,h1->GetBinContent(i+1)*0.02); // Tony said to overwrite errors to 2%

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
	      //cout << "before norm errors: " << quad[arm][pt][i] << endl;
	      // cout << "after norm errors: " << h_sum[arm][pt]->GetBinError(i+1) << endl;
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
	      //h_pt_int[arm]->Add(h_sum[arm][pt]);
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
  int pt_hi = 0;

  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      cout << "a: " << a_array[pt] << ", b: " << b_array[pt] << endl;
      if(a_array[pt] == 0)
	pt_lo = pt;
      if(b_array[pt] == 2) 
	pt_hi =  pt;
    }

   cout << "pt lo: " << pt_lo << ", and pt hi: " << pt_hi << endl;

   ///////////////////////  MAKE HISTOGRAM PLOTS ///////////////////////////

  char unique_l[800];
  TLegend *leg_h[2];
 
  for(int arm = 0; arm < 2; arm++)
    {
      leg_h[arm] = new TLegend(0.51, 0.6, 0.8, 0.9);  
      leg_h[arm]->SetFillColor(0); 
      leg_h[arm]->SetTextSize(0.035);
	  	 
      for(int pt = pt_lo; pt < pt_hi+1; pt++)
	{
	  cout << "a: " << a_array[pt] << ", b: " << b_array[pt] << endl;
	  sprintf(unique_l,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_h[arm]->AddEntry(h_sum[arm][pt],unique_l, "l");
	}
    }

  TCanvas *c1 = new TCanvas("c1","Histogram South",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gPad->SetLogy();
  gPad->SetGrid();
  
  for(int arm = 0; arm < 2; arm++)
    {
      h_sum[arm][pt_lo]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_sum[arm][pt_lo]->GetXaxis()->SetLabelSize(0.06);
      h_sum[arm][pt_lo]->GetXaxis()->SetTitleSize(0.05);
      h_sum[arm][pt_lo]->GetXaxis()->SetTitleOffset(1.00);
      h_sum[arm][pt_lo]->SetAxisRange(2,5); 
    }

  h_sum[0][pt_lo]->SetTitle("Yue Hang Corr BG (cc+bb_UL+dy-bb_LS+corr_had_UL-corr_had_LS), South");
  
  double Counter = 1;

  for(int pt = pt_lo; pt < pt_hi+1; pt++)  
    {  
      h_sum[0][pt]->SetLineColor(Counter); 
      if(Counter == 10)
	h_sum[0][pt]->SetLineColor(40); 
      h_sum[0][pt]->SetLineWidth(3.0);
      if(pt == pt_lo)
	h_sum[0][pt]->Draw();
      else
	h_sum[0][pt]->Draw("SAME");
      h_sum[0][pt_lo]->SetAxisRange(2,5); 
      Counter++;
    }
   
  leg_h[0]->Draw();

  /////////////////////////////////// Make Mass Ratio TGraph /////////////////////////////////

  TGraph *gr[2][38];
     
  TCanvas *c3 = new TCanvas("c3","TGraph Mass Ratio",200,10,700,500);
  gPad->SetGrid();
  // gPad->SetLogy();
 
  Counter = 1;

  for(int arm = 0; arm < 2; arm++) 
    {
      for(int pt = pt_lo; pt < pt_hi+1; pt++)
	{
	  gr[arm][pt] = new TGraph(bins,x_array,mass_ratio[arm][pt]); 
	  gr[arm][pt]->SetLineColor(Counter);
	  if(Counter == 10)
	    gr[arm][pt]->SetLineColor(40);
	  gr[arm][pt]->SetLineWidth(3.0);
	       
	  if(pt == pt_lo)
	    {	     
	      gr[arm][pt]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
	      gr[arm][pt]->GetXaxis()->SetLabelSize(0.04);
	      gr[arm][pt]->GetXaxis()->SetTitleSize(0.04);
	      gr[arm][pt]->GetXaxis()->SetTitleOffset(0.9);
	      gr[arm][pt]->GetXaxis()->SetLimits(2,5);
	      gr[arm][pt]->GetYaxis()->SetLabelSize(0.04);
	      gr[arm][pt]->GetYaxis()->SetTitleSize(0.1); 
	      gr[arm][pt]->GetYaxis()->SetTitleOffset(0.52);
	      gr[arm][pt]->SetMaximum(2);
	      gr[arm][pt]->SetMinimum(0);
	    }
	  Counter++;
	}
    }

  gr[0][pt_lo]->SetTitle(" #frac{pT}{pT Int} Mass Ratio, South");   
 
  for(int pt = pt_lo; pt < pt_hi; pt++)
    {
      if(pt == pt_lo)
	gr[0][pt]->Draw("AL");
      else
	gr[0][pt]->Draw("L");
    }
  
  char unique_gr[800];
  TLegend *leg_gr[2];
  
  for(int arm = 0; arm < 2; arm++)
    {
      leg_gr[arm] = new TLegend(0.51, 0.6, 0.8, 0.9);  
      leg_gr[arm]->SetFillColor(0); 
      leg_gr[arm]->SetTextSize(0.035);
      
      for(int pt = pt_lo; pt < pt_hi+1; pt++)
	{
	  sprintf(unique_gr,"Corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_gr[arm]->AddEntry(gr[arm][pt],unique_gr, "l");
	}
    }
      
  leg_gr[0]->Draw();


  
  // // pT Integrated Histogram Results
  // TCanvas *c2 = new TCanvas("c2","pT Integrated Histogram, South",200,10,700,500);
  // gPad->SetLeftMargin(0.15); 
  // gPad->SetBottomMargin(0.15);
  // gStyle->SetOptStat(0);
  // gPad->SetLogy();
  // gPad->SetGrid();
  // h_pt_int[0]->SetTitle("Yue Hang Corr BG (pT Int), South");
   
 
  // h_pt_int[0]->SetLineColor(5);
  // h_pt_int[0]->SetLineWidth(3.0);
  // h_pt_int[0]->Draw();
  
     

  // FOR LIKE SIGN BACKGROUND
  TGraphErrors *gr_corr_pt = new TGraphErrors(bins,x_array,corr_bg[0][pt_l],x_errors[0],corr_bg_err[0][pt_l]);

  // TGraphErrors *gr_corr_pt = new TGraphErrors(bins,x_array,corr_total[0],x_errors[0],corr_bg_err[0][pt_l]);
   // TGraphErrors *gr_corr_pt = new TGraphErrors(bins,x_array,corr_total[0],x_errors[0],quad[0][pt_l]);
  

  TCanvas *c6 = new TCanvas("c6","Corr BG pT Int",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetGrid();

  //gr_corr_pt->SetTitle("pT Int. Corr BG (cc+bb_UL+dy-bb_LS), North");
  gr_corr_pt->SetTitle("pT Int. Corr BG (cc+bb(UL)+dy-bb(LS)+corr_had(UL)-corr_had(LS)), South");
  gr_corr_pt->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
  gr_corr_pt->GetXaxis()->SetLabelSize(0.04);
  gr_corr_pt->GetXaxis()->SetTitleSize(0.04);
  gr_corr_pt->GetXaxis()->SetTitleOffset(0.9);
  gr_corr_pt->GetXaxis()->SetLimits(0,6);
  gr_corr_pt->GetYaxis()->SetLabelSize(0.04);
  gr_corr_pt->GetYaxis()->SetTitleSize(0.1); 
  gr_corr_pt->GetYaxis()->SetTitleOffset(0.52);
 
  gr_corr_pt->SetMarkerColor(kCyan+1);  
  gr_corr_pt->SetMarkerSize(0.7);
  gr_corr_pt->SetMarkerStyle(20);

  gr_corr_pt->Draw("AP");

  //////////////////////////// Fit histograms //////////////////////////
  
  char name1[800];
  char name2[800];

  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;
  int bins_fail[bins];

  TF1 *mix_ul_fit;
  TF1 *mix_ul_fit_norm;
  double fit_area[2][38];
  double sum_area[2] = {0.0,0.0};
  double sum_norm_area[2] = {0.0,0.0};
  double sum_fit_area[2] = {0.0,0.0};
  double f5,f6,f10,f11,f12;
   
  mix_ul_fit = new TF1("mix_ul_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);
  //mix_ul_fit = new TF1("mix_ul_fit"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",lim_a,lim_b);

  for(int arm = 0; arm < 1; arm++)   /////////SELECT ARM TO FIT
    {
      //for(int pt = 0; pt < num_bins_on_plots[pt_binwidth]; pt++)
      //	for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
      for(int pt = pt_l; pt < pt_l+1; pt++)
      //for(int pt = 0; pt < 1; pt++)
      {
	  cout << "Arm " << arm << ", Bin " << pt+1 << endl;
	
	  // For fitting bg as function of pT and then plotting parameter values as function of pT
	  //**************************************************************************************
	  f5 = 1;
	  f6 = 0.1;
	  f10 = 100;
	  f11 = 1;
	  f12 = 1;

	  // f5 = 1;
	  // f6 = 0.1;
	  // f10 = 100;
	  // f11 = 1; 
	  // f12 = 1;

	  // if(pt_l == 2)
	  //   {
	  //     f5 = 1;
	  //     f6 = 1;
	  //     f10 = 100;
	  //     f11 = 1;
	  //     f12 = 1;
	  //   }

	  // if(pt_l == 0)
	  //   {
	  //     f5 = 1;
	  //     f6 = 0.01;
	  //     f10 = 100;
	  //     f11 = 1;
	  //     f12 = 0.01;
	  //   }
	  
	  if((pt_l == 3))
	    {
	      f5 = 1;
	      f6 = 0.01;
	      f10 = 100;
	      f11 = 1;
	      f12 = 10;
	    }

	  if((pt_l == 4) || (pt_l == 5) || (pt_l == 6) || (pt_l == 7) || (pt_l == 8) || (pt_l == 9) || (pt_l == 0) || (pt_l == 1) || (pt_l == 2) || (pt_l == 3))
	    {
	      f5 = 1;
	      f6 = 0.01;
	      f10 = 100;
	      f11 = 10;
	      f12 = 10;
	    }
	  
	  if((pt_l == 10))
	    {
	      f5 = 1;
	      f6 = 0.01;
	      f10 = 1;
	      f11 = 10;
	      f12 = 10;
	    }
	
	  mix_ul_fit->SetParameter(0, f5);
	  mix_ul_fit->SetParameter(1, f6);
	  mix_ul_fit->SetParameter(2, f10);
	  mix_ul_fit->SetParameter(3, f11);
	  mix_ul_fit->SetParameter(4, f12);
	 
	  // for calculating the parameter from a function plugging in the pt value
	  //************************************************************************
	  /*
	  f10 = 100;

	  if((pt_l == 8) || (pt_l == 9))
	    f10 = 10;

	  double par_a[10];
	  double par_b[10];
	  double par_d[10];
	  double par_e[10];

	  double par_a_fx;
	  double par_b_fx;
	  double par_d_fx;
	  double par_e_fx;

	  double x;

	  for(int par = 0; par < 10; par++)
	    {	 
	      par_a[par] = coeff_par[0][par];
	      par_b[par] = coeff_par[1][par];
	      par_d[par] = coeff_par[2][par];
	      par_e[par] = coeff_par[3][par];
	      // cout << coeff_par[0][par] << endl;
	      // cout << coeff_par[1][par] << endl;
	      // cout << coeff_par[2][par] << endl;
	      // cout << coeff_par[3][par] << endl;
	      cout << par_b[par] << endl;
	      
	    }

	  x = pt_center[2][pt_l];
	  cout << "x = pT center: " << x << endl;
	  
	  par_a_fx = par_a[0] + par_a[1]*x + par_a[2]*x*x + par_a[3]*x*x*x + par_a[4]*x*x*x*x + par_a[5]*x*x*x*x*x + par_a[6]*x*x*x*x*x*x + par_a[7]*x*x*x*x*x*x*x + par_a[8]*x*x*x*x*x*x*x*x + par_a[9]*x*x*x*x*x*x*x*x*x;
	  par_b_fx = par_b[0] + par_b[1]*x + par_b[2]*x*x + par_b[3]*x*x*x + par_b[4]*x*x*x*x + par_b[5]*x*x*x*x*x + par_b[6]*x*x*x*x*x*x + par_b[7]*x*x*x*x*x*x*x + par_b[8]*x*x*x*x*x*x*x*x + par_b[9]*x*x*x*x*x*x*x*x*x;
	  par_d_fx = par_d[0] + par_d[1]*x + par_d[2]*x*x + par_d[3]*x*x*x + par_d[4]*x*x*x*x + par_d[5]*x*x*x*x*x + par_d[6]*x*x*x*x*x*x + par_d[7]*x*x*x*x*x*x*x + par_d[8]*x*x*x*x*x*x*x*x + par_d[9]*x*x*x*x*x*x*x*x*x;
	  par_e_fx = par_e[0] + par_e[1]*x + par_e[2]*x*x + par_e[3]*x*x*x + par_e[4]*x*x*x*x + par_e[5]*x*x*x*x*x + par_e[6]*x*x*x*x*x*x + par_e[7]*x*x*x*x*x*x*x + par_e[8]*x*x*x*x*x*x*x*x + par_e[9]*x*x*x*x*x*x*x*x*x;
	  
	  cout << "a0: " << par_a[0] << ", a1: " << par_a[1] << ", a2: " << par_a[2] << "a3: " << par_a[3] << ", a4: " << par_a[4] << ", a5: " << par_a[5] << ", a6: " << par_a[6] << ", a7: " << par_a[7] << ", a8: " << par_a[8] << ", a9: " << par_a[9] << endl;

	  cout << "b0: " << par_b[0] << ", b1: " << par_b[1] << ", b2: " << par_b[2] << "b3: " << par_b[3] << ", b4: " << par_b[4] << ", b5: " << par_b[5] << ", b6: " << par_b[6] << ", b7: " << par_b[7] << ", b8: " << par_b[8] << ", b9: " << par_b[9] << endl;

	  cout << "d0: " << par_d[0] << ", d1: " << par_d[1] << ", d2: " << par_d[2] << "d3: " << par_d[3] << ", d4: " << par_d[4] << ", d5: " << par_d[5] << ", d6: " << par_d[6] << ", d7: " << par_d[7] << ", d8: " << par_d[8] << ", d9: " << par_d[9] << endl;

	  cout << "e0: " << par_e[0] << ", e1: " << par_e[1] << ", e2: " << par_e[2] << "e3: " << par_e[3] << ", e4: " << par_e[4] << ", e5: " << par_e[5] << ", e6: " << par_e[6] << ", e7: " << par_e[7] << ", e8: " << par_e[8] << ", e9: " << par_e[9] << endl;
	  
	  mix_ul_fit->FixParameter(0, par_a_fx);
	  mix_ul_fit->FixParameter(1, par_b_fx);
	  mix_ul_fit->SetParameter(2, f10);
	  mix_ul_fit->FixParameter(3, par_d_fx);
	  mix_ul_fit->FixParameter(4, par_e_fx);
	  */
	  //****************************************************************************

     	  mix_ul_fit->SetLineColor(kRed);
	  mix_ul_fit->SetLineStyle(5);
	  mix_ul_fit->SetLineWidth(2);
	 
	  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
	     
	  // TFitResultPtr r = h_sum[0][pt]->Fit(mix_ul_fit,"S R");
	  // h_sum[arm][pt_lo]->SetAxisRange(2,5); 	    
	 
	  //TFitResultPtr r = gr_corr_pt[0].Fit(mix_ul_fit,"S R");

	  TFitResultPtr r = gr_corr_pt->Fit(mix_ul_fit,"S R");	 
	 	  	  
	  // fit_area[arm][pt] = mix_ul_fit->Integral(h_sum[0][0]->FindBin(2.001), h_sum[0][0]->FindBin(4.999) );
	 	 
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
	    
	  Int_t fitStatus = r;  
	  //r->Print("V");
	  cout << "status: " << fitStatus << " for pt bin " << pt_l << endl;

	  if(fitStatus == 0) // 0 indicates it converged, and 4 indicates it failed 

	    // https://root.cern.ch/doc/master/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd

	    {
	      if(write)	     
		{
		  if(arm == 1)
		    {
		      sprintf(name1,"yuehang_bestfit_500mev/bestfit_parameters_N_%i.dat",pt_l); 
		      std::fstream bestfit_parameters_1(name1,std::ofstream::out); 
		      bestfit_parameters_1 <<  p0  <<  " " << e0 << " " << p1 << " " <<  e1 << " " << p2 <<  " " << e2 << " " << p3 << " " << e3 << " " << p4 << " "  << e4 << " " << endl;
		      bestfit_parameters_1.close();
		    }
		  else
		    {
		      sprintf(name2,"yuehang_bestfit_500mev/bestfit_parameters_S_%i.dat",pt_l); 
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
  
  
      // sum_area[0] = 0.0;
      // sum_area[1] = 0.0;
  
      // for(int arm = 0; arm < 2; arm++)
      //   {
      //     for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
  // 	{
  // 	  sum_area[arm] += area[arm][pt];
  // 	  sum_fit_area[arm] += fit_area[arm][pt];

  // 	  cout << "For pt slice " << pt+1 << ", TH1D un-norm area is: " << area[arm][pt] << ", TH1D norm area: " << norm_area[arm][pt] << ", sum of un-norm area: " << sum_area[arm] << endl;
  // 	  cout << "Fit area: " << fit_area[arm][pt] << ", and fit area sum: " << sum_fit_area[arm] << endl;

  // 	  //cout << "bin low: " << h_sum[arm][pt]->FindBin(2.001) << ", and high: " << h_sum[arm][pt]->FindBin(4.999) << endl;
  // 	} 
  //   }

 
  

  h_sum[1][0]->SetMinimum(pow(10,-3));
  h_sum[1][0]->SetMaximum(pow(10,0));  
  h_sum[0][pt_lo]->SetMinimum(pow(10,-3));
  h_sum[0][pt_lo]->SetMaximum(pow(10,0));
  
  gr_corr_pt[0].SetMaximum(pow(10,4));
  gr_corr_pt[0].SetMinimum(pow(10,1));

  //  h_pt_int[0]->SetAxisRange(2,5); 
  


  	  

}// void end macro

  
  
  
  
