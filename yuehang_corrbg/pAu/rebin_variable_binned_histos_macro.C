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

double Interpolate(double m,TH1D *h)
{
  double RpA = 0;
  int bin = h->FindBin(m);
   
  if( (bin - 1 < 1) || (bin + 1 > h->GetNbinsX() - 1))
    return RpA;

  // double m_1 = h->GetBinCenter(bin - 1);
  // double m_2 = h->GetBinCenter(bin);
  // double m_3 = h->GetBinCenter(bin+1);

  // double RpA_1 = h->GetBinContent(bin - 1);
  // double RpA_2 =  h->GetBinContent(bin);
  // double RpA_3 =  h->GetBinContent(bin + 1);
  
  //int bin = bin_2;
  
  double sum = 0;
  double weight = 0;
  for(int i = bin - 1; i <= bin + 1; i++)
    {
      if(h->GetBinContent(i) < 0)
	{
	  double temp = 0;
	  for(int j = 0; j < 5; j++)
	    {
	      if( h->GetBinContent(i- j) > 0)
		{
		  temp += h->GetBinContent(i - j);
		  break;
		}
	    }
	  for(int j = 0; j < 5; j++)
	    {
	      if( h->GetBinContent(i+ j) > 0)
		{
		  temp += h->GetBinContent(i+ j);
		  break;
		}
	    }
	  temp /= 2.0;
	  h->SetBinContent(i,temp);
	}
    }
  
  if(0)
    {     
      double slope = (h->GetBinContent(bin+1) - h->GetBinContent(bin-1)) / (h->GetBinCenter(bin+1) - h->GetBinCenter(bin-1));
      
      RpA = h->GetBinContent(bin-1) + slope*(m - h->GetBinCenter(bin-1));
    }

  RpA = (h->GetBinContent(bin-1) + h->GetBinContent(bin)  + h->GetBinContent(bin+1) )/3;

  //cout << "m1: " << m_1 << ", m_2: " << m_2 << ", m_3: " << m_3 << ", m: " << m << ", RpA: " << RpA << ", slope: " << slope << ", RpA_1: " << RpA_1 << ", RpA_2: " << RpA_2 << ", RpA_3: " << RpA_3 << endl;
  
  return RpA;
  
}

void rebin_variable_binned_histos_macro()
{
  
  bool write = true;
  
  TH1D *temp_had[2][2][4];
  
  std::string obj_filename[2][2][4] = {
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
	      file_corrhad->GetObject(obj_filename[arm][i_histo][pt].c_str(),temp_had[arm][i_histo][pt]);
	    }
	}	  
    }

 for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      cout << "obj_filename[" << arm << "][" << i_histo << "][" << pt << "]: " << obj_filename[arm][i_histo][pt].c_str() << endl;
	    }
	}	  
    }
  
  TFile *fout; 

  std::string rebinned_filename[2] = {"yuehang_Run15pAu_S_rebinned_corr_had.root","yuehang_Run15pAu_N_rebinned_corr_had.root"};
 
  char unique1[500];
  char unique2[500];
 
 ///////////////////////////////////////////// Generate mass bin array for desired rebinning ///////////////////////////
 
  int mass_slices = 100;  
  double mass_width[100];
  double mass_center[100];
  double bin_edges[100] = {0};
  double binwidth = 0.1;
  int last_element= mass_slices - 1;

  bin_edges[0] = 0;

  for(int i = 1; i < mass_slices; i++)
    {
      bin_edges[i]  = i*binwidth;
      mass_width[i] = binwidth;
      mass_center[i-1] = (bin_edges[i] + bin_edges[i-1])/2;
      if(i == (mass_slices - 1) )
	mass_center[last_element] = bin_edges[i] + binwidth/2;
    }
 
  for(int i = 0; i < 100; i++)
    {
      cout << mass_center[i] << endl;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  int a[4] = {0,1,2,3};
  int b[4] = {1,2,3,8};
 
  int numbinsx[2][2][4] = {0};
  TH1D *corrhad[2][2][4] = {0};
  int numBins;

  TAxis *xaxis[2][2][4]; 

  double bincenter_0[2][2][4][100];
  double binwidth_0[2][2][4][100];

  double corr_hadrons_mod[2][2][4][100];
  double RpA[2][2][4][100] = {0};

  double x1,x2,x3,y1,y2,y3,aa,f_aa,slope,mass,rpa;
  int pt_low, pt_high;
  
  char unique_ch[2][2][4][500];
 
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {
	      if(pt < 3)
		{
		  sprintf(unique_ch[arm][0][pt],"corr_had_UL_%i_%i",pt,pt+1);
		  sprintf(unique_ch[arm][1][pt],"corr_had_LS_%i_%i",pt,pt+1);
		}
	      else
		{
		  sprintf(unique_ch[arm][0][pt],"corr_had_UL_%i_%i",pt,pt+5);
		  sprintf(unique_ch[arm][1][pt],"corr_had_LS_%i_%i",pt,pt+5);
		}
	      cout << unique_ch[arm][i_histo][pt] << " " << endl;
	    }
	}
    }
  
  TH1D *h1 = new TH1D("h1","h1 distribution",100,0,10); 
   
  // Generate 16 different TH1D histograms
   for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {  
	      corrhad[arm][i_histo][pt] =  (TH1D *) h1->Clone(unique_ch[arm][i_histo][pt]);
	    }
	}
    }

   // Calculate and fill TH1D array with the interpolated RpA values
  for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {  
	      for(int k = 0; k < 100; k++)
		{
		  mass = mass_center[k];
		  double rpa = Interpolate(mass,temp_had[arm][i_histo][pt]);
		  corrhad[arm][i_histo][pt]->SetBinContent(k+1, rpa);
		  corrhad[arm][i_histo][pt]->SetBinError(k+1, rpa*0.005);

		}
	    }	  
	}
    }
   
  // Check the RpA interpolated values
 for(int arm = 0; arm < 2; arm++)
    {
      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {  
	      for(int k = 0; k < 100; k++)
		{
		  // cout << "corrhad[" << arm << "][" << i_histo << "][" << pt << "][" << k << "] = " << corrhad[arm][i_histo][pt]->GetBinContent(k+1) << endl;
		}
	    }
	}
    }

 // Write the RpA interpolated values to unique TH1D histograms in same root file
 for(int arm = 0; arm < 2; arm++)
    {
      fout = new TFile(rebinned_filename[arm].c_str(), "RECREATE");

      for(int i_histo = 0; i_histo < 2; i_histo++)
	{
	  for(int pt = 0; pt < 4; pt++)
	    {  
	      if(write == true)
		{
		  corrhad[arm][i_histo][pt]->Write(); 
		}
	    }
	}
    }

 fout->Close();

}// void end macro
