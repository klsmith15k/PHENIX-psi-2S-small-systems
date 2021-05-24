
// for special binning (to match any PPG width)
// Need to run this twice for each rebinning (north and south arms)

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


void rebinning_macro_ptint_cent()
{
  
  bool North = true;


  char name4[500];
  TFile *fout;
  TFile *file_data;	  


  if(North == true)
    {
      sprintf(name4,"psi2S_root_files/fvtx_batch105108_nofvtxchi2_fid188N_badprob_test_N.root");
      fout = new TFile(name4, "RECREATE");	  
       file_data = TFile::Open("matt_centrality_macros/CONDOR/fvtx_batch105108_nofvtxchi2_fid188N_badprob_test/result.root");
     
    }
else
    {
      sprintf(name4,"psi2S_root_files/fvtx_batch105108_nofvtxchi2_fid188N_badprob_test_S.root");
      fout = new TFile(name4, "RECREATE");	  
        file_data = TFile::Open("matt_centrality_macros/CONDOR/fvtx_batch105108_nofvtxchi2_fid188N_badprob_test/result.root");
    
 }

// if(North == true)
//     {
//       sprintf(name4,"psi2S_root_files/yes_sngdbl_test_N.root");
//       fout = new TFile(name4, "RECREATE");	  
//        file_data = TFile::Open("psi2S_root_files/sngdbl_test.root");
     
//     }
// else
//     {
//       sprintf(name4,"psi2S_root_files/yes_sng_test_S.root");
//       fout = new TFile(name4, "RECREATE");	  
//       file_data = TFile::Open("psi2S_root_files/sng_test.root");
    
//  }
 
  TH1D *ul[7];  // only works if arrays are constant static integers
  TH1D *pp[7];
  TH1D *mm[7];
  TH1D *mix_ul[7];
  
  TH2D *t[7];
  TH2D *u[7];
  TH2D *v[7];
  TH2D *w[7];
   
  int j = 0;
  int k = 0;
  char unique[500];
  double sum = 0.0;
  
  
  if(North == true)
    {
      file_data->GetObject("mass_fvtx_pt_N_ul",t[0]);
      file_data->GetObject("mass_fvtx_pt_N_pp",u[0]);
      file_data->GetObject("mass_fvtx_pt_N_mm",v[0]);
       file_data->GetObject("mass_fvtx_pt_N_mix",w[0]);

       file_data->GetObject("mass_fvtx_pt_020_N_ul",t[1]);
      file_data->GetObject("mass_fvtx_pt_020_N_pp",u[1]);
      file_data->GetObject("mass_fvtx_pt_020_N_mm",v[1]);
      file_data->GetObject("mass_fvtx_pt_020_N_mix",w[1]);

      file_data->GetObject("mass_fvtx_pt_2040_N_ul",t[2]);
      file_data->GetObject("mass_fvtx_pt_2040_N_pp",u[2]);
      file_data->GetObject("mass_fvtx_pt_2040_N_mm",v[2]);
       file_data->GetObject("mass_fvtx_pt_2040_N_mix",w[2]);

      file_data->GetObject("mass_fvtx_pt_4060_N_ul",t[3]);
      file_data->GetObject("mass_fvtx_pt_4060_N_pp",u[3]);
      file_data->GetObject("mass_fvtx_pt_4060_N_mm",v[3]);
        file_data->GetObject("mass_fvtx_pt_4060_N_mix",w[3]);

      file_data->GetObject("mass_fvtx_pt_6084_N_ul",t[4]);
      file_data->GetObject("mass_fvtx_pt_6084_N_pp",u[4]);
      file_data->GetObject("mass_fvtx_pt_6084_N_mm",v[4]);
       file_data->GetObject("mass_fvtx_pt_6084_N_mix",w[4]);

      file_data->GetObject("mass_fvtx_pt_4084_N_ul",t[5]);
      file_data->GetObject("mass_fvtx_pt_4084_N_pp",u[5]);
      file_data->GetObject("mass_fvtx_pt_4084_N_mm",v[5]);
       file_data->GetObject("mass_fvtx_pt_4084_N_mix",w[5]);

      file_data->GetObject("mass_fvtx_pt_2084_N_ul",t[6]);
      file_data->GetObject("mass_fvtx_pt_2084_N_pp",u[6]);
      file_data->GetObject("mass_fvtx_pt_2084_N_mm",v[6]);
        file_data->GetObject("mass_fvtx_pt_2084_N_mix",w[6]);

      // file_data->GetObject("mass_fvtx_pt_6084_N_ul",t[7]);
      // file_data->GetObject("mass_fvtx_pt_6084_N_pp",u[7]);
      // file_data->GetObject("mass_fvtx_pt_6084_N_mm",v[7]);
      //  file_data->GetObject("mass_fvtx_pt_6084_N_mix",w[7]);

      // file_data->GetObject("mass_fvtx_pt_015_N_ul",t[8]);
      // file_data->GetObject("mass_fvtx_pt_015_N_pp",u[8]);
      // file_data->GetObject("mass_fvtx_pt_015_N_mm",v[8]);
      // file_data->GetObject("mass_fvtx_pt_015_N_mix",w[8]);

      // file_data->GetObject("mass_fvtx_pt_1530_N_ul",t[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_N_pp",u[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_N_mm",v[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_N_mix",w[9]);

      // file_data->GetObject("mass_fvtx_pt_3050_N_ul",t[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_N_pp",u[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_N_mm",v[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_N_mix",w[10]);

      // file_data->GetObject("mass_fvtx_pt_5084_N_ul",t[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_N_pp",u[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_N_mm",v[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_N_mix",w[11]);

      // file_data->GetObject("mass_fvtx_pt_515_N_ul",t[12]);
      // file_data->GetObject("mass_fvtx_pt_515_N_pp",u[12]);
      // file_data->GetObject("mass_fvtx_pt_515_N_mm",v[12]);
      // file_data->GetObject("mass_fvtx_pt_515_N_mix",w[12]);

      // file_data->GetObject("mass_fvtx_pt_1030_N_ul",t[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_N_pp",u[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_N_mm",v[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_N_mix",w[13]);

    }
  else
    {
      file_data->GetObject("mass_fvtx_pt_S_ul",t[0]);
      file_data->GetObject("mass_fvtx_pt_S_pp",u[0]);
      file_data->GetObject("mass_fvtx_pt_S_mm",v[0]);
      file_data->GetObject("mass_fvtx_pt_S_mix",w[0]);

      file_data->GetObject("mass_fvtx_pt_020_S_ul",t[1]);
      file_data->GetObject("mass_fvtx_pt_020_S_pp",u[1]);
      file_data->GetObject("mass_fvtx_pt_020_S_mm",v[1]);
      file_data->GetObject("mass_fvtx_pt_020_S_mix",w[1]);

      file_data->GetObject("mass_fvtx_pt_2040_S_ul",t[2]);
      file_data->GetObject("mass_fvtx_pt_2040_S_pp",u[2]);
      file_data->GetObject("mass_fvtx_pt_2040_S_mm",v[2]);
       file_data->GetObject("mass_fvtx_pt_2040_S_mix",w[2]);

      file_data->GetObject("mass_fvtx_pt_4060_S_ul",t[3]);
      file_data->GetObject("mass_fvtx_pt_4060_S_pp",u[3]);
      file_data->GetObject("mass_fvtx_pt_4060_S_mm",v[3]);
      file_data->GetObject("mass_fvtx_pt_4060_S_mix",w[3]);

      file_data->GetObject("mass_fvtx_pt_6084_S_ul",t[4]);
      file_data->GetObject("mass_fvtx_pt_6084_S_pp",u[4]);
      file_data->GetObject("mass_fvtx_pt_6084_S_mm",v[4]);
      file_data->GetObject("mass_fvtx_pt_6084_S_mix",w[4]);

      file_data->GetObject("mass_fvtx_pt_4084_S_ul",t[5]);
      file_data->GetObject("mass_fvtx_pt_4084_S_pp",u[5]);
      file_data->GetObject("mass_fvtx_pt_4084_S_mm",v[5]);
       file_data->GetObject("mass_fvtx_pt_4084_S_mix",w[5]);

      file_data->GetObject("mass_fvtx_pt_2084_S_ul",t[6]);
      file_data->GetObject("mass_fvtx_pt_2084_S_pp",u[6]);
      file_data->GetObject("mass_fvtx_pt_2084_S_mm",v[6]);
      file_data->GetObject("mass_fvtx_pt_2084_S_mix",w[6]);

      // file_data->GetObject("mass_fvtx_pt_6084_S_ul",t[7]);
      // file_data->GetObject("mass_fvtx_pt_6084_S_pp",u[7]);
      // file_data->GetObject("mass_fvtx_pt_6084_S_mm",v[7]);
      //  file_data->GetObject("mass_fvtx_pt_6084_S_mix",w[7]);

      // file_data->GetObject("mass_fvtx_pt_015_S_ul",t[8]);
      // file_data->GetObject("mass_fvtx_pt_015_S_pp",u[8]);
      // file_data->GetObject("mass_fvtx_pt_015_S_mm",v[8]);
      // file_data->GetObject("mass_fvtx_pt_015_S_mix",w[8]);

      // file_data->GetObject("mass_fvtx_pt_1530_S_ul",t[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_S_pp",u[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_S_mm",v[9]);
      // file_data->GetObject("mass_fvtx_pt_1530_S_mix",w[9]);

      // file_data->GetObject("mass_fvtx_pt_3050_S_ul",t[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_S_pp",u[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_S_mm",v[10]);
      // file_data->GetObject("mass_fvtx_pt_3050_S_mix",w[10]);

      // file_data->GetObject("mass_fvtx_pt_5084_S_ul",t[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_S_pp",u[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_S_mm",v[11]);
      // file_data->GetObject("mass_fvtx_pt_5084_S_mix",w[11]);

      // file_data->GetObject("mass_fvtx_pt_515_S_ul",t[12]);
      // file_data->GetObject("mass_fvtx_pt_515_S_pp",u[12]);
      // file_data->GetObject("mass_fvtx_pt_515_S_mm",v[12]);
      // file_data->GetObject("mass_fvtx_pt_515_S_mix",w[12]);

      // file_data->GetObject("mass_fvtx_pt_1030_S_ul",t[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_S_pp",u[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_S_mm",v[13]);
      // file_data->GetObject("mass_fvtx_pt_1030_S_mix",w[13]);
    }
  

  for(int i = 0; i < 7; i++)
    { 
      sprintf(unique,"ul_%d",i+1); 
      ul[i] = t[i]->ProjectionX(unique); 
      cout << "ul entries " << ul[i]->GetEntries() << " for bin " << i+1 << endl;
      fout->cd();  
      ul[i]->Write(); 
      if(i==0)
	sum+= ul[i]->GetEntries();
	     
    }

  for(int i = 0; i < 7; i++)
    { 
      sprintf(unique,"pp_%d",i+1); 
      pp[i] = u[i]->ProjectionX(unique); 
      cout << "ul entries " << pp[i]->GetEntries() << " for bin " << i+1 << endl;
      fout->cd();  
      pp[i]->Write(); 
     	     
    }

  for(int i = 0; i < 7; i++)
    { 
      sprintf(unique,"mm_%d",i+1); 
      mm[i] = v[i]->ProjectionX(unique); 
      cout << "ul entries " << mm[i]->GetEntries() << " for bin " << i+1 << endl;
      fout->cd();  
      mm[i]->Write(); 
    	     
    }

  for(int i = 0; i < 7; i++)
    { 
      sprintf(unique,"mix_ul_%d",i+1); 
      mix_ul[i] = w[i]->ProjectionX(unique); 
      cout << "ul entries " << mix_ul[i]->GetEntries() << " for bin " << i+1 << endl;
      fout->cd();  
      mix_ul[i]->Write(); 
     	     
    }

	 
  fout->Close();
  //cout << "ul total entries " << sum << endl;

}// void macro





















