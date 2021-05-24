
#include <Math/MinimizerOptions.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TEfficiency.h>
#include <iostream>
#include <TProfile2D.h>

using namespace std;

void mass_sim_rebinning_macro_cent_newZ()
{


   bool north_arm = true;


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
    
  char name4[500];
  char name5[500];
  char name6[500];
  char name7[500];


  bool embed = true;

  TFile *file_data;	 
 
  
  if(north_arm == true && !embed)
    {
      file_data = TFile::Open("HEPMC_info_noembed_300.root"); 
      cout << "north arm no embed" << endl;
    }
  if(north_arm == false && !embed)
    {
      file_data = TFile::Open("HEPMC_info_noembed_300.root"); 
      cout << "south arm no embed" << endl;
    }
  if(north_arm == true && embed)
    {
       file_data = TFile::Open("sanghoon_simulations/fvtx_jpsi_newlib_trig_FID188N_notrigZ_all.root");
       //  file_data = TFile::Open("sanghoon_simulations/fvtx_psi2s_newlib_trig_FID188N_notrigZ_all.root");
      // file_data = TFile::Open("sanghoon_simulations/fvtx_jpsi_newlib_acceff_FID188N_nochi2fvtx_trigZ_all.root");
      //   file_data = TFile::Open("sanghoon_simulations/fvtx_psi2s_newlib_acceff_num_NOFID_all.root");
      cout << "north arm embed" << endl;
    }
  if(north_arm == false && embed)
    {
       file_data = TFile::Open("sanghoon_simulations/fvtx_jpsi_newlib_trig_FID188N_notrigZ_all.root");
       //  file_data = TFile::Open("sanghoon_simulations/fvtx_psi2s_newlib_trig_FID188N_notrigZ_all.root");
      //  file_data = TFile::Open("sanghoon_simulations/fvtx_jpsi_newlib_acceff_FID188N_nochi2fvtx_trigZ_all.root");
      //  file_data = TFile::Open("sanghoon_simulations/fvtx_psi2s_newlib_acceff_num_NOFID_all.root");
       cout << "south arm embed" << endl;
    }
  

  int pt_nbins = 27;
  int pt_nbins2 = 5;
  int pt_nbins3 = 1;
  
  double pt_edges[28] = {0,0.25,0.5,0.75,1.0,1.25,1.50,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.50,4.75,5.0,5.25,5.50,5.75,6.0,6.50,7.0,8.0};
  double pt_1GeV_edges[6] = {0,1.0,2.0,3.0,4.0,6.0};
  double pt_int_edges[2] = {0,12.0};

  ///////////////////////////////////////

  const int arm = 2;
  
  TH3F *t0 = new TH3F("t0","",300,0,15,48,0,12,60,-30,30);
  TH3F *t0_pp = new TH3F("t0_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t0_mm = new TH3F("t0_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t1 = new TH3F("t1","",300,0,15,48,0,12,60,-30,30);
  TH3F *t1_pp = new TH3F("t1_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t1_mm = new TH3F("t1_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t2 = new TH3F("t2","",300,0,15,48,0,12,60,-30,30);
  TH3F *t2_pp = new TH3F("t2_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t2_mm = new TH3F("t2_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t3 = new TH3F("t3","",300,0,15,48,0,12,60,-30,30);
  TH3F *t3_pp = new TH3F("t3_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t3_mm = new TH3F("t3_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t4 = new TH3F("t4","",300,0,15,48,0,12,60,-30,30);
  TH3F *t4_pp = new TH3F("t4_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t4_mm = new TH3F("t4_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t5 = new TH3F("t5","",300,0,15,48,0,12,60,-30,30);
  TH3F *t5_pp = new TH3F("t5_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t5_mm = new TH3F("t5_mm","",300,0,15,48,0,12,60,-30,30);

  TH3F *t6 = new TH3F("t6","",300,0,15,48,0,12,60,-30,30);
  TH3F *t6_pp = new TH3F("t6_pp","",300,0,15,48,0,12,60,-30,30);
  TH3F *t6_mm = new TH3F("t6_mm","",300,0,15,48,0,12,60,-30,30);
   

  /// copied th FillAnaTree module after running Run15pAu, and left the centrality formatting.  Can remove the centrality histograms 
    if(north_arm)
      {
	file_data->GetObject("hmass_N_fvtx_010",t0);  
	file_data->GetObject("hmass_N_fvtx_pp_010",t0_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_010",t0_mm);  
	
	file_data->GetObject("hmass_N_fvtx_1030",t1);  
	file_data->GetObject("hmass_N_fvtx_pp_1030",t1_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_1030",t1_mm);  
	
	file_data->GetObject("hmass_N_fvtx_020",t2);  
	file_data->GetObject("hmass_N_fvtx_pp_020",t2_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_020",t2_mm);  
	
	file_data->GetObject("hmass_N_fvtx_2040",t3);  
	file_data->GetObject("hmass_N_fvtx_pp_2040",t3_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_2040",t3_mm);  
	
	file_data->GetObject("hmass_N_fvtx_4060",t4);  
	file_data->GetObject("hmass_N_fvtx_pp_4060",t4_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_4060",t4_mm);  
	
	file_data->GetObject("hmass_N_fvtx_6084",t5);  
	file_data->GetObject("hmass_N_fvtx_pp_6084",t5_pp);  
	file_data->GetObject("hmass_N_fvtx_mm_6084",t5_mm);  
	
	file_data->GetObject("hmass_N_fvtx",t6);  
	file_data->GetObject("hmass_N_fvtx_pp",t6_pp);  
	file_data->GetObject("hmass_N_fvtx_mm",t6_mm);  
      }
    else
      {
	file_data->GetObject("hmass_S_fvtx_010",t0);  
	file_data->GetObject("hmass_S_fvtx_pp_010",t0_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_010",t0_mm);  
	
	file_data->GetObject("hmass_S_fvtx_1030",t1);  
	file_data->GetObject("hmass_S_fvtx_pp_1030",t1_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_1030",t1_mm);  
	
	file_data->GetObject("hmass_S_fvtx_020",t2);  
	file_data->GetObject("hmass_S_fvtx_pp_020",t2_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_020",t2_mm);  
	
	file_data->GetObject("hmass_S_fvtx_2040",t3);  
	file_data->GetObject("hmass_S_fvtx_pp_2040",t3_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_2040",t3_mm);  
	
	file_data->GetObject("hmass_S_fvtx_4060",t4);  
	file_data->GetObject("hmass_S_fvtx_pp_4060",t4_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_4060",t4_mm);  
	
	file_data->GetObject("hmass_S_fvtx_6084",t5);  
	file_data->GetObject("hmass_S_fvtx_pp_6084",t5_pp);  
	file_data->GetObject("hmass_S_fvtx_mm_6084",t5_mm);  
	
	file_data->GetObject("hmass_S_fvtx",t6);  
	file_data->GetObject("hmass_S_fvtx_pp",t6_pp);  
	file_data->GetObject("hmass_S_fvtx_mm",t6_mm);
      }

  t0->GetXaxis()->SetRangeUser(2,5);
  t0_pp->GetXaxis()->SetRangeUser(2,5);
  t0_mm->GetXaxis()->SetRangeUser(2,5);

  t1->GetXaxis()->SetRangeUser(2,5);
  t1_pp->GetXaxis()->SetRangeUser(2,5);
  t1_mm->GetXaxis()->SetRangeUser(2,5);

  t2->GetXaxis()->SetRangeUser(2,5);
  t2_pp->GetXaxis()->SetRangeUser(2,5);
  t2_mm->GetXaxis()->SetRangeUser(2,5);

  t3->GetXaxis()->SetRangeUser(2,5);
  t3_pp->GetXaxis()->SetRangeUser(2,5);
  t3_mm->GetXaxis()->SetRangeUser(2,5);

  t4->GetXaxis()->SetRangeUser(2,5);
  t4_pp->GetXaxis()->SetRangeUser(2,5);
  t4_mm->GetXaxis()->SetRangeUser(2,5);

  t5->GetXaxis()->SetRangeUser(2,5);
  t5_pp->GetXaxis()->SetRangeUser(2,5);
  t5_mm->GetXaxis()->SetRangeUser(2,5);

  t6->GetXaxis()->SetRangeUser(2,5);
  t6_pp->GetXaxis()->SetRangeUser(2,5);
  t6_mm->GetXaxis()->SetRangeUser(2,5);
 

  // lo and hi bins will be the same in all histograms since identical binning
  int hi_xbin = t0->GetXaxis()->FindBin(8);
  int lo_xbin = t0->GetXaxis()->FindBin(0.01);
   
  cout << lo_xbin << endl;
  cout << hi_xbin << endl;
  
  
  double pt_t1[400] = {0};
  double pt_t2[400] = {0};
  ////// ul
  TH1D *th0_N = new TH1D("th0_N","th0_N",48,0,12); 
  TH1D *th1_N = new TH1D("th1_N","th1_N",48,0,12); 
  TH1D *th2_N = new TH1D("th2_N","th2_N",48,0,12);  
  TH1D *th3_N = new TH1D("th3_N","th3_N",48,0,12);  
  TH1D *th4_N = new TH1D("th4_N","th4_N",48,0,12);  
  TH1D *th5_N = new TH1D("th5_N","th5_N",48,0,12);  
  TH1D *th6_N = new TH1D("th6_N","th6_N",48,0,12);  
  
  TH1D *th0_S = new TH1D("th0_S","th0_S",48,0,12); 
  TH1D *th1_S = new TH1D("th1_S","th1_S",48,0,12);  
  TH1D *th2_S = new TH1D("th2_S","th2_S",48,0,12);  
  TH1D *th3_S = new TH1D("th3_S","th3_S",48,0,12);  
  TH1D *th4_S = new TH1D("th4_S","th4_S",48,0,12);  
  TH1D *th5_S = new TH1D("th5_S","th5_S",48,0,12);  
  TH1D *th6_S = new TH1D("th6_S","th6_S",48,0,12);  
  ////// mm
  TH1D *th0_pp_N = new TH1D("th0_pp_N","th0_pp_N",48,0,12);  
  TH1D *th1_pp_N = new TH1D("th1_pp_N","th1_pp_N",48,0,12);  
  TH1D *th2_pp_N = new TH1D("th2_pp_N","th2_pp_N",48,0,12);  
  TH1D *th3_pp_N = new TH1D("th3_pp_N","th3_pp_N",48,0,12);  
  TH1D *th4_pp_N = new TH1D("th4_pp_N","th4_pp_N",48,0,12); 
  TH1D *th5_pp_N = new TH1D("th5_pp_N","th5_pp_N",48,0,12); 
  TH1D *th6_pp_N = new TH1D("th6_pp_N","th6_pp_N",48,0,12); 
  
  TH1D *th0_pp_S = new TH1D("th0_pp_S","th0_pp_S",48,0,12); 
  TH1D *th1_pp_S = new TH1D("th1_pp_S","th1_pp_S",48,0,12); 
  TH1D *th2_pp_S = new TH1D("th2_pp_S","th2_pp_S",48,0,12); 
  TH1D *th3_pp_S = new TH1D("th3_pp_S","th3_pp_S",48,0,12);  
  TH1D *th4_pp_S = new TH1D("th4_pp_S","th4_pp_S",48,0,12);  
  TH1D *th5_pp_S = new TH1D("th5_pp_S","th5_pp_S",48,0,12); 
  TH1D *th6_pp_S = new TH1D("th6_pp_S","th6_pp_S",48,0,12);  
  //////// pp
  TH1D *th0_mm_N = new TH1D("th0_mm_N","th0_mm_N",48,0,12); 
  TH1D *th1_mm_N = new TH1D("th1_mm_N","th1_mm_N",48,0,12); 
  TH1D *th2_mm_N = new TH1D("th2_mm_N","th2_mm_N",48,0,12);  
  TH1D *th3_mm_N = new TH1D("th3_mm_N","th3_mm_N",48,0,12);  
  TH1D *th4_mm_N = new TH1D("th4_mm_N","th4_mm_N",48,0,12); 
  TH1D *th5_mm_N = new TH1D("th5_mm_N","th5_mm_N",48,0,12);  
  TH1D *th6_mm_N = new TH1D("th6_mm_N","th6_mm_N",48,0,12);  
  
  TH1D *th0_mm_S = new TH1D("th0_mm_S","th0_mm_S",48,0,12);  
  TH1D *th1_mm_S = new TH1D("th1_mm_S","th1_mm_S",48,0,12);  
  TH1D *th2_mm_S = new TH1D("th2_mm_S","th2_mm_S",48,0,12);  
  TH1D *th3_mm_S = new TH1D("th3_mm_S","th3_mm_S",48,0,12);  
  TH1D *th4_mm_S = new TH1D("th4_mm_S","th4_mm_S",48,0,12);  
  TH1D *th5_mm_S = new TH1D("th5_mm_S","th5_mm_S",48,0,12);  
  TH1D *th6_mm_S = new TH1D("th6_mm_S","th6_mm_S",48,0,12);  

    
  // for pT int in N and S 
  if(north_arm)
    {
      th0_N = (TH1D *)t0->ProjectionX("ul_1",1,48,1,60);
      th1_N = (TH1D *)t1->ProjectionX("ul_2",1,48,1,60);
      th2_N = (TH1D *)t2->ProjectionX("ul_3",1,48,1,60);
      th3_N = (TH1D *)t3->ProjectionX("ul_4",1,48,1,60);
      th4_N = (TH1D *)t4->ProjectionX("ul_5",1,48,1,60);
      th5_N = (TH1D *)t5->ProjectionX("ul_6",1,48,1,60);
      th6_N = (TH1D *)t6->ProjectionX("ul_7",1,48,1,60);
    
      th0_pp_N = (TH1D *)t0_pp->ProjectionX("pp_1",1,48,1,60);
      th1_pp_N = (TH1D *)t1_pp->ProjectionX("pp_2",1,48,1,60);
      th2_pp_N = (TH1D *)t2_pp->ProjectionX("pp_3",1,48,1,60);
      th3_pp_N = (TH1D *)t3_pp->ProjectionX("pp_4",1,48,1,60);
      th4_pp_N = (TH1D *)t4_pp->ProjectionX("pp_5",1,48,1,60);
      th5_pp_N = (TH1D *)t5_pp->ProjectionX("pp_6",1,48,1,60);
      th6_pp_N = (TH1D *)t6_pp->ProjectionX("pp_7",1,48,1,60);
    
      th0_mm_N = (TH1D *)t0_mm->ProjectionX("mm_1",1,48,1,60);
      th1_mm_N = (TH1D *)t1_mm->ProjectionX("mm_2",1,48,1,60);
      th2_mm_N = (TH1D *)t2_mm->ProjectionX("mm_3",1,48,1,60);
      th3_mm_N = (TH1D *)t3_mm->ProjectionX("mm_4",1,48,1,60);
      th4_mm_N = (TH1D *)t4_mm->ProjectionX("mm_5",1,48,1,60);
      th5_mm_N = (TH1D *)t5_mm->ProjectionX("mm_6",1,48,1,60);
      th6_mm_N = (TH1D *)t6_mm->ProjectionX("mm_7",1,48,1,60);
    
       
    }
  else
    {
      th0_S = (TH1D *)t0->ProjectionX("ul_1",1,48,1,60);
      th1_S = (TH1D *)t1->ProjectionX("ul_2",1,48,1,60);
      th2_S = (TH1D *)t2->ProjectionX("ul_3",1,48,1,60);
      th3_S = (TH1D *)t3->ProjectionX("ul_4",1,48,1,60);
      th4_S = (TH1D *)t4->ProjectionX("ul_5",1,48,1,60);
      th5_S = (TH1D *)t5->ProjectionX("ul_6",1,48,1,60);
      th6_S = (TH1D *)t6->ProjectionX("ul_7",1,48,1,60);
   
      th0_pp_S = (TH1D *)t0_pp->ProjectionX("pp_1",1,48,1,60);
      th1_pp_S = (TH1D *)t1_pp->ProjectionX("pp_2",1,48,1,60);
      th2_pp_S = (TH1D *)t2_pp->ProjectionX("pp_3",1,48,1,60);
      th3_pp_S = (TH1D *)t3_pp->ProjectionX("pp_4",1,48,1,60);
      th4_pp_S = (TH1D *)t4_pp->ProjectionX("pp_5",1,48,1,60);
      th5_pp_S = (TH1D *)t5_pp->ProjectionX("pp_6",1,48,1,60);
      th6_pp_S = (TH1D *)t6_pp->ProjectionX("pp_7",1,48,1,60);
   
      th0_mm_S = (TH1D *)t0_mm->ProjectionX("mm_1",1,48,1,60);
      th1_mm_S = (TH1D *)t1_mm->ProjectionX("mm_2",1,48,1,60);
      th2_mm_S = (TH1D *)t2_mm->ProjectionX("mm_3",1,48,1,60);
      th3_mm_S = (TH1D *)t3_mm->ProjectionX("mm_4",1,48,1,60);
      th4_mm_S = (TH1D *)t4_mm->ProjectionX("mm_5",1,48,1,60);
      th5_mm_S = (TH1D *)t5_mm->ProjectionX("mm_6",1,48,1,60);
      th6_mm_S = (TH1D *)t6_mm->ProjectionX("mm_7",1,48,1,60);
    
     
    }
  
    char const *outfile;
    if(north_arm && !embed)
      outfile = "./jpsi_root_files/Run15pp_N_NOEMBEDC_int.root";  
    if(!north_arm && !embed)
      outfile = "./jpsi_root_files/Run15pp_S_NOEMBEDC_int.root";  
    if(north_arm && embed)
      outfile = "sanghoon_simulations/fvtx_jpsi_newlib_trig_FID188N_notrigZ_N.root";
      // outfile = "sanghoon_simulations/fvtx_psi2s_newlib_trig_FID188N_notrigZ_N.root";
      // outfile = "sanghoon_simulations/mass_fvtx_jpsi_newlib_acceff_FID188N_nochi2fvtx_trigZ_N.root";
      //  outfile = "sanghoon_simulations/mass_fvtx_psi2s_newlib_acceff_num_NOFID_nochi2_newZ_N.root";
    if(!north_arm && embed)
      outfile = "sanghoon_simulations/fvtx_jpsi_newlib_trig_FID188N_notrigZ_S.root";
      //  outfile = "sanghoon_simulations/fvtx_psi2s_newlib_trig_FID188N_notrigZ_S.root";
      //  outfile = "sanghoon_simulations/mass_fvtx_jpsi_newlib_acceff_FID188N_nochi2fvtx_trigZ_S.root";
      //  outfile = "sanghoon_simulations/mass_fvtx_psi2s_newlib_acceff_num_NOFID_nochi2_newZ_nS.root";
    
  TFile *h = new TFile(outfile, "RECREATE");
  h->cd();

  if(north_arm)
    {
      th0_N->Write();
      th1_N->Write();
      th2_N->Write();
      th3_N->Write();
      th4_N->Write();
      th5_N->Write();
      th6_N->Write();

      th0_pp_N->Write();
      th1_pp_N->Write();
      th2_pp_N->Write();
      th3_pp_N->Write();
      th4_pp_N->Write();
      th5_pp_N->Write();
      th6_pp_N->Write();

      th0_mm_N->Write();
      th1_mm_N->Write();
      th2_mm_N->Write();
      th3_mm_N->Write();
      th4_mm_N->Write();
      th5_mm_N->Write();
      th6_mm_N->Write();
    }
  else
    {
      th0_S->Write();
      th1_S->Write();
      th2_S->Write();
      th3_S->Write();
      th4_S->Write();
      th5_S->Write();
      th6_S->Write();

      th0_pp_S->Write();
      th1_pp_S->Write();
      th2_pp_S->Write();
      th3_pp_S->Write();
      th4_pp_S->Write();
      th5_pp_S->Write();
      th6_pp_S->Write();

      th0_mm_S->Write();
      th1_mm_S->Write();
      th2_mm_S->Write();
      th3_mm_S->Write();
      th4_mm_S->Write();
      th5_mm_S->Write();
      th6_mm_S->Write();
    }
  

    h->Close();
}

























