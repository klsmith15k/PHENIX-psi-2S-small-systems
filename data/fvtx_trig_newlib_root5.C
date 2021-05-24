
#include "TEfficiency.h"
#include<TTreeReader.h>
#include<TTreeReaderValue.h>
#include<TTreeReaderArray.h>
#include <TH2D.h>
#include <TFile.h>
#include<TSystem.h>
#include<TROOT.h>
#include<TLorentzVector.h>
#include<TH1D.h>
#include<TGraphAsymmErrors.h>
#include<sstream>
#include<iostream>
#include<string>

using namespace std;

void fvtx_trig_newlib_root5(char const *infile, char const *outfile)  // for running on CONDOR
//void fvtx_trig_newlib_root5() // for running locally over one fiel to test before submitting to CONDOR
 {

   // only for testing code locally before submitting to condor
   // char const *infile = "/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/jobfiles_PPG188/Run15pp_NOFVTX_227.txt";   /// it says NOFVTX but just ignore
   // const char *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/pDST_105intersect108_mixed/total.root";
   //const char *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/batch_pro108_newlib/Run15pp_64.root";

  ifstream FILES (infile);
  if(!FILES)
    {
      cout<<"The computer can't open the files for some reason."<<endl;
    }

  char filepath [1000];
  FILES.getline(filepath,1000);
  cout<< "filepath is: " <<filepath<<endl; 
  FILES.close();
  
// if using CONDOR (opens text file)
  TFile *f = new TFile(filepath);

 // for using a root file locally 
 // TFile *f = TFile::Open(infile);

  /////////////////////
  TTree *T = (TTree*) f->Get("T");
  /////////////////////

  
  TH2D *mass_fvtx_pt_N_ul = new TH2D("mass_fvtx_pt_N_ul","mass_fvtx_pt_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_N_ul->Sumw2(); mass_fvtx_pt_N_ul->SetLineWidth(2);
  mass_fvtx_pt_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_N_pp = new TH2D("mass_fvtx_pt_N_pp","mass_fvtx_pt_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_N_pp->Sumw2(); mass_fvtx_pt_N_pp->SetLineWidth(2);
  mass_fvtx_pt_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_N_mm = new TH2D("mass_fvtx_pt_N_mm","mass_fvtx_pt_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_N_mm->Sumw2(); mass_fvtx_pt_N_mm->SetLineWidth(2);
  mass_fvtx_pt_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_N_bg = new TH2D("mass_fvtx_pt_N_bg","mass_fvtx_pt_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_N_bg->Sumw2(); mass_fvtx_pt_N_bg->SetLineWidth(2);
  mass_fvtx_pt_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_N_mix = new TH2D("mass_fvtx_pt_N_mix","mass_fvtx_pt_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_N_mix->Sumw2(); mass_fvtx_pt_N_mix->SetLineWidth(2);
  mass_fvtx_pt_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_N_ul->SetLineWidth(2); mass_fvtx_pt_N_ul->SetMarkerStyle(20); mass_fvtx_pt_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_N_bg->SetLineColor(kRed);  mass_fvtx_pt_N_bg->SetLineWidth(2);
  mass_fvtx_pt_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_N_mix->SetLineWidth(2);
  mass_fvtx_pt_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_S_ul = new TH2D("mass_fvtx_pt_S_ul","mass_fvtx_pt_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_S_ul->Sumw2(); mass_fvtx_pt_S_ul->SetLineWidth(2);
  mass_fvtx_pt_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_S_pp = new TH2D("mass_fvtx_pt_S_pp","mass_fvtx_pt_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_S_pp->Sumw2(); mass_fvtx_pt_S_pp->SetLineWidth(2);
  mass_fvtx_pt_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_S_mm = new TH2D("mass_fvtx_pt_S_mm","mass_fvtx_pt_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_S_mm->Sumw2(); mass_fvtx_pt_S_mm->SetLineWidth(2);
  mass_fvtx_pt_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_S_bg = new TH2D("mass_fvtx_pt_S_bg","mass_fvtx_pt_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_S_bg->Sumw2(); mass_fvtx_pt_S_bg->SetLineWidth(2);
  mass_fvtx_pt_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_S_mix = new TH2D("mass_fvtx_pt_S_mix","mass_fvtx_pt_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_S_mix->Sumw2(); mass_fvtx_pt_S_mix->SetLineWidth(2);
  mass_fvtx_pt_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
   mass_fvtx_pt_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_S_ul->SetLineWidth(2); mass_fvtx_pt_S_ul->SetMarkerStyle(20); mass_fvtx_pt_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_S_bg->SetLineColor(kRed);  mass_fvtx_pt_S_bg->SetLineWidth(2);
  mass_fvtx_pt_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_S_mix->SetLineWidth(2);
  mass_fvtx_pt_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_S_pp->GetXaxis()->SetRangeUser(2,5);


 TH2D *mass_fvtx_pt_020_N_ul = new TH2D("mass_fvtx_pt_020_N_ul","mass_fvtx_pt_020_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_N_ul->Sumw2(); mass_fvtx_pt_020_N_ul->SetLineWidth(2);
  mass_fvtx_pt_020_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_N_pp = new TH2D("mass_fvtx_pt_020_N_pp","mass_fvtx_pt_020_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_N_pp->Sumw2(); mass_fvtx_pt_020_N_pp->SetLineWidth(2);
  mass_fvtx_pt_020_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_N_mm = new TH2D("mass_fvtx_pt_020_N_mm","mass_fvtx_pt_020_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_N_mm->Sumw2(); mass_fvtx_pt_020_N_mm->SetLineWidth(2);
  mass_fvtx_pt_020_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_N_bg = new TH2D("mass_fvtx_pt_020_N_bg","mass_fvtx_pt_020_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_N_bg->Sumw2(); mass_fvtx_pt_020_N_bg->SetLineWidth(2);
  mass_fvtx_pt_020_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_N_mix = new TH2D("mass_fvtx_pt_020_N_mix","mass_fvtx_pt_020_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_N_mix->Sumw2(); mass_fvtx_pt_020_N_mix->SetLineWidth(2);
  mass_fvtx_pt_020_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_020_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_020_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_020_N_ul->SetLineWidth(2); mass_fvtx_pt_020_N_ul->SetMarkerStyle(20); mass_fvtx_pt_020_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_020_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_020_N_bg->SetLineColor(kRed);  mass_fvtx_pt_020_N_bg->SetLineWidth(2);
  mass_fvtx_pt_020_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_020_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_020_N_mix->SetLineWidth(2);
  mass_fvtx_pt_020_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_020_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_020_S_ul = new TH2D("mass_fvtx_pt_020_S_ul","mass_fvtx_pt_020_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_S_ul->Sumw2(); mass_fvtx_pt_020_S_ul->SetLineWidth(2);
  mass_fvtx_pt_020_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_S_pp = new TH2D("mass_fvtx_pt_020_S_pp","mass_fvtx_pt_020_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_S_pp->Sumw2(); mass_fvtx_pt_020_S_pp->SetLineWidth(2);
  mass_fvtx_pt_020_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_S_mm = new TH2D("mass_fvtx_pt_020_S_mm","mass_fvtx_pt_020_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_S_mm->Sumw2(); mass_fvtx_pt_020_S_mm->SetLineWidth(2);
  mass_fvtx_pt_020_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_S_bg = new TH2D("mass_fvtx_pt_020_S_bg","mass_fvtx_pt_020_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_S_bg->Sumw2(); mass_fvtx_pt_020_S_bg->SetLineWidth(2);
  mass_fvtx_pt_020_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_020_S_mix = new TH2D("mass_fvtx_pt_020_S_mix","mass_fvtx_pt_020_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_020_S_mix->Sumw2(); mass_fvtx_pt_020_S_mix->SetLineWidth(2);
  mass_fvtx_pt_020_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_020_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_020_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_020_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_020_S_ul->SetLineWidth(2); mass_fvtx_pt_020_S_ul->SetMarkerStyle(20); mass_fvtx_pt_020_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_020_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_020_S_bg->SetLineColor(kRed);  mass_fvtx_pt_020_S_bg->SetLineWidth(2);
  mass_fvtx_pt_020_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_020_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_020_S_mix->SetLineWidth(2);
  mass_fvtx_pt_020_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_020_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_020_S_pp->GetXaxis()->SetRangeUser(2,5);
  
  //////////

 TH2D *mass_fvtx_pt_2040_N_ul = new TH2D("mass_fvtx_pt_2040_N_ul","mass_fvtx_pt_2040_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_N_ul->Sumw2(); mass_fvtx_pt_2040_N_ul->SetLineWidth(2);
  mass_fvtx_pt_2040_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_N_pp = new TH2D("mass_fvtx_pt_2040_N_pp","mass_fvtx_pt_2040_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_N_pp->Sumw2(); mass_fvtx_pt_2040_N_pp->SetLineWidth(2);
  mass_fvtx_pt_2040_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_N_mm = new TH2D("mass_fvtx_pt_2040_N_mm","mass_fvtx_pt_2040_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_N_mm->Sumw2(); mass_fvtx_pt_2040_N_mm->SetLineWidth(2);
  mass_fvtx_pt_2040_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_N_bg = new TH2D("mass_fvtx_pt_2040_N_bg","mass_fvtx_pt_2040_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_N_bg->Sumw2(); mass_fvtx_pt_2040_N_bg->SetLineWidth(2);
  mass_fvtx_pt_2040_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_N_mix = new TH2D("mass_fvtx_pt_2040_N_mix","mass_fvtx_pt_2040_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_N_mix->Sumw2(); mass_fvtx_pt_2040_N_mix->SetLineWidth(2);
  mass_fvtx_pt_2040_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
 mass_fvtx_pt_2040_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_2040_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_2040_N_ul->SetLineWidth(2); mass_fvtx_pt_2040_N_ul->SetMarkerStyle(20); mass_fvtx_pt_2040_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_2040_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_2040_N_bg->SetLineColor(kRed);  mass_fvtx_pt_2040_N_bg->SetLineWidth(2);
  mass_fvtx_pt_2040_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_2040_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_2040_N_mix->SetLineWidth(2);
  mass_fvtx_pt_2040_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_2040_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_2040_S_ul = new TH2D("mass_fvtx_pt_2040_S_ul","mass_fvtx_pt_2040_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_S_ul->Sumw2(); mass_fvtx_pt_2040_S_ul->SetLineWidth(2);
  mass_fvtx_pt_2040_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_S_pp = new TH2D("mass_fvtx_pt_2040_S_pp","mass_fvtx_pt_2040_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_S_pp->Sumw2(); mass_fvtx_pt_2040_S_pp->SetLineWidth(2);
  mass_fvtx_pt_2040_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_S_mm = new TH2D("mass_fvtx_pt_2040_S_mm","mass_fvtx_pt_2040_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_S_mm->Sumw2(); mass_fvtx_pt_2040_S_mm->SetLineWidth(2);
  mass_fvtx_pt_2040_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_S_bg = new TH2D("mass_fvtx_pt_2040_S_bg","mass_fvtx_pt_2040_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_S_bg->Sumw2(); mass_fvtx_pt_2040_S_bg->SetLineWidth(2);
  mass_fvtx_pt_2040_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2040_S_mix = new TH2D("mass_fvtx_pt_2040_S_mix","mass_fvtx_pt_2040_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2040_S_mix->Sumw2(); mass_fvtx_pt_2040_S_mix->SetLineWidth(2);
  mass_fvtx_pt_2040_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2040_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
   mass_fvtx_pt_2040_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_2040_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_2040_S_ul->SetLineWidth(2); mass_fvtx_pt_2040_S_ul->SetMarkerStyle(20); mass_fvtx_pt_2040_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_2040_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_2040_S_bg->SetLineColor(kRed);  mass_fvtx_pt_2040_S_bg->SetLineWidth(2);
  mass_fvtx_pt_2040_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_2040_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_2040_S_mix->SetLineWidth(2);
  mass_fvtx_pt_2040_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2040_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_2040_S_pp->GetXaxis()->SetRangeUser(2,5);
  
  //////////////

 TH2D *mass_fvtx_pt_4060_N_ul = new TH2D("mass_fvtx_pt_4060_N_ul","mass_fvtx_pt_4060_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_N_ul->Sumw2(); mass_fvtx_pt_4060_N_ul->SetLineWidth(2);
  mass_fvtx_pt_4060_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_N_pp = new TH2D("mass_fvtx_pt_4060_N_pp","mass_fvtx_pt_4060_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_N_pp->Sumw2(); mass_fvtx_pt_4060_N_pp->SetLineWidth(2);
  mass_fvtx_pt_4060_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_N_mm = new TH2D("mass_fvtx_pt_4060_N_mm","mass_fvtx_pt_4060_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_N_mm->Sumw2(); mass_fvtx_pt_4060_N_mm->SetLineWidth(2);
  mass_fvtx_pt_4060_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_N_bg = new TH2D("mass_fvtx_pt_4060_N_bg","mass_fvtx_pt_4060_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_N_bg->Sumw2(); mass_fvtx_pt_4060_N_bg->SetLineWidth(2);
  mass_fvtx_pt_4060_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_N_mix = new TH2D("mass_fvtx_pt_4060_N_mix","mass_fvtx_pt_4060_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_N_mix->Sumw2(); mass_fvtx_pt_4060_N_mix->SetLineWidth(2);
  mass_fvtx_pt_4060_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_4060_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_4060_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_4060_N_ul->SetLineWidth(2); mass_fvtx_pt_4060_N_ul->SetMarkerStyle(20); mass_fvtx_pt_4060_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_4060_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_4060_N_bg->SetLineColor(kRed);  mass_fvtx_pt_4060_N_bg->SetLineWidth(2);
  mass_fvtx_pt_4060_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_4060_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_4060_N_mix->SetLineWidth(2);
  mass_fvtx_pt_4060_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_4060_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_4060_S_ul = new TH2D("mass_fvtx_pt_4060_S_ul","mass_fvtx_pt_4060_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_S_ul->Sumw2(); mass_fvtx_pt_4060_S_ul->SetLineWidth(2);
  mass_fvtx_pt_4060_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_S_pp = new TH2D("mass_fvtx_pt_4060_S_pp","mass_fvtx_pt_4060_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_S_pp->Sumw2(); mass_fvtx_pt_4060_S_pp->SetLineWidth(2);
  mass_fvtx_pt_4060_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_S_mm = new TH2D("mass_fvtx_pt_4060_S_mm","mass_fvtx_pt_4060_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_S_mm->Sumw2(); mass_fvtx_pt_4060_S_mm->SetLineWidth(2);
  mass_fvtx_pt_4060_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_S_bg = new TH2D("mass_fvtx_pt_4060_S_bg","mass_fvtx_pt_4060_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_S_bg->Sumw2(); mass_fvtx_pt_4060_S_bg->SetLineWidth(2);
  mass_fvtx_pt_4060_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4060_S_mix = new TH2D("mass_fvtx_pt_4060_S_mix","mass_fvtx_pt_4060_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4060_S_mix->Sumw2(); mass_fvtx_pt_4060_S_mix->SetLineWidth(2);
  mass_fvtx_pt_4060_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4060_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
 mass_fvtx_pt_4060_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_4060_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_4060_S_ul->SetLineWidth(2); mass_fvtx_pt_4060_S_ul->SetMarkerStyle(20); mass_fvtx_pt_4060_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_4060_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_4060_S_bg->SetLineColor(kRed);  mass_fvtx_pt_4060_S_bg->SetLineWidth(2);
  mass_fvtx_pt_4060_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_4060_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_4060_S_mix->SetLineWidth(2);
  mass_fvtx_pt_4060_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4060_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_4060_S_pp->GetXaxis()->SetRangeUser(2,5);
  
  //////////////

 TH2D *mass_fvtx_pt_6084_N_ul = new TH2D("mass_fvtx_pt_6084_N_ul","mass_fvtx_pt_6084_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_N_ul->Sumw2(); mass_fvtx_pt_6084_N_ul->SetLineWidth(2);
  mass_fvtx_pt_6084_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_N_pp = new TH2D("mass_fvtx_pt_6084_N_pp","mass_fvtx_pt_6084_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_N_pp->Sumw2(); mass_fvtx_pt_6084_N_pp->SetLineWidth(2);
  mass_fvtx_pt_6084_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_N_mm = new TH2D("mass_fvtx_pt_6084_N_mm","mass_fvtx_pt_6084_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_N_mm->Sumw2(); mass_fvtx_pt_6084_N_mm->SetLineWidth(2);
  mass_fvtx_pt_6084_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_N_bg = new TH2D("mass_fvtx_pt_6084_N_bg","mass_fvtx_pt_6084_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_N_bg->Sumw2(); mass_fvtx_pt_6084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_6084_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_N_mix = new TH2D("mass_fvtx_pt_6084_N_mix","mass_fvtx_pt_6084_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_N_mix->Sumw2(); mass_fvtx_pt_6084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_6084_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
   mass_fvtx_pt_6084_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_6084_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_6084_N_ul->SetLineWidth(2); mass_fvtx_pt_6084_N_ul->SetMarkerStyle(20); mass_fvtx_pt_6084_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_6084_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_6084_N_bg->SetLineColor(kRed);  mass_fvtx_pt_6084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_6084_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_6084_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_6084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_6084_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_6084_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_6084_S_ul = new TH2D("mass_fvtx_pt_6084_S_ul","mass_fvtx_pt_6084_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_S_ul->Sumw2(); mass_fvtx_pt_6084_S_ul->SetLineWidth(2);
  mass_fvtx_pt_6084_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_S_pp = new TH2D("mass_fvtx_pt_6084_S_pp","mass_fvtx_pt_6084_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_S_pp->Sumw2(); mass_fvtx_pt_6084_S_pp->SetLineWidth(2);
  mass_fvtx_pt_6084_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_S_mm = new TH2D("mass_fvtx_pt_6084_S_mm","mass_fvtx_pt_6084_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_S_mm->Sumw2(); mass_fvtx_pt_6084_S_mm->SetLineWidth(2);
  mass_fvtx_pt_6084_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_S_bg = new TH2D("mass_fvtx_pt_6084_S_bg","mass_fvtx_pt_6084_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_S_bg->Sumw2(); mass_fvtx_pt_6084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_6084_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_6084_S_mix = new TH2D("mass_fvtx_pt_6084_S_mix","mass_fvtx_pt_6084_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_6084_S_mix->Sumw2(); mass_fvtx_pt_6084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_6084_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_6084_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_6084_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_6084_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_6084_S_ul->SetLineWidth(2); mass_fvtx_pt_6084_S_ul->SetMarkerStyle(20); mass_fvtx_pt_6084_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_6084_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_6084_S_bg->SetLineColor(kRed);  mass_fvtx_pt_6084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_6084_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_6084_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_6084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_6084_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_6084_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_6084_S_pp->GetXaxis()->SetRangeUser(2,5);
  
  //////////////

  TH2D *mass_good_N = new TH2D("mass_good_N","mass_good_N", 300, 0, 15, 120, 0, 12); mass_good_N->Sumw2(); mass_good_N->SetLineWidth(2);
  mass_good_N->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_good_N->GetYaxis()->SetTitle("p_{T} (GeV/c)");
 
  TH2D *mass_bad_N = new TH2D("mass_bad_N","mass_bad_N", 300, 0, 15, 120, 0, 12); mass_bad_N->Sumw2(); mass_bad_N->SetLineWidth(2);
  mass_bad_N->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_bad_N->GetYaxis()->SetTitle("p_{T} (GeV/c)");

 TH2D *mass_good_S = new TH2D("mass_good_S","mass_good_S", 300, 0, 15, 120, 0, 12); mass_good_S->Sumw2(); mass_good_S->SetLineWidth(2);
  mass_good_S->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_good_S->GetYaxis()->SetTitle("p_{T} (GeV/c)");
 
  TH2D *mass_bad_S = new TH2D("mass_bad_S","mass_bad_S", 300, 0, 15, 120, 0, 12); mass_bad_S->Sumw2(); mass_bad_S->SetLineWidth(2);
  mass_bad_S->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_bad_S->GetYaxis()->SetTitle("p_{T} (GeV/c)");

  ///////////////////////////////////////////////////////////////

 //////////////////////////////////////////////
  TTreeReader myReader("T1", f);
  /////////////////////////////////////////////////
  TTreeReaderValue<int> run_number(myReader, "RunNumber");

  unsigned int run_num = 0;
  int badrun_S = 0;
  int badrun_N = 0;
  int counter = 0;

 std::string  filename[2] = {"/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/matt_centrality_macros/sanghoon_badruns_list_S_v2.txt","/phenix/hhj/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/matt_centrality_macros/sanghoon_badruns_list_N_v2.txt"};

 while( myReader.Next())
    {
      run_num = *run_number;
    } // while

 ifstream run_files(filename[0].c_str() );
 if(run_files)
   {
     do 
       {
  	 run_files >> badrun_S;
	 //	 cout << "badrun S: " << badrun_S << endl;
	 //	 cout << "run number: " << run_num << endl;

 	 counter++;

	 if(badrun_S == run_num)
	   break;  // continue here does not work

       }while(counter< 171);

   } // end if
 else
   cout << "test read file does not exist" << endl;
 
 counter = 0;

ifstream run_files2(filename[1].c_str() );
 if(run_files2)
   {
     do 
       {
 	 run_files2 >> badrun_N;
	 //	 cout << "badrun N: " << badrun_N << endl;
	 // cout << "run number: " << run_num << endl;

 	 counter++;

	 if(badrun_N == run_num)
	   break;  // continue here does not work

       }while(counter< 250);

   } // end if
 else
   cout << "test read file does not exist" << endl;

 //cout << "out of ifstream" << endl;
 if(badrun_S == run_num)
   cout << "South arm bad run" << endl;
 if(badrun_N == run_num)
   cout << "North arm bad run" << endl;
 
 /////////////////////////////////////////////////////////////////////////////////////
 
 string mass_cut = "mass_fvtxmutr > 1.5";  // single tracks are stored in fvtxmutr variable
  string evt_cut_N = "(Evt_bbcZ<30) && (Evt_bbcZ>-30)";  // already applied in dimuon container make_dimu.C
  string evt_cut_S = "(Evt_bbcZ>-30) && (Evt_bbcZ<30)";  // already applied in dimuon container make_dimu.C
  string trig_cut_N = "lvl1_trigscaled&0x00100000";
  string trig_cut_S = "lvl1_trigscaled&0x00200000";
  string nidhits_cut = "Tr0_idhits>14 && Tr1_idhits>14";
  string lastgap_cut = "(Tr0_lastgap>=3 && Tr1_lastgap>=3)";

  // cuts from all runs
  // string DG0_cut_N = "Tr0_DG0<(8.49898+100.736/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DG0<(8.49898+100.736/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  // string DG0_cut_S = "Tr0_DG0<(8.79546+376.293/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DG0<(8.79546+376.293/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  // string DDG0_cut_N = "Tr0_DDG0<(4.51376+137.931/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DDG0<(4.51376+137.931/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  // string DDG0_cut_S = "Tr0_DDG0<(4.71542+119.735/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DDG0<(4.71542+119.735/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";

  // cuts from good runs only
  string DG0_cut_N = "Tr0_DG0<(8.14581+106.858/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DG0<(8.14581+106.858/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string DG0_cut_S = "Tr0_DG0<(8.93709+366.856/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DG0<(8.93709+366.856/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string DDG0_cut_N = "Tr0_DDG0<(4.3602+142.651/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DDG0<(4.3602+142.651/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string DDG0_cut_S = "Tr0_DDG0<(4.66658+120.627/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_DDG0<(4.66658+120.627/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  
  string tkchi2_cut = "Tr0_trchi2<23 && Tr1_trchi2<23";
  string pz_fvtxmutr_cut_N = "abs(Tr0_pz_fvtxmutr)>2 && abs(Tr1_pz_fvtxmutr)>2 && Tr0_pz_fvtxmutr>0 && Tr1_pz_fvtxmutr>0 ";
  string pz_fvtxmutr_cut_S = "abs(Tr0_pz_fvtxmutr)>2 && abs(Tr1_pz_fvtxmutr)>2 &&  Tr0_pz_fvtxmutr<0 && Tr1_pz_fvtxmutr<0 ";
  string asym_cut = "abs( (sqrt(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr ) - sqrt(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr ))/ (sqrt(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr ) + sqrt(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr )) ) < 1 ";
  string octant_cut ="abs(  ( ((int((atan2(Tr0_xst1,Tr0_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) - ( ((int((atan2(Tr1_xst1,Tr1_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) ) >0";
 
  
  // PPG188 fiducial cuts
  string fid_cut_N = " !(  atan2(Tr0_yst3,Tr0_xst3)>0.35 && atan2(Tr0_yst3,Tr0_xst3)<0.77  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>0.77 && atan2(Tr0_yst3,Tr0_xst3)<1.15 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>160 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>-2.75 && atan2(Tr0_yst3,Tr0_xst3)<-1.92 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>237 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst2,Tr0_xst2)>-2.75 && atan2(Tr0_yst2,Tr0_xst2)<-1.92 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)>140 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)<163  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>-2.75 && atan2(Tr0_yst1,Tr0_xst1)<-2 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>76 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<80  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>0.52 && atan2(Tr0_yst1,Tr0_xst1)<1.11 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>60 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<75 )  && !(  atan2(Tr1_yst3,Tr1_xst3)>0.35 && atan2(Tr1_yst3,Tr1_xst3)<0.77  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>0.77 && atan2(Tr1_yst3,Tr1_xst3)<1.15 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>160 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>-2.75 && atan2(Tr1_yst3,Tr1_xst3)<-1.92 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>237 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst2,Tr1_xst2)>-2.75 && atan2(Tr1_yst2,Tr1_xst2)<-1.92 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)>140 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)<163  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>-2.75 && atan2(Tr1_yst1,Tr1_xst1)<-2 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>76 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<80  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>0.52 && atan2(Tr1_yst1,Tr1_xst1)<1.11 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>60 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<75 )";
  string fid_cut_S = "Evt_Cent!=-12";
  // string fid_cut_N = "Evt_Cent!=-12";

//////////////////////////////////////////////////////////////  PPG228 fiducial cuts January 17, 2021
//string fid_cut_N = "!(  atan2(Tr0_yst3,Tr0_xst3)>0.35 && atan2(Tr0_yst3,Tr0_xst3)<0.77  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>0.77 && atan2(Tr0_yst3,Tr0_xst3)<1.15 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>160 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>-2.75 && atan2(Tr0_yst3,Tr0_xst3)<-1.92 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>237 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst2,Tr0_xst2)>-2.75 && atan2(Tr0_yst2,Tr0_xst2)<-1.92 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)>140 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)<163  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>-2.75 && atan2(Tr0_yst1,Tr0_xst1)<-2 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>76 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<80  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>0.52 && atan2(Tr0_yst1,Tr0_xst1)<1.11 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>60 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<75 )  && !(  atan2(Tr1_yst3,Tr1_xst3)>0.35 && atan2(Tr1_yst3,Tr1_xst3)<0.77  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>0.77 && atan2(Tr1_yst3,Tr1_xst3)<1.15 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>160 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>-2.75 && atan2(Tr1_yst3,Tr1_xst3)<-1.92 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>237 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst2,Tr1_xst2)>-2.75 && atan2(Tr1_yst2,Tr1_xst2)<-1.92 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)>140 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)<163  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>-2.75 && atan2(Tr1_yst1,Tr1_xst1)<-2 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>76 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<80  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>0.52 && atan2(Tr1_yst1,Tr1_xst1)<1.11 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>60 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<75 ) && !( atan2(Tr0_yst2,Tr0_xst2)>-1.20 && atan2(Tr0_yst2,Tr0_xst2)<-0.90 )&& !( atan2(Tr1_yst2,Tr1_xst2)>-1.20 && atan2(Tr1_yst2,Tr1_xst2)<-0.90 )";
//string fid_cut_S = "!( atan2(Tr0_yst1,Tr0_xst1)>0 && atan2(Tr0_yst1,Tr0_xst1)<0.4 ) && !( atan2(Tr0_yst2,Tr0_xst2)>0 && atan2(Tr0_yst2,Tr0_xst2)<0.4 ) && !( atan2(Tr0_yst3,Tr0_xst3)>0 && atan2(Tr0_yst3,Tr0_xst3)<0.4 ) && !( atan2(Tr0_yst1,Tr0_xst1)>1.55 && atan2(Tr0_yst1,Tr0_xst1)<1.9 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>1.55 && atan2(Tr0_yst2,Tr0_xst2)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>1.55 && atan2(Tr0_yst3,Tr0_xst3)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>-0.60 && atan2(Tr0_yst3,Tr0_xst3)<-0.40 ) && !( atan2(Tr1_yst1,Tr1_xst1)>0 && atan2(Tr1_yst1,Tr1_xst1)<0.4 ) && !( atan2(Tr1_yst2,Tr1_xst2)>0 && atan2(Tr1_yst2,Tr1_xst2)<0.4 ) && !( atan2(Tr1_yst3,Tr1_xst3)>0 && atan2(Tr1_yst3,Tr1_xst3)<0.4 ) && !( atan2(Tr1_yst1,Tr1_xst1)>1.55 && atan2(Tr1_yst1,Tr1_xst1)<1.9 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>1.55 && atan2(Tr1_yst2,Tr1_xst2)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>1.55 && atan2(Tr1_yst3,Tr1_xst3)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>-0.60 && atan2(Tr1_yst3,Tr1_xst3)<-0.40 ) && !( atan2(Tr0_yst1,Tr0_xst1)>-2.35 && atan2(Tr0_yst1,Tr0_xst1)<-1.95 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>-2.35 && atan2(Tr0_yst2,Tr0_xst2)<-1.95 ) &&  !( atan2(Tr0_yst3,Tr0_xst3)>-2.35 && atan2(Tr0_yst3,Tr0_xst3)<-1.95 ) && !( atan2(Tr1_yst1,Tr1_xst1)>-2.35 && atan2(Tr1_yst1,Tr1_xst1)<-1.95 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>-2.35 && atan2(Tr1_yst2,Tr1_xst2)<-1.95 ) &&  !( atan2(Tr1_yst3,Tr1_xst3)>-2.35 && atan2(Tr1_yst3,Tr1_xst3)<-1.95 )"; 
 ////////////////////////////////////////////////////////////

  string run_cut_S = "runnumber!=-999";   
  string run_cut_N = "runnumber!=-999";

  string muon_arm_cut_N =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_N + "&&" + DDG0_cut_N + "&&" + mass_cut + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_N + "&&" + run_cut_N;
  string muon_arm_cut_S =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_S + "&&" + DDG0_cut_S + "&&" + mass_cut + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_S + "&&" + run_cut_S;

  string chi2_fvtx_cut = "Tr0_chi2_fvtx<10 && Tr0_chi2_fvtx>0 && Tr1_chi2_fvtx<10 && Tr1_chi2_fvtx>0";
  string chi2_fvtxmutr_cut = "Tr0_chi2_fvtxmutr<10 && Tr0_chi2_fvtxmutr>0 && Tr1_chi2_fvtxmutr<10 && Tr1_chi2_fvtxmutr>0";

  // cuts from all runs
//  string dr_cut_N = "Tr0_dr_fvtx<(1.69833+24.7377/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dr_fvtx<(1.69833+24.7377/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
// string dr_cut_S = "Tr0_dr_fvtx<(1.54454+27.5314/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dr_fvtx<(1.54454+27.5314/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
// string dtheta_cut_N = "Tr0_dtheta_fvtx<(0.0482243+0.732457/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dtheta_fvtx<(0.0482243+0.732457/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
// string dtheta_cut_S = "Tr0_dtheta_fvtx<(0.0559102+0.580277/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dtheta_fvtx<(0.0559102+0.580277/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
// string dphi_cut_N = "Tr0_dphi_fvtx<(0.117266+1.53911/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dphi_fvtx<(0.117266+1.53911/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
// string dphi_cut_S = "Tr0_dphi_fvtx<(0.113169+1.6207/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dphi_fvtx<(0.113169+1.6207/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";

// cuts from only good runs
  string dr_cut_N = "Tr0_dr_fvtx<(1.76048+23.9016/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dr_fvtx<(1.76048+23.9016/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string dr_cut_S = "Tr0_dr_fvtx<(1.39017+29.279/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dr_fvtx<(1.39017+29.279/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string dtheta_cut_N = "Tr0_dtheta_fvtx<(0.0499702+0.694368/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dtheta_fvtx<(0.0499702+0.694368/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string dtheta_cut_S = "Tr0_dtheta_fvtx<(0.0555553+0.557225/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dtheta_fvtx<(0.0555553+0.557225/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string dphi_cut_N = "Tr0_dphi_fvtx<(0.11936+1.44372/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dphi_fvtx<(0.11936+1.44372/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";
  string dphi_cut_S = "Tr0_dphi_fvtx<(0.116105+1.51064/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)) && Tr1_dphi_fvtx<(0.116105+1.51064/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))";


  string sng_fvtx_cut = chi2_fvtx_cut + "&&" + chi2_fvtxmutr_cut;
  string rap_cut_S = "rapidity_fvtxmutr<=-1.2 && rapidity_fvtxmutr>=-2.2";
  string rap_cut_N = "rapidity_fvtxmutr>=1.2 && rapidity_fvtxmutr<=2.2";

  string fvtx_trk_vtx = "( (Tr0_nhits_fvtx > 0)  && (Tr1_nhits_fvtx > 0) )";
  //string sng_trk_vtx = "( ( (Tr0_nhits_fvtx > 0)  && (Tr1_nhits_fvtx < 0) ) || ( (Tr1_nhits_fvtx > 0)  && (Tr0_nhits_fvtx < 0) ) )";

  //string cuts_N = "&&" + rap_cut_N  + "&&"  + pz_fvtxmutr_cut_N + "&&"+ muon_arm_cut_N + "&&" + sng_fvtx_cut + "&&" + evt_cut_N + "&&" + dr_cut_N + "&&" + dtheta_cut_N + "&&" + dphi_cut_N + "&&" + fvtx_trk_vtx ;
  //string cuts_S = "&&" + rap_cut_S  + "&&"  + pz_fvtxmutr_cut_S + "&&" + muon_arm_cut_S + "&&" + sng_fvtx_cut + "&&" + evt_cut_S + "&&" + dr_cut_S + "&&" + dtheta_cut_S + "&&" + dphi_cut_S + "&&" + fvtx_trk_vtx;

string cuts_N = "&&" + rap_cut_N  + "&&"  + pz_fvtxmutr_cut_N + "&&"+ muon_arm_cut_N + "&&" + evt_cut_N + "&&" + dr_cut_N + "&&" + dtheta_cut_N + "&&" + dphi_cut_N + "&&" + fvtx_trk_vtx ;
  string cuts_S = "&&" + rap_cut_S  + "&&"  + pz_fvtxmutr_cut_S + "&&" + muon_arm_cut_S + "&&" + evt_cut_S + "&&" + dr_cut_S + "&&" + dtheta_cut_S + "&&" + dphi_cut_S + "&&" + fvtx_trk_vtx;
 
  ////////////////////////////


  int lo_bin =  mass_fvtx_pt_N_pp->GetXaxis()->FindBin(2);
  int hi_bin =  mass_fvtx_pt_N_pp->GetXaxis()->FindBin(5.5);
 

 ///project North arm histograms
  if(badrun_N != run_num)
  //  if( (fvtx_prob0 > 0.05) && (fvtx_prob1 > 0.05) )
   {
  T->Project("mass_fvtx_pt_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1" + cuts_N).c_str() );
    
  T->Project("mass_fvtx_pt_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1" + cuts_N).c_str() ); 
    
  T->Project("mass_fvtx_pt_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1" + cuts_N).c_str() );
  
  T->Project("mass_fvtx_pt_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_N).c_str() );

/////////////////////////////////////
  ///project North arm histograms, 0-20% centrality
  
  T->Project("mass_fvtx_pt_020_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N).c_str() );
    
  T->Project("mass_fvtx_pt_020_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N).c_str() ); 
  
  T->Project("mass_fvtx_pt_020_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N).c_str() );
 
  T->Project("mass_fvtx_pt_020_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_N).c_str() );

 /////////////////////////////////////
  ///project North arm histograms, 20-40% centrality
  
  T->Project("mass_fvtx_pt_2040_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N).c_str() );
    
  T->Project("mass_fvtx_pt_2040_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N).c_str() ); 
  
  T->Project("mass_fvtx_pt_2040_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N).c_str() );
  
  T->Project("mass_fvtx_pt_2040_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_N).c_str() );

 /////////////////////////////////////
 ///project North arm histograms, 40-60% centrality
 
  T->Project("mass_fvtx_pt_4060_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N).c_str() );
    
  T->Project("mass_fvtx_pt_4060_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N).c_str() ); 
  
  T->Project("mass_fvtx_pt_4060_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N).c_str() );
 
  T->Project("mass_fvtx_pt_4060_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_N).c_str() );

 /////////////////////////////////////
 ///project North arm histograms, 60-84% centrality
 
  T->Project("mass_fvtx_pt_6084_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N).c_str() );
  
  T->Project("mass_fvtx_pt_6084_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N).c_str() ); 
  
  T->Project("mass_fvtx_pt_6084_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N).c_str() );
  
  T->Project("mass_fvtx_pt_6084_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_N).c_str() );

   }

  if(run_num!=badrun_S)
   //if( (fvtx_prob0 > 0.05) && (fvtx_prob1 > 0.05) )
    {
   ///project South arm histograms
  T->Project("mass_fvtx_pt_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1" + cuts_S).c_str() );
   
  T->Project("mass_fvtx_pt_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1" + cuts_S).c_str() ); 
   
  T->Project("mass_fvtx_pt_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1" + cuts_S).c_str() );
 
  T->Project("mass_fvtx_pt_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_S).c_str() );

   ///project South arm histograms
  T->Project("mass_fvtx_pt_020_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_020_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S).c_str() ); 
   
  T->Project("mass_fvtx_pt_020_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_020_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_S).c_str() );
 
   ///project South arm histograms
  T->Project("mass_fvtx_pt_2040_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S).c_str() );
 
  T->Project("mass_fvtx_pt_2040_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S).c_str() ); 
    
  T->Project("mass_fvtx_pt_2040_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_2040_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_S).c_str() );
  
  /////////////////////////////////////

   ///project South arm histograms
  T->Project("mass_fvtx_pt_4060_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_4060_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S).c_str() ); 
   
  T->Project("mass_fvtx_pt_4060_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_4060_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_S).c_str() );

   ///project South arm histograms
  T->Project("mass_fvtx_pt_6084_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S).c_str() );
  
  T->Project("mass_fvtx_pt_6084_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S).c_str() ); 
   
  T->Project("mass_fvtx_pt_6084_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S).c_str() );
 
  T->Project("mass_fvtx_pt_6084_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_S).c_str() );
    } // end S
  //////////////////////////////////////

  ///look at best and worst tracks to determine systematic on second gaussian
  T->Project("mass_good_S","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==1&&Tr0_ntrhits>=14&&Tr1_ntrhits>=14" + cuts_S).c_str() );
  T->Project("mass_bad_S","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==1&&Tr0_ntrhits<14&&Tr1_ntrhits<14" + cuts_S).c_str() );
 
 T->Project("mass_good_N","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==1&&Tr0_ntrhits>=14&&Tr1_ntrhits>=14" + cuts_N).c_str() );
  T->Project("mass_bad_N","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==1&&Tr0_ntrhits<14&&Tr1_ntrhits<14" + cuts_N).c_str() );
 
  string muon2_arm_cut_N = lastgap_cut + "&&" + DG0_cut_N + "&&" + DDG0_cut_N + "&&" + tkchi2_cut + "&&" + fid_cut_N; 
  string muon2_arm_cut_S = lastgap_cut + "&&" + DG0_cut_S + "&&" + DDG0_cut_S + "&&" + tkchi2_cut + "&&" + fid_cut_S;

 TH1D *runs = new TH1D("runs","runs",450000.5-422000.5,422000.5,450000.5);  
 T->Project("runs","runnumber","abs(Evt_bbcZ)<30");

  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  
  c1->Divide(2,1);
  c1->cd(1); c1->cd(1)->SetLogy();
  // mass_fvtx_pt_S_mix->ProjectionX()->Draw();
  mass_fvtx_pt_S_ul->ProjectionX()->Draw();
  c1->cd(2);  c1->cd(2)->SetLogy();
  // mass_fvtx_pt_N_mix->ProjectionX()->Draw();
  mass_fvtx_pt_N_ul->ProjectionX()->Draw();
  
  //Plot the variables:
  string nidhits_cut2 = "Tr0_idhits>8 && Tr1_idhits>8";
  string lastgap_cut2 = "(Tr0_lastgap>=3 && Tr1_lastgap>=3)";
  string DG0_cut2 = "Tr0_DG0<30 && Tr1_DG0<30";
  string DDG0_cut2 = "Tr0_DDG0<30 && Tr1_DDG0<30";
  string tkchi2_cut2 = "Tr0_trchi2<23 && Tr1_trchi2<23";
  string octant_cut2 ="abs(  ( ((int((atan2(Tr0_xst1,Tr0_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) - ( ((int((atan2(Tr1_xst1,Tr1_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) ) >0";
  string mass_cut2 = " mass_fvtx >2.8 && mass_fvtx <3.4 ";
  string pz_fvtxmutr_cut2_N = " Tr0_pz_fvtxmutr>0 && Tr1_pz_fvtxmutr>0 ";//NORTH ARM ONLY
  string pz_fvtxmutr_cut2_S = " Tr0_pz_fvtxmutr<0 && Tr1_pz_fvtxmutr<0 ";//SOUTH ARM ONLY
  string ul_cut2 = " charge == 0 && same_event==1";  
  
  string ul_N_data = evt_cut_N + "&&" + nidhits_cut2 + "&&" + lastgap_cut2 + "&&" + DG0_cut2 + "&&" + DDG0_cut2 + "&&" + tkchi2_cut2 + "&&" + octant_cut2 + "&&" + mass_cut2 + "&&" + pz_fvtxmutr_cut2_N + "&&" + ul_cut2;  
  string ul_S_data = evt_cut_S + "&&" + nidhits_cut2 + "&&" + lastgap_cut2 + "&&" + DG0_cut2 + "&&" + DDG0_cut2 + "&&" + tkchi2_cut2 + "&&" + octant_cut2 + "&&" + mass_cut2 + "&&" + pz_fvtxmutr_cut2_S + "&&" + ul_cut2;  
  //lets start with dr
  TH2D *dr_N_data = new TH2D("dr_N_data","dr_N_data",150,0,15,50,0,10);  dr_N_data->Sumw2(); dr_N_data->GetXaxis()->SetTitle("track p (GeV/c)");  dr_N_data->GetYaxis()->SetTitle("dr_fvtx");
  TH2D *dr_S_data = new TH2D("dr_S_data","dr_S_data",150,0,15,50,0,10);  dr_S_data->Sumw2(); dr_S_data->GetXaxis()->SetTitle("track p (GeV/c)");  dr_S_data->GetYaxis()->SetTitle("dr_fvtx");
  
  if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("dr_N_data","Tr0_dr_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_N_data).c_str());
      T->Project("dr_S_data","Tr0_dr_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_S_data).c_str());
    }

  //now how about track matching dphi
  TH2D *dphi_N_data = new TH2D("dphi_N_data","dphi_N_data",150,0,15,100,-0.5,0.5);  dphi_N_data->Sumw2(); dphi_N_data->GetXaxis()->SetTitle("track p (GeV/c)");  dphi_N_data->GetYaxis()->SetTitle("dphi_fvtx");
  TH2D *dphi_S_data = new TH2D("dphi_S_data","dphi_S_data",150,0,15,100,-0.5,0.5);  dphi_S_data->Sumw2(); dphi_S_data->GetXaxis()->SetTitle("track p (GeV/c)");  dphi_S_data->GetYaxis()->SetTitle("dphi_fvtx");
  
 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("dphi_N_data","Tr0_dphi_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_N_data).c_str());
      T->Project("dphi_S_data","Tr0_dphi_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_S_data).c_str());
    }

  //and now track matching dtheta
  TH2D *dtheta_N_data = new TH2D("dtheta_N_data","dtheta_N_data",150,0,15,50,0,0.5);  dtheta_N_data->Sumw2(); dtheta_N_data->GetXaxis()->SetTitle("track p (GeV/c)");  dtheta_N_data->GetYaxis()->SetTitle("dtheta_fvtx");
  TH2D *dtheta_S_data = new TH2D("dtheta_S_data","dtheta_S_data",150,0,15,50,0,0.5);  dtheta_S_data->Sumw2(); dtheta_S_data->GetXaxis()->SetTitle("track p (GeV/c)");  dtheta_S_data->GetYaxis()->SetTitle("dtheta_fvtx");
  
 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("dtheta_N_data","Tr0_dtheta_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_N_data).c_str());
      T->Project("dtheta_S_data","Tr0_dtheta_fvtx:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_S_data).c_str());
    }

  //and now track matching DG0
  TH2D *DG0_N_data = new TH2D("DG0_N_data","DG0_N_data",150,0,15,50,0,30);  DG0_N_data->Sumw2(); DG0_N_data->GetXaxis()->SetTitle("track p (GeV/c)");  DG0_N_data->GetYaxis()->SetTitle("DG0");
  TH2D *DG0_S_data = new TH2D("DG0_S_data","DG0_S_data",150,0,15,50,0,30);  DG0_S_data->Sumw2(); DG0_S_data->GetXaxis()->SetTitle("track p (GeV/c)");  DG0_S_data->GetYaxis()->SetTitle("DG0");
  
 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("DG0_N_data","Tr0_DG0:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_N_data).c_str());
      T->Project("DG0_S_data","Tr0_DG0:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_S_data).c_str());
    }

  //and now track matching DDG0
  TH2D *DDG0_N_data = new TH2D("DDG0_N_data","DDG0_N_data",150,0,15,50,0,30);  DDG0_N_data->Sumw2(); DDG0_N_data->GetXaxis()->SetTitle("track p (GeV/c)");  DDG0_N_data->GetYaxis()->SetTitle("DDG0");
  TH2D *DDG0_S_data = new TH2D("DDG0_S_data","DDG0_S_data",150,0,15,50,0,30);  DDG0_S_data->Sumw2(); DDG0_S_data->GetXaxis()->SetTitle("track p (GeV/c)");  DDG0_S_data->GetYaxis()->SetTitle("DDG0");
  
 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("DDG0_N_data","Tr0_DDG0:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_N_data).c_str());
      T->Project("DDG0_S_data","Tr0_DDG0:sqrt(Tr0_px_fvtxmutr*Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr)",(ul_S_data).c_str());
    }
  
  ////draw the Tr0_trchi2 distribution
  TH1D *h_trchi2_N = new TH1D("h_trchi2_N","h_trchi2_N",200,0,20);
  TH1D *h_trchi2_S = new TH1D("h_trchi2_S","h_trchi2_S",200,0,20);

 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("h_trchi2_N","Tr0_trchi2",(ul_N_data).c_str());
      T->Project("h_trchi2_S","Tr0_trchi2",(ul_S_data).c_str());
    }

  ////draw the Tr0_chi2_fvtx distribution
  TH1D *h_chi2_fvtx_N = new TH1D("h_chi2_fvtx_N","h_chi2_fvtx_N",200,0,20);
  TH1D *h_chi2_fvtx_S = new TH1D("h_chi2_fvtx_S","h_chi2_fvtx_S",200,0,20);

 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("h_chi2_fvtx_N","Tr0_chi2_fvtx",(ul_N_data).c_str());
      T->Project("h_chi2_fvtx_S","Tr0_chi2_fvtx",(ul_S_data).c_str());
    }
  
  ////draw the Tr0_chi2_fvtxmutr distribution
  TH1D *h_chi2_fvtxmutr_N = new TH1D("h_chi2_fvtxmutr_N","h_chi2_fvtxmutr_N",200,0,20);
  TH1D *h_chi2_fvtxmutr_S = new TH1D("h_chi2_fvtxmutr_S","h_chi2_fvtxmutr_S",200,0,20);

 if( (run_num!=badrun_N) && (run_num!=badrun_S) )
    {
      T->Project("h_chi2_fvtxmutr_N","Tr0_chi2_fvtxmutr",(ul_N_data).c_str());
      T->Project("h_chi2_fvtxmutr_S","Tr0_chi2_fvtxmutr",(ul_S_data).c_str());
    }

  // only for testing locally
 // char const *outfile = "test.root";  

  TFile *h = new TFile(outfile, "RECREATE");
  h->cd();
    
  
  mass_fvtx_pt_N_ul->Write();
  mass_fvtx_pt_N_pp->Write();
  mass_fvtx_pt_N_mm->Write();
  mass_fvtx_pt_N_bg->Write();
  mass_fvtx_pt_N_mix->Write();
  
  mass_fvtx_pt_S_ul->Write();
  mass_fvtx_pt_S_pp->Write();
  mass_fvtx_pt_S_mm->Write();
  mass_fvtx_pt_S_bg->Write();
  mass_fvtx_pt_S_mix->Write();
 
  mass_fvtx_pt_020_N_ul->Write();
  mass_fvtx_pt_020_N_pp->Write();
  mass_fvtx_pt_020_N_mm->Write();
  mass_fvtx_pt_020_N_bg->Write();
  mass_fvtx_pt_020_N_mix->Write();
  
  mass_fvtx_pt_020_S_ul->Write();
  mass_fvtx_pt_020_S_pp->Write();
  mass_fvtx_pt_020_S_mm->Write();
  mass_fvtx_pt_020_S_bg->Write();
  mass_fvtx_pt_020_S_mix->Write();
 
  mass_fvtx_pt_2040_N_ul->Write();
  mass_fvtx_pt_2040_N_pp->Write();
  mass_fvtx_pt_2040_N_mm->Write();
  mass_fvtx_pt_2040_N_bg->Write();
  mass_fvtx_pt_2040_N_mix->Write();
  
  mass_fvtx_pt_2040_S_ul->Write();
  mass_fvtx_pt_2040_S_pp->Write();
  mass_fvtx_pt_2040_S_mm->Write();
  mass_fvtx_pt_2040_S_bg->Write();
  mass_fvtx_pt_2040_S_mix->Write();
 
  mass_fvtx_pt_4060_N_ul->Write();
  mass_fvtx_pt_4060_N_pp->Write();
  mass_fvtx_pt_4060_N_mm->Write();
  mass_fvtx_pt_4060_N_bg->Write();
  mass_fvtx_pt_4060_N_mix->Write();
  
  mass_fvtx_pt_4060_S_ul->Write();
  mass_fvtx_pt_4060_S_pp->Write();
  mass_fvtx_pt_4060_S_mm->Write();
  mass_fvtx_pt_4060_S_bg->Write();
  mass_fvtx_pt_4060_S_mix->Write();
 
  mass_fvtx_pt_6084_N_ul->Write();
  mass_fvtx_pt_6084_N_pp->Write();
  mass_fvtx_pt_6084_N_mm->Write();
  mass_fvtx_pt_6084_N_bg->Write();
  mass_fvtx_pt_6084_N_mix->Write();
  
  mass_fvtx_pt_6084_S_ul->Write();
  mass_fvtx_pt_6084_S_pp->Write();
  mass_fvtx_pt_6084_S_mm->Write();
  mass_fvtx_pt_6084_S_bg->Write();
  mass_fvtx_pt_6084_S_mix->Write();
 
  mass_good_N->Write();
  mass_bad_N->Write();
  mass_good_S->Write();
  mass_bad_S->Write();

  dr_N_data->Write();
  dr_S_data->Write();
  
  dphi_N_data->Write();
  dphi_S_data->Write();
  
  dtheta_N_data->Write();
  dtheta_S_data->Write();
  
  DG0_N_data->Write();
  DG0_S_data->Write();
  
  DDG0_N_data->Write();
  DDG0_S_data->Write();
  
  runs->Write();
  
  h_trchi2_N->Write();
  h_trchi2_S->Write();
  
  h_chi2_fvtx_N->Write();
  h_chi2_fvtx_S->Write();
  
  h_chi2_fvtxmutr_N->Write();
  h_chi2_fvtxmutr_S->Write(); 
  
}
