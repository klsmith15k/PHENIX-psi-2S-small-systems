
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
#include <iostream>

using namespace std;

void sng_trig_newlib_root5(char const *infile, char const *outfile)  // for running on CONDOR
//void sng_trig_newlib_root5() // for running locally over one fiel to test before submitting to CONDOR
 {

   // only for testing code locally before submitting to condor
   // char const *infile = "/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/jobfiles_PPG188/Run15pAu_NOFVTX_227.txt";   /// it says NOFVTX but just ignore
   //const char *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/jobfiles_local_sanghoon_version/Run15pAu_sanghoon_local_verison_1022.root";
   //const char *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/matt_centrality_macros/CONDOR/batch_pro105_newlib/Run15pAu_64.root";

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
  //  TFile *f = TFile::Open(infile);

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

 ////////////////////////////////////////////
  // MB North
  TH2D *h1 = new TH2D("h1","h1", 300, 0, 15, 120, 0, 12); h1->Sumw2(); h1->SetLineWidth(2);
  h1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *h2 = new TH2D("h2","h2", 300, 0, 15, 120, 0, 12); h2->Sumw2(); h2->SetLineWidth(2);
  h2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *h3 = new TH2D("h3","h3", 300, 0, 15, 120, 0, 12); h3->Sumw2(); h3->SetLineWidth(2);
  h3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *h4 = new TH2D("h4","h4", 300, 0, 15, 120, 0, 12); h4->Sumw2(); h4->SetLineWidth(2);
  h4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // MB South
  TH2D *h5 = new TH2D("h5","h5", 300, 0, 15, 120, 0, 12); h5->Sumw2(); h5->SetLineWidth(2);
  h5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *h6 = new TH2D("h6","h6", 300, 0, 15, 120, 0, 12); h6->Sumw2(); h6->SetLineWidth(2);
  h6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *h7 = new TH2D("h7","h7", 300, 0, 15, 120, 0, 12); h7->Sumw2(); h7->SetLineWidth(2);
  h7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *h8 = new TH2D("h8","h8", 300, 0, 15, 120, 0, 12); h8->Sumw2(); h8->SetLineWidth(2);
  h8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    h8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////
 // 0-20 North
  TH2D *t1 = new TH2D("t1","t1", 300, 0, 15, 120, 0, 12); t1->Sumw2(); t1->SetLineWidth(2);
  t1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *t2 = new TH2D("t2","t2", 300, 0, 15, 120, 0, 12); t2->Sumw2(); t2->SetLineWidth(2);
  t2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *t3 = new TH2D("t3","t3", 300, 0, 15, 120, 0, 12); t3->Sumw2(); t3->SetLineWidth(2);
  t3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *t4 = new TH2D("t4","t4", 300, 0, 15, 120, 0, 12); t4->Sumw2(); t4->SetLineWidth(2);
  t4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 0-20 South
  TH2D *t5 = new TH2D("t5","t5", 300, 0, 15, 120, 0, 12); t5->Sumw2(); t5->SetLineWidth(2);
  t5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *t6 = new TH2D("t6","t6", 300, 0, 15, 120, 0, 12); t6->Sumw2(); t6->SetLineWidth(2);
  t6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *t7 = new TH2D("t7","t7", 300, 0, 15, 120, 0, 12); t7->Sumw2(); t7->SetLineWidth(2);
  t7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *t8 = new TH2D("t8","t8", 300, 0, 15, 120, 0, 12); t8->Sumw2(); t8->SetLineWidth(2);
  t8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    t8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////
 // 20-40 North
  TH2D *s1 = new TH2D("s1","s1", 300, 0, 15, 120, 0, 12); s1->Sumw2(); s1->SetLineWidth(2);
  s1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *s2 = new TH2D("s2","s2", 300, 0, 15, 120, 0, 12); s2->Sumw2(); s2->SetLineWidth(2);
  s2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *s3 = new TH2D("s3","s3", 300, 0, 15, 120, 0, 12); s3->Sumw2(); s3->SetLineWidth(2);
  s3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *s4 = new TH2D("s4","s4", 300, 0, 15, 120, 0, 12); s4->Sumw2(); s4->SetLineWidth(2);
  s4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 20-40 South
  TH2D *s5 = new TH2D("s5","s5", 300, 0, 15, 120, 0, 12); s5->Sumw2(); s5->SetLineWidth(2);
  s5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *s6 = new TH2D("s6","s6", 300, 0, 15, 120, 0, 12); s6->Sumw2(); s6->SetLineWidth(2);
  s6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *s7 = new TH2D("s7","s7", 300, 0, 15, 120, 0, 12); s7->Sumw2(); s7->SetLineWidth(2);
  s7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *s8 = new TH2D("s8","s8", 300, 0, 15, 120, 0, 12); s8->Sumw2(); s8->SetLineWidth(2);
  s8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    s8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////
 // 40-60 North
  TH2D *w1 = new TH2D("w1","w1", 300, 0, 15, 120, 0, 12); w1->Sumw2(); w1->SetLineWidth(2);
  w1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *w2 = new TH2D("w2","w2", 300, 0, 15, 120, 0, 12); w2->Sumw2(); w2->SetLineWidth(2);
  w2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *w3 = new TH2D("w3","w3", 300, 0, 15, 120, 0, 12); w3->Sumw2(); w3->SetLineWidth(2);
  w3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *w4 = new TH2D("w4","w4", 300, 0, 15, 120, 0, 12); w4->Sumw2(); w4->SetLineWidth(2);
  w4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 40-60 South
  TH2D *w5 = new TH2D("w5","w5", 300, 0, 15, 120, 0, 12); w5->Sumw2(); w5->SetLineWidth(2);
  w5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *w6 = new TH2D("w6","w6", 300, 0, 15, 120, 0, 12); w6->Sumw2(); w6->SetLineWidth(2);
  w6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *w7 = new TH2D("w7","w7", 300, 0, 15, 120, 0, 12); w7->Sumw2(); w7->SetLineWidth(2);
  w7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *w8 = new TH2D("w8","w8", 300, 0, 15, 120, 0, 12); w8->Sumw2(); w8->SetLineWidth(2);
  w8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    w8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////
 // 60-84 North
  TH2D *z1 = new TH2D("z1","z1", 300, 0, 15, 120, 0, 12); z1->Sumw2(); z1->SetLineWidth(2);
  z1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *z2 = new TH2D("z2","z2", 300, 0, 15, 120, 0, 12); z2->Sumw2(); z2->SetLineWidth(2);
  z2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *z3 = new TH2D("z3","z3", 300, 0, 15, 120, 0, 12); z3->Sumw2(); z3->SetLineWidth(2);
  z3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *z4 = new TH2D("z4","z4", 300, 0, 15, 120, 0, 12); z4->Sumw2(); z4->SetLineWidth(2);
  z4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 60-84 South
  TH2D *z5 = new TH2D("z5","z5", 300, 0, 15, 120, 0, 12); z5->Sumw2(); z5->SetLineWidth(2);
  z5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *z6 = new TH2D("z6","z6", 300, 0, 15, 120, 0, 12); z6->Sumw2(); z6->SetLineWidth(2);
  z6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *z7 = new TH2D("z7","z7", 300, 0, 15, 120, 0, 12); z7->Sumw2(); z7->SetLineWidth(2);
  z7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *z8 = new TH2D("z8","z8", 300, 0, 15, 120, 0, 12); z8->Sumw2(); z8->SetLineWidth(2);
  z8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    z8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////

// 40-84 North
  TH2D *b1 = new TH2D("b1","b1", 300, 0, 15, 120, 0, 12); b1->Sumw2(); b1->SetLineWidth(2);
  b1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *b2 = new TH2D("b2","b2", 300, 0, 15, 120, 0, 12); b2->Sumw2(); b2->SetLineWidth(2);
  b2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *b3 = new TH2D("b3","b3", 300, 0, 15, 120, 0, 12); b3->Sumw2(); b3->SetLineWidth(2);
  b3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *b4 = new TH2D("b4","b4", 300, 0, 15, 120, 0, 12); b4->Sumw2(); b4->SetLineWidth(2);
  b4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 40-84 South
  TH2D *b5 = new TH2D("b5","b5", 300, 0, 15, 120, 0, 12); b5->Sumw2(); b5->SetLineWidth(2);
  b5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *b6 = new TH2D("b6","b6", 300, 0, 15, 120, 0, 12); b6->Sumw2(); b6->SetLineWidth(2);
  b6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *b7 = new TH2D("b7","b7", 300, 0, 15, 120, 0, 12); b7->Sumw2(); b7->SetLineWidth(2);
  b7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *b8 = new TH2D("b8","b8", 300, 0, 15, 120, 0, 12); b8->Sumw2(); b8->SetLineWidth(2);
  b8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    b8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////

// 20-84 North
  TH2D *d1 = new TH2D("d1","d1", 300, 0, 15, 120, 0, 12); d1->Sumw2(); d1->SetLineWidth(2);
  d1->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d1->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *d2 = new TH2D("d2","d2", 300, 0, 15, 120, 0, 12); d2->Sumw2(); d2->SetLineWidth(2);
  d2->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d2->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *d3 = new TH2D("d3","d3", 300, 0, 15, 120, 0, 12); d3->Sumw2(); d3->SetLineWidth(2);
  d3->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *d4 = new TH2D("d4","d4", 300, 0, 15, 120, 0, 12); d4->Sumw2(); d4->SetLineWidth(2);
  d4->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  // 20-84 South
  TH2D *d5 = new TH2D("d5","d5", 300, 0, 15, 120, 0, 12); d5->Sumw2(); d5->SetLineWidth(2);
  d5->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d5->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *d6 = new TH2D("d6","d6", 300, 0, 15, 120, 0, 12); d6->Sumw2(); d6->SetLineWidth(2);
  d6->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d6->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *d7 = new TH2D("d7","d7", 300, 0, 15, 120, 0, 12); d7->Sumw2(); d7->SetLineWidth(2);
  d7->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d7->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   
  TH2D *d8 = new TH2D("d8","d8", 300, 0, 15, 120, 0, 12); d8->Sumw2(); d8->SetLineWidth(2);
  d8->GetXaxis()->SetTitle("mass (GeV/c^{2})");    d8->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  ///////////////////////////////////////////////////////////////////////////////////


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
/////////// 40 - 84


 TH2D *mass_fvtx_pt_4084_N_ul = new TH2D("mass_fvtx_pt_4084_N_ul","mass_fvtx_pt_4084_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_N_ul->Sumw2(); mass_fvtx_pt_4084_N_ul->SetLineWidth(2);
  mass_fvtx_pt_4084_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_N_pp = new TH2D("mass_fvtx_pt_4084_N_pp","mass_fvtx_pt_4084_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_N_pp->Sumw2(); mass_fvtx_pt_4084_N_pp->SetLineWidth(2);
  mass_fvtx_pt_4084_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_N_mm = new TH2D("mass_fvtx_pt_4084_N_mm","mass_fvtx_pt_4084_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_N_mm->Sumw2(); mass_fvtx_pt_4084_N_mm->SetLineWidth(2);
  mass_fvtx_pt_4084_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_N_bg = new TH2D("mass_fvtx_pt_4084_N_bg","mass_fvtx_pt_4084_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_N_bg->Sumw2(); mass_fvtx_pt_4084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_4084_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_N_mix = new TH2D("mass_fvtx_pt_4084_N_mix","mass_fvtx_pt_4084_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_N_mix->Sumw2(); mass_fvtx_pt_4084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_4084_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
   mass_fvtx_pt_4084_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_4084_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_4084_N_ul->SetLineWidth(2); mass_fvtx_pt_4084_N_ul->SetMarkerStyle(20); mass_fvtx_pt_4084_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_4084_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_4084_N_bg->SetLineColor(kRed);  mass_fvtx_pt_4084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_4084_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_4084_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_4084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_4084_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_4084_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_4084_S_ul = new TH2D("mass_fvtx_pt_4084_S_ul","mass_fvtx_pt_4084_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_S_ul->Sumw2(); mass_fvtx_pt_4084_S_ul->SetLineWidth(2);
  mass_fvtx_pt_4084_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_S_pp = new TH2D("mass_fvtx_pt_4084_S_pp","mass_fvtx_pt_4084_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_S_pp->Sumw2(); mass_fvtx_pt_4084_S_pp->SetLineWidth(2);
  mass_fvtx_pt_4084_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_S_mm = new TH2D("mass_fvtx_pt_4084_S_mm","mass_fvtx_pt_4084_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_S_mm->Sumw2(); mass_fvtx_pt_4084_S_mm->SetLineWidth(2);
  mass_fvtx_pt_4084_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_S_bg = new TH2D("mass_fvtx_pt_4084_S_bg","mass_fvtx_pt_4084_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_S_bg->Sumw2(); mass_fvtx_pt_4084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_4084_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_4084_S_mix = new TH2D("mass_fvtx_pt_4084_S_mix","mass_fvtx_pt_4084_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_4084_S_mix->Sumw2(); mass_fvtx_pt_4084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_4084_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_4084_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_4084_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_4084_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_4084_S_ul->SetLineWidth(2); mass_fvtx_pt_4084_S_ul->SetMarkerStyle(20); mass_fvtx_pt_4084_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_4084_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_4084_S_bg->SetLineColor(kRed);  mass_fvtx_pt_4084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_4084_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_4084_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_4084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_4084_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_4084_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_4084_S_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///////////  20 - 84


 TH2D *mass_fvtx_pt_2084_N_ul = new TH2D("mass_fvtx_pt_2084_N_ul","mass_fvtx_pt_2084_N_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_N_ul->Sumw2(); mass_fvtx_pt_2084_N_ul->SetLineWidth(2);
  mass_fvtx_pt_2084_N_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_N_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_N_pp = new TH2D("mass_fvtx_pt_2084_N_pp","mass_fvtx_pt_2084_N_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_N_pp->Sumw2(); mass_fvtx_pt_2084_N_pp->SetLineWidth(2);
  mass_fvtx_pt_2084_N_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_N_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_N_mm = new TH2D("mass_fvtx_pt_2084_N_mm","mass_fvtx_pt_2084_N_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_N_mm->Sumw2(); mass_fvtx_pt_2084_N_mm->SetLineWidth(2);
  mass_fvtx_pt_2084_N_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_N_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_N_bg = new TH2D("mass_fvtx_pt_2084_N_bg","mass_fvtx_pt_2084_N_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_N_bg->Sumw2(); mass_fvtx_pt_2084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_2084_N_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_N_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_N_mix = new TH2D("mass_fvtx_pt_2084_N_mix","mass_fvtx_pt_2084_N_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_N_mix->Sumw2(); mass_fvtx_pt_2084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_2084_N_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_N_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
   mass_fvtx_pt_2084_N_ul->SetMarkerColor(kBlue); mass_fvtx_pt_2084_N_ul->SetLineColor(kBlue);  mass_fvtx_pt_2084_N_ul->SetLineWidth(2); mass_fvtx_pt_2084_N_ul->SetMarkerStyle(20); mass_fvtx_pt_2084_N_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_2084_N_bg->SetMarkerColor(kRed);  mass_fvtx_pt_2084_N_bg->SetLineColor(kRed);  mass_fvtx_pt_2084_N_bg->SetLineWidth(2);
  mass_fvtx_pt_2084_N_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_2084_N_mix->SetLineColor(kBlack);  mass_fvtx_pt_2084_N_mix->SetLineWidth(2);
  mass_fvtx_pt_2084_N_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_N_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_N_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_N_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_2084_N_pp->GetXaxis()->SetRangeUser(2,5);
  
  ///Now S arm histos
  TH2D *mass_fvtx_pt_2084_S_ul = new TH2D("mass_fvtx_pt_2084_S_ul","mass_fvtx_pt_2084_S_ul", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_S_ul->Sumw2(); mass_fvtx_pt_2084_S_ul->SetLineWidth(2);
  mass_fvtx_pt_2084_S_ul->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_S_ul->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_S_pp = new TH2D("mass_fvtx_pt_2084_S_pp","mass_fvtx_pt_2084_S_pp", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_S_pp->Sumw2(); mass_fvtx_pt_2084_S_pp->SetLineWidth(2);
  mass_fvtx_pt_2084_S_pp->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_S_pp->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_S_mm = new TH2D("mass_fvtx_pt_2084_S_mm","mass_fvtx_pt_2084_S_mm", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_S_mm->Sumw2(); mass_fvtx_pt_2084_S_mm->SetLineWidth(2);
  mass_fvtx_pt_2084_S_mm->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_S_mm->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_S_bg = new TH2D("mass_fvtx_pt_2084_S_bg","mass_fvtx_pt_2084_S_bg", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_S_bg->Sumw2(); mass_fvtx_pt_2084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_2084_S_bg->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_S_bg->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  TH2D *mass_fvtx_pt_2084_S_mix = new TH2D("mass_fvtx_pt_2084_S_mix","mass_fvtx_pt_2084_S_mix", 300, 0, 15, 120, 0, 12); mass_fvtx_pt_2084_S_mix->Sumw2(); mass_fvtx_pt_2084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_2084_S_mix->GetXaxis()->SetTitle("mass (GeV/c^{2})");    mass_fvtx_pt_2084_S_mix->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  
  //pretty stuff
  mass_fvtx_pt_2084_S_ul->SetMarkerColor(kBlue); mass_fvtx_pt_2084_S_ul->SetLineColor(kBlue);  mass_fvtx_pt_2084_S_ul->SetLineWidth(2); mass_fvtx_pt_2084_S_ul->SetMarkerStyle(20); mass_fvtx_pt_2084_S_ul->SetMarkerSize(0.75);
  mass_fvtx_pt_2084_S_bg->SetMarkerColor(kRed);  mass_fvtx_pt_2084_S_bg->SetLineColor(kRed);  mass_fvtx_pt_2084_S_bg->SetLineWidth(2);
  mass_fvtx_pt_2084_S_mix->SetMarkerColor(kBlack);  mass_fvtx_pt_2084_S_mix->SetLineColor(kBlack);  mass_fvtx_pt_2084_S_mix->SetLineWidth(2);
  mass_fvtx_pt_2084_S_ul->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_S_bg->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_S_mix->GetXaxis()->SetRangeUser(2,5);  mass_fvtx_pt_2084_S_mm->GetXaxis()->SetRangeUser(2,5); mass_fvtx_pt_2084_S_pp->GetXaxis()->SetRangeUser(2,5);
  
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

 // //////////////////////////////////////////////
     TTreeReader myReader("T1", f);
// //   /////////////////////////////////////////////////
  TTreeReaderValue<int> run_number(myReader, "RunNumber");

  unsigned int run_num = 0;
  int badrun_S = 0;
  int badrun_N = 0;
  int counter = 0;

 std::string  filename[2] = {"/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/matt_centrality_macros/badruns_S.txt","/phenix/hhj/klsmith/Jpsis/Run15pAu/fit_NOFVTX/krista_fits/matt_centrality_macros/badruns_N.txt"};

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
	 /// cout << "badrun S: " << badrun_S << endl;
	 // cout << "run number: " << run_num << endl;

 	 counter++;

	 if(badrun_S == run_num)
	   break;  // continue here does not work

       }while(counter< 10);

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
	 // cout << "badrun N: " << badrun_N << endl;
	 //	 cout << "run number: " << run_num << endl;

 	 counter++;

	 if(badrun_N == run_num)
	   break;  // continue here does not work

       }while(counter< 30);

   } // end if
 else
   cout << "test read file does not exist" << endl;

 cout << "out of ifstream" << endl;
 if(badrun_S == run_num)
   cout << "South arm bad run" << endl;
 if(badrun_N == run_num)
   cout << "North arm bad run" << endl;
 
 //////////////////////////////////////////////////////////////////////////////////////
  string run_cut_S = "runnumber!=-999";   
  string run_cut_N = "runnumber!=-999";

  string mass_cut = "mass_fvtxmutr > 1.5";  // single tracks are stored in fvtxmutr variable
  string evt_cut = "(Evt_bbcZ>-30) && (Evt_bbcZ<30)";  // already applied in dimuon container make_dimu.C

  string trig_cut_N = "lvl1_trigscaled&0x00100000";
  string trig_cut_S = "lvl1_trigscaled&0x00200000";
  string nidhits_cut = "Tr0_idhits>14 && Tr1_idhits>14";
  string lastgap_cut = "(Tr0_lastgap>=3 && Tr1_lastgap>=3)";


  string DG0_cut_N1 = " ((Tr0_nhits_fvtx > 0) && (Tr0_DG0<(8.30848+100.812/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz))))";
  string DG0_cut_S1 = " ((Tr0_nhits_fvtx > 0) && (Tr0_DG0<(10.6218+337.951/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz))))";
  string DDG0_cut_N1 = " ((Tr0_nhits_fvtx > 0) && (Tr0_DDG0<(4.84265+127.384/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz))))";
  string DDG0_cut_S1 = " ((Tr0_nhits_fvtx > 0) && (Tr0_DDG0<(5.71576+111.407/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz))))";

  string DG0_cut_N2 = " ((Tr1_nhits_fvtx > 0) && (Tr1_DG0<(8.30848+100.812/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz))))";
  string DG0_cut_S2 = "((Tr1_nhits_fvtx > 0) && (Tr1_DG0<(10.6218+337.951/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz))))";
  string DDG0_cut_N2 = "((Tr1_nhits_fvtx > 0) && (Tr1_DDG0<(4.84265+127.384/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz))))";
  string DDG0_cut_S2 = " ((Tr1_nhits_fvtx > 0) && (Tr1_DDG0<(5.71576+111.407/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz))))";

  string tkchi2_cut = "Tr0_trchi2<23 && Tr1_trchi2<23";
  string pz_cut_N = "abs(Tr0_pz)>2 && abs(Tr1_pz)>2&&Tr0_pz>0&&Tr1_pz>0";
  string pz_cut_S = "abs(Tr0_pz)>2 && abs(Tr1_pz)>2&&Tr0_pz<0&&Tr1_pz<0";
  string asym_cut = "abs( (sqrt(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz ) - sqrt(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz ))/ (sqrt(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz ) + sqrt(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz )) ) < 1 ";
  string octant_cut ="abs(  ( ((int((atan2(Tr0_xst1,Tr0_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) - ( ((int((atan2(Tr1_xst1,Tr1_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) ) >0";
 
  // PPG188 fiducial cuts
//string fid_cut_N = " !(  atan2(Tr0_yst3,Tr0_xst3)>0.35 && atan2(Tr0_yst3,Tr0_xst3)<0.77  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>0.77 && atan2(Tr0_yst3,Tr0_xst3)<1.15 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>160 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>-2.75 && atan2(Tr0_yst3,Tr0_xst3)<-1.92 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>237 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst2,Tr0_xst2)>-2.75 && atan2(Tr0_yst2,Tr0_xst2)<-1.92 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)>140 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)<163  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>-2.75 && atan2(Tr0_yst1,Tr0_xst1)<-2 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>76 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<80  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>0.52 && atan2(Tr0_yst1,Tr0_xst1)<1.11 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>60 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<75 )  && !(  atan2(Tr1_yst3,Tr1_xst3)>0.35 && atan2(Tr1_yst3,Tr1_xst3)<0.77  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>0.77 && atan2(Tr1_yst3,Tr1_xst3)<1.15 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>160 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>-2.75 && atan2(Tr1_yst3,Tr1_xst3)<-1.92 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>237 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst2,Tr1_xst2)>-2.75 && atan2(Tr1_yst2,Tr1_xst2)<-1.92 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)>140 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)<163  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>-2.75 && atan2(Tr1_yst1,Tr1_xst1)<-2 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>76 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<80  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>0.52 && atan2(Tr1_yst1,Tr1_xst1)<1.11 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>60 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<75 )";
  //string fid_cut_N = "Evt_Cent!=-12";
  //string fid_cut_S = "Evt_Cent!=-12";

//////////////////////////////////////////////////////////////  PPG228 fiducial cuts January 17, 2021
//string fid_cut_N = "!(  atan2(Tr0_yst3,Tr0_xst3)>0.35 && atan2(Tr0_yst3,Tr0_xst3)<0.77  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>0.77 && atan2(Tr0_yst3,Tr0_xst3)<1.15 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>160 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>-2.75 && atan2(Tr0_yst3,Tr0_xst3)<-1.92 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>237 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst2,Tr0_xst2)>-2.75 && atan2(Tr0_yst2,Tr0_xst2)<-1.92 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)>140 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)<163  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>-2.75 && atan2(Tr0_yst1,Tr0_xst1)<-2 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>76 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<80  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>0.52 && atan2(Tr0_yst1,Tr0_xst1)<1.11 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>60 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<75 )  && !(  atan2(Tr1_yst3,Tr1_xst3)>0.35 && atan2(Tr1_yst3,Tr1_xst3)<0.77  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>0.77 && atan2(Tr1_yst3,Tr1_xst3)<1.15 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>160 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>-2.75 && atan2(Tr1_yst3,Tr1_xst3)<-1.92 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>237 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst2,Tr1_xst2)>-2.75 && atan2(Tr1_yst2,Tr1_xst2)<-1.92 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)>140 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)<163  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>-2.75 && atan2(Tr1_yst1,Tr1_xst1)<-2 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>76 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<80  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>0.52 && atan2(Tr1_yst1,Tr1_xst1)<1.11 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>60 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<75 ) && !( atan2(Tr0_yst2,Tr0_xst2)>-1.20 && atan2(Tr0_yst2,Tr0_xst2)<-0.90 )&& !( atan2(Tr1_yst2,Tr1_xst2)>-1.20 && atan2(Tr1_yst2,Tr1_xst2)<-0.90 )";
  string fid_cut_S = "!( atan2(Tr0_yst1,Tr0_xst1)>0 && atan2(Tr0_yst1,Tr0_xst1)<0.4 ) && !( atan2(Tr0_yst2,Tr0_xst2)>0 && atan2(Tr0_yst2,Tr0_xst2)<0.4 ) && !( atan2(Tr0_yst3,Tr0_xst3)>0 && atan2(Tr0_yst3,Tr0_xst3)<0.4 ) && !( atan2(Tr0_yst1,Tr0_xst1)>1.55 && atan2(Tr0_yst1,Tr0_xst1)<1.9 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>1.55 && atan2(Tr0_yst2,Tr0_xst2)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>1.55 && atan2(Tr0_yst3,Tr0_xst3)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>-0.60 && atan2(Tr0_yst3,Tr0_xst3)<-0.40 ) && !( atan2(Tr1_yst1,Tr1_xst1)>0 && atan2(Tr1_yst1,Tr1_xst1)<0.4 ) && !( atan2(Tr1_yst2,Tr1_xst2)>0 && atan2(Tr1_yst2,Tr1_xst2)<0.4 ) && !( atan2(Tr1_yst3,Tr1_xst3)>0 && atan2(Tr1_yst3,Tr1_xst3)<0.4 ) && !( atan2(Tr1_yst1,Tr1_xst1)>1.55 && atan2(Tr1_yst1,Tr1_xst1)<1.9 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>1.55 && atan2(Tr1_yst2,Tr1_xst2)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>1.55 && atan2(Tr1_yst3,Tr1_xst3)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>-0.60 && atan2(Tr1_yst3,Tr1_xst3)<-0.40 ) && !( atan2(Tr0_yst1,Tr0_xst1)>-2.35 && atan2(Tr0_yst1,Tr0_xst1)<-1.95 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>-2.35 && atan2(Tr0_yst2,Tr0_xst2)<-1.95 ) &&  !( atan2(Tr0_yst3,Tr0_xst3)>-2.35 && atan2(Tr0_yst3,Tr0_xst3)<-1.95 ) && !( atan2(Tr1_yst1,Tr1_xst1)>-2.35 && atan2(Tr1_yst1,Tr1_xst1)<-1.95 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>-2.35 && atan2(Tr1_yst2,Tr1_xst2)<-1.95 ) &&  !( atan2(Tr1_yst3,Tr1_xst3)>-2.35 && atan2(Tr1_yst3,Tr1_xst3)<-1.95 )"; 
 ////////////////////////////////////////////////////////////
  // string fid_cut_S = "Evt_Cent!=-12";
 string fid_cut_N = "Evt_Cent!=-12";

  string muon_arm_cut_N1 =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_N1 + "&&" + DDG0_cut_N1 + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_N + "&&" + run_cut_N;
  string muon_arm_cut_S1 =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_S1 + "&&" + DDG0_cut_S1 + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_S + "&&" + run_cut_S;

  string muon_arm_cut_N2 =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_N2 + "&&" + DDG0_cut_N2 + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_N + "&&" + run_cut_N;
  string muon_arm_cut_S2 =  lastgap_cut + "&&" + nidhits_cut + "&&" + DG0_cut_S2 + "&&" + DDG0_cut_S2 + "&&" + tkchi2_cut + "&&" + octant_cut + "&&" + asym_cut + "&&" + fid_cut_S + "&&" + run_cut_S;

  //string chi2_fvtx_cut = "Tr0_chi2_fvtx<10 && Tr0_chi2_fvtx>0 && Tr1_chi2_fvtx<10 && Tr1_chi2_fvtx>0";
  // string chi2_cut = "Tr0_chi2<10 && Tr0_chi2>0 && Tr1_chi2<10 && Tr1_chi2>0";

  string dr_cut_N1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dr_fvtx<(1.38061+28.89/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";
  string dr_cut_S1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dr_fvtx<(1.67962+26.8308/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";
  string dtheta_cut_N1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dtheta_fvtx<(0.0473169+0.665664/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";
  string dtheta_cut_S1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dtheta_fvtx<(0.062556+0.580508/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";
  string dphi_cut_N1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dphi_fvtx<(0.121509+1.55712/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";
  string dphi_cut_S1 = "(Tr0_nhits_fvtx > 0) && (Tr0_dphi_fvtx<(0.121325+1.79439/(Tr0_px * Tr0_px + Tr0_py*Tr0_py + Tr0_pz*Tr0_pz)))";

  string dr_cut_N2 = " (Tr1_nhits_fvtx > 0) && (Tr1_dr_fvtx<(1.38061+28.89/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";
  string dr_cut_S2 = "(Tr1_nhits_fvtx > 0) && (Tr1_dr_fvtx<(1.67962+26.8308/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";
  string dtheta_cut_N2 = " (Tr1_nhits_fvtx > 0) && (Tr1_dtheta_fvtx<(0.0473169+0.665664/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";
  string dtheta_cut_S2 = " (Tr1_nhits_fvtx > 0) && (Tr1_dtheta_fvtx<(0.062556+0.580508/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";
  string dphi_cut_N2 = "(Tr1_nhits_fvtx > 0) && (Tr1_dphi_fvtx<(0.121509+1.55712/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";
  string dphi_cut_S2 = " (Tr1_nhits_fvtx > 0) && (Tr1_dphi_fvtx<(0.121325+1.79439/(Tr1_px * Tr1_px + Tr1_py*Tr1_py + Tr1_pz*Tr1_pz)))";

  //  string sng_fvtx_cut = chi2_fvtx_cut + "&&" + chi2_cut;
  string rap_cut_S = "rapidity<=-1.2 && rapidity>=-2.2";
  string rap_cut_N = "rapidity>=1.2 && rapidity<=2.2";
  
  string sng_trk_vtx1 = "Tr0_nhits_fvtx > 0 && Tr1_nhits_fvtx == 0";
  string sng_trk_vtx2 = "Tr1_nhits_fvtx > 0 && Tr0_nhits_fvtx == 0";

  string cuts_N1 = "&&" + rap_cut_N  + "&&" + pz_cut_N + "&&"+ mass_cut + "&&" + evt_cut + "&&" + muon_arm_cut_N1 + "&&" + dr_cut_N1 + "&&" + dtheta_cut_N1 + "&&" + dphi_cut_N1 + "&&" + sng_trk_vtx1;
  string cuts_S1 = "&&" + rap_cut_S  + "&&" + pz_cut_S + "&&" + mass_cut + "&&" + evt_cut + "&&" + muon_arm_cut_S1 + "&&" + dr_cut_S1 + "&&" + dtheta_cut_S1 + "&&" + dphi_cut_S1 + "&&" + sng_trk_vtx1;

  string cuts_N2 = "&&" + rap_cut_N  + "&&" + pz_cut_N + "&&"+ mass_cut + "&&" + evt_cut + "&&" + muon_arm_cut_N2 + "&&" + dr_cut_N2 + "&&" + dtheta_cut_N2 + "&&" + dphi_cut_N2 + "&&" + sng_trk_vtx2;
  string cuts_S2 = "&&" + rap_cut_S  + "&&"  + pz_cut_S + "&&" + mass_cut + "&&" + evt_cut + "&&"+ muon_arm_cut_S2 + "&&" + dr_cut_S2 + "&&" + dtheta_cut_S2 + "&&" + dphi_cut_S2 + "&&" + sng_trk_vtx2;

  //string sng12 = " (Tr0_nhits_fvtx > 0)  && (Tr1_nhits_fvtx == 0) ";
  ////////////////////////////

  int lo_bin =  mass_fvtx_pt_N_pp->GetXaxis()->FindBin(2);
  int hi_bin =  mass_fvtx_pt_N_pp->GetXaxis()->FindBin(5.5);
 
 
  
  if(badrun_N != run_num)
    {
      T->Project("mass_fvtx_pt_N_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1" + cuts_N1).c_str() );
    
      T->Project("mass_fvtx_pt_N_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1" + cuts_N1).c_str() ); 
    
      T->Project("mass_fvtx_pt_N_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1" + cuts_N1).c_str() );
  
      T->Project("mass_fvtx_pt_N_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_N1).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 0-20% centrality
  
      T->Project("mass_fvtx_pt_020_N_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N1).c_str() );
    
      T->Project("mass_fvtx_pt_020_N_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N1).c_str() ); 
  
      T->Project("mass_fvtx_pt_020_N_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N1).c_str() );
 
      T->Project("mass_fvtx_pt_020_N_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_N1).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 20-40% centrality
  
      T->Project("mass_fvtx_pt_2040_N_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N1).c_str() );
    
      T->Project("mass_fvtx_pt_2040_N_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N1).c_str() ); 
  
      T->Project("mass_fvtx_pt_2040_N_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N1).c_str() );
  
      T->Project("mass_fvtx_pt_2040_N_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_N1).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 40-60% centrality
 
      T->Project("mass_fvtx_pt_4060_N_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N1).c_str() );
    
      T->Project("mass_fvtx_pt_4060_N_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N1).c_str() ); 
  
      T->Project("mass_fvtx_pt_4060_N_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N1).c_str() );
 
      T->Project("mass_fvtx_pt_4060_N_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_N1).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 60-84% centrality
 
      T->Project("mass_fvtx_pt_6084_N_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N1).c_str() );
  
      T->Project("mass_fvtx_pt_6084_N_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N1).c_str() ); 
  
      T->Project("mass_fvtx_pt_6084_N_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N1).c_str() );
  
      T->Project("mass_fvtx_pt_6084_N_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_N1).c_str() );


  /// 40-84
  T->Project("mass_fvtx_pt_4084_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N1).c_str() );
  
  T->Project("mass_fvtx_pt_4084_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N1).c_str() ); 
   
  T->Project("mass_fvtx_pt_4084_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N1).c_str() );
 
  T->Project("mass_fvtx_pt_4084_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=84" + cuts_N1).c_str() );

  ///project 20 - 84
  T->Project("mass_fvtx_pt_2084_N_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N1).c_str() );
  
  T->Project("mass_fvtx_pt_2084_N_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N1).c_str() ); 
   
  T->Project("mass_fvtx_pt_2084_N_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N1).c_str() );
 
  T->Project("mass_fvtx_pt_2084_N_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=84" + cuts_N1).c_str() );

    } // end N

  //=====================================================

 
  if(badrun_N != run_num)
    {

      T->Project("h1","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1" + cuts_N2).c_str() );
    
      T->Project("h2","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1" + cuts_N2).c_str() ); 
    
      T->Project("h3","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1" + cuts_N2).c_str() );
  
      T->Project("h4","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_N2).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 0-20% centrality
  
      T->Project("t1","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N2).c_str() );
    
      T->Project("t2","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N2).c_str() ); 
  
      T->Project("t3","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_N2).c_str() );
 
      T->Project("t4","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_N2).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 20-40% centrality
  
      T->Project("s1","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N2).c_str() );
    
      T->Project("s2","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N2).c_str() ); 
  
      T->Project("s3","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_N2).c_str() );
  
      T->Project("s4","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_N2).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 40-60% centrality
 
      T->Project("w1","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N2).c_str() );
    
      T->Project("w2","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N2).c_str() ); 
  
      T->Project("w3","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_N2).c_str() );
 
      T->Project("w4","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_N2).c_str() );

      /////////////////////////////////////
      ///project North arm histograms, 60-84% centrality
 
      T->Project("z1","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N2).c_str() );
  
      T->Project("z2","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N2).c_str() ); 
  
      T->Project("z3","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_N2).c_str() );
  
      T->Project("z4","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_N2).c_str() );


      /// 40-84
      T->Project("b1","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N2).c_str() );
  
      T->Project("b2","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N2).c_str() ); 
   
      T->Project("b3","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>40&&centrality<=84" + cuts_N2).c_str() );
 
      T->Project("b4","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=84" + cuts_N2).c_str() );

      ///project 20 - 84
      T->Project("d1","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge==0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N2).c_str() );
  
      T->Project("d2","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge>0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N2).c_str() ); 
   
      T->Project("d3","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00100000&&charge<0&&same_event==1&&centrality>20&&centrality<=84" + cuts_N2).c_str() );
 
      T->Project("d4","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=84" + cuts_N2).c_str() );



      mass_fvtx_pt_N_ul->Add(h1,1);  // add 2nd track combination to the first track combination ( otherwise 1st histogram is overwritten)
      mass_fvtx_pt_N_pp->Add(h2, 1);
      mass_fvtx_pt_N_mm->Add(h3, 1);
      mass_fvtx_pt_N_mix->Add(h4, 1);

      mass_fvtx_pt_020_N_ul->Add(t1,1); 
      mass_fvtx_pt_020_N_pp->Add(t2, 1);
      mass_fvtx_pt_020_N_mm->Add(t3, 1);
      mass_fvtx_pt_020_N_mix->Add(t4, 1);

      mass_fvtx_pt_2040_N_ul->Add(s1,1); 
      mass_fvtx_pt_2040_N_pp->Add(s2, 1);
      mass_fvtx_pt_2040_N_mm->Add(s3, 1);
      mass_fvtx_pt_2040_N_mix->Add(s4, 1);

      mass_fvtx_pt_4060_N_ul->Add(w1,1); 
      mass_fvtx_pt_4060_N_pp->Add(w2, 1);
      mass_fvtx_pt_4060_N_mm->Add(w3, 1);
      mass_fvtx_pt_4060_N_mix->Add(w4, 1);

      mass_fvtx_pt_6084_N_ul->Add(z1,1);  
      mass_fvtx_pt_6084_N_pp->Add(z2, 1);
      mass_fvtx_pt_6084_N_mm->Add(z3, 1);
      mass_fvtx_pt_6084_N_mix->Add(z4, 1);

      mass_fvtx_pt_4084_N_ul->Add(b1,1);  
      mass_fvtx_pt_4084_N_pp->Add(b2, 1);
      mass_fvtx_pt_4084_N_mm->Add(b3, 1);
      mass_fvtx_pt_4084_N_mix->Add(b4, 1);

      mass_fvtx_pt_2084_N_ul->Add(d1,1);  
      mass_fvtx_pt_2084_N_pp->Add(d2, 1);
      mass_fvtx_pt_2084_N_mm->Add(d3, 1);
      mass_fvtx_pt_2084_N_mix->Add(d4, 1);


    } // end N

  //=========================================================================

 
  if(run_num!=badrun_S)
    {
      ///project South arm histograms
      T->Project("mass_fvtx_pt_S_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1" + cuts_S1).c_str() );
   
      T->Project("mass_fvtx_pt_S_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_S_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1" + cuts_S1).c_str() );
 
      T->Project("mass_fvtx_pt_S_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_S1).c_str() );

      ///project South arm histograms
      T->Project("mass_fvtx_pt_020_S_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_020_S_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_020_S_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_020_S_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_S1).c_str() );
 
      ///project South arm histograms
      T->Project("mass_fvtx_pt_2040_S_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S1).c_str() );
 
      T->Project("mass_fvtx_pt_2040_S_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S1).c_str() ); 
    
      T->Project("mass_fvtx_pt_2040_S_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_2040_S_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_S1).c_str() );
  
      /////////////////////////////////////

      ///project South arm histograms
      T->Project("mass_fvtx_pt_4060_S_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_4060_S_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_4060_S_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_4060_S_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_S1).c_str() );

      ///project South arm histograms
      T->Project("mass_fvtx_pt_6084_S_ul","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_6084_S_pp","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_6084_S_mm","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S1).c_str() );
 
      T->Project("mass_fvtx_pt_6084_S_mix","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_S1).c_str() );

 /// 40-84
      T->Project("mass_fvtx_pt_4084_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_4084_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_4084_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S1).c_str() );
 
      T->Project("mass_fvtx_pt_4084_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=84" + cuts_S1).c_str() );

      ///project 20 - 84
      T->Project("mass_fvtx_pt_2084_S_ul","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S1).c_str() );
  
      T->Project("mass_fvtx_pt_2084_S_pp","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S1).c_str() ); 
   
      T->Project("mass_fvtx_pt_2084_S_mm","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S1).c_str() );
 
      T->Project("mass_fvtx_pt_2084_S_mix","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=84" + cuts_S1).c_str() );


    } // end S

  //=========================================================================

 
  if(run_num!=badrun_S)
    {
      ///project South arm histograms
      T->Project("h5","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1" + cuts_S2).c_str() );
   
      T->Project("h6","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1" + cuts_S2).c_str() ); 
   
      T->Project("h7","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1" + cuts_S2).c_str() );
 
      T->Project("h8","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0" + cuts_S2).c_str() );

      ///project South arm histograms
      T->Project("t5","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S2).c_str() );
  
      T->Project("t6","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S2).c_str() ); 
   
      T->Project("t7","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>=0&&centrality<=20" + cuts_S2).c_str() );
  
      T->Project("t8","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>=0&&centrality<=20" + cuts_S2).c_str() );
 
      ///project South arm histograms
      T->Project("s5","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S2).c_str() );
 
      T->Project("s6","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S2).c_str() ); 
    
      T->Project("s7","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>20&&centrality<=40" + cuts_S2).c_str() );
  
      T->Project("s8","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=40" + cuts_S2).c_str() );
  
      /////////////////////////////////////

      ///project South arm histograms
      T->Project("w5","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S2).c_str() );
  
      T->Project("w6","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S2).c_str() ); 
   
      T->Project("w7","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>40&&centrality<=60" + cuts_S2).c_str() );
  
      T->Project("w8","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=60" + cuts_S2).c_str() );

      ///project South arm histograms
      T->Project("z5","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S2).c_str() );
  
      T->Project("z6","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S2).c_str() ); 
   
      T->Project("z7","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>60&&centrality<=84" + cuts_S2).c_str() );
 
      T->Project("z8","sqrt(  ((Tr0_px + Tr1_px)*(Tr0_px + Tr1_px)) +  ((Tr0_py + Tr1_py)*(Tr0_py + Tr1_py))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>60&&centrality<=84" + cuts_S2).c_str() );


 /// 40-84
      T->Project("b5","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S2).c_str() );
  
      T->Project("b6","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S2).c_str() ); 
   
      T->Project("b7","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>40&&centrality<=84" + cuts_S2).c_str() );
 
      T->Project("b8","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>40&&centrality<=84" + cuts_S2).c_str() );

      ///project 20 - 84
      T->Project("d5","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge==0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S2).c_str() );
  
      T->Project("d6","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge>0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S2).c_str() ); 
   
      T->Project("d7","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("lvl1_trigscaled&0x00200000&&charge<0&&same_event==1&&centrality>20&&centrality<=84" + cuts_S2).c_str() );
 
      T->Project("d8","sqrt(  ((Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)*(Tr0_px_fvtxmutr + Tr1_px_fvtxmutr)) +  ((Tr0_py_fvtxmutr + Tr1_py_fvtxmutr)*(Tr0_py_fvtxmutr + Tr1_py_fvtxmutr))  ):mass_fvtxmutr", ("charge==0&&same_event==0&&centrality>20&&centrality<=84" + cuts_S2).c_str() );


      mass_fvtx_pt_S_ul->Add(h5,1);  // add 6nd track combination to the first track combination ( otherwise 1st histogram is overwritten)
      mass_fvtx_pt_S_pp->Add(h6, 1);
      mass_fvtx_pt_S_mm->Add(h7, 1);
      mass_fvtx_pt_S_mix->Add(h8, 1);

      mass_fvtx_pt_020_S_ul->Add(t5,1); 
      mass_fvtx_pt_020_S_pp->Add(t6, 1);
      mass_fvtx_pt_020_S_mm->Add(t7, 1);
      mass_fvtx_pt_020_S_mix->Add(t8, 1);

      mass_fvtx_pt_2040_S_ul->Add(s5,1); 
      mass_fvtx_pt_2040_S_pp->Add(s6, 1);
      mass_fvtx_pt_2040_S_mm->Add(s7, 1);
      mass_fvtx_pt_2040_S_mix->Add(s8, 1);

      mass_fvtx_pt_4060_S_ul->Add(w5,1); 
      mass_fvtx_pt_4060_S_pp->Add(w6, 1);
      mass_fvtx_pt_4060_S_mm->Add(w7, 1);
      mass_fvtx_pt_4060_S_mix->Add(w8, 1);

      mass_fvtx_pt_6084_S_ul->Add(z5,1);  
      mass_fvtx_pt_6084_S_pp->Add(z6, 1);
      mass_fvtx_pt_6084_S_mm->Add(z7, 1);
      mass_fvtx_pt_6084_S_mix->Add(z8, 1);

      mass_fvtx_pt_4084_S_ul->Add(b5,1);  
      mass_fvtx_pt_4084_S_pp->Add(b6, 1);
      mass_fvtx_pt_4084_S_mm->Add(b7, 1);
      mass_fvtx_pt_4084_S_mix->Add(b8, 1);

      mass_fvtx_pt_2084_S_ul->Add(d5,1);  
      mass_fvtx_pt_2084_S_pp->Add(d6, 1);
      mass_fvtx_pt_2084_S_mm->Add(d7, 1);
      mass_fvtx_pt_2084_S_mix->Add(d8, 1);

    } // end S


  //////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  
  c1->Divide(2,1);
  c1->cd(1); c1->cd(1)->SetLogy();
  // mass_fvtx_pt_S_mix->ProjectionX()->Draw();
  mass_fvtx_pt_S_ul->ProjectionX()->Draw();
  c1->cd(2);  c1->cd(2)->SetLogy();
  // mass_fvtx_pt_N_mix->ProjectionX()->Draw();
  mass_fvtx_pt_N_ul->ProjectionX()->Draw();

  // only for testing locally
  //char const *outfile = "test.root";  

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


  mass_fvtx_pt_2084_N_ul->Write();
  mass_fvtx_pt_2084_N_pp->Write();
  mass_fvtx_pt_2084_N_mm->Write();
  mass_fvtx_pt_2084_N_bg->Write();
  mass_fvtx_pt_2084_N_mix->Write();
  
  mass_fvtx_pt_2084_S_ul->Write();
  mass_fvtx_pt_2084_S_pp->Write();
  mass_fvtx_pt_2084_S_mm->Write();
  mass_fvtx_pt_2084_S_bg->Write();
  mass_fvtx_pt_2084_S_mix->Write();


  mass_fvtx_pt_4084_N_ul->Write();
  mass_fvtx_pt_4084_N_pp->Write();
  mass_fvtx_pt_4084_N_mm->Write();
  mass_fvtx_pt_4084_N_bg->Write();
  mass_fvtx_pt_4084_N_mix->Write();
  
  mass_fvtx_pt_4084_S_ul->Write();
  mass_fvtx_pt_4084_S_pp->Write();
  mass_fvtx_pt_4084_S_mm->Write();
  mass_fvtx_pt_4084_S_bg->Write();
  mass_fvtx_pt_4084_S_mix->Write();

 
}
