// run on condor to fill histograms first
// run this macro (Matt Durham's) to fit histograms using the method root fitslicesy
// add the cuts listed here to the analysis code and run on condor again


#include <THistPainter.h>
#include <TH2.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>


void sigmalizer_macro_durham()
{

  bool TTreeReader = false;


  //TFile *file_data = new TFile("CONDOR/red_16662/result.root");  
  //TFile *file_data = new TFile("Red_v0_local.root");
  //TFile *file_data = new TFile("Red_v0_16662_sanghoon_local_mixed_muid1D_cut.root");  
  //TFile *file_data = new TFile("CONDOR/fvtx_batch105108_sanghoon_N/result.root");
   //TFile *file_data = new TFile("CONDOR/YHcuts_sanghoon_reasstn_105108_job5/result.root");
  TFile *file_data = new TFile("CONDOR/fvtx_batch105108_sanghoon/result.root");
  //TFile *file_data = new TFile("CONDOR/red_16662/result.root");
  //  TFile *file_data = new TFile("Run15pAu_mass_pt_cent_5sigma_PPG188.root");

  TH2D *dphi_nhits_fvtx_Tr0 = new TH2D("dphi_nhits_fvtx_Tr0","dphi_fvtx_Tr0",150,0,15,2200,-110,110);  
  TH2D *dphi_nhits_fvtx_Tr1 = new TH2D("dphi_nhits_fvtx_Tr1","dphi_fvtx_Tr1",150,0,15,2200,-110,110);  

  // for rebinned histograms 

  TH2D *rDG0_N = new TH2D("rDG0_N","rDG0_N",37,0,15,50,0,30);  
  TH2D *rDG0_S = new TH2D("rDG0_S","rDG0_S",37,0,15,50,0,30);  

  TH2D *rDDG0_N = new TH2D("rDDG0_N","rDDG0_N",37,0,15,50,0,30);  
  TH2D *rDDG0_S = new TH2D("rDDG0_S","rDDG0_S",37,0,15,50,0,30);   

  TH2D *rdr_N = new TH2D("rdr_N","rdr_N",37,0,15,50,0,10);  
  TH2D *rdr_S = new TH2D("rdr_S","rdr_S",37,0,15,50,0,10);  

  TH2D *rdphi_N = new TH2D("rdphi_N","rdphi_N",37,0,15,100,-0.5,0.5);  
  TH2D *rdphi_S = new TH2D("rdphi_S","rdphi_S",37,0,15,100,-0.5,0.5);  

  TH2D *rdtheta_N = new TH2D("rdtheta_N","rdtheta_N",37,0,15,50,0,0.5);  
  TH2D *rdtheta_S = new TH2D("rdtheta_S","rdtheta_S",37,0,15,50,0,0.5);  

  //  to hold root file histograms 

  TH2D *dr_0;
  TH2D *dr_1;

  TH2D *dtheta_0;
  TH2D *dtheta_1;

  TH2D *dphi_0;
  TH2D *dphi_1;

  TH2D *DG0_0;
  TH2D *DG0_1;

  TH2D *DDG0_0;
  TH2D *DDG0_1;

  // histograms to be rebinned for fitting
  file_data->GetObject("dr_S_data",rdr_S);
  file_data->GetObject("dr_N_data",rdr_N);
     
  file_data->GetObject("dtheta_S_data",rdtheta_S);
  file_data->GetObject("dtheta_N_data",rdtheta_N);
     
  file_data->GetObject("dphi_S_data",rdphi_S);
  file_data->GetObject("dphi_N_data",rdphi_N);
     
  file_data->GetObject("DG0_S_data",rDG0_S);
  file_data->GetObject("DG0_N_data",rDG0_N);
       
  file_data->GetObject("DDG0_S_data",rDDG0_S);
  file_data->GetObject("DDG0_N_data",rDDG0_N);

 // histograms to remain unmodified 

  file_data->GetObject("dr_S_data",dr_0);
  file_data->GetObject("dr_N_data",dr_1);
     
  file_data->GetObject("dtheta_S_data",dtheta_0);
  file_data->GetObject("dtheta_N_data",dtheta_1);
     
  file_data->GetObject("dphi_S_data",dphi_0);
  file_data->GetObject("dphi_N_data",dphi_1);
     
  file_data->GetObject("DG0_S_data",DG0_0);
  file_data->GetObject("DG0_N_data",DG0_1);
       
  file_data->GetObject("DDG0_S_data",DDG0_0);
  file_data->GetObject("DDG0_N_data",DDG0_1);
   

  //////////////////////

  int Nsigma = 5;

 /////////////////////////////////// dr south
  rdr_S->Rebin2D(4,1);
  rdr_S->FitSlicesY(0,0,150);
  
  TH1D *dr_S_data_1 = (TH1D *)gDirectory->Get("dr_S_data_1"); 
  TH1D *dr_S_data_2 = (TH1D *)gDirectory->Get("dr_S_data_2");
  
  if(dr_S_data_1!=0 && dr_S_data_2!=0) 
    {
      dr_S_data_1->Draw();
      dr_S_data_2->Draw();
    }
  else 
    {
      cout << "Retrieval of histograms from current directory was unsuccessful!" << endl;
    }

  TH1D *hdr_S = (TH1D *)dr_S_data_1->Clone("hdr_S"); 
  hdr_S->Add(dr_S_data_2,Nsigma);
 
  TF1 *dr_S_fit = new TF1("dr_S_fit","[0] + [1]/(x*x)",0,15);
  hdr_S->Fit(dr_S_fit);

 /////////////////////////////////// dr north
  rdr_N->Rebin2D(4,1);
  rdr_N->FitSlicesY(0,0,150);
  
  TH1D *dr_N_data_1 = (TH1D *)gDirectory->Get("dr_N_data_1"); 
  TH1D *dr_N_data_2 = (TH1D *)gDirectory->Get("dr_N_data_2");
  
  TH1D *hdr_N = (TH1D *)dr_N_data_1->Clone("hdr_N"); 
  hdr_N->Add(dr_N_data_2,Nsigma);
 
  TF1 *dr_N_fit = new TF1("dr_N_fit","[0] + [1]/(x*x)",0,15);
  hdr_N->Fit(dr_N_fit);

 /////////////////////////////////// dtheta south
  rdtheta_S->Rebin2D(4,1);
  rdtheta_S->FitSlicesY(0,0,150);
  
  TH1D *dtheta_S_data_1 = (TH1D *)gDirectory->Get("dtheta_S_data_1"); 
  TH1D *dtheta_S_data_2 = (TH1D *)gDirectory->Get("dtheta_S_data_2");
  
  TH1D *hdtheta_S = (TH1D *)dtheta_S_data_1->Clone("hdtheta_S"); 
  hdtheta_S->Add(dtheta_S_data_2,Nsigma);
 
  TF1 *dtheta_S_fit = new TF1("dtheta_S_fit","[0] + [1]/(x*x)",0,15);
  hdtheta_S->Fit(dtheta_S_fit);

 /////////////////////////////////// dtheta north
  rdtheta_N->Rebin2D(4,1);
  rdtheta_N->FitSlicesY(0,0,150);
  
  TH1D *dtheta_N_data_1 = (TH1D *)gDirectory->Get("dtheta_N_data_1"); 
  TH1D *dtheta_N_data_2 = (TH1D *)gDirectory->Get("dtheta_N_data_2");
  
  TH1D *hdtheta_N = (TH1D *)dtheta_N_data_1->Clone("hdtheta_N"); 
  hdtheta_N->Add(dtheta_N_data_2,Nsigma);
 
  TF1 *dtheta_N_fit = new TF1("dtheta_N_fit","[0] + [1]/(x*x)",0,15);
  hdtheta_N->Fit(dtheta_N_fit);

  /////////////////////////////////// dphi south
  rdphi_S->Rebin2D(4,1);
  rdphi_S->FitSlicesY(0,0,150);
  
  TH1D *dphi_S_data_1 = (TH1D *)gDirectory->Get("dphi_S_data_1"); 
  TH1D *dphi_S_data_2 = (TH1D *)gDirectory->Get("dphi_S_data_2");
  
  TH1D *hdphi_S = (TH1D *)dphi_S_data_1->Clone("hdphi_S"); 
  hdphi_S->Add(dphi_S_data_2,Nsigma);
 
  TF1 *dphi_S_fit = new TF1("dphi_S_fit","[0] + [1]/(x*x)",0,15);
  hdphi_S->Fit(dphi_S_fit);

 /////////////////////////////////// dphi north
  rdphi_N->Rebin2D(4,1);
  rdphi_N->FitSlicesY(0,0,150);
  
  TH1D *dphi_N_data_1 = (TH1D *)gDirectory->Get("dphi_N_data_1"); 
  TH1D *dphi_N_data_2 = (TH1D *)gDirectory->Get("dphi_N_data_2");
  
  TH1D *hdphi_N = (TH1D *)dphi_N_data_1->Clone("hdphi_N"); 
  hdphi_N->Add(dphi_N_data_2,Nsigma);
 
  TF1 *dphi_N_fit = new TF1("dphi_N_fit","[0] + [1]/(x*x)",0,15);
  hdphi_N->Fit(dphi_N_fit);

  /////////////////////////////////// DG0 south
  rDG0_S->Rebin2D(4,1);
  rDG0_S->FitSlicesY(0,0,150);
  
  TH1D *DG0_S_data_1 = (TH1D *)gDirectory->Get("DG0_S_data_1"); 
  TH1D *DG0_S_data_2 = (TH1D *)gDirectory->Get("DG0_S_data_2");

  TH1D *hDG0_S = (TH1D *)DG0_S_data_1->Clone("hDG0_S"); 
  hDG0_S->Add(DG0_S_data_2,Nsigma);
 
  TF1 *DG0_S_fit = new TF1("DG0_S_fit","[0] + [1]/(x*x)",0,15);
  hDG0_S->Fit(DG0_S_fit);

  /////////////////////////////////// DG0 north
  rDG0_N->Rebin2D(4,1);
  rDG0_N->FitSlicesY(0,0,150);
  
  TH1D *DG0_N_data_1 = (TH1D *)gDirectory->Get("DG0_N_data_1"); 
  TH1D *DG0_N_data_2 = (TH1D *)gDirectory->Get("DG0_N_data_2");

  TH1D *hDG0_N = (TH1D *)DG0_N_data_1->Clone("hDG0_N"); 
  hDG0_N->Add(DG0_N_data_2,Nsigma);
 
  TF1 *DG0_N_fit = new TF1("DG0_N_fit","[0] + [1]/(x*x)",0,15);
  hDG0_N->Fit(DG0_N_fit);

  /////////////////////////////////// DDG0 south
  rDDG0_S->Rebin2D(4,1);
  rDDG0_S->FitSlicesY(0,0,150);
  
  TH1D *DDG0_S_data_1 = (TH1D *)gDirectory->Get("DDG0_S_data_1"); 
  TH1D *DDG0_S_data_2 = (TH1D *)gDirectory->Get("DDG0_S_data_2");

  TH1D *hDDG0_S = (TH1D *)DDG0_S_data_1->Clone("hDDG0_S"); 
  hDDG0_S->Add(DDG0_S_data_2,Nsigma);
 
  TF1 *DDG0_S_fit = new TF1("DDG0_S_fit","[0] + [1]/(x*x)",0,15);
  hDDG0_S->Fit(DDG0_S_fit);
  
  /////////////////////////////////// DG0 north
  rDDG0_N->Rebin2D(4,1);
  rDDG0_N->FitSlicesY(0,0,150);
  
  TH1D *DDG0_N_data_1 = (TH1D *)gDirectory->Get("DDG0_N_data_1"); 
  TH1D *DDG0_N_data_2 = (TH1D *)gDirectory->Get("DDG0_N_data_2");

  TH1D *hDDG0_N = (TH1D *)DDG0_N_data_1->Clone("hDDG0_N"); 
  hDDG0_N->Add(DDG0_N_data_2,Nsigma);
 
  TF1 *DDG0_N_fit = new TF1("DDG0_N_fit","[0] + [1]/(x*x)",0,15);
  hDDG0_N->Fit(DDG0_N_fit);
   
  //////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *dr_S = new TCanvas("dr_S","dr_S",1500,500);
  dr_S->Divide(3,1); 
  dr_S->cd(1);  dr_S_data_1->GetYaxis()->SetRangeUser(0,5);  
  dr_S_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dr_S_data_1->GetYaxis()->SetTitle("dr [cm]");
  dr_S_data_1->GetYaxis()->SetTitleOffset(1.4);
  dr_S_data_1->Draw();   

  dr_S->cd(2);  
  hdr_S->SetTitle("mean plus 5 sigma");  
  hdr_S->GetYaxis()->SetRangeUser(0,5);  
  hdr_S->GetXaxis()->SetTitle("p [GeV/c]");
  hdr_S->GetYaxis()->SetTitle("dr [cm]");
  hdr_S->GetYaxis()->SetTitleOffset(1.4);

  hdr_S->Draw();  
  dr_S->cd(3);  
  dr_0->Draw("colz"); 
  dr_0->GetXaxis()->SetTitle("p [GeV/c]");
  dr_0->GetYaxis()->SetTitle("dr [cm]");
  dr_0->GetYaxis()->SetTitleOffset(1.4);
  dr_S_fit->Draw("same"); 
  
  
  TCanvas *dr_N = new TCanvas("dr_N","dr_N",1500,500);
  dr_N->Divide(3,1); 
  dr_N->cd(1);  
  dr_N_data_1->GetYaxis()->SetRangeUser(0,5);  
  dr_N_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dr_N_data_1->GetYaxis()->SetTitle("dr [cm]");
  dr_N_data_1->GetYaxis()->SetTitleOffset(1.4);
  dr_N_data_1->Draw();

  dr_N->cd(2);  
  hdr_N->SetTitle("mean plus 5 sigma");  
  hdr_N->GetYaxis()->SetRangeUser(0,5);  
  hdr_N->GetXaxis()->SetTitle("p [GeV/c]");
  hdr_N->GetYaxis()->SetTitle("dr [cm]");
  hdr_N->GetYaxis()->SetTitleOffset(1.4);
  hdr_N->Draw();

  dr_N->cd(3);  
  dr_1->Draw("colz"); 
  dr_1->GetXaxis()->SetTitle("p [GeV/c]");
  dr_1->GetYaxis()->SetTitle("dr [cm]");
  dr_1->GetYaxis()->SetTitleOffset(1.4);
  dr_N_fit->Draw("same");
  
  TCanvas *dtheta_S = new TCanvas("dtheta_S","dtheta_S",1500,500);
  dtheta_S->Divide(3,1); 
  dtheta_S->cd(1);  
  dtheta_S_data_1->GetYaxis()->SetRangeUser(0,0.2);  
  dtheta_S_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dtheta_S_data_1->GetYaxis()->SetTitle("d#theta [deg]");
  dtheta_S_data_1->GetYaxis()->SetTitleOffset(1.4);
  dtheta_S_data_1->GetYaxis()->SetLabelSize(0.03);
  dtheta_S_data_1->Draw();

  dtheta_S->cd(2);  
  hdtheta_S->SetTitle("mean plus 5 sigma");  
  hdtheta_S->GetYaxis()->SetRangeUser(0,0.2); 
  hdtheta_S->GetXaxis()->SetTitle("p [GeV/c]");
  hdtheta_S->GetYaxis()->SetTitle("d#theta [deg]");
  hdtheta_S->GetYaxis()->SetTitleOffset(1.4);
  hdtheta_S->GetYaxis()->SetLabelSize(0.03);
  hdtheta_S->Draw();

  dtheta_S->cd(3);   
  dtheta_0->GetYaxis()->SetRangeUser(0,0.2); 
  dtheta_0->Draw("colz"); 
  dtheta_0->GetXaxis()->SetTitle("p [GeV/c]");
  dtheta_0->GetYaxis()->SetTitle("d#theta [deg]");
  dtheta_0->GetYaxis()->SetTitleOffset(1.4);
  dtheta_0->GetYaxis()->SetLabelSize(0.03);
  dtheta_S_fit->Draw("same");
  
  TCanvas *dtheta_N = new TCanvas("dtheta_N","dtheta_N",1500,500);
  dtheta_N->Divide(3,1); 
  dtheta_N->cd(1);  
  dtheta_N_data_1->GetYaxis()->SetRangeUser(0,0.2); 
  dtheta_N_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dtheta_N_data_1->GetYaxis()->SetTitle("d#theta [deg]");
  dtheta_N_data_1->GetYaxis()->SetTitleOffset(1.4);  
  dtheta_N_data_1->GetYaxis()->SetLabelSize(0.03);
  dtheta_N_data_1->Draw();

  dtheta_N->cd(2);  
  hdtheta_N->SetTitle("mean plus 5 sigma");  
  hdtheta_N->GetYaxis()->SetRangeUser(0,0.2);  
  hdtheta_N->GetXaxis()->SetTitle("p [GeV/c]");
  hdtheta_N->GetYaxis()->SetTitle("d#theta [deg]");
  hdtheta_N->GetYaxis()->SetTitleOffset(1.4);
  hdtheta_N->GetYaxis()->SetLabelSize(0.03);
  hdtheta_N->Draw();

  dtheta_N->cd(3);  
  dtheta_1->GetYaxis()->SetRangeUser(0,0.2); 
  dtheta_1->Draw("colz"); 
  dtheta_1->GetXaxis()->SetTitle("p [GeV/c]");
  dtheta_1->GetYaxis()->SetTitle("d#theta [deg]");
  dtheta_1->GetYaxis()->SetTitleOffset(1.4);
  dtheta_1->GetYaxis()->SetLabelSize(0.03);
  dtheta_N_fit->Draw("same");
  
  TCanvas *dphi_S = new TCanvas("dphi_S","dphi_S",1500,500);
  dphi_S->Divide(3,1); 
  dphi_S->cd(1);  
  dphi_S_data_1->GetYaxis()->SetRangeUser(-0.2,0.2); 
  dphi_S_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dphi_S_data_1->GetYaxis()->SetTitle("d#phi [deg]");
  dphi_S_data_1->GetYaxis()->SetTitleOffset(1.4);
  dphi_S_data_1->GetYaxis()->SetLabelSize(0.03);
  dphi_S_data_1->Draw();

  dphi_S->cd(2);  
  hdphi_S->SetTitle("mean plus 5 sigma");  
  hdphi_S->GetYaxis()->SetRangeUser(-0.2,0.2); 
  hdphi_S->GetXaxis()->SetTitle("p [GeV/c]");
  hdphi_S->GetYaxis()->SetTitle("d#phi [deg]");
  hdphi_S->GetYaxis()->SetTitleOffset(1.4); 
  hdphi_S->GetYaxis()->SetLabelSize(0.03);
  hdphi_S->Draw();

  dphi_S->cd(3); 
  dphi_0->GetYaxis()->SetRangeUser(-0.2,0.2);  
  dphi_0->Draw("colz"); 
  dphi_0->GetXaxis()->SetTitle("p [GeV/c]");
  dphi_0->GetYaxis()->SetTitle("d#phi [deg]");
  dphi_0->GetYaxis()->SetTitleOffset(1.4);
  dphi_0->GetYaxis()->SetLabelSize(0.03);
  dphi_S_fit->Draw("same");
  
  TCanvas *dphi_N = new TCanvas("dphi_N","dphi_N",1500,500);
  dphi_N->Divide(3,1); 
  dphi_N->cd(1);  
  dphi_N_data_1->GetYaxis()->SetRangeUser(-0.2,0.2); 
  dphi_N_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  dphi_N_data_1->GetYaxis()->SetTitle("d#phi [deg]");
  dphi_N_data_1->GetYaxis()->SetTitleOffset(1.4);
  dphi_N_data_1->GetYaxis()->SetLabelSize(0.03);
  dphi_N_data_1->Draw();

  dphi_N->cd(2);  
  hdphi_N->SetTitle("mean plus 5 sigma");  
  hdphi_N->GetYaxis()->SetRangeUser(-0.2,0.2);  
  hdphi_N->GetXaxis()->SetTitle("p [GeV/c]");
  hdphi_N->GetYaxis()->SetTitle("d#phi [deg]");
  hdphi_N->GetYaxis()->SetTitleOffset(1.4); 
  hdphi_N->GetYaxis()->SetLabelSize(0.03);
  hdphi_N->Draw();

  dphi_N->cd(3);  
  dphi_1->GetYaxis()->SetRangeUser(-0.2,0.2);  
  dphi_1->Draw("colz"); 
  dphi_1->GetXaxis()->SetTitle("p [GeV/c]");
  dphi_1->GetYaxis()->SetTitle("d#phi [deg]");
  dphi_1->GetYaxis()->SetTitleOffset(1.4);
  dphi_1->GetYaxis()->SetLabelSize(0.03);
  dphi_N_fit->Draw("same");
  
  TCanvas *DG0_S = new TCanvas("DG0_S","DG0_S",1500,500);
  DG0_S->Divide(3,1); 
  DG0_S->cd(1);  
  DG0_S_data_1->GetYaxis()->SetRangeUser(0,30); 
  DG0_S_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  DG0_S_data_1->GetYaxis()->SetTitle("DG0 [cm]");
  DG0_S_data_1->GetYaxis()->SetTitleOffset(1.4);
  DG0_S_data_1->Draw();

  DG0_S->cd(2);  
  hDG0_S->SetTitle("mean plus 5 sigma");  
  hDG0_S->GetYaxis()->SetRangeUser(0,30); 
  hDG0_S->GetXaxis()->SetTitle("p [GeV/c]");
  hDG0_S->GetYaxis()->SetTitle("DG0 [cm]");
  hDG0_S->GetYaxis()->SetTitleOffset(1.4);
  hDG0_S->Draw();

  DG0_S->cd(3);  
  DG0_0->GetYaxis()->SetRangeUser(0,30);  
  DG0_0->Draw("colz"); 
  DG0_0->GetXaxis()->SetTitle("p [GeV/c]");
  DG0_0->GetYaxis()->SetTitle("DG0 [cm]");
  DG0_0->GetYaxis()->SetTitleOffset(1.4);
  DG0_S_fit->Draw("same");
  
  TCanvas *DG0_N = new TCanvas("DG0_N","DG0_N",1500,500);
  DG0_N->Divide(3,1); 
  DG0_N->cd(1);  
  DG0_N_data_1->GetYaxis()->SetRangeUser(0,30);  
  DG0_N_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  DG0_N_data_1->GetYaxis()->SetTitle("DG0 [cm]");
  DG0_N_data_1->GetYaxis()->SetTitleOffset(1.4);
  DG0_N_data_1->Draw();

  DG0_N->cd(2);  
  hDG0_N->SetTitle("mean plus 5 sigma");  
  hDG0_N->GetYaxis()->SetRangeUser(0,30); 
  hDG0_N->GetXaxis()->SetTitle("p [GeV/c]");
  hDG0_N->GetYaxis()->SetTitle("DG0 [cm]");
  hDG0_N->GetYaxis()->SetTitleOffset(1.4); 
  hDG0_N->Draw();

  DG0_N->cd(3);  
  DG0_1->GetYaxis()->SetRangeUser(0,30);  
  DG0_1->Draw("colz"); 
  DG0_1->GetXaxis()->SetTitle("p [GeV/c]");
  DG0_1->GetYaxis()->SetTitle("DG0 [cm]");
  DG0_1->GetYaxis()->SetTitleOffset(1.4);
  DG0_N_fit->Draw("same");
  
  TCanvas *DDG0_S = new TCanvas("DDG0_S","DDG0_S",1500,500);
  DDG0_S->Divide(3,1); 
  DDG0_S->cd(1); 
  DDG0_S_data_1->GetYaxis()->SetRangeUser(0,30);  
  DDG0_S_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  DDG0_S_data_1->GetYaxis()->SetTitle("DDG0 [deg]");
  DDG0_S_data_1->GetYaxis()->SetTitleOffset(1.4);
  DDG0_S_data_1->Draw();

  DDG0_S->cd(2);  
  hDDG0_S->SetTitle("mean plus 5 sigma");  
  hDDG0_S->GetYaxis()->SetRangeUser(0,30);  
  hDDG0_S->GetXaxis()->SetTitle("p [GeV/c]");
  hDDG0_S->GetYaxis()->SetTitle("DDG0 [deg]");
  hDDG0_S->GetYaxis()->SetTitleOffset(1.4);
  hDDG0_S->Draw();

  DDG0_S->cd(3);  
  DDG0_0->GetYaxis()->SetRangeUser(0,30);  
  DDG0_0->Draw("colz"); 
  DDG0_0->GetXaxis()->SetTitle("p [GeV/c]");
  DDG0_0->GetYaxis()->SetTitle("DDG0 [deg]");
  DDG0_0->GetYaxis()->SetTitleOffset(1.4);
  DDG0_S_fit->Draw("same");
  
  TCanvas *DDG0_N = new TCanvas("DDG0_N","DDG0_N",1500,500);
  DDG0_N->Divide(3,1); 
  DDG0_N->cd(1);  
  DDG0_N_data_1->GetYaxis()->SetRangeUser(0,30); 
  DDG0_N_data_1->GetXaxis()->SetTitle("p [GeV/c]");
  DDG0_N_data_1->GetYaxis()->SetTitle("DDG0 [deg]");
  DDG0_N_data_1->GetYaxis()->SetTitleOffset(1.4);
  DDG0_N_data_1->Draw();

  DDG0_N->cd(2);  
  hDDG0_N->SetTitle("mean plus 5 sigma");  
  hDDG0_N->GetYaxis()->SetRangeUser(0,30);  
  hDDG0_N->GetXaxis()->SetTitle("p [GeV/c]");
  hDDG0_N->GetYaxis()->SetTitle("DDG0 [deg]");
  hDDG0_N->GetYaxis()->SetTitleOffset(1.4);
  hDDG0_N->Draw();

  DDG0_N->cd(3);  
  DDG0_1->GetYaxis()->SetRangeUser(0,30);  
  DDG0_1->Draw("colz"); 
  DDG0_1->GetXaxis()->SetTitle("p [GeV/c]");
  DDG0_1->GetYaxis()->SetTitle("DDG0 [deg]");
  DDG0_1->GetYaxis()->SetTitleOffset(1.4);
  DDG0_N_fit->Draw("same");
  


  if(!TTreeReader)
    {
  cout<<"string dr_cut_N = \"Tr0_dr_fvtx<("<<dr_N_fit->GetParameter(0)<<"+"<<dr_N_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dr_fvtx<("<<dr_N_fit->GetParameter(0)<<"+"<<dr_N_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

  cout<<"string dr_cut_S = \"Tr0_dr_fvtx<("<<dr_S_fit->GetParameter(0)<<"+"<<dr_S_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dr_fvtx<("<<dr_S_fit->GetParameter(0)<<"+"<<dr_S_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

 cout<<"string dtheta_cut_N = \"Tr0_dtheta_fvtx<("<<dtheta_N_fit->GetParameter(0)<<"+"<<dtheta_N_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dtheta_fvtx<("<<dtheta_N_fit->GetParameter(0)<<"+"<<dtheta_N_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

cout<<"string dtheta_cut_S = \"Tr0_dtheta_fvtx<("<<dtheta_S_fit->GetParameter(0)<<"+"<<dtheta_S_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dtheta_fvtx<("<<dtheta_S_fit->GetParameter(0)<<"+"<<dtheta_S_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

 cout<<"string dphi_cut_N = \"Tr0_dphi_fvtx<("<<dphi_N_fit->GetParameter(0)<<"+"<<dphi_N_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dphi_fvtx<("<<dphi_N_fit->GetParameter(0)<<"+"<<dphi_N_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

cout<<"string dphi_cut_S = \"Tr0_dphi_fvtx<("<<dphi_S_fit->GetParameter(0)<<"+"<<dphi_S_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_dphi_fvtx<("<<dphi_S_fit->GetParameter(0)<<"+"<<dphi_S_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;


  

cout<<"string DG0_cut_N = \"Tr0_DG0<("<<DG0_N_fit->GetParameter(0)<<"+"<<DG0_N_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_DG0<("<<DG0_N_fit->GetParameter(0)<<"+"<<DG0_N_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

cout<<"string DG0_cut_S = \"Tr0_DG0<("<<DG0_S_fit->GetParameter(0)<<"+"<<DG0_S_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_DG0<("<<DG0_S_fit->GetParameter(0)<<"+"<<DG0_S_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

cout<<"string DDG0_cut_N = \"Tr0_DDG0<("<<DDG0_N_fit->GetParameter(0)<<"+"<<DDG0_N_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_DDG0<("<<DDG0_N_fit->GetParameter(0)<<"+"<<DDG0_N_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

cout<<"string DDG0_cut_S = \"Tr0_DDG0<("<<DDG0_S_fit->GetParameter(0)<<"+"<<DDG0_S_fit->GetParameter(1)<<"/(Tr0_px_fvtxmutr * Tr0_px_fvtxmutr + Tr0_py_fvtxmutr*Tr0_py_fvtxmutr + Tr0_pz_fvtxmutr*Tr0_pz_fvtxmutr))\ && Tr1_DDG0<("<<DDG0_S_fit->GetParameter(0)<<"+"<<DDG0_S_fit->GetParameter(1)<<"/(Tr1_px_fvtxmutr * Tr1_px_fvtxmutr + Tr1_py_fvtxmutr*Tr1_py_fvtxmutr + Tr1_pz_fvtxmutr*Tr1_pz_fvtxmutr))\";"<<endl;

    }
  ///////////////////////////////////////////////////////////////////////////////////////////
  if(TTreeReader)
    {
cout<<" dr_cut_N = Tr0_dr_fvtx[i]<("<<dr_N_fit->GetParameter(0)<<"+"<<dr_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dr_fvtx[i]<("<<dr_N_fit->GetParameter(0)<<"+"<<dr_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

  cout<<" dr_cut_S = Tr0_dr_fvtx[i]<("<<dr_S_fit->GetParameter(0)<<"+"<<dr_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dr_fvtx[i]<("<<dr_S_fit->GetParameter(0)<<"+"<<dr_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

 cout<<" dtheta_cut_N = Tr0_dtheta_fvtx[i]<("<<dtheta_N_fit->GetParameter(0)<<"+"<<dtheta_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dtheta_fvtx[i]<("<<dtheta_N_fit->GetParameter(0)<<"+"<<dtheta_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

cout<<" dtheta_cut_S = Tr0_dtheta_fvtx[i]<("<<dtheta_S_fit->GetParameter(0)<<"+"<<dtheta_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dtheta_fvtx[i]<("<<dtheta_S_fit->GetParameter(0)<<"+"<<dtheta_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

 cout<<" dphi_cut_N = Tr0_dphi_fvtx[i]<("<<dphi_N_fit->GetParameter(0)<<"+"<<dphi_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dphi_fvtx[i]<("<<dphi_N_fit->GetParameter(0)<<"+"<<dphi_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

cout<<" dphi_cut_S = Tr0_dphi_fvtx[i]<("<<dphi_S_fit->GetParameter(0)<<"+"<<dphi_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && Tr1_dphi_fvtx[i]<("<<dphi_S_fit->GetParameter(0)<<"+"<<dphi_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;


  

cout<<"DG0_N = DG0_Tr0[i]<("<<DG0_N_fit->GetParameter(0)<<"+"<<DG0_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && DG0_Tr1[i]<("<<DG0_N_fit->GetParameter(0)<<"+"<<DG0_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

 cout<<" DG0_S = DG0_Tr0[i]<("<<DG0_S_fit->GetParameter(0)<<"+"<<DG0_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && DG0_Tr1[i]<("<<DG0_S_fit->GetParameter(0)<<"+"<<DG0_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

cout<<" DDG0_N = DDG0_Tr0[i]<("<<DDG0_N_fit->GetParameter(0)<<"+"<<DDG0_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && DDG0_Tr1[i]<("<<DDG0_N_fit->GetParameter(0)<<"+"<<DDG0_N_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

cout<<" DDG0_S = DDG0_Tr0[i]<("<<DDG0_S_fit->GetParameter(0)<<"+"<<DDG0_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr0[i] * px_fvtxmutr_Tr0[i] + py_fvtxmutr_Tr0[i]*py_fvtxmutr_Tr0[i] + pz_fvtxmutr_Tr0[i]*pz_fvtxmutr_Tr0[i])) && DDG0_Tr1[i]<("<<DDG0_S_fit->GetParameter(0)<<"+"<<DDG0_S_fit->GetParameter(1)<<"/(px_fvtxmutr_Tr1[i] * px_fvtxmutr_Tr1[i] + py_fvtxmutr_Tr1[i]*py_fvtxmutr_Tr1[i] + pz_fvtxmutr_Tr1[i]*pz_fvtxmutr_Tr1[i]));"<<endl;

    }



  
 

}
