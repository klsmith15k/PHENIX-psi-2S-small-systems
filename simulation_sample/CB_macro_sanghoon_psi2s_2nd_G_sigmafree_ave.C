// LS bgnorm is applied (line 259)


#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <Math/MinimizerOptions.h>
#include <fstream>
#include <iostream>

using namespace std;

TF1 *combg_fit;

Double_t CBcalc_LL(Double_t *xx, Double_t *par) // CB sigma scaled by 0.15
{
 
  


  double f;
  double x = xx[0];
  
  int s=81;

  const double sigma_step = 0.001;
  double fixed = 0.100 + s*sigma_step;


  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:
  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double N_2nd = par[5];
  double X_mean_2nd = x_mean;
  double Sigma_2nd = par[6];
  
  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

    double gauss2 = N_2nd*exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(Sigma_2nd,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2;
    }
 //Like sign combinatoric background fit
 // double bg12 = combg_fit->Eval(x);

  // Add total fit with like sign background fit
  //f += bg12;  // comb bg fit results are added to psi2s CB, psi(2s)CB and corr bg 



   if(x == 2.60)
    {
      cout << "==================== " << endl;
      cout << "J/psi 2nd gaussian parameters:" << endl;
      cout << "==================== " << endl;
      cout <<" N_2nd : " << N_2nd << endl;
      cout << " X_mean_2nd : " << X_mean_2nd << endl;
      cout << " Sigma 2nd : " << sigma << endl;
      cout << " N CB : " << N << endl;
     
      if(Sigma_2nd == fixed)
	cout << " --> 1. gaussian width fixed to " << fixed << " " << endl;

      if(X_mean_2nd == x_mean)
	cout << " --> 2. 2nd gaussian mean equal to Psi2s CB mean " << endl;

    }

  return f;
}

Double_t CBcalc_psi2s(Double_t *xx, Double_t *par)  // CB calc only
{
  double f;
  double x = xx[0];

  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:
  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];
  double N_2nd = par[5];

  double X_mean_2nd = par[6];
  double Sigma_2nd = par[7];
 
  int s=81;

  const double sigma_step = 0.001;
  double fixed = 0.100 + s*sigma_step;

  double gauss2 = N_2nd*exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(Sigma_2nd,2)));

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2;
    }

  return f;
}


Double_t Gauss2(Double_t *xx, Double_t *par)  
{
  double f;
  double x = xx[0];

  double N_2nd = par[0];
  double x_mean_2nd = par[1];
  double sigma_2nd = par[2];
  
  int s=81;

  const double sigma_step = 0.001;
  double fixed = 0.100 + s*sigma_step;


  f = N_2nd*exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(sigma_2nd,2)));

  return f;

}

void CB_macro_sanghoon_psi2s_2nd_G_sigmafree_ave()
{
  //////////////////////////////////////////////
  bool north_arm = true;
  //bool north_arm = false;
  ///////////////////////////////////////////////

  TH1 *sim_mass;
  TH1 *recomassBG;
  TH1 *recomassBG1;
  TH1 *recomassBG2;
  bool sanghoon_files = true;
  bool embed = true;
  TFile *file1S;  
  bool chisquare_write = false;
  int bin;
 

  int s=81;

  const double sigma_step = 0.001;
  double fixed = 0.100 + s*sigma_step;

  //cout << "Enter the bin number " << endl;//  (1,2,3,4,5) = rap1,rap2,rap3,rap4,MB
  //cin >> bin;
  bin = 7;

  int com_bg;
  // cout << "Enter combinatoric BG initial parameter set (0-13) " << endl;
  // cin >> com_bg;

  if(north_arm)
    com_bg = 0;
  else
    com_bg = 6;

  if(north_arm && embed)
    {
       if(sanghoon_files)
	 file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_psi2s_newlib_trig_fid188228_nochi2fvtx_trigZ_N.root");
      else
      	file1S = TFile::Open("my_test_simulations/mass_sngdbl_jpsi_embedC_N.root");
      
      cout << " north arm embed " << endl;
    }
  if(!north_arm && embed)
    {
    if(sanghoon_files)
      file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_psi2s_newlib_trig_fid188228_nochi2fvtx_trigZ_S.root");
    else
     file1S = TFile::Open("my_test_simulations/mass_sngdbl_jpsi_dimu_02_embedC_S.root");
      
      cout << " south arm embed " << endl;
    }

  TCanvas *psi2s_sim = new TCanvas("psi2s_sim","psi2s_sim",5,5,800,600);
  psi2s_sim->SetLogy();
 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  

  char ulname[500];
  char mmname[500];
  char ppname[500];
    
  sprintf(ulname,"ul_%i",bin);
  cout << "Open " << ulname << endl;
  sim_mass = (TH1 *)file1S->Get(ulname);
  if(!sim_mass)
    cout << "Failed to get " << ulname << endl; 
 
  double temp[60] = {0};

  for(int i = 0; i < 60; i++)
    {
      //temp[i] = sim_mass->GetBinContent(i+1)/12800;
      temp[i] = sim_mass->GetBinContent(i+1)/1;
      sim_mass->SetBinContent(i+1, temp[i]);

    }
  double f0,f1,f2,f3,f4,f5,f6,f7 = 0;
  double Npsi2s = 0;
  double bin_width = sim_mass->GetBinWidth(1);
  int binl = sim_mass->FindBin(2.6); 
  int binh = sim_mass->FindBin(3.6); 
  double renorm = 1.0/bin_width; 
   
  sprintf(mmname,"mm_%i",bin);
  cout << "Open " << mmname << endl;
  recomassBG1 = (TH1 *)file1S->Get(mmname);

  if(!recomassBG1)
    cout << "Failed to get " << mmname << endl; 

  sprintf(ppname,"pp_%i",bin);
  cout << "Open " << ppname << endl;
  recomassBG2 = (TH1 *)file1S->Get(ppname);
  
  double dNpp, dNmm;
  double Npp = recomassBG2->IntegralAndError(recomassBG2->GetXaxis()->FindBin(2), recomassBG2->GetXaxis()->FindBin(5),dNpp,"");
  double Nmm = recomassBG1->IntegralAndError(recomassBG1->GetXaxis()->FindBin(2), recomassBG1->GetXaxis()->FindBin(5),dNmm,"");
  double bgnorm = 2*sqrt(Npp*Nmm)/(Npp+Nmm);
  cout << "bgnorm " << bgnorm << " Npp " << Npp << " Nmm " << Nmm << endl;

  recomassBG = (TH1D*)recomassBG1->Clone("recomassBG");  
  recomassBG->SetName("recomassBG"); 
  recomassBG->SetTitle("recomassBG"); 
  recomassBG->SetLineColor(kBlue+4);
  for(int i = 0; i < recomassBG->GetNbinsX(); i++)
    {
      recomassBG->SetBinContent(i+1,(recomassBG2->GetBinContent(i+1) + recomassBG1->GetBinContent(i+1))*bgnorm ); 
      recomassBG->SetBinError(i+1,(recomassBG2->GetBinError(i+1) + recomassBG1->GetBinError(i+1))*bgnorm);
    }
  // cout <<  recomassBG->GetXaxis()->GetNbins() << endl;



  double coeff;
  
  ////////////////////////////////

double m0,m1,m2,m3,m4;

  if(com_bg == 0)
    {
      m0 = 0.1;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 10;
    }
  if(com_bg == 1)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 10;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 2)   
    {
      m0 = 0.01;
      m1 = 0.001;
      m2 = 100;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 3)   
    {
      m0 = 0.1;
      m1 = 0.1;      
      m2 = 10;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 4)    
    {
      m0 = 0.1;
      m1 = 0.01;      
      m2 = 0.01;
      m3 = 10;
      m4 = 10;
    }
   if(com_bg == 5)  
    {
      m0 = -0.001;
      m1 = -0.01;
      m2 = 10; 
      m3 = 1; 
      m4 = 1;
    }
  
  if(com_bg == 6)
    {
      m0 = 0.01;
      m1 = 0.1;
      m2 = 1; 
      m3 = 10; 
      m4 = 10;
    }
 
 if(com_bg == 7)
    {
      m0 = -0.01;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 1;
    }

 if(com_bg == 8)
    {
      m0 = 0.01;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 1;
    }
  if(com_bg == 9)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 100;
      m3 = 10;
      m4 = 10;
    }
  if(com_bg == 10)
    {
      m0 = 0.1;
      m1 = -0.01;      
      m2 = 1;
      m3 = 10;
      m4 = 10;
    }
 if(com_bg == 11)
    {
      m0 = 10;
      m1 = 10;
      m2 = 10;
      m3 = 100; 
      m4 = 0.01;
    }
 if(com_bg == 12)
    {
      m0 = 0.01;
      m1 = 10;
      m2 = 10;
      m3 = 100; 
      m4 = 0.01;
    }
  if(com_bg == 13)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 1;
      m3 = 10;
      m4 = 10;
    }

ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000); 
 
  combg_fit = new TF1("combg_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);
  combg_fit->SetParameter(0, m0);  
  combg_fit->SetParameter(1, m1);
  combg_fit->SetParameter(2, m2); 
  combg_fit->SetParameter(3, m3);
  combg_fit->SetParameter(4, m4);
  combg_fit->SetLineColor(kViolet+3);
 
  combg_fit->SetLineStyle(2);

  
 
  //recomassBG->Fit(combg_fit,"LL","",2,5); 
 
  
  ////////////////////////////////

  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 

  double N_start = 0;
 
 

  /////////////////////////////////////////////////////////////////////////////////
  coeff = 1.5;

  N_start = coeff*sim_mass->Integral(binl,binh) / renorm; 

  if(north_arm && embed)
    {
     
      f0 = 1;
      f1 = 5;
      f2 = 3.7;
      f3 = 0.11; 
      f4 = N_start; 
      f5 = 10; // N_2nd
      f6 = 0.215;
     
    }
  /////////////////////////////////////////////////////////////////////////////////

  coeff = 1.5;
  // no prob cut
  if(!north_arm && embed)
    {
       f0 = 1; 
      f1 = 10;
      f2 = 3.7;
      f3 = 0.12;
      f4 = N_start;
      f5 = 10; 
      f6 = 0.215; 
    }



 // TF1 *total_fit = new TF1("total_fit",CBcalc_LL,2.0,4.2,7);
 TF1 *total_fit = new TF1("total_fit",CBcalc_LL,2.0,5.0,7);
  total_fit->SetLineColor(kBlue+1);
  

  // ///////////////////////////////////////////////////////////////////// previosu results
     
   if(north_arm)
     {
       total_fit->SetParameter(0, f0);
       total_fit->SetParLimits(0, 0, 2);
       total_fit->SetParameter(1, f1);
       total_fit->SetParLimits(1, 0, 10);
     }
   else
     {
       total_fit->SetParameter(0, f0);
       total_fit->SetParLimits(0, 0, 2);
       total_fit->SetParameter(1, f1);
       total_fit->SetParLimits(1, 0, 10);
     }
   
  //  /////////////////////////////////////////////////////////////////////

   // test with no limits on tail --> maybe only necessary in data

 total_fit->SetParameter(0, f0);
 total_fit->SetParameter(1, f1);
   total_fit->SetParameter(2, f2);
  total_fit->SetParameter(3, f3);
   total_fit->SetParameter(4, f4);
   total_fit->SetParameter(5, f5);
    total_fit->FixParameter(6, f6);
   
 
   if(north_arm)
     sim_mass->Fit(total_fit,"","",2.0,4.4); 
   else
     sim_mass->Fit(total_fit,"","",2.0,4.4); 

  cout << "par 0 error: " << (  total_fit->GetParError(0) / total_fit->GetParameter(0) ) * 100 << endl;
  cout << "par 1 error: " << (  total_fit->GetParError(1) / total_fit->GetParameter(1) ) * 100 << endl;
  cout << "par 2 error: " << (  total_fit->GetParError(2) / total_fit->GetParameter(2) ) * 100 << endl;
  cout << "par 3 error: " << (  total_fit->GetParError(3) / total_fit->GetParameter(3) ) * 100 << endl;
  cout << "par 4 error: " << (  total_fit->GetParError(4) / total_fit->GetParameter(4) ) * 100 << endl;
  cout << "par 5 error: " << (  total_fit->GetParError(5) / total_fit->GetParameter(5) ) * 100 << endl;
  cout << "par 5 error: " << (  total_fit->GetParError(6) / total_fit->GetParameter(6) ) * 100 << endl;
  




  sim_mass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
  sim_mass->GetYaxis()->SetTitle("raw counts/(50 MeV/c^{2})");
  sim_mass->GetXaxis()->SetTitleOffset(1.25);
  sim_mass->Draw();
  sim_mass->GetXaxis()->SetRangeUser(2,5);
  // sim_mass->Fit(total_fit, "LS");
  total_fit->Draw("same");
  // combg_fit->Draw("same");

  
 

 double X_mean_2nd;

  if(north_arm)
    X_mean_2nd = total_fit->GetParameter(2);
  else
    X_mean_2nd = total_fit->GetParameter(2);

  
  TF1 *psi2s_fit;

  if(north_arm)
    psi2s_fit  = new TF1("psi2s_fit",CBcalc_psi2s,2.0,4.4, 8);  
  else
    psi2s_fit  = new TF1("psi2s_fit",CBcalc_psi2s,2.0,4.4, 8);  

  psi2s_fit->SetLineColor(kRed+1);
  psi2s_fit->SetParameter(0, total_fit->GetParameter(0));
  psi2s_fit->SetParError(0, total_fit->GetParError(0));
  psi2s_fit->SetParameter(1, total_fit->GetParameter(1));
  psi2s_fit->SetParError(1, total_fit->GetParError(1));
  psi2s_fit->SetParameter(2, total_fit->GetParameter(2));
  psi2s_fit->SetParError(2, total_fit->GetParError(2));
  psi2s_fit->SetParameter(3, total_fit->GetParameter(3));
  psi2s_fit->SetParError(3, total_fit->GetParError(3));
  psi2s_fit->SetParameter(4, total_fit->GetParameter(4));
  psi2s_fit->SetParError(4, total_fit->GetParError(4));
  psi2s_fit->SetParameter(5, total_fit->GetParameter(5));
  psi2s_fit->SetParError(5, total_fit->GetParError(5));
  psi2s_fit->FixParameter(6,X_mean_2nd);
  psi2s_fit->SetParameter(7, total_fit->GetParameter(6));
  psi2s_fit->SetParError(7, total_fit->GetParError(6));

    

  TF1 *gauss2_fit = new TF1("gauss2_fit",Gauss2,2.0,5.0,3); // gauss2_psi2s
  gauss2_fit->SetParameter(0,total_fit->GetParameter(5));
  gauss2_fit->SetParError(0,total_fit->GetParError(5));
  gauss2_fit->SetParameter(1,total_fit->GetParameter(2));
  gauss2_fit->SetParError(1,total_fit->GetParError(2));
  gauss2_fit->SetParameter(2,total_fit->GetParameter(6));
  gauss2_fit->SetParError(2, total_fit->GetParError(6));

  gauss2_fit->SetLineColor(kViolet+2);
  gauss2_fit->SetLineStyle(10);
   
   
  // if(north_arm)
  //   {
  //     psi2s_fit->SetParLimits(0,0,1);
  //     psi2s_fit->SetParLimits(1,0,100);
  //   }
  // else
  //   {
  //     psi2s_fit->SetParLimits(0,0,100);
  //     psi2s_fit->SetParLimits(1,0,1000);
  //   }

      
  // Double_t psi2spar[6] = {0};
      
  // for(int i = 0;i < 6; i++)
  //   {
  //     psi2spar[i] = total_fit->GetParameter(i); 
  //   }
  
  
  double err_Npsi2s = 0;
  Double_t width = sim_mass->GetBinWidth(1);
  //  cout << "binwidth: " << width << endl;
  Npsi2s = psi2s_fit->Integral(2.0, 5.0)/width;  //  psi2s counts are calcualted from both the 2nd gaussian and the CB fit
  //err_Npsi2s = psi2s_fit->IntegralError(2.0,5,psi2spar,psi2s_array)/width;

  double chisquare;
  double ndf;

  gauss2_fit->Draw("same");
  psi2s_fit->Draw("same");
  //recomassBG->Draw("same");

  chisquare = total_fit->GetChisquare();
  ndf =  total_fit->GetNDF();
 
  
  TLatex l0;
  l0.SetTextSize(0.06);
  l0.SetTextAlign(13);
  l0.SetTextColor(4);

  char text4[100];
  if(north_arm && embed)
     sprintf(text4,"p + p North #psi(2S), Embedded");
  if(!north_arm && embed)
    sprintf(text4,"p + p South #psi(2S), Embedded");
 if(north_arm && !embed)
     sprintf(text4,"p + p North #psi(2S), No Embed");
  if(!north_arm && !embed)
    sprintf(text4,"p + p South #psi(2S), No Embed");

  l0.SetTextAlign(12);
  l0.DrawLatexNDC(0.4, 0.93, text4); //4.4,150
 
  Char_t message[200];
  Char_t message2[200];
 
  sprintf(message,"#chi^{2}/NDF = %.1f / %.d",total_fit->GetChisquare(),total_fit->GetNDF());
  sprintf(message2,"Psi2s Counts =  %.0f +/- %.2f", Npsi2s, sqrt(Npsi2s));
  
  TPaveText *mytext = new TPaveText(0.5,0.8,0.88,0.7,"NDC"); // x0,y0,x1,y1
  mytext->SetTextSize(0.035);
  mytext->SetFillColor(0);
  mytext->SetTextAlign(12);
  mytext->AddText(message);
  mytext->AddText(message2);
  mytext->Draw();


 TLatex l2;
  l2.SetTextSize(0.05);
  l2.SetTextAlign(13);
  l2.SetTextColor(1);

  
  char text3[100];
  if(bin == 1)
    {
      if(north_arm == true)
	sprintf(text3,"1.2 < y < 1.45");
      if(north_arm == false)
	sprintf(text3,"-1.45 < y < -1.2");
        }
     
  if(bin == 2)
    {
      if(north_arm == true)
	sprintf(text3,"1.45 < y < 1.7");
      if(north_arm == false)
	sprintf(text3,"-1.7 < y < -1.45");
    }
  if(bin == 3)
    {
      if(north_arm == true)
	sprintf(text3,"1.7 < y < 1.95");
      if(north_arm == false)
	sprintf(text3,"-1.95 < y < -1.7");
    }
 if(bin == 4)
    {
      if(north_arm == true)
	sprintf(text3,"1.95 < y < 2.2");
      if(north_arm == false)
	sprintf(text3,"-2.2 < y < -1.95");
    }
 if(bin == 7)
    {
      if(north_arm == true)
	sprintf(text3,"p_{T} Integrated, Sng+Dbl Tracks");
      if(north_arm == false)
	sprintf(text3,"p_{T} Integrated, Sng+Dbl Tracks");
    };
       
  l2.SetTextAlign(12);
  l2.SetTextSize(0.04);
  l2.DrawLatexNDC(0.63, 0.84, text3);
  
  char text14[100];
  sprintf(text14,"No Prob Cut");
  l2.SetTextAlign(12);
  l2.SetTextSize(0.035);
  // l2.DrawLatexNDC(0.13, 0.69, text14);
  
  char name200[500]; 

  if(chisquare_write)
    {
    
      if(north_arm)
	sprintf(name200,"chisquare_2nd_gaussian/sigmafree/north/N_gaussian_%i.txt", s);	 
      else
	sprintf(name200,"chisquare_2nd_gaussian/sigmafree/south/S_gaussian_%i.txt", s);	  	
     
      // whether north arm or south arm, write to unique filename  
      std::fstream Run15pp_width(name200,std::ofstream::out); 
     Run15pp_width <<  Npsi2s << " " <<  Npsi2s << " " << Npsi2s << " " << err_Npsi2s << " " << Npsi2s << " " << chisquare << " " << chisquare/ndf << " " << fixed << " " << endl;
      Run15pp_width.close();
   

    } // if width_write
  ///////////////////////////////////////

double free=total_fit->GetParameter(6);
double sigma =  total_fit->GetParameter(3);
double mean = total_fit->GetParameter(2);

 cout << "N_2nd/N: " <<  total_fit->GetParameter(5)/total_fit->GetParameter(4) << endl;
 cout << "===========" << endl;
 cout <<  "N_2nd: " << total_fit->GetParameter(5) << endl;  /// ratio numerator
 cout << "N: " << total_fit->GetParameter(4) << endl;   // ratio denominator


  TLatex l3;
  l3.SetTextSize(0.06);
  l3.SetTextAlign(13);
  l3.SetTextColor(4);

 

  char text5[100];
  char text6[100];
  char text7[100];
  char text8[100];
  l3.SetTextColor(1);
  l3.SetTextSize(0.035);
  sprintf(text5,"#sigma_{2nd} =  %.3f ",free);
  sprintf(text6,"CB #bar{x} =  %.5f ", mean);
  sprintf(text7,"Double Tracks");
  sprintf(text8,"CB #sigma =  %.5f ", sigma);
  l3.SetTextAlign(12);
  l3.DrawLatexNDC(0.13, 0.84, text5); //4.4,150
  // l3.DrawLatexNDC(0.13, 0.74, text7); //4.4,150
  l3.DrawLatexNDC(0.13, 0.79, text6); //4.4,150
  l3.DrawLatexNDC(0.13, 0.74, text8); //4.4,150

  if(north_arm)
    cout << "North int gauss / psi2s : " << gauss2_fit->Integral(2,5) / psi2s_fit->Integral(2,5) << ", for sigma_2nd = " << fixed << endl;
  else
    cout << "South int gauss / psi2s : " << gauss2_fit->Integral(2,5) / psi2s_fit->Integral(2,5) << ", for sigma_2nd = " << fixed << endl;
}


