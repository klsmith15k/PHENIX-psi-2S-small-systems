
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
#include <TRandom3.h>

using namespace std;

 TF1 *combg_fit = new TF1("combg_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);

Double_t CBcalc_LL(Double_t *xx, Double_t *par) // CB sigma scaled by 0.20
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
  double X_mean_2nd = x_mean;
  double fixed = par[6];

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
    
  double gauss2_psip = N_2nd*exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psip;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psip;
    }

 //Like sign combnatoric background fit
  double bg12 = combg_fit->Eval(x);

  // Add total fit with like sign background fit
  f += bg12;  // comb bg fit results are added to jpsi CB, psi(2s)CB and corr bg 

  return f;
}

Double_t CBcalc_jpsi(Double_t *xx, Double_t *par)  // CB calc only
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
  double fixed = par[7];
 
  double gauss2_psip = N_2nd*exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psip;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psip;
    }

  return f;
}

Double_t Gauss(Double_t *xx, Double_t *par)   // Gaussian sigma scaled by 0.23
{
  double f;
  double x = xx[0];

  double N_2nd = par[0];
  double x_mean_2nd = par[1];
  double fixed = par[2];
 

  f = N_2nd*exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  return f;

}

void CB_macro_sanghoon_for_sim_check()
{

 ///////////////////////////////////////////////////////
 bool north_arm = false;
 bool jpsi = true;

 bool fvtx = true;
 bool sngtrk = false;
 bool fvtxsngtrk = false;
 bool nofvtx = false;
 ///////////////////////////////////////////////////////
 
 TH1 *recomassBG;
 TH1 *recomassBG1;
 TH1 *recomassBG2;
 TF1 *combg_fit;
 TH1 *sim_mass;
 bool embed = true;
  bool sanghoon_files = true;

  bool chisquare_write = true;

  int bin; 


  int com_bg;
  // cout << "Enter combinatoric BG initial parameter set (0-13) " << endl;
  // cin >> com_bg;

  if(north_arm)
    com_bg = 0;
  else
    com_bg = 0;
 
   //cout << "Enter the bin number " << endl;//  (1,2,3,4,5) = rap1,rap2,rap3,rap4,MB
   //cin >> bin;
   bin = 7;
  
   TFile *file1S; 

   if(north_arm && embed)
     {
       if(jpsi)
	 file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_jpsi_newlib_trig_fid188228_nochi2fvtx_trigZ_N.root");
       else
	 file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_psi2s_newlib_trig_fid188228_nochi2fvtx_trigZ_N.root");
     }
   if(!north_arm && embed)
     {
       if(jpsi)
	 file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_jpsi_newlib_trig_fid188228_nochi2fvtx_trigZ_S.root");
       else
	 file1S = TFile::Open("sanghoon_simulations/yes_sngdbl_psi2s_newlib_trig_fid188228_nochi2fvtx_trigZ_S.root");
     }


  TCanvas *jpsi_sim = new TCanvas("jpsi_sim","jpsi_sim",5,5,800,600);
  jpsi_sim->SetLogy();
 
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
  double Njpsi = 0;
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


 double fixed = 0.20;
  double factor = 0.1;

 const double sigma_step = 0.001;
 
 //for(int i = 0; i < 300; i++)
   {
    
     // fixed =  0.10 + i*sigma_step;
    

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

  
 
      recomassBG->Fit(combg_fit,"LL","",2,5); 
 
  
      ////////////////////////////////

      // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 

      double N_start = 0;
 

      /////////////////////////////////////////////////////////////////////////////////
      coeff = 1.0;

      N_start = coeff*sim_mass->Integral(binl,binh) / renorm; 

      if(north_arm && jpsi)
	{
	  // coeff = 1.25; 
	  f0 = 1;
	  f1 = 10;
	  f2 = 3.1;
	  f3 = 0.07;  // 0.17
	  f4 = N_start; 
	  f5 = factor*N_start;
	   f6 = fixed;
	}
      if(north_arm && !jpsi) // psi2s
	{
	  f0 = 1;  // 3
	  f1 = 10; // 0
	  f2 = 3.7;
	  f3 = 0.12;  // 0.17
	  //f4 = N_start; 
	  f5 = factor*N_start;
	  f4 = 10*N_start;
	  // f5 = 0.12;
	  f6 = fixed;
	}
      /////////////////////////////////////////////////////////////////////////////////

      coeff = 1.0;

      if(!north_arm && jpsi)
	{
	  // coeff = 1.25; 
	  f0 = 1;
	  f1 = 10;
	  f2 = 3.1;
	  f3 = 0.10;  // 0.10
	  f4 = N_start;
	  f5 = factor*N_start;
	  f6 = fixed;
	}
      if(!north_arm && !jpsi)  // psi2s
	{
	  f0 = 1;
	  f1 = 10;
	  f2 = 3.7;
	  f3 = 0.08;  // 0.17
	  //f4 = N_start; 
	  f5 = factor*N_start; 
	  f4 = 10*N_start;
	  // f5 = 0.12;
	  f6 = fixed;
	}


      char const *outfile;
      if(jpsi && north_arm)
      	outfile = "jpsi_mass_dist_for_MC_random_N.root"; 
      if(!jpsi && north_arm)
      	outfile = "psi2s_mass_dist_for_MC_random_N.root"; 
      if(jpsi && !north_arm)
      	outfile = "jpsi_mass_dist_for_MC_random_S.root"; 
      if(!jpsi && !north_arm)
      	outfile = "psi2s_mass_dist_for_MC_random_S.root"; 

      // char const *outfile;
      // if(jpsi && north_arm)
      // 	outfile = "jpsi_mass_BG_dist_for_random_N.root"; 
      // if(!jpsi && north_arm)
      // 	outfile = "psi2s_mass_BG_dist_for_random_N.root"; 
      // if(jpsi && !north_arm)
      // 	outfile = "jpsi_mass_BG_dist_for_random_S.root"; 
      // if(!jpsi && !north_arm)
      // 	outfile = "psi2s_mass_BG_dist_for_random_S.root"; 

      TFile *h = new TFile(outfile, "RECREATE");
      h->cd();

      jpsi_sim->SetName("h1");
      jpsi_sim->Write();
      sim_mass->SetName("h2");
      sim_mass->Write();
      recomassBG->SetName("h3");
      recomassBG->Write();
      h->Close();


      TF1 *total_fit = new TF1("total_fit",CBcalc_LL,2.0,5.0,7);
      total_fit->SetLineColor(kBlue+1);
      total_fit->SetParameter(0, f0);
      // total_fit->SetParLimits(0, 2, 4);
      total_fit->SetParLimits(0, 0, 2);
      total_fit->SetParameter(1, f1);
      // total_fit->SetParLimits(1, 0, 1);
      total_fit->SetParLimits(1, 5, 15);
      total_fit->SetParameter(2, f2);
      total_fit->SetParameter(3, f3);
      total_fit->SetParameter(4, f4);
      total_fit->SetParameter(5, f5);
      total_fit->FixParameter(6, f6);
  
   
      sim_mass->Fit(total_fit,"LL","",2,5); 
      //sim_mass->Fit(total_fit,"","",2,5); 

      
	// TF1 *fresult= 0;
  
	// TFitResultPtr r = sim_mass->Fit(total_fit,"S");////// does the fit
	// fresult = sim_mass->GetFunction("total_fit");
      
	// TMatrixDSym cov = r->GetCovarianceMatrix();
	// r->Print("V"); //verbose setting
	// Double_t *fullmat = 0;
	// fullmat = cov.GetMatrixArray();

	// double jpsi_array[36] = {0};

	// for(Int_t i = 0;i < 6; i++)
	// {
	// for(Int_t j = 0;j < 6; j++)
	// {
	// jpsi_array[6*i+j] = fullmat[6*i+j];
	// //cout << "jpsi array " << i << " " << j << " " << jpsi_array[i+j] << endl;
	// }  
	// }
  
	// sim_mass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
	// sim_mass->Draw();
	// sim_mass->GetXaxis()->SetRangeUser(2,5);
	// // sim_mass->Fit(total_fit);
       


      sim_mass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
      sim_mass->Draw();
      sim_mass->GetXaxis()->SetRangeUser(2,5);
      // sim_mass->Fit(total_fit, "LS");
      total_fit->Draw("same");
      combg_fit->Draw("same");

  
      TF1 *jpsi_fit = new TF1("jpsi_fit",CBcalc_jpsi,2.0, 5.0, 8);  
      jpsi_fit->SetLineColor(kRed+1);
      jpsi_fit->SetParameter(0, total_fit->GetParameter(0));
      jpsi_fit->SetParError(0, total_fit->GetParError(0));
      jpsi_fit->SetParameter(1, total_fit->GetParameter(1));
      jpsi_fit->SetParError(1, total_fit->GetParError(1));
      jpsi_fit->SetParameter(2, total_fit->GetParameter(2));
      jpsi_fit->SetParError(2, total_fit->GetParError(2));
      jpsi_fit->SetParameter(3, total_fit->GetParameter(3));
      jpsi_fit->SetParError(3, total_fit->GetParError(3));
      jpsi_fit->SetParameter(4, total_fit->GetParameter(4));
      jpsi_fit->SetParError(4, total_fit->GetParError(4));
      jpsi_fit->SetParameter(5, total_fit->GetParameter(5));
      jpsi_fit->SetParError(5, total_fit->GetParError(5));
      jpsi_fit->FixParameter(6, total_fit->GetParameter(2));
      jpsi_fit->SetParameter(7, fixed);
    


      TF1 *gauss_fit = new TF1("gauss_fit",Gauss,2.0,5.0,3); // gauss_jpsi
      gauss_fit->SetParameter(0,total_fit->GetParameter(5));
      gauss_fit->SetParError(0,total_fit->GetParError(5));
      gauss_fit->SetParameter(1,total_fit->GetParameter(2));
      gauss_fit->SetParError(1,total_fit->GetParError(2));
      gauss_fit->SetParameter(2,fixed);
    
      gauss_fit->SetLineColor(kViolet+2);
      gauss_fit->SetLineStyle(10);
   
   
      // if(north_arm)
      //   {
      //     jpsi_fit->SetParLimits(0,0,1);
      //     jpsi_fit->SetParLimits(1,0,100);
      //   }
      // else
      //   {
      //     jpsi_fit->SetParLimits(0,0,100);
      //     jpsi_fit->SetParLimits(1,0,1000);
      //   }

      
      // Double_t jpsipar[6] = {0};
      
      // for(int i = 0;i < 6; i++)
      //   {
      //     jpsipar[i] = total_fit->GetParameter(i); 
      //   }
  
      double Njpsi = 0;
      double err_Njpsi = 0;
      Double_t width = sim_mass->GetBinWidth(1);
      //  cout << "binwidth: " << width << endl;
      Njpsi = jpsi_fit->Integral(2.0, 5.0)/width;  
      //err_Njpsi = jpsi_fit->IntegralError(2.0,5,jpsipar,jpsi_array)/width;
      err_Njpsi = sqrt(Njpsi);

      double chisquare;
      double ndf;

      gauss_fit->Draw("same");
      jpsi_fit->Draw("same");
      recomassBG->Draw("same");

      chisquare = total_fit->GetChisquare();
      ndf =  total_fit->GetNDF();
 
  
      TLatex l3;
      l3.SetTextSize(0.06);
      l3.SetTextAlign(13);
      l3.SetTextColor(4);

      char text4[100];
      if(north_arm && embed)
	sprintf(text4,"p + p North J/#psi, Embedded");
      if(!north_arm && embed)
	sprintf(text4,"p + p South J/#psi, Embedded");
      if(north_arm && !jpsi)
	sprintf(text4,"p + p North #psi(2S), Embedded");
      if(!north_arm && !jpsi)
	sprintf(text4,"p + p South #psi(2S), Embedded");

      l3.SetTextAlign(12);
      l3.DrawLatexNDC(0.4, 0.93, text4); //4.4,150
 
      Char_t message[200];
      Char_t message2[200];
 
      sprintf(message,"#chi^{2}/NDF = %.1f / %.d",total_fit->GetChisquare(),total_fit->GetNDF());
      sprintf(message2,"Jpsi Counts =  %.0f +/- %.0f", Njpsi, sqrt(Njpsi));
      TPaveText *mytext = new TPaveText(0.12,0.85,0.42,0.75,"NDC"); // x0,y0,x1,y1
      mytext->SetTextSize(0.03);
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
      if(bin == 5)
	{
	  if(north_arm == true)
	    sprintf(text3,"p_{T} Integrated");
	  if(north_arm == false)
	    sprintf(text3,"p_{T} Integrated");
	};
       
      l2.SetTextAlign(12);
      l2.SetTextSize(0.04);
      l2.DrawLatexNDC(0.7, 0.84, text3);
  
      char text5[100];
      char text6[100];
      l3.SetTextColor(1);
      l3.SetTextSize(0.035);
      sprintf(text5,"#sigma_{2nd} =  %.3f ",fixed);
      sprintf(text6,"Double Tracks");
      l3.SetTextAlign(12);
      l3.DrawLatexNDC(0.7, 0.78, text5); //4.4,150
      l3.DrawLatexNDC(0.7, 0.72, text6); //4.4,150


      double X_mean =  total_fit->GetParameter(2);
      double CB_mean =  total_fit->GetParameter(2);
      double alpha =  total_fit->GetParameter(0);
      double n =  total_fit->GetParameter(1);
      double sigma =  total_fit->GetParameter(3);
      double N =  total_fit->GetParameter(4);
      double N_2nd = total_fit->GetParameter(5);

      ///////////////////////////////////////
      char name700[500]; 


      // cout << "Njpsi: " << Njpsi << " +/- " << sqrt(Njpsi) << endl;
      // cout << "Chi2/NDF: " << total_fit->GetChisquare() << "/" << total_fit->GetNDF() << endl;
      // cout << "N_2nd/N: " <<  total_fit->GetParameter(5)/total_fit->GetParameter(4) << endl;
      // cout <<  total_fit->GetParameter(5) << endl;
      // cout << total_fit->GetParameter(4) << endl;

      double ratio_int;
      double num, denom;
      // double num,denom,ratio,total_int,combg_int,gauss_int,jpsi_int;
      num = gauss_fit->Integral(2.0,5.0);
      denom = jpsi_fit->Integral(2.0, 5.0); 
      ratio_int = num/denom;
      // total_int = total_fit->Integral(2.0,5.0);
      // gauss_int = gauss_fit->Integral(2.0,5.0);
      // combg_int = combg_fit->Integral(2.0,5.0);
      // jpsi_int = jpsi_fit->Integral(2.0,5.0);

      // cout << ratio << endl;

      // cout << total_fit->Integral(2.0,5.0) << endl;
      // cout << gauss_fit->Integral(2.0,5.0) << endl;
      // cout << combg_fit->Integral(2.0,5.0) << endl;
      // cout << jpsi_fit->Integral(2.0,5.0) << endl;

      //cout << N_start << endl;

      // if(chisquare_write)
      // 	{
    
      // 	  if(north_arm)
      // 	    {
      // 	      if(fvtx)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_fvtx_v2_north_loop/N_gaussian_%i.txt", i);
      // 	      if(fvtxsngtrk)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_fvtxsngtrk_south_loop/N_gaussian_%i.txt", i);
      // 	      if(nofvtx)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_nofvtx_south_loop/N_gaussian_%i.txt", i);
      // 	      if(sngtrk)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_sngtrk_south_loop/N_gaussian_%i.txt", i);
      // 	    }	 
      // 	  if(!north_arm)
      // 	    {
      // 	      if(fvtx)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_fvtx_v2_south_loop/S_gaussian_%i.txt", i);
      // 	      if(fvtxsngtrk)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_fvtxsngtrk_south_loop/S_gaussian_%i.txt", i);
      // 	      if(nofvtx)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_nofvtx_south_loop/S_gaussian_%i.txt", i);
      // 	      if(sngtrk)
      // 		sprintf(name700,"CB_parameter_scan_vs_sigma_2nd/jpsi_sngtrk_south_loop/S_gaussian_%i.txt", i);
      // 	    }	  
	      
      // 	  // whether north arm or south arm, write to unique filename  
      // 	    std::fstream Run15pp_width(name700,std::ofstream::out); 
      // 	  Run15pp_width << Njpsi << " " <<  chisquare/ndf << " " << fixed << " " << ratio_int << " " << CB_mean << " " << N_2nd << " " << alpha << " " << n << " " << sigma << " " << N << " " <<  endl;
      // 	  Run15pp_width.close();
   

      // 	} // if write chisquare



   } // end for loop


// // create a random number generator
//         gRandom = new TRandom3();

//         // create a histogram 
//         TH1D * hist = new TH1D("data", ";x;N", 20, 0.0, 100.0);

//         // fill in the histogram
//         for (int i = 0; i < 100; ++i)
// 	  // hist->Fill(gRandom->Gaus(65.0, 5.0));
// 	  hist->Fill(gRandom->CBcalc_LL(f0,f1,f2,f3,f4,f5,f6));

//         TCanvas * c1= new TCanvas("c1", "random",5,5,800,600);
//         hist->Draw();
//         hist->SaveAs("random.eps");
















}// end void macro



