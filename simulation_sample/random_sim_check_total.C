

// For p + p

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
#include <TAttText.h>
#include <TRandom3.h>

using namespace std;



bool north_arm = false;

Double_t CBcalc_jpsi(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
{
  double f;
  double x = xx[0];

  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:
 

 double ratio = 0;
  double fraction = 0;

  // if(north_arm)
  //   ratio = 976.375/8661.16;
  // else
  //   ratio = 1751.11/9888.28;
  // if(north_arm)
  //   fraction = 0.160/0.08068;
  // else
  //   fraction = 0.160/0.08655;

 //  if(north_arm)
 //    ratio = 998.054/9504.71;
 //  else
 //    ratio = 1902.86/9616.81;
  
 //  if(north_arm)
 //    fraction = 0.19/0.09309;
 //  else
 //    fraction = 0.19/0.10030;


  // if(north_arm)
  //   ratio = 622.847/9821.05;
  // else
  //   ratio = 1168.2/10267.1;
  
  // if(north_arm)
  //   fraction = 0.219/0.09607;
  // else
  //   fraction = 0.219/0.10537;

 // if(north_arm)
 //   ratio = 735.23/9920.44;
 //  else
 //    ratio = 1525.86/11328.2;
  
 //  if(north_arm)
 //    fraction = 0.212/0.09495;
 //  else
 //    fraction = 0.212/0.10406;



  // if(north_arm)
  //   ratio = 830.814/10856.3;
  // else
  //   ratio = 1350.98/10816.9;
  
  // if(north_arm)
  //   fraction = 0.216/0.09466;
  // else
  //   fraction = 0.216/0.10299;


  if(north_arm)
    ratio = 843.22/10846.3;
  else
    ratio = 1372.08/10799.2;
  
  if(north_arm)
    fraction = 0.215/0.09456;
  else
    fraction = 0.215/0.10282;



  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];
  double N_2nd = ratio*N;
  double X_mean_2nd = x_mean;
 
  double fixed = fraction*sigma;

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd  * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
 
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
f = N  * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }

  return f;
}

Double_t CBcalc_psi2s(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
{
  double f;
  double x = xx[0];
 
  

 double ratio = 0;
 double fraction = 0;



  if(north_arm)
    ratio = 843.22/10846.3;
  else
    ratio = 1372.08/10799.2;
  
  if(north_arm)
    fraction = 0.215/0.09456;
  else
    fraction = 0.215/0.10282;






  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double n_2nd = ratio*N;
  double x_mean_2nd = x_mean;
  double fixed = fraction*sigma;

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psi2s;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psi2s;
    }

  return f;
}


Double_t CBcalc_jpsi_G(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
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

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
 
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
      f = N  * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }

  return f;
}


Double_t CBcalc_psi2s_G(Double_t *xx, Double_t *par)  // CB fucntion for psi 2s
{
  double f;
  double x = xx[0];
 
  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:

  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double n_2nd = par[5];
  double x_mean_2nd = par[6];
  double fixed = par[7];
 
  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
 
  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psi2s;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psi2s;
    }

  return f;
}

///the second gaussian alone to draw fit for J/psi and psi2s
Double_t Gauss(Double_t *xx, Double_t *par) 
{
  double f;
  double x = xx[0];

  double N_2nd = par[0];
  double x_mean_2nd = par[1];
  double fixed = par[2];


  f = N_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));  // multiply by tot_fit->GetParameter(f4) for J/psi and tot_fit->GetParameter(f10) for psi2S because this is only for drawing purposes

  return f;

}

Double_t CBcalc_LL(Double_t *xx, Double_t *par) // this is total fit function.  Jpsi CB ([0], [1], [2], [3], [4]) + jpsi(2s) CB ([7], [8], [9]) + correlated bg ([5], [6], [10], [11], [12] and then to this the comb fit is added at line 177.  Added also second gaussian under J/psi peak to better resemble data in the mass region between the J/psi and psi(2S) peaks [13], [14], and [12]
{
  double f;
  double x = xx[0];

    

 double ratio = 0;
 double fraction = 0;

 
  
  if(north_arm)
    ratio = 843.22/10846.3;
  else
    ratio = 1372.08/10799.2;
  
  if(north_arm)
    fraction = 0.215/0.09456;
  else
    fraction = 0.215/0.10282;





   // J/psi CB parameters
  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];
  
  // the second gaussian for J/psi peak
  double N_2nd = ratio * N;
  double  X_mean_2nd = x_mean;
 
  // psi2s CB parameters (alpha and n equal J/psi values)
  double N_psi2s = par[5];
  double mass_offset = par[6];
  double width_factor = par[7];
  
  // the second gaussian for psi2s peak
  double n_2nd = ratio * N_psi2s;
  double x_mean_2nd = x_mean + mass_offset;
 
  // The other psi2s parameters will be constrained by the jpsi parameters, only the normalization will vary
  double x_mean_psi2s = x_mean + mass_offset;
  double sigma_psi2s = sigma*width_factor;

 double fixed = fraction*sigma;  
  double fixed_psi2s = fraction*sigma_psi2s;

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed_psi2s,2)));

  // for J/psi CB function
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }
  // for psi(2S) CB function
  if( (x-x_mean_psi2s)/sigma_psi2s > -alpha)
    {
      f += N_psi2s * exp( -pow(x-x_mean_psi2s,2) / (2.0*pow(sigma_psi2s,2))) + gauss2_psi2s;
    }
  else
    {
      f += N_psi2s * A * pow(B - (x-x_mean_psi2s)/sigma_psi2s, -n) + gauss2_psi2s;
    }
  // add correlated backgroudn to CB functions
  // f += par_c /( pow( exp(-par_a*x - par_b*x*x) + x/par_d, par_e)); // NEW corrbg formula

  //Like sign combinatoric background fit
  // double bg12 = combg_fit->Eval(x);

  // Add total fit with like sign background fit
  //f += bg12;  // comb bg fit results are added to jpsi CB, psi(2s)CB and corr bg 

  if(x == 2.60)
    {
      cout << "==================== " << endl;
      cout << "J/psi 2nd gaussian parameters:" << endl;
      cout << "==================== " << endl;
      cout <<" N_2nd : " << N_2nd << endl;
      cout << " X_mean_2nd : " << X_mean_2nd << endl;
      cout << " Sigma 2nd : " << fixed << endl;
      cout << " N CB : " << N << endl;
      cout << " Sigma_CB: " << sigma << endl;
      cout << "==================== " << endl;
      cout << "psi2S 2nd gaussian parameters:" << endl;
      cout << "==================== " << endl;
      cout <<" n_2nd : " << n_2nd << endl;
      cout << " x_mean_2nd : " << x_mean_2nd << endl;
      cout << " sigma 2nd : " << fixed << endl;
      cout << " n CB : " << N_psi2s << endl;
      cout << " sigma_CB: " << sigma_psi2s << endl;
      cout << " ratio : " << n_2nd / N_psi2s << endl;
      cout << " ratio: " << ratio << endl;
      cout << "==================== " << endl;

      if( N_2nd/N == n_2nd/N_psi2s)
   	cout << " --> 1. normalization ratios are equal " << endl;
      if(  N_2nd/N != n_2nd/N_psi2s)
   	cout << " --> 1. f: " << N_2nd / N << " and f': " << n_2nd / N_psi2s << " " << endl;

      if(fixed)
   	cout << " --> 2. gaussian widths both fixed to " << fixed << " " << endl;

      if(X_mean_2nd == x_mean)
   	cout << " --> 3. 2nd gaussian mean equal to Jpsi CB mean " << endl;

      if(x_mean_2nd == x_mean_psi2s)
   	cout << " --> 4. 2nd gaussian mean equal to psi2S CB mean " << endl;

    }

  return f;
}



void random_sim_check_total()
{




  double ave = 0;
  double weight = 0;
  double sigma = 0;

  ///////////////////////////////////////////////////////
  //bool north_arm = true;
  bool north_arm = false;

  bool fvtx = true;
  bool sngtrk = false;
  bool fvtxsngtrk = false;
  bool nofvtx = false;

  bool pt_low = false;
  bool pt_high = false;
  ///////////////////////////////////////////////////////

  bool embed = true;
  bool write = true;

  TH1D *sim_mass = new TH1D("sim_mass","",60,2,5);
  TH1D *sim_psi2s_mass = new TH1D("sim_psi2s_mass","",60,2,5);
  
  TFile *file_data1; 
  TFile *file_data2; 

  if(north_arm)
    {
    file_data1 = TFile::Open("jpsi_mass_dist_for_MC_random_N.root");
    file_data2 = TFile::Open("psi2s_mass_dist_for_MC_random_N.root");
    }
  else
    {
      file_data1 = TFile::Open("jpsi_mass_dist_for_MC_random_S.root");
      file_data2 = TFile::Open("psi2s_mass_dist_for_MC_random_S.root");
    }
  // if(north_arm)
  //   {
  //   file_data1 = TFile::Open("jpsi_mass_dist_for_random_N_above_2pT.root");
  //   file_data2 = TFile::Open("psi2s_mass_dist_for_random_N_above_2pT.root");
  //   }
  // else
  //   {
  //     file_data1 = TFile::Open("jpsi_mass_dist_for_random_S_above_2pT.root");
  //     file_data2 = TFile::Open("psi2s_mass_dist_for_random_S_above_2pT.root");
  //   }

  file_data1->GetObject("h2",sim_mass);  
  file_data2->GetObject("h2",sim_psi2s_mass);  

  
  int step;
  int binl = sim_mass->FindBin(2.6); 
  int binh = sim_mass->FindBin(3.6); 
  double bin_width = sim_mass->GetBinWidth(1);
  double renorm = 1.0/bin_width;

  TCanvas * c1= new TCanvas("c1", "random",5,5,800,600);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLogy();
  c1->cd();
  c1->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.04);
  gPad->SetBottomMargin(0.13);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
 
 TH1D *hmass = new TH1D("hmass", "", 60,2,5);  // outside 
 TH1D *hpsi2s_mass = new TH1D("hpsi2s_mass", "", 60,2,5);  // outside loop

 TH1D *hbias_jpsi = new TH1D("hbias_jpsi","",100,0,100);
 TH1D *hbias_psi2s = new TH1D("hbias_psi2s","",100,0,100);
	  


 //Double_t gRandom = sim_mass->GetRandom();
 
  if(gRandom)
    {
      delete gRandom;
      gRandom = new TRandom3(0);
    }
  
   for(int f=0; f < 150; f++) 
    {
      hmass->Reset(); //clear contents of previous h1 (We are generating n4 number of h1 histograms)
      gRandom->SetSeed(0);
      hmass->FillRandom(sim_mass,5000);


      hpsi2s_mass->Reset(); //clear contents of previous h1 (We are generating n4 number of h1 histograms)
      gRandom->SetSeed(0);
      hpsi2s_mass->FillRandom(sim_psi2s_mass,150);

      hmass->Add(hpsi2s_mass,1);
      /////////////////////////////////////////////////////////////////////////////////////////////

      double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,coeff;
      double N_start = 0;
      double fixed = 0.165;

      /////////////////////////////////////////////////////////////////////////////////
      coeff = 1.5;
      step = f;

      N_start = coeff*sim_mass->Integral(binl,binh) / renorm; 

     
      int fit = 1;
    // combinatoric BG initital parameters for pT dependent fitting - USE SAME SET FOR ALL SYSTEMS AND ALL DEPENDENCES

      //  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 
 
      
      if(north_arm)
	{
	  f0 = 1;
	  f1 = 10;
	  f2 = 3.1; // Jpsi x_mean
	  f3 = 0.07;	 
	  f4 = N_start;  // Jpsi norm
	  f5 = 0.05*N_start; // norm psi2s
	  f6 = 0.589;  //  x_mean ppsi2s
	  f7 = 1.15;   // width factor psi2s  ------- > used 1.15 in pdf sent to sanghoon
	    
	}
 
      else
	{
	  f0 = 2;  // alpha 1
	  f1 = 5; // n 10
	  f2 = 3.1; // Jpsi x_mean
	   f3 = 0.08;  // Jpsi sigma 
	  f4 = N_start;  // Jpsi norm
	  f5 = 0.05*N_start; // norm psi2s
	  f6 = 0.589;  //  x_mean ppsi2s
	  f7 = 1.15; // calcualted from sim no comb bg
   
	}

      //////////////////////////////////////////////////////////////
 
 TF1 *total_fit = new TF1("total_fit",CBcalc_LL,2.0,5.0,8);

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
 
      total_fit->SetParameter(2, f2); // Jpsi CB mean
      total_fit->SetParameter(3, f3);  // Jpsi CB sigma
      total_fit->SetParameter(4, f4);    // N J/psi
      total_fit->SetParameter(5, f5);    //norm 2nd G Jpsi
      total_fit->FixParameter(6, f6);    // mass offset psi2s
      total_fit->FixParameter(7, f7);    // width factor psi2s
              

      total_fit->SetLineColor(kBlue+1);
      total_fit->SetMarkerStyle(42);
      total_fit->SetLineWidth(5);
      total_fit->SetLineStyle(10);
      //hmass->Fit(total_fit,"LL","",2,4.2); 
      hmass->Fit(total_fit,"","",2,4.4); 
  
      hmass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
      hmass->Draw();
      hmass->GetXaxis()->SetRangeUser(2,5);
      total_fit->Draw("same");
      

  TF1 *fresult= 0;

 
  /////////////////////////////////////////////////////////////BEGIN ORIGINAL CODE
  


  
      // log likelihood fit with integer data option is "L"
      TFitResultPtr r = hmass->Fit(total_fit,"S");////// does the fit
      fresult = hmass->GetFunction("total_fitLL");
      
      TMatrixDSym cov = r->GetCovarianceMatrix();
      r->Print("V"); //verbose setting
      Double_t *fullmat = 0;
      fullmat = cov.GetMatrixArray();

      TF1 *jpsi_fit = new TF1("jpsi_fit",CBcalc_jpsi,2.0,4.0,5); // use CBcalc_jpsi instead
      for(int i=0;i<5;i++)
	{
	  jpsi_fit->SetParameter(i,total_fit->GetParameter(i));
	  jpsi_fit->SetParError(i,total_fit->GetParError(i));
	}
 
      double Jpsi_array[25] = {0};
          // can only use this is J/psi CB parameters are sequential in the total fit function
         
     for(Int_t i = 0;i < 5; i++)
	{
	  for(Int_t j = 0;j < 5; j++)
	    {
	      // Jpsi_array[5*i+j] = fullmat[10*i+j];
	      Jpsi_array[5*i+j] = fullmat[8*i+j];
	      // cout << "Jpsi array " << i << " " << j << " " << Jpsi_array[8*i+j] << endl;
	    }  
	}
   
        
    Double_t Jpsipar[5] = {0};
      
      for(int i = 0;i < 5; i++)
	{
	  Jpsipar[i] = total_fit->GetParameter(i); 
	}
   
      double Njpsi = 0;
      double err_NJpsi = 0;
      Double_t width = sim_mass->GetBinWidth(1);
      Njpsi = jpsi_fit->Integral(2.0, 5.0)/width;  
      //err_Njpsi = sqrt(Njpsi);
      err_NJpsi = jpsi_fit->IntegralError(2.0,5.0,Jpsipar,Jpsi_array)/width;

      double chisquare;
      double ndf;
       
      double psi2s_array[25] = {0};

      // CB calc   --- >    fit total fucntion parameters.  psi2s array must be in order of parameters( rows [16],[17],[18] must go last)

      // par 0 row (alpha)
      psi2s_array[0] = fullmat[0]; // alpha
      psi2s_array[1] = fullmat[1]; // n
      psi2s_array[2] = fullmat[2]; // x_mean
      psi2s_array[3] = fullmat[3]; // sigma
      psi2s_array[4] = fullmat[5]; // N psi2S
               
    // row 1 n
      psi2s_array[5] = fullmat[8]; 
      psi2s_array[6] = fullmat[9]; 
      psi2s_array[7] = fullmat[10]; 
      psi2s_array[8] = fullmat[11]; 
      psi2s_array[9] = fullmat[13]; 
          
      // row 2 x_mean
      psi2s_array[10] = fullmat[16]; 
      psi2s_array[11] = fullmat[17]; 
      psi2s_array[12] = fullmat[18]; 
      psi2s_array[13] = fullmat[19]; 
      psi2s_array[14] = fullmat[21]; 
         
      // row 3 sigma
      psi2s_array[15] = fullmat[24];  // N psi2S
      psi2s_array[16] = fullmat[25]; 
      psi2s_array[17] = fullmat[26]; 
      psi2s_array[18] = fullmat[27]; 
      psi2s_array[19] = fullmat[29]; 
         
      // row 8 n_CB
      psi2s_array[20] = fullmat[40]; 
      psi2s_array[21] = fullmat[41];  // N psi2S
      psi2s_array[22] = fullmat[42]; 
      psi2s_array[23] = fullmat[43]; 
      psi2s_array[24] = fullmat[45]; 
        
      //////////////////////////////////////

      double frac = 1.15;
   
      TF1 *psi2s_fit = new TF1("psi2s_fit",CBcalc_psi2s,2.0, 4.4, 5);  // use CBcalc_psi2s instead
      psi2s_fit->SetParameter(0,total_fit->GetParameter(0));
      psi2s_fit->SetParError(0,total_fit->GetParError(0));
      psi2s_fit->SetParameter(1,total_fit->GetParameter(1));
      psi2s_fit->SetParError(1,total_fit->GetParError(1));
      psi2s_fit->SetParameter(2,total_fit->GetParameter(2) + total_fit->GetParameter(6));
      psi2s_fit->SetParError(2,total_fit->GetParError(2));  // par 6 is fixed
      psi2s_fit->SetParameter(3,total_fit->GetParameter(3)*total_fit->GetParameter(7));
      psi2s_fit->SetParError(3,total_fit->GetParError(3)*total_fit->GetParameter(7)); 
      psi2s_fit->SetParameter(4,total_fit->GetParameter(5));
      psi2s_fit->SetParError(4,total_fit->GetParError(5));
                      
      Double_t psi2s_par[5] = {0}; // 5 CB parameters 

      psi2s_par[0] =  total_fit->GetParameter(0);     
      psi2s_par[1] =  total_fit->GetParameter(1);     
      psi2s_par[2] =  total_fit->GetParameter(2) + total_fit->GetParameter(6);         
      psi2s_par[3] =  total_fit->GetParameter(3)*total_fit->GetParameter(7);       
      psi2s_par[4] =  total_fit->GetParameter(5);          
              
      double Npsi2s = 0;
      double err_Npsi2s = 0;
      Npsi2s = psi2s_fit->Integral(2.0, 5.0)/width;  
      //err_Npsi2s = sqrt(Npsi2s);
      err_Npsi2s = psi2s_fit->IntegralError(2.0,5.0,psi2s_par,psi2s_array)/width;

      hmass->SetTitle("");
      hmass->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass(GeV/c^{2})");
      hmass->GetXaxis()->SetLabelSize(0.05);
      hmass->GetXaxis()->SetNdivisions(5,5,0);
      hmass->GetXaxis()->SetTitleOffset(1.11);  
      // hmass->GetYaxis()->SetLabelSize(0.05);
      hmass->GetYaxis()->SetTitleOffset(1.3);  
      hmass->GetYaxis()->SetTitle("raw counts/(50 MeV/c^{2})");
      hmass->GetYaxis()->SetTitleFont(102);
      hmass->GetYaxis()->SetLabelFont(102);
      hmass->GetYaxis()->SetTitleSize(0.045);
      hmass->GetYaxis()->SetLabelSize(0.045);
      hmass->GetXaxis()->SetTitleSize(0.045);
      hmass->GetXaxis()->SetTitleOffset(0.85);
      hmass->GetXaxis()->SetLabelSize(0.045);
      hmass->GetXaxis()->SetTitleFont(102);
      hmass->GetXaxis()->SetLabelFont(102);
      hmass->DrawCopy();

  
      chisquare = total_fit->GetChisquare();
      ndf =  total_fit->GetNDF();
      double chi2ndf = chisquare/ndf;
  
    
     
      ///////////////////////////////////////////////////////////////////////////////END ORIGINAL CODE

      // Plot the fit results
  
 double ratio = 0;
 double fraction = 0;


  if(north_arm)
    ratio = 843.22/10846.3;
  else
    ratio = 1372.08/10799.2;
  
  if(north_arm)
    fraction = 0.215/0.09456;
  else
    fraction = 0.215/0.10282;




      double Fixed = fraction*total_fit->GetParameter(3);
      double Fixed_psi2s = fraction*total_fit->GetParameter(3)*total_fit->GetParameter(7);

      double norm_2nd = ratio* total_fit->GetParameter(4);
      double norm_psi2s = ratio* total_fit->GetParameter(5);

      double x_mean_2nd = total_fit->GetParameter(2) + 0.589;
      double X_mean_2nd = total_fit->GetParameter(2);

      TF1 *jpsi_fit2 = new TF1("jpsi_fit2",CBcalc_jpsi_G,2.0,5.0,8); // only total fit fucntion has correct normalzation for 2nd gaussian
      for(int i=0;i<5;i++)
	{
	  jpsi_fit2->SetParameter(i,total_fit->GetParameter(i));
	  jpsi_fit2->SetParError(i,total_fit->GetParError(i));
	}
      jpsi_fit2->SetParameter(5,norm_2nd);
      jpsi_fit2->FixParameter(6,X_mean_2nd);
      jpsi_fit2->FixParameter(7,Fixed);
       
      jpsi_fit2->SetLineColor(12);
      jpsi_fit2->SetLineWidth(3);
      jpsi_fit2->SetLineStyle(1);
      jpsi_fit2->Draw("same");


      TF1 *gauss_fit = new TF1("gauss_fit",Gauss,2.0,5.0,3); // gauss2  J/psi
      gauss_fit->SetParameter(0,norm_2nd);
      gauss_fit->SetParError(0,norm_2nd);
      gauss_fit->FixParameter(1,X_mean_2nd);
      gauss_fit->SetParameter(2,Fixed);
      gauss_fit->SetParError(2,Fixed);
  
      gauss_fit->SetLineColor(kRed+1);
      gauss_fit->SetLineStyle(10);
      gauss_fit->Draw("same");

      TF1 *psi2s_fit2 = new TF1("psi2s_fit2",CBcalc_psi2s_G,2.0,5.0,8);  // CBcalc_psi2s
      psi2s_fit2->SetParameter(0,total_fit->GetParameter(0));
      psi2s_fit2->SetParError(0,total_fit->GetParError(0));
      psi2s_fit2->SetParameter(1,total_fit->GetParameter(1));
      psi2s_fit2->SetParError(1,total_fit->GetParError(1));
      psi2s_fit2->SetParameter(2,total_fit->GetParameter(2) + total_fit->GetParameter(6));
      psi2s_fit2->SetParError(2,total_fit->GetParError(2));  // par 6 is Fixed
      psi2s_fit2->SetParameter(3,total_fit->GetParameter(3)*total_fit->GetParameter(7));
      psi2s_fit2->SetParError(3,total_fit->GetParError(3)*total_fit->GetParameter(7)); 
      psi2s_fit2->SetParameter(4,total_fit->GetParameter(5));
      psi2s_fit2->SetParError(4,total_fit->GetParError(5));
     
      psi2s_fit2->SetParameter(5,norm_psi2s);
      psi2s_fit2->SetParError(5,norm_psi2s);
      psi2s_fit2->FixParameter(6,x_mean_2nd);
      psi2s_fit2->FixParameter(7,Fixed_psi2s);
     
      psi2s_fit2->SetLineColor(kBlack);
      psi2s_fit2->Draw("same");
      psi2s_fit2->SetLineColor(12);
      psi2s_fit2->SetLineWidth(3);
      psi2s_fit2->SetLineStyle(1);


      TF1 *gauss2_fit = new TF1("gauss2_fit",Gauss,2.0,5.0,3); // gauss2_psi2s
      gauss2_fit->SetParameter(0,norm_psi2s);
      gauss2_fit->SetParError(0,norm_psi2s);
      gauss2_fit->FixParameter(1,x_mean_2nd);
      gauss2_fit->SetParameter(2,Fixed_psi2s);
      gauss2_fit->SetParError(2,Fixed_psi2s);

      gauss2_fit->SetLineColor(kRed+1);
      gauss2_fit->SetLineStyle(10);
      gauss2_fit->Draw("same");


      
      TLatex l3;
      l3.SetTextSize(0.06);
      l3.SetTextAlign(13);
      l3.SetTextColor(kBlue+1);

      char text4[100];

      if(north_arm && embed)
	sprintf(text4,"p + p North, Embedded");
      if(!north_arm && embed)
	sprintf(text4,"p + p South, Embedded");
      if(north_arm && !embed)
	sprintf(text4,"p + p North J/#psi, No Embed");
      if(!north_arm && !embed)
	sprintf(text4,"p + p South J/#psi, No Embed");

      l3.SetTextAlign(12);
      l3.DrawLatexNDC(0.48, 0.89, text4); //4.4,150
 
      Char_t message[200];
      Char_t message2[200];
      Char_t message3[200];
      Char_t message4[200];
 
      sprintf(message,"#chi^{2}/NDF = %.1f / %.d",total_fit->GetChisquare(),total_fit->GetNDF());
      sprintf(message2,"J/#psi Counts =  %.0f +/- %.0f", Njpsi, err_NJpsi);
      sprintf(message3,"#psi(2s) Counts =  %.0f +/- %.0f", Npsi2s, err_Npsi2s);
      if(pt_low)
	sprintf(message4,"     p_{T} < 2 GeV/c");
      if(pt_high)
	sprintf(message4,"     p_{T} > 2 GeV/c");
      TPaveText *mytext = new TPaveText(0.59,0.83,0.89,0.64,"NDC"); // x0,y0,x1,y1
      mytext->SetTextSize(0.03);
      mytext->SetFillColor(0);
      mytext->SetTextAlign(12);
      mytext->AddText(message);
      mytext->AddText(message2);
      mytext->AddText(message3);
      if(pt_high || pt_low)
	mytext->AddText(message4);
      mytext->Draw();


      TLatex l2;
      l2.SetTextSize(0.05);
      l2.SetTextAlign(13);
      l2.SetTextColor(1);

      char text5[100];
      char text6[100];
      char text7[100];
      l3.SetTextColor(1);
      l3.SetTextSize(0.04);
      sprintf(text5,"#sigma_{2nd} =  %.3f ",Fixed);
      sprintf(text6,"Sng+Dbl Tracks");
      sprintf(text7,"No Prob Cut");
      l3.SetTextAlign(12);
      l3.DrawLatexNDC(0.16, 0.87, text5); //4.4,150
      l3.DrawLatexNDC(0.16, 0.82, text6); //4.4,150
      // l3.DrawLatexNDC(0.16, 0.92, text7); //4.4,150

      ///////////////////////////////////////
      char name700[500]; 


      cout << "Njpsi: " << Njpsi << " +/- " << sqrt(Njpsi) << endl;
      cout << "Chi2/NDF: " << total_fit->GetChisquare() << "/" << total_fit->GetNDF() << endl;
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
      cout << gauss_fit->Integral(2.0,5.0) << endl;
      // cout << combg_fit->Integral(2.0,5.0) << endl;
      // cout << jpsi_fit->Integral(2.0,5.0) << endl;

      //cout << N_start << endl;


      // double Ngauss = gauss_fit->Integral(2.0, 5.0)/width;  
      // double Ngauss2 = gauss2_fit->Integral(2.0, 5.0)/width;  

      // cout << Ngauss << endl;
      // cout << Ngauss2 << endl;

      // cout << "subtracted Njpsi: " << Njpsi - Ngauss << endl;
      // cout << "subtracted Npsi2s: " << Npsi2s - Ngauss2 << endl;
      
	 
      //  hbias_psi2s->Fill(i+1,Npsi2s);
      //  hbias_jpsi->Fill(i+1, Njpsi);



      Int_t fitStatus = r;  
           
      cout << "status: " << fitStatus << endl;
      
      if(fitStatus != 4) 
	{
	  if(write)
	    {
	      
	      if(north_arm)
		{
		  sprintf(name700,"sim_check_ROOT/newlib_MC/N_sim_%i.txt", step);
		}	 
	      if(!north_arm)
		{
		  sprintf(name700,"sim_check_ROOT/newlib_MC/S_sim_%i.txt", step);
		}	  	
	      
	      // whether north arm or south arm, write to unique filename  
	      std::fstream Run15pp(name700,std::ofstream::out); 
	      Run15pp<< Njpsi << " " <<  err_NJpsi << " " << fixed << " " << Npsi2s << " " <<  err_Npsi2s << " " << endl;
	      Run15pp.close();
	      
	    } // if write
	}
      
  cout << "par 0 error: " << (  total_fit->GetParError(0) / total_fit->GetParameter(0) ) * 100 << endl;
  cout << "par 1 error: " << (  total_fit->GetParError(1) / total_fit->GetParameter(1) ) * 100 << endl;
  cout << "par 2 error: " << (  total_fit->GetParError(2) / total_fit->GetParameter(2) ) * 100 << endl;
  cout << "par 3 error: " << (  total_fit->GetParError(3) / total_fit->GetParameter(3) ) * 100 << endl;
  cout << "par 4 error: " << (  total_fit->GetParError(4) / total_fit->GetParameter(4) ) * 100 << endl;
  cout << "par 5 error: " << (  total_fit->GetParError(5) / total_fit->GetParameter(5) ) * 100 << endl;
  cout << "par 6 error: " << (  total_fit->GetParError(6) / total_fit->GetParameter(6) ) * 100 << endl;
  cout << "par 7 error: " << (  total_fit->GetParError(7) / total_fit->GetParameter(7) ) * 100 << endl;
  cout << "par 8 error: " << (  total_fit->GetParError(8) / total_fit->GetParameter(8) ) * 100 << endl;
  cout << "par 9 error: " << (  total_fit->GetParError(9) / total_fit->GetParameter(9) ) * 100 << endl;
  

   sigma = total_fit->GetParameter(3);

  char name200[500]; 

  if(write)
    {
    
      if(north_arm)
	sprintf(name200,"chisquare_2nd_gaussian/sigmafree/north/N_gaussian_%i.txt", step);	 
      else
	sprintf(name200,"chisquare_2nd_gaussian/sigmafree/south/S_gaussian_%i.txt", step);	  	
     
      // whether north arm or south arm, write to unique filename  
      std::fstream Run15pp_width(name200,std::ofstream::out); 
     Run15pp_width <<  Njpsi << " " <<  Njpsi << " " << Njpsi << " " << sigma << " " << Njpsi << " " << chisquare << " " << chisquare/ndf << " " << fixed << " " << endl;
      Run15pp_width.close();
   
      ave+= sigma;
      weight++;

    } // if width_write
  ///////////////////////////////////////


    }// end for loop jpsi


   ave /= weight;
   
   char const *outfile;
      
   outfile = "hbias.root"; 
    

   TFile *h = new TFile(outfile, "RECREATE");
   h->cd();

     
   hbias_jpsi->SetName("hbias_jpsi");
   hbias_jpsi->Write();

   hbias_psi2s->SetName("hbias_psi2s");
   hbias_psi2s->Write();
     
   h->Close();
   
   cout << " Sigma Average : " << ave << endl;

} // void end macro




  
