// SEE LINE 111


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

void plot_param_macro_newfit_low_pt()
{
  bool write = true;
  int Series;
 
  cout << "Enter the series you wish to plot ( '0' for 100 MeV binwidth, '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth)" << endl;
  cin >> Series;

  int pt_slices[4] = {10,5,12,10}; 
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;

  //double results_matrix[2][4][38][5]; // [number of series][arms][largest number of bins in one series][number of parameters]
  double results_matrix[4][2][5][pt_slices[Series]]; // [number of series][arms][largest number of bins in one series][number of parameters]

  double pt_center[4][38] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  char basename[8][500] = {"yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_","yuehang_bestfit_200mev/bestfit_parameters_low_pt_S_","yuehang_bestfit_500mev/bestfit_parameters_low_pt_S_","yuehang_bestfit_1gev/bestfit_parameters_low_pt_S_",
			   
			   "yuehang_bestfit_100mev/bestfit_parameters_extra_low_pt_S_","yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_","yuehang_bestfit_500mev/bestfit_parameters_low_pt_N_","yuehang_bestfit_1gev/bestfit_parameters_low_pt_N_"};
  
  vector < vector < vector <std::string> > > bestfit_filename;
  
  for(int j_arm = 0; j_arm < 1; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm;
      
      for(int j_series = 0; j_series < 4; j_series++)
	{
	  vector < std::string >  filename_series;
	  
	  int j = j_arm*4 + j_series;
	  cout << "j: " << j << endl;
	  
	  for(int i = 0; i < pt_slices[j_series]; i++)
	    {
	      char name[500];
	      sprintf(name,"%s%i%s",basename[j],i+1,".dat");	     
	      std::string filename_tmp = name;
	      filename_series.push_back(filename_tmp);
	    } 
	  
	  filename_arm.push_back(filename_series);
	}
      
      bestfit_filename.push_back(filename_arm);
      
    } 
  
  for(int arm = 0; arm < 1; arm++)
    {
      for(int series = 0; series < 4; series++)
	{
	  for(int pt = 0; pt < pt_slices[series]; pt++)
	    {
	      cout << bestfit_filename[arm][series][pt].c_str() << endl; 
	    }
	}
    }
  
  char name[800];
  double par_array[5];// just a place holder, so 1D is fine
  double par_error[2][5][pt_slices[Series]];
  double pt_error[pt_slices[Series]]; 
  // int upp_limit = 10; 
  // int upp_limit = 8; 
  
  for(int arm = 0; arm < 1; arm++)
    {
      for(int series = 1; series < 2; series++) // for series you want
	{
	  for(int pt = 0; pt < pt_slices[Series]; pt++)  
	    {
	      for(int i = 0; i < 5; i++)
		{
		  sprintf(name,bestfit_filename[arm][series][pt].c_str(),pt+1); 
		  
		  //cout << bestfit_filename[arm][series][pt].c_str() << endl;
		  
		  std::fstream bestfit_parameters(name,std::ifstream::in); 

		  if(bestfit_parameters)
		    {
		      bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      cout << "p3 " << p3 << endl;
		      par_array[4] = p4;

		      par_error[arm][0][pt] = e0;
		      par_error[arm][1][pt] = e1;
		      par_error[arm][2][pt] = e2;
		      par_error[arm][3][pt] = e3;
		      par_error[arm][4][pt] = e4;
		      
		      //cout << par_array[4] << endl;

		      bestfit_parameters.close();
		    }
		  else
		    {
		      //par_array[i] = 0;
		      //par_error[arm][i][pt] = 0.0;
		    }
	
		  results_matrix[series][arm][i][pt] = par_array[i];
		  pt_error[pt] = 0;
		  // cout << results_matrix[series][arm][i][pt] << ", for arm: " << arm << ", series:" << series << ", pt: " << pt << ", parameter: " << i << endl;
		}
	    }
	} // for loop arm
      
    } // for series loop
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  TGraphErrors *gr[2][14];
  double points;
  
  int a = 1;
  int b = 2;
    
  for(int arm = 0; arm < 1; arm++)  // only plotting south arm now
    {
      for(int series = a; series < b; series++)
	{
	  for(int i = 0; i < 5; i++) 
	    {
	      points = pt_slices[series];
	      //points = upp_limit;
	      gr[arm][i] = new TGraphErrors(points,pt_center[series],results_matrix[series][arm][i],pt_error,par_error[arm][i]);
	      gr[arm][i]->SetMarkerStyle(20);	      
	      gr[arm][i]->SetMarkerSize(0.25);	     
	      gr[arm][i]->SetMarkerColor(i+1);
	     
	      gr[arm][i]->GetXaxis()->SetTitle("p_{T} GeV/c");
	      gr[arm][i]->GetXaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetXaxis()->SetTitleSize(0.05);
	      gr[arm][i]->GetXaxis()->SetTitleOffset(0.9);
	      gr[arm][i]->GetYaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetYaxis()->SetTitleSize(0.1); 
	      gr[arm][i]->GetYaxis()->SetTitleOffset(0.52);
	      // gr[arm][i]->GetXaxis()->SetLimits(0,5.0);
	      gr[arm][i]->GetXaxis()->SetLimits(0,1);
	      gr[arm][i]->SetMinimum(-10);
	      gr[arm][i]->SetMaximum(10);
	      
	    }
	}
    }
 
  char unique[800];
  TLegend *leg[2];  

  TCanvas *c2 = new TCanvas("c2","South",200,10,700,500);
  gPad->SetGrid();
  
  for(int i = 0; i < 5; i++)
    {
      if(i == 0)
	{
	  gr[0][i]->SetTitle("Corr BG Bestfit Parameters, South");
	  gr[0][i]->Draw("AP");
	}
      if(i != 0 && i != 2)
	gr[0][i]->Draw("P");
    }
  
  char unique2[800];
    
  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = a; series < b; series++)
  	{
  	  leg[0] = new TLegend(0.51, 0.6, 0.8, 0.9);  
  	  leg[0]->SetFillColor(0); 
  	  leg[0]->SetTextSize(0.035);
	  
  	  for(int i = 0; i < 5; i++)
  	    {
  	      sprintf(unique2,"parameter %d",i);
	      if(i != 2)
		leg[0]->AddEntry(gr[0][i],unique2, "p");
	    }
	}
    }

  leg[0]->Draw(); 

  TLegend *leg3 = new TLegend(0.51, 0.6, 0.8, 0.9); 
  char unique3[800];
 
  TCanvas *c3 = new TCanvas("c3","South par a",200,10,700,500);
  gPad->SetGrid();
  gr[0][0]->SetTitle("Parameter a");
  gr[0][0]->SetMarkerSize(2.0); 
  gr[0][0]->SetMarkerColor(kRed);
  gr[0][0]->Draw("AP");

  TCanvas *c4 = new TCanvas("c4","South par b",200,10,700,500);
  gPad->SetGrid();
  gr[0][1]->SetTitle("Parameter b");
  gr[0][1]->SetMarkerSize(2.0); 
  gr[0][1]->SetMarkerColor(kYellow);
  gr[0][1]->Draw("AP"); 

  // TCanvas *c5 = new TCanvas("c5","South par c",200,10,700,500);
  // gPad->SetGrid();
  // gr[0][2]->SetMarkerSize(2.0); 
  // gr[0][2]->SetMarkerColor(kGreen);
  // gr[0][2]->Draw("AP"); //par c

  TCanvas *c6 = new TCanvas("c6","South par d",200,10,700,500);
  gPad->SetGrid();
  gr[0][3]->SetTitle("Parameter d");
  gr[0][3]->SetMarkerSize(2.0); 
  gr[0][3]->SetMarkerColor(kAzure);
  gr[0][3]->Draw("AP"); 

  TCanvas *c7 = new TCanvas("c7","South par e",200,10,700,500);
  gPad->SetGrid();
  gr[0][4]->SetTitle("Parameter e");
  gr[0][4]->SetMarkerSize(2.0); 
  gr[0][4]->SetMarkerColor(kViolet);
  gr[0][4]->Draw("AP"); 

  TF1 *par_fit[4];
  TFitResultPtr r[4];
  double par_list[4][10] = {1, 10, -100, 10, 0, 0, 0, 0, 0, 0,
			    -10, 200, -300, 200, 0, 0, 0, 0, 0, 0,
			    1*pow(10,4),  2*pow(10,5), -6*pow(10,5), 4*pow(10,5),0 , 0, 0, 0, 0,0, 
			    1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
  
  par_fit[0] = new TF1("par_fit_0"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",0,1);
  par_fit[1] = new TF1("par_fit_1"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",0,1);
  par_fit[2] = new TF1("par_fit_2"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",0,1);
  par_fit[3] = new TF1("par_fit_3"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",0,1);

  for(int i = 0; i < 4; i++)
    {
      for(int par = 0; par < 10; par++)
	{
	  par_fit[i]->SetParameter(par, par_list[i][par]);
	  if(par_list[i][par] == 0)
	    par_fit[i]->FixParameter(par, par_list[i][par]);
      
	  par_fit[0]->SetLineColor(kRed);
	  par_fit[1]->SetLineColor(kYellow);
	  par_fit[2]->SetLineColor(kBlue);
	  par_fit[3]->SetLineColor(kViolet);

	  par_fit[i]->SetLineStyle(5);
	  par_fit[i]->SetLineWidth(2);

	}
    }
  
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
  
  r[0] = gr[0][0]->Fit(par_fit[0],"S R");
  r[1] = gr[0][1]->Fit(par_fit[1],"S R");
  r[2] = gr[0][3]->Fit(par_fit[2],"S R");
  r[3] = gr[0][4]->Fit(par_fit[3],"S R");	 
  
  cout << "Fit a Chi Square/NDF: " << par_fit[0]->GetChisquare() << "/" << par_fit[0]->GetNDF() << "" << endl;
  cout << "Fit b Chi Square/NDF: " << par_fit[1]->GetChisquare() << "/" << par_fit[1]->GetNDF() << "" << endl;
  cout << "Fit d Chi Square/NDF: " << par_fit[2]->GetChisquare() << "/" << par_fit[2]->GetNDF() << "" << endl;
  cout << "Fit e Chi Square/NDF: " << par_fit[3]->GetChisquare() << "/" << par_fit[3]->GetNDF() << "" << endl;
  

  //////////////////////////////////////////////// write out the bestfit results for the coefficients of each function
   char unique_par1[500];
   char unique_par0[500];
   double coeff_array[4][12]; // corr bg parameters,arm,coefficient values of f(x)
   
   // Fill an array with the coefficient values for each parameter function

   for(int i = 0; i < 4; i++)
     {
       for(int par = 0; par < 10; par++)
	 {
	   coeff_array[i][par] = r[i]->Parameter(par);
	 }
     }

   std::fstream coeff_array_0("yuehang_fit_coeff_500/fit_coeff_low_pt_S_array.C", std::ofstream::out);
     
   coeff_array_0 << "double coeff_par[4][10] = { " ;

   if(write) // only looking at south arm now
     {
       for(int i = 0; i < 4; i++)
	 {
	   for(int par = 0; par < 10; par++)
	     {
	       if((par < 10) && (i < 3))
		 coeff_array_0 <<  coeff_array[i][par]  <<  ", ";
	       if((i == 3) && (par < 9))
		 coeff_array_0 <<  coeff_array[i][par]  <<  ", ";
	       if((i == 3) && (par == 9))
		 coeff_array_0 <<  coeff_array[i][par] << "};" << endl;
	     }
	   coeff_array_0 << endl;
	 }
     }

}// end void macro







      // Int_t fitStatus_a = r[0];  
      // Int_t fitStatus_b = r[1];  
      // Int_t fitStatus_c = r[2];  
      // Int_t fitStatus_d = r[3];  

      //r->Print("V");


