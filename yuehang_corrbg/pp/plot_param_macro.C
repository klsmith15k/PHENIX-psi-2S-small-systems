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

void plot_param_macro()
{

  int Series;

  cout << "Enter the series you wish to plot ( '0' for 200 MeV binwidth, enter '1' for 500 MeV binwidth, '2' for 1 GeV binwidth, '3' for 2 GeV binwidth)" << endl;
  cin >> Series;

  int pt_slices[4] = {38,20,10,5}; 
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;

  //double results_matrix[2][4][38][5]; // [number of series][arms][largest number of bins in one series][number of parameters]
  double results_matrix[4][2][5][pt_slices[Series]]; // [number of series][arms][largest number of bins in one series][number of parameters]

  double pt_center[4][38] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     1.0,3.0,5.0,7.0,9.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  char basename[8][500] = {"yuehang_bestfit_200mev/bestfit_parameters_S_","yuehang_bestfit_500mev/bestfit_parameters_S_","yuehang_bestfit_1gev/bestfit_parameters_S_","yuehang_bestfit_2gev/bestfit_parameters_S_",

			   "yuehang_bestfit_200mev/bestfit_parameters_N_","yuehang_bestfit_500mev/bestfit_parameters_N_","yuehang_bestfit_1gev/bestfit_parameters_N_","yuehang_bestfit_2gev/bestfit_parameters_N_"};
  
  vector < vector < vector <std::string> > > bestfit_filename;
    
  for(int j_arm = 0; j_arm < 2; j_arm++)
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
  
  for(int arm = 0; arm < 2; arm++)
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
  int upp_limit = 7;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = 2; series < 3; series++) // for series you want
	{
	  for(int pt = 1; pt < pt_slices[series]; pt++) // to plot all bins in series
	    // for(int pt = 0; pt < upp_limit; pt++) // to plot up to a certain cutoff
	    {
	      for(int i = 0; i < 5; i++)
		{
		  sprintf(name,bestfit_filename[arm][series][pt].c_str(),pt+1); 
		  std::fstream bestfit_parameters(name,std::ifstream::in); 
		  if(bestfit_parameters)
		    {
		      bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      par_array[4] = p4;

		      par_error[arm][0][pt] = e0;
		      par_error[arm][1][pt] = e1;
		      par_error[arm][2][pt] = e2;
		      par_error[arm][3][pt] = e3;
		      par_error[arm][4][pt] = e4;
		      
		      bestfit_parameters.close();
		    }
		  else
		    {
		      par_array[i] = 0;
		      par_error[arm][i][pt] = 0.0;
		    }
	
		  results_matrix[series][arm][i][pt] = par_array[i];
		  pt_error[pt] = 0;
		}
	    }
	} // for loop arm
      
    } // for series loop
  
  
  
  TGraphErrors *gr[2][5];
  double points;
  
  int a = 2; // series 2 (pt_binwidth = 1 GeV)
  int b = 3;
    
  TCanvas *c1 = new TCanvas("c1","North",200,10,700,500);
  gPad->SetGrid();
  cout << "HERE" << endl;
  for(int arm = 0; arm < 2; arm++) 
    {
      for(int series = a; series < b; series++)
	{
	  for(int i = 0; i < 5; i++)
	    {
	      //points = pt_slices[series];
	      points = 5;
	      gr[arm][i] = new TGraphErrors(points,pt_center[series],results_matrix[series][arm][i],pt_error,par_error[arm][i]);
	      gr[arm][i]->SetLineColor(i+1);
	      gr[arm][i]->SetLineWidth(2.5);
	      gr[arm][i]->GetXaxis()->SetTitle("p_{T} GeV/c");
	      gr[arm][i]->GetXaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetXaxis()->SetTitleSize(0.05);
	      gr[arm][i]->GetXaxis()->SetTitleOffset(0.9);
	      gr[arm][i]->GetYaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetYaxis()->SetTitleSize(0.1); 
	      gr[arm][i]->GetYaxis()->SetTitleOffset(0.52);
	      gr[arm][i]->GetXaxis()->SetLimits(1.0,5.0);
	      gr[0][i]->SetMinimum(-10);
	      gr[0][i]->SetMaximum(80);
	      //gr[1][i]->SetMinimum(-200);
	      //gr[1][i]->SetMaximum(200);
	      
	    }
	}
    }
    
gr[1][0]->SetMinimum(-20);
gr[1][0]->SetMaximum(200);

  for(int i = 0; i < 5; i++)
    {
      if(i == 0)
	{
	  gr[1][i]->SetTitle("Corr BG Bestfit Parameters, North");
	  gr[1][i]->Draw("AL");
	}
      else
	gr[1][i]->Draw("L");
    }
  
  char unique[800];
  TLegend *leg[2]; 
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = a; series < b; series++)
  	{
  	  leg[1] = new TLegend(0.51, 0.6, 0.8, 0.9);  
  	  leg[1]->SetFillColor(0); 
  	  leg[1]->SetTextSize(0.035);
	  
  	  for(int i = 0; i < 5; i++)
  	    {
  	      sprintf(unique,"parameter %d",i+1);
	      leg[1]->AddEntry(gr[1][i],unique, "l");
	    }
	}
    }

  leg[1]->Draw();  

  TCanvas *c2 = new TCanvas("c2","South",200,10,700,500);
  gPad->SetGrid();
  
  for(int i = 0; i < 5; i++)
    {
      if(i == 0)
	{
	  gr[0][i]->SetTitle("Corr BG Bestfit Parameters, South");
	  gr[0][i]->Draw("AL");
	}
      else
	gr[0][i]->Draw("L");
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
  	      sprintf(unique2,"parameter %d",i+1);
	      leg[0]->AddEntry(gr[0][i],unique2, "l");
	    }
	}
    }

  leg[0]->Draw(); 

  TLegend *leg3 = new TLegend(0.51, 0.6, 0.8, 0.9); 
  char unique3[800];
 
  TCanvas *c3 = new TCanvas("c3","South 2",200,10,700,500);
  gPad->SetGrid();
  
  gr[0][0]->Draw("AL");
  gr[0][1]->Draw("L");
  gr[0][2]->Draw("L");
  //gr[0][0]->SetMaximum(20);
  //gr[0][0]->SetMinimum(-20);

  leg3->SetFillColor(0); 
  leg3->SetTextSize(0.035);
  
  for(int i = 0; i < 3; i++)
    {
      sprintf(unique3,"parameter %d",i+1);
      leg3->AddEntry(gr[0][i],unique3, "l");
    }
  
  leg3->Draw();

  double counter = 0.0;
  for(int arm = 0; arm < 2; arm++) 
    {
      for(int series = a; series < b; series++)
	{
	  for(int i = 0; i < 5; i++)
	    {
	      for(int pt = 0; pt < pt_slices[series]; pt++)
		{
		  //cout << "Arm: " << arm << ", series: " << series << ", parameter: " << i << ", pt bin: " << pt+1 << endl;
		  // cout << results_matrix[series][arm][i][pt] << endl;
		  // counter++;
		  // cout << par_error[arm][i][pt] << endl;
		  if((results_matrix[series][arm][i][pt] > 100) || (par_error[arm][i][pt] > 100))
		    cout << "bad fit in arm " << arm << ", bin " << pt+1 << endl; 
		}
	    }
	}
    }
  //cout << "counter:" << counter << endl;

}// end void macro




  // char unique[800];
  // TLegend *leg[2]; 
  
  // for(int arm = 0; arm < 2; arm++)
  //   {
  //     for(int series = a; series < b; series++)
  // 	{
  // 	  leg[1] = new TLegend(0.51, 0.6, 0.8, 0.9);  
  // 	  leg[1]->SetFillColor(0); 
  // 	  leg[1]->SetTextSize(0.035);
	  
  // 	  for(int pt = 0; pt < pt_slices[series] - 1; pt++)
  // 	    {
  // 	      sprintf(unique,"bin center: %.1f GeV/c",pt_center[series][pt]);
  // 	      // leg[1]->AddEntry(gr[1][series],unique, "p");
  // 	    }
  // 	}
