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

void plot_param_macro_south_300()
{
  bool write = true;
  int Series = 0;
 
  //cout << "Enter the series you wish to plot ( '0' for 300 MeV binwidth, enter '1' for 500 MeV binwidth, '2' for 1 GeV binwidth, '3' for 2 GeV binwidth)" << endl;
  //cin >> Series;

    int pt_slices[4] = {17,12,10,3};  
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;

  double results_matrix[4][5][38] = {0};  // [number of series][number of parameters][largest number of bins in one series]

  double pt_center[4][38] = {0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85,3.15,3.45,3.75,4.05,4.35,4.65,4.95,5.5,6.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     1.0,2.0,4.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};




  for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 38; j++)
	{
	  //cout << " pt center: " << pt_center[i][j] << ", for i = " <<  i << ", j = " << j << endl;
	}
    }


  
  char basename[8][500] = {"yuehang_bestfit_300mev/bestfit_parameters_S_","yuehang_bestfit_500mev/bestfit_parameters_S_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_","yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_",

			   "yuehang_bestfit_300mev/bestfit_parameters_N_","yuehang_bestfit_500mev/bestfit_parameters_N_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_N_","yuehang_bestfit_2gev/bestfit_parameters_high_pt_N_"};
  
  vector < vector < vector <std::string> > > bestfit_filename;
  
  for(int j_arm = 0; j_arm < 2; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm;
      
      for(int j_series = 0; j_series < 4; j_series++)
	{
	  vector < std::string >  filename_series;
	  
	  int j = j_arm*4 + j_series;
	  //cout << "j: " << j << endl;
	  
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
	      // cout << "bestfit_filename[" << arm << "][" << series << "][" << pt << "] = " << bestfit_filename[arm][series][pt].c_str() << endl; 
	    }
	}
    }
  
  char name[800];
  char name1000[800];
  double par_array[5];
  double par_error[2][5][38] = {0};
  double pt_error[38] = {0};
  
  int a;  
  int b;
  
  for(int i = 0; i < 20; i++)
    {
      //cout << pt_center[0][i] << endl;
    }

  if(Series == 0)
    {
      a = 0;
      b = 1;
    }

  if(Series == 1)
    {
      a = 1;
      b = 2;
    }
  
  
  if(Series == 2)
    {
      a = 2;
      b = 3;
    }
 
  for(int series = a; series < b; series++) 
    {
      for(int pt = 0; pt < 17; pt++)  
	{
	  for(int i = 0; i < 5; i++)
	   {
	     ifstream bestfit_parameters(bestfit_filename[0][series][pt].c_str());  // [1] is north arm
	     
	     if(bestfit_parameters)
	       {
		 bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		 par_array[0] = p0;
		 par_array[1] = p1;
		 par_array[2] = p2;
		 par_array[3] = p3;
		 par_array[4] = p4;
	
		 par_error[0][0][pt] = e0; // par error [1]  is arm
		 par_error[0][1][pt] = e1;
		 par_error[0][2][pt] = e2;
		 par_error[0][3][pt] = e3;
		 par_error[0][4][pt] = e4;
		 
		 bestfit_parameters.close();
	       }
	     else
	       {
		 par_array[i] = 0;
		 par_error[0][i][pt] = 0.0;
	       }

	     results_matrix[0][i][pt] = par_array[i]; // [0] is series array
	     pt_error[pt] = 0;
	     
	   }

	  if(pt > 3)
	    {
	      //cout << p0 << endl;
	      // cout << p1 << endl;
	      // cout << p3 << endl;
	      cout << p4 << endl;
	    }
	} 
     
   } // for series loop
 
  
  for(int series = a+2; series < b+2; series++) // for 1 GeV
    {
      for(int pt = 5; pt < 7; pt++) 
	{
	  for(int i = 0; i < 5; i++)
	   {
	     ifstream bestfit_parameters2(bestfit_filename[0][series][pt].c_str()); 
	     
	     if(bestfit_parameters2)
	       {
		 bestfit_parameters2 >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		 par_array[0] = p0;
		 par_array[1] = p1;
		 par_array[2] = p2;
		 par_array[3] = p3;
		 par_array[4] = p4;
	
		 //cout << p0 << endl;
		 //cout << p1 << endl;
		 //  cout << p3 << endl;
		  cout << p4 << endl;
		 
		 if(pt == 5)
		   {
		     par_error[0][0][17] = e0;
		     par_error[0][1][17] = e1;
		     par_error[0][2][17] = e2;
		     par_error[0][3][17] = e3;
		     par_error[0][4][17] = e4;
		     results_matrix[0][i][17] = par_array[i]; // series [0] location
		   }
		 if(pt == 6)
		   {
		     par_error[0][0][18] = e0;
		     par_error[0][1][18] = e1;
		     par_error[0][2][18] = e2;
		     par_error[0][3][18] = e3;
		     par_error[0][4][18] = e4;
		     results_matrix[0][i][18] = par_array[i]; // series [0] location
		   }

		 bestfit_parameters2.close();
	       }
	     else
	       {
		 par_array[i] = 0;
		 par_error[0][i][pt] = 0.0;
	       }
	    
	     pt_error[pt] = 0;
	     
	   }
	} // for loop pt
      
    } // for series loop
  

 for(int series = a; series < b; series++)  
   {
     for(int i = 0; i < 5; i++) 
       {
	 for(int pt = 4; pt < pt_slices[Series]; pt++) 
	   {
	     // if(pt == 2)
	     //     cout << "par " << i << endl;
	    
	     // cout << "results_matrix[" << series << "][" << i << "][" << pt << "] = " << results_matrix[series][i][pt] << endl;
	     // cout  << results_matrix[series][i][pt] << endl;
	   }
       }
   }

 
 
 //////////////////////////////////////////////////////////////////////////////////////////////////
 TGraphErrors *gr[2][14] ={0};
 double points;
 
 for(int series = a; series < b; series++)  
   {
     for(int arm = 0; arm < 1; arm++)  
	{
	  for(int i = 0; i < 5; i++) 
	    {
	      
	      points = 19;
	      gr[arm][i] = new TGraphErrors(points,pt_center[series],results_matrix[series][i],pt_error,par_error[arm][i]);
	      // cout <<"points: " << points << ", pt center: " << pt_center[series][5] << ", results matrix: " << results_matrix[series][i][5] << ", pt_error: " << pt_error[5] << ", par_error: " << par_error[arm][i][5] << ", arm: " << arm << ", series: " << series << endl;
	      
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
	      gr[arm][i]->GetXaxis()->SetLimits(1.1,7);
	      gr[arm][i]->SetMinimum(-15);
	      gr[arm][i]->SetMaximum(15);
	      
	    }
	}
    }
 
  char unique[800];
  TLegend *leg[2];  

  TCanvas *c2 = new TCanvas("c2","South",200,10,700,500);
  gPad->SetGrid();
  

  TCanvas *c7;
  TCanvas *c3;
  TCanvas *c4;
 TCanvas *c6;
 TLegend *leg3 = new TLegend(0.51, 0.6, 0.8, 0.9); 

 TFitResultPtr r[4] = {0};

 //////////////////////////////////////////////////

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
           
      for(int series = a; series < b; series++)
	{
	  leg[0] = new TLegend(0.51, 0.6, 0.8, 0.9);  
	  leg[0]->SetFillColor(0); 
	  leg[0]->SetTextSize(0.035);
	  
  	  for(int i = 0; i < 5; i++)
  	    {
  	      sprintf(unique2,"parameter %d",i+1);
	      if(i != 2)
		leg[0]->AddEntry(gr[0][i],unique2, "p");
	    }
	}
      
      
      leg[0]->Draw(); 
      
      
      char unique3[800];
      c3 = new TCanvas("c3","South par a",200,10,700,500);
        
      gPad->SetGrid();
      gr[0][0]->SetTitle("Parameter a");
      gr[0][0]->SetMarkerSize(2.0); 
      gr[0][0]->SetMarkerColor(kRed);
      gr[0][0]->Draw("AP");
      
      c4  = new TCanvas("c4","South par b",200,10,700,500);
           
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
      
      c6  = new TCanvas("c6","South par d",200,10,700,500);
      gPad->SetGrid();
      gr[0][3]->SetTitle("Parameter d");
      gr[0][3]->SetMarkerSize(2.0); 
      gr[0][3]->SetMarkerColor(kAzure);
      gr[0][3]->Draw("AP"); 
      
      c7 = new TCanvas("c7","South par e",200,10,700,500);
        
      gPad->SetGrid();
      gr[0][4]->SetTitle("Parameter e");
      gr[0][4]->SetMarkerSize(2.0); 
      gr[0][4]->SetMarkerColor(kViolet);
      gr[0][4]->Draw("AP"); 
      
      TF1 *par_fit[4];
      


  // for 0.25 - 5 fit range
      double par_list[4][10] = { //-0.1, 0.1, -0.1, 100, -0.01, 10, -0.0001, 0, 0, 0,   // 10.8/8
	-0.01, 0.01, -0.01, 1, -0.01, 10, -0.0001, 0, 0, 0,
      	-10, 1, -1, 10, -0.001, -0.01, -0.001, 0, 0, 0,
       	-100, 100, -100, 100, -0.001, 0.001, -0.001, 0, 0, 0, 
      	-1, 1000, -10, 10000, -0.1, 0.1, 0, 0, 0, 0};
   


      // for 0.5 - 5 fit range
      // double par_list[4][10] = {       -100, 100, -100, 100, -0.01, 0.01, -0.0001, 0.001, -0.01, 10,
      // 	-10, 1, -1, 10, -0.001, -0.01, -0.001, 0.0001, 0, 0,
      //  	-100, 100, -100, 100, -0.001, 0.001, -0.001, 0, 0, 0, 
      // 	-1, 1, -100, 100, -0.001, 0.001, -0.001, 0.001, -0.001, 0};
                  
      par_fit[0] = new TF1("par_fit_0"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1.1,7);
      par_fit[1] = new TF1("par_fit_1"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1.1,7);
      par_fit[2] = new TF1("par_fit_2"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1.1,7);
      par_fit[3] = new TF1("par_fit_3"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1.1,7);
       
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
   double coeff_array[4][15]; // corr bg parameters,arm,coefficient values of f(x)

   // Fill an array with the coefficient values for each parameter function

   for(int i = 0; i < 4; i++)
     {
       for(int par = 0; par < 10; par++)
	 {
	   coeff_array[i][par] = r[i]->Parameter(par);
	 }
     }
   
   
   std::fstream coeff_array_0("yuehang_fit_coeff_300/fit_coeff_S_array_nineteen_point.C", std::ofstream::out);

   coeff_array_0 << "double coeff_par_S[4][10] = { " ;

   if(write)
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




