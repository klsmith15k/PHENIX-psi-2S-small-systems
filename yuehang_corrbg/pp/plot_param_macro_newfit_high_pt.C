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

void plot_param_macro_newfit_high_pt()
{
  // bool south_arm = false;
  bool write = true;
  int Series;
 
  cout << "Enter the series you wish to plot ( '0' for 200 MeV binwidth, enter '1' for 500 MeV binwidth, '2' for 1 GeV binwidth, '3' for 2 GeV binwidth)" << endl;
  cin >> Series;

  // int pt_slices[4] = {38,14,10,5}; 
  int pt_slices[4] = {38,12,10,5};  
  double p0,p1,p2,p3,p4,e0,e1,e2,e3,e4;

  //double results_matrix[2][4][38][5]; // [number of series][arms][largest number of bins in one series][number of parameters]
  double results_matrix[4][5][38] = {0};  // [number of series][largest number of bins in one series][number of parameters]

  double pt_center[4][38] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     // 0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.5,6.0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     1.0,3.0,5.0,7.0,9.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  char basename[8][500] = {"yuehang_bestfit_200mev/bestfit_parameters_S_","yuehang_bestfit_500mev/bestfit_parameters_S_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_S_","yuehang_bestfit_2gev/bestfit_parameters_high_pt_S_",

			   "yuehang_bestfit_200mev/bestfit_parameters_N_","yuehang_bestfit_500mev/bestfit_parameters_N_","yuehang_bestfit_1gev/bestfit_parameters_high_pt_N_","yuehang_bestfit_2gev/bestfit_parameters_high_pt_N_"};
  
  vector < vector < vector <std::string> > > bestfit_filename;
    
  for(int j_arm = 0; j_arm < 2; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm;
      
      for(int j_series = 0; j_series < 4; j_series++)
	{
	  vector < std::string >  filename_series;

	  int j = j_arm*4 + j_series;
	  //  cout << "j: " << j << endl;
	  
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
	      // cout << bestfit_filename[arm][series][pt].c_str() << ", for arm: " << arm << ", series: " << series << ", pt: " << pt << endl; 
	    }
       }
   }

  char name[800];
  char name1000[800];
  double par_array[5];// just a place holder, so 1D is fine
  double par_error[2][5][38] = {0};
  double pt_error[38] = {0};
  
  int a;  // series 2 (pt_binwidth = 1 GeV), series 1 pt_binwidth = 500 MeV
  int b;

  if(Series == 1)
    {
      a = 1;
      b = 2;
    }

  // for ten points
  /*
  for(int series = a; series < b; series++) 
    {
      // for(int pt = 0; pt < 10; pt++)  
       	for(int pt = 0; pt < 8; pt++)  
	{
	  for(int i = 0; i < 5; i++)
	   {
	     ifstream bestfit_parameters(bestfit_filename[0][series][pt].c_str()); 
	     
	     if(bestfit_parameters)
	       {
		 bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		 par_array[0] = p0;
		 par_array[1] = p1;
		 par_array[2] = p2;
		 par_array[3] = p3;
		 par_array[4] = p4;
		 
		 //cout << p0 <<  " " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << endl;

		 par_error[0][0][pt] = e0;
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

	     results_matrix[1][i][pt] = par_array[i];
	     pt_error[pt] = 0;
	     
	   }

	  cout << " " << p0 << endl;
	  // cout << " " << p1 << endl;
	  // cout << " " << p3 << endl;
	  // cout << " " << p4 << endl;

	} 
     
   } // for series loop
  */

  // for twelve points
 for(int series = a; series < b; series++) // for series you want
    {
     for(int arm = 0; arm < 1; arm++)
	{
	  for(int pt = 0; pt < 10; pt++) 
	    {
	      for(int i = 0; i < 5; i++)
		{
		 
		  ifstream bestfit_parameters(bestfit_filename[0][series][pt].c_str()); 

		  if(bestfit_parameters)
		    {
		      bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      par_array[4] = p4;

		      // cout << "par b: " << p1 << ", for bin: " << pt+1 << endl;

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
	
		  results_matrix[1][i][pt] = par_array[i];
		  pt_error[pt] = 0;
		}
	    }

	  //cout << " " << p0 << endl;
	  // cout << " " << p1 << endl;
	  // cout << " " << p3 << endl;
	   cout << " " << p4 << endl;

	} // for loop arm
      
    } // for series loop
  
 // for twleve points
for(int series = a+1; series < b+1; series++) // for series you want
    {
     for(int arm = 0; arm < 1; arm++)
	{
	  for(int pt = 5; pt < 7; pt++) 
	    {
	      for(int i = 0; i < 5; i++)
		{
		 
		  ifstream bestfit_parameters(bestfit_filename[arm][series][pt].c_str()); 

		  if(bestfit_parameters)
		    {
		      bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		      par_array[0] = p0;
		      par_array[1] = p1;
		      par_array[2] = p2;
		      par_array[3] = p3;
		      par_array[4] = p4;

		      cout << "par b: " << p1 << ", for bin: " << pt+1 << endl;

		      if(pt == 5)
			{
			  par_error[arm][0][10] = e0;
			  par_error[arm][1][10] = e1;
			  par_error[arm][2][10] = e2;
			  par_error[arm][3][10] = e3;
			  par_error[arm][4][10] = e4;
			}
		      if(pt == 6)
			{
			  par_error[arm][0][11] = e0;
			  par_error[arm][1][11] = e1;
			  par_error[arm][2][11] = e2;
			  par_error[arm][3][11] = e3;
			  par_error[arm][4][11] = e4;
			}


		      bestfit_parameters.close();
		    }
		  else
		    {
		      par_array[i] = 0;
		      par_error[arm][i][pt] = 0.0;
		    }
		  if(pt == 5)
		  results_matrix[1][i][10] = par_array[i];

		  if(pt == 6)
		  results_matrix[1][i][11] = par_array[i];

		  pt_error[pt] = 0;

		}
	    }
	} // for loop arm
      
    } // for series loop


// for ten points
/*
  for(int series = a+1; series < b+1; series++) // for 1 GeV 
    {
      for(int pt = 4; pt < 5; pt++) 
	{
	  for(int i = 0; i < 5; i++)
	   {
	     ifstream bestfit_parameters(bestfit_filename[0][series][pt].c_str()); 
	     
	     if(bestfit_parameters)
	       {
		 bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		 par_array[0] = p0;
		 par_array[1] = p1;
		 par_array[2] = p2;
		 par_array[3] = p3;
		 par_array[4] = p4;
		 
		 if(pt == 4)
		   {
		     par_error[1][0][8] = e0;
		     par_error[1][1][8] = e1;
		     par_error[1][2][8] = e2;
		     par_error[1][3][8] = e3;
		     par_error[1][4][8] = e4;
		     results_matrix[1][i][8] = par_array[i];
		   }
	
		 bestfit_parameters.close();
	       }
	     else
	       {
		 par_array[i] = 0;
		 par_error[0][i][pt] = 0.0;
	       }
	  
	     pt_error[pt] = 0;
	     
	   }

	  cout << " " << p0 << endl;
	  //cout << " " << p1 << endl;
	  // cout << " " << p3 << endl;
	  //cout << " " << p4 << endl;

       }
     
   } // for series loop


  for(int series = a+2; series < b+2; series++) // for 2 GeV 
    {
      for(int pt = 3; pt < 4; pt++) 
	{
	  for(int i = 0; i < 5; i++)
	   {
	     ifstream bestfit_parameters(bestfit_filename[0][series][pt].c_str()); 
	     
	     if(bestfit_parameters)
	       {
		 bestfit_parameters >> p0 >> e0 >> p1 >> e1 >> p2 >> e2 >> p3 >> e3 >> p4 >> e4;
		 par_array[0] = p0;
		 par_array[1] = p1;
		 par_array[2] = p2;
		 par_array[3] = p3;
		 par_array[4] = p4;
		 
		 if(pt == 3)
		   {
		     par_error[0][0][9] = e0;
		     par_error[0][1][9] = e1;
		     par_error[0][2][9] = e2;
		     par_error[0][3][9] = e3;
		     par_error[0][4][9] = e4;
		     results_matrix[1][i][9] = par_array[i];
		   }
	
		 bestfit_parameters.close();
	       }
	     else
	       {
		 par_array[i] = 0;
		 par_error[0][i][pt] = 0.0;
	       }
	  
	     pt_error[pt] = 0;
	     
	   }

	  cout << " " << p0 << endl;
	  //cout << " " << p1 << endl;
	  //cout << " " << p3 << endl;
	  //cout << " " << p4 << endl;

       }
     
   } // for series loop
*/

  for(int series = 1; series < 2; series++)  
    {
      for(int arm = 0; arm < 1; arm++)  
	{
	  for(int i = 0; i < 5; i++) 
	    {
	      for(int pt = 2; pt < pt_slices[Series]; pt++) 
		{
		  //cout << "results_matrix[" << series << "][" <<  i << "][" << pt << "] = " << results_matrix[series][i][pt] << endl;
		  // if(pt == 2)
		  //   cout << "par " << i << endl;
		  //  cout  << results_matrix[series][i][pt] << endl;
		}
	    }
	}
    }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  TGraphErrors *gr[2][14];
  double points;

  for(int series = a; series < b; series++)  
    {
      for(int arm = 0; arm < 1; arm++)  
	{
	  for(int i = 0; i < 5; i++) 
	    {
	      //points = pt_slices[series];
	      points = 12;
	      //points = 10;
	      gr[arm][i] = new TGraphErrors(points,pt_center[series],results_matrix[series][i],pt_error,par_error[arm][i]);
	      // cout <<"points: " << points << ", pt center: " << pt_center[series][6] << ", results matrix: " << results_matrix[series][i][6] << ", pt_error: " << pt_error[6] << ", par_error: " << par_error[arm][i][6] << ", arm: " << arm << ", series: " << series << endl;
	      
	      gr[arm][i]->SetMarkerStyle(20);	      
	      gr[arm][i]->SetMarkerSize(0.25);	     
	      gr[arm][i]->SetMarkerColor(i+1);
	      // if(i == 5)
	      //	 gr[arm][i]->SetMarkerColor(14);
	      gr[arm][i]->GetXaxis()->SetTitle("p_{T} GeV/c");
	      gr[arm][i]->GetXaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetXaxis()->SetTitleSize(0.05);
	      gr[arm][i]->GetXaxis()->SetTitleOffset(0.9);
	      gr[arm][i]->GetYaxis()->SetLabelSize(0.06);
	      gr[arm][i]->GetYaxis()->SetTitleSize(0.1); 
	      gr[arm][i]->GetYaxis()->SetTitleOffset(0.52);
	      // gr[arm][i]->GetXaxis()->SetLimits(0,5.0);
	      gr[arm][i]->GetXaxis()->SetLimits(1,7);
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
  for(int arm = 0; arm < 1; arm++)
    {
      for(int i = 0; i < 5; i++)
	{
	  if(i == 0)
	    {
	      gr[arm][i]->SetTitle("Corr BG Bestfit Parameters, South");
	      gr[arm][i]->Draw("AP");
	    }
	  if(i != 0 && i != 2)
	    gr[arm][i]->Draw("P");
	}
      
      char unique2[800];
           
      for(int series = a; series < b; series++)
	{
	  leg[arm] = new TLegend(0.51, 0.6, 0.8, 0.9);  
	  leg[arm]->SetFillColor(0); 
	  leg[arm]->SetTextSize(0.035);
	  
  	  for(int i = 0; i < 5; i++)
  	    {
  	      sprintf(unique2,"parameter %d",i+1);
	      if(i != 2)
		leg[arm]->AddEntry(gr[arm][i],unique2, "p");
	    }
	}
      
      
      leg[arm]->Draw(); 
      
      
      char unique3[800];
      c3 = new TCanvas("c3","South par a",200,10,700,500);
        
      gPad->SetGrid();
      gr[arm][0]->SetTitle("Parameter a");
      gr[arm][0]->SetMarkerSize(2.0); 
      gr[arm][0]->SetMarkerColor(kRed);
      gr[arm][0]->Draw("AP");
      
      c4  = new TCanvas("c4","South par b",200,10,700,500);
           
      gPad->SetGrid();
      gr[arm][1]->SetTitle("Parameter b");
      gr[arm][1]->SetMarkerSize(2.0); 
      gr[arm][1]->SetMarkerColor(kYellow);
      gr[arm][1]->Draw("AP"); 
      
      // TCanvas *c5 = new TCanvas("c5","South par c",200,10,700,500);
      // gPad->SetGrid();
      // gr[arm][2]->SetMarkerSize(2.0); 
      // gr[arm][2]->SetMarkerColor(kGreen);
      // gr[arm][2]->Draw("AP"); //par c
      
      c6  = new TCanvas("c6","South par d",200,10,700,500);
      gPad->SetGrid();
      gr[arm][3]->SetTitle("Parameter d");
      gr[arm][3]->SetMarkerSize(2.0); 
      gr[arm][3]->SetMarkerColor(kAzure);
      gr[arm][3]->Draw("AP"); 
      
      c7 = new TCanvas("c7","South par e",200,10,700,500);
        
      gPad->SetGrid();
      gr[arm][4]->SetTitle("Parameter e");
      gr[arm][4]->SetMarkerSize(2.0); 
      gr[arm][4]->SetMarkerColor(kViolet);
      gr[arm][4]->Draw("AP"); 
      
      TF1 *par_fit[4];
      


      // double par_list[4][10] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
      // 				1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
      // 				500, -1500, 1500, -700, 160, -20, 0.9, 0, 0, 0,
      // 				-60, 150, -120, 50, -10, 1, -1, 0, 0, 0}; 

      double par_list[4][10] = {25, -50, 40, -15, 3, -0.3, 0.01, 1, 0, 0,
      				3, 6, 5, -1, 0.5, -0.05, 0.003, 0, 0, 0,
      				500, -1500, 1500, -700, 160, -20, 0.9, 0, 0, 0,
      				-60, 150, -120, 50, -10, 0, 0, 0, 0, 0}; 
    
          // double par_list[4][10] = {10, -25, 20, -5, 0, 0, 0, 0, 0, 0,
      	  // 			//	0.5, 1, -10, 0, 0, 0, 0, 0, 0, 0,
      	  // 			0.01, -0.05, 0.01, -0.01, 0.001, -0.001, 0, 0, 0, 0,   // 0.72/1 chisquare
      	  // 			50, -50, 25, -5, 0.1, 0, 0, 0, 0, 0,
      	  // 			10, -5, 1, -0.1, 0.01, 0, 0, 0, 0, 0};


    

  
      par_fit[0] = new TF1("par_fit_0"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1,7);
      par_fit[1] = new TF1("par_fit_1"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1,7);
      par_fit[2] = new TF1("par_fit_2"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1,7);
      par_fit[3] = new TF1("par_fit_3"," [0] + [1]*x + [2]*x*x + [3]*x*x*x +[4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x*x",1,7);
      
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

      r[0] = gr[arm][0]->Fit(par_fit[0],"S R");
      r[1] = gr[arm][1]->Fit(par_fit[1],"S R");
      r[2] = gr[arm][3]->Fit(par_fit[2],"S R");
      r[3] = gr[arm][4]->Fit(par_fit[3],"S R");	 
      
      cout << "Fit a Chi Square/NDF: " << par_fit[0]->GetChisquare() << "/" << par_fit[0]->GetNDF() << "" << endl;
      cout << "Fit b Chi Square/NDF: " << par_fit[1]->GetChisquare() << "/" << par_fit[1]->GetNDF() << "" << endl;
      cout << "Fit d Chi Square/NDF: " << par_fit[2]->GetChisquare() << "/" << par_fit[2]->GetNDF() << "" << endl;
      cout << "Fit e Chi Square/NDF: " << par_fit[3]->GetChisquare() << "/" << par_fit[3]->GetNDF() << "" << endl;
      
    }


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

   // std::fstream coeff_array_0("yuehang_fit_coeff_500/fit_coeff_S_array.C", std::ofstream::out);

   // coeff_array_0 << "double coeff_par_S[4][10] = { " ;

   // if(write) // only looking at south arm now
   //   {
   //     for(int i = 0; i < 4; i++)
   // 	 {
   // 	   for(int par = 0; par < 10; par++)
   // 	     {
   // 	       if((par < 10) && (i < 3))
   // 		 coeff_array_0 <<  coeff_array[i][par]  <<  ", ";
   // 	       if((i == 3) && (par < 9))
   // 		 coeff_array_0 <<  coeff_array[i][par]  <<  ", ";
   // 	       if((i == 3) && (par == 9))
   // 		 coeff_array_0 <<  coeff_array[i][par] << "};" << endl;
   // 	     }
   // 	   coeff_array_0 << endl;
   // 	 }
   //   }

}// end void macro







      // Int_t fitStatus_a = r[0];  
      // Int_t fitStatus_b = r[1];  
      // Int_t fitStatus_c = r[2];  
      // Int_t fitStatus_d = r[3];  

      //r->Print("V");


