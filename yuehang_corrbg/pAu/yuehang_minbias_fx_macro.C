// A(pt) = this macros reads in the area under the curve of each initial fit of yue hang's corr BG 
// f(x|pt) = it also reads in the arrays of each fitted Yue hang corr bg evaluated at each mass bin (arrays have 100 entries)
// the formula is the sum of f(x|pt)*A(pt) divided by the sum of all A(pt) , per arm
// for 200 mev, f(x|pt) is the initital fit not the refit, since no parameter functions were used


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
#include <TPaveText.h>
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

void yuehang_minbias_fx_macro()
{
  bool write_minbias = true;
 
  double x_array[100] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,5.05,5.15,5.25,5.35,5.45,5.55,5.65,5.75,5.85,5.95,6.05,6.15,6.25,6.35,6.45,6.55,6.65,6.75,6.85,6.95,7.05,7.15,7.25,7.35,7.45,7.55,7.65,7.75,7.85,7.95,8.05,8.15,8.25,8.35,8.45,8.55,8.65,8.75,8.85,8.95,9.05,9.15,9.25,9.35,9.45,9.55,9.65,9.75,9.85,9.95};


//for each fit evaluated at each mass bin
#include "pt_int/fit_array_200mev_N.C" 
#include "pt_int/fit_array_200mev_S.C" 
#include "pt_int/fit_array_300mev_N.C" 
#include "pt_int/fit_array_300mev_S.C" 
#include "pt_int/fit_array_1gev_N.C" 
#include "pt_int/fit_array_1gev_S.C" 
  
  bool print_screen = true;

  //bool south_arm = false;
  bool south_arm = true;
  
  TF1 *int_fit_pt;
  TFitResultPtr r;

  double pt_slices[3] = {5,17,7};

  // read in area under curve for 200, 300 and 1,000 mev binwidths
  char basename2[6][500] = {"pt_int/weighted_fits_200mev_S_","pt_int/weighted_fits_300mev_S_","pt_int/weighted_fits_1gev_S_",
			    "pt_int/weighted_fits_200mev_N_","pt_int/weighted_fits_300mev_N_","pt_int/weighted_fits_1gev_N_"};
		  

  int arm_low;

  if(south_arm)
    arm_low = 0;
  else
    arm_low = 1;

 
  vector < vector < vector <std::string> > > bestfit_filename;
    
  for(int j_arm = 0; j_arm < 2; j_arm++)
    {
      vector  < vector <std::string> >  filename_arm;
  
      for(int j_series = 0; j_series < 3; j_series++)
	{
	  vector < std::string >  filename_series;
	  int j = j_arm*3 + j_series;
	  // cout << "j: " << j << endl;
	  
	  for(int i = 0; i < pt_slices[j_series]; i++)
	    {
	      char name[500];
	      sprintf(name,"%s%i%s",basename2[j],i+1,".dat");	     
	      std::string filename_tmp = name;
	      filename_series.push_back(filename_tmp);
	    } 

	  filename_arm.push_back(filename_series);
	}

      bestfit_filename.push_back(filename_arm);

    } 

  for(int arm = 0; arm < 2; arm++)
    {
      for(int series = 0; series < 3; series++)
	{
	  for(int pt = 0; pt < pt_slices[series]; pt++)
	    {
	      cout << bestfit_filename[arm][series][pt].c_str() << endl; 
	    }
	}
    }
   
  double area_corrbg[2][19] ={0};
  double sum_area[2] = {0};
  int counter1 = 0;
  int counter2 = 0;
 
  
  for(int arm = 0; arm < 2; arm++)
    {
      counter1 = 4; 
      counter2 = 17;
     
      for(int series = 0; series < 3; series++) 
	{
	  for(int pt = 0; pt < 19; pt++)  
	    {
	     		 
	      ifstream bestfit_parameters(bestfit_filename[arm][series][pt].c_str()); 
	      
	      if(bestfit_parameters)
		{
		  double a0;
		  bestfit_parameters >> a0;
		  if(series == 0 && pt < 5) 
		    {
		      if(pt < 2)
			area_corrbg[arm][pt] = a0;
		      if( (pt == 3) || (pt == 4) )
			area_corrbg[arm][pt-1] = a0; 
		    }
		  if( (series == 1) && (pt < 17) )
		    {
		      if( (pt > 3) && (pt < 17) ) 
			{
			  area_corrbg[arm][counter1] = a0;
			  counter1++;
			}
		    }
		  if((series == 2) && ( (pt == 5) ||  (pt == 6) ) )
		    {
		      area_corrbg[arm][counter2] = a0;
		      counter2++;
		    } 

		  bestfit_parameters.close();
		}
	      else
		continue;
	     	     
	    }// pt for loop

	} // for series arm
      
      counter1 = 0;
      counter2 = 0;
    } // for arm loop
  

  // counter1 = 0;
  // counter2 = 0;

  // // load south arm seperatley from north arm
  // for(int arm = 1; arm < 2; arm++)
  //   {
  //     counter1 = 4; 
  //     counter2 = 17;
     
  //     for(int series = 0; series < 3; series++) 
  // 	{
  // 	  for(int pt = 0; pt < 19; pt++)  
  // 	    {
	     		 
  // 	      ifstream bestfit_parameters2(bestfit_filename[arm][series][pt].c_str()); 
	      
  // 	      if(bestfit_parameters2)
  // 		{
  // 		  double a0;
  // 		  bestfit_parameters2 >> a0;
  // 		  if(series == 0 && pt < 5) 
  // 		    {
  // 		      if(pt < 2)
  // 			area_corrbg[arm][pt] = a0;
  // 		      if( (pt == 3) || (pt == 4) )
  // 			area_corrbg[arm][pt-1] = a0; 
  // 		    }
  // 		  if( (series == 1) && (pt < 17) )
  // 		    {
  // 		      if( (pt > 3) && (pt < 17) ) 
  // 			{
  // 			  area_corrbg[arm][counter1] = a0;
  // 			  counter1++;
  // 			}
  // 		    }
  // 		  if((series == 2) && ( (pt == 5) ||  (pt == 6) ) )
  // 		    {
  // 		      area_corrbg[arm][counter2] = a0;
  // 		      counter2++;
  // 		    } 

  // 		  bestfit_parameters2.close();
  // 		}
  // 	      else
  // 		continue;
	     	     
  // 	    }
  // 	} // for loop arm
      
  //   } // for series loop
  

  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < 19; k++)
	{
	  cout << "arm " << arm << ", bin: " << k+1 << endl;
	  cout << area_corrbg[arm][k] << endl;
	}
    }
 
  double x_errors[2][100] = {0};  
  double y_errors[2][100] = {0};
  
  double num[2][100] = {0};
  double fit_array[19][2][100] = {0};   
  double Run15pAu_minbias[2][100] = {0};

 int bins = 100;
  
    for(int arm = 0; arm < 2; arm++)
   {
     for(int k = 0; k < 19; k++)
       {
	 sum_area[arm]+= area_corrbg[arm][k];
       }
   }
    
 cout << "sum_area S: " << sum_area[0] << endl;
 cout << "sum_area N: " << sum_area[1] << endl;

 double scale[2] = {0.0025,0.0075};
 
 for(int k = 0; k < 19; k++)
   {
     for(int arm = 0; arm < 2; arm++)
       {
	 for(int bin = 0; bin < 100; bin++)
	   {
	     if(k < 4) // elements 0,1,2,3 ---> match elements 0,1,2,3 in 200 mev array (skipped 400-600 when writing)  
	       {
		 fit_array[k][0][bin] = fit_array_200mev_S[k][bin];
		 fit_array[k][1][bin] = fit_array_200mev_N[k][bin];
	       }
	     if((k < 17) && (k > 3))  // match elements 4-16---> 4-16 in 300 mev array
	       {
		 fit_array[k][0][bin] = fit_array_300mev_S[k][bin];
		 fit_array[k][1][bin] = fit_array_300mev_N[k][bin];
	       }
	     if(k > 16)  // elements 17,18---> 5-6 in 1000 mev array
	       {
		 fit_array[k][0][bin] = fit_array_1gev_S[k-12][bin];
		 fit_array[k][1][bin] = fit_array_1gev_N[k-12][bin];
	       }
	     
	     //   cout << "Sum of area: " << sum_area[arm] << endl;
	     num[arm][bin] += fit_array[k][arm][bin] * area_corrbg[arm][k];
	     //  cout << "numerator: " << num[arm][bin] << endl;
	     Run15pAu_minbias[arm][bin] =   num[arm][bin] /sum_area[arm];
	    	     
	      // cout << "North Arm Ratio: " << Run15pAu_minbias[arm_low][bin] << endl;

	     y_errors[arm][bin] = sqrt(scale[arm]*Run15pAu_minbias[arm][bin]);
	   }
       }
   }
 
 for(int arm = 0; arm < 2; arm++)
   {
     for(int k = 0; k < 19; k++)
       {
	 cout << fit_array[k][arm][1] << endl;
       }
   }
 
 //PLOT MINBIAS AND FIT
 
   TCanvas *c20 = new TCanvas("c20","Corr BG pT Integrated",200,10,700,500);
   gPad->SetGrid();
     
   TGraphErrors *corrbg_int = new TGraphErrors(bins,x_array,Run15pAu_minbias[arm_low],x_errors[arm_low],y_errors[arm_low]);
   
   if(south_arm)
     corrbg_int->SetTitle("Corr BG pT Integrated, South");  
   else
     corrbg_int->SetTitle("Corr BG pT Integrated, North");  
   corrbg_int->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
   corrbg_int->GetXaxis()->SetLabelSize(0.04);
   corrbg_int->GetXaxis()->SetTitleSize(0.04);
   corrbg_int->GetXaxis()->SetTitleOffset(0.9);
   corrbg_int->GetXaxis()->SetLimits(0,5);
   corrbg_int->GetYaxis()->SetLabelSize(0.04);
   corrbg_int->GetYaxis()->SetTitleSize(0.1); 
   corrbg_int->GetYaxis()->SetTitleOffset(0.52);
   //corrbg_int->SetMaximum(0.3);
       // corrbg_int->SetMinimum(pow());
    
   /// assign different colors and draw each bg
   if(south_arm)
     corrbg_int->SetMarkerColor(kViolet);  
   else
     corrbg_int->SetMarkerColor(kAzure);  
   corrbg_int->SetMarkerSize(0.75);
   corrbg_int->SetMarkerStyle(20);
   corrbg_int->Draw("AP");
   

   
   ////////
   double f0,f1,f2,f3,f4;

   if(south_arm == false)
     {
       f0 = 1;
       f1 = 0.01;
       f2 = 100;
       f3 = 10;
       f4 = 10;


       // f0 = -1;
       // f1 = -0.01;
       // f2 = 0.0001;
       // f3 = -100;
       // f4 = 10;
      
     }
   else
     {
       f0 = 1;
       f1 = 0.01;
       f2 = 10;
       f3 = 10;
       f4 = 10;
     }

   double pt_int_array[2][100] = {0};
 
   int_fit_pt = new TF1("int_fit_pt"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",1.75,5);

   int_fit_pt->SetParameter(0, f0);
   int_fit_pt->SetParameter(1, f1);
   int_fit_pt->SetParameter(2, f2);
   int_fit_pt->SetParameter(3, f3);
   int_fit_pt->SetParameter(4, f4);
   
   int_fit_pt->SetLineColor(kRed);
   int_fit_pt->SetLineStyle(5);
   int_fit_pt->SetLineWidth(2);
 

for(int arm = 0; arm < 2; arm++)
     {
       for(int bin = 0; bin < 100; bin++)
	 {
	   //num_ratio[arm][bin] = Run15pp_minbias[arm][bin];
	   pt_int_array[arm][bin] = int_fit_pt->Eval(x_array[bin]);
	  
	 }
     }

   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
 
   r = corrbg_int->Fit(int_fit_pt,"S R");
   
   cout << "Chi Square/NDF: " << int_fit_pt->GetChisquare() << "/" << int_fit_pt->GetNDF() << "" << endl;

   // r->Parameter(0);
   // r->Parameter(1);
   // r->Parameter(2);
   // p3 = r->Parameter(3);
   // p4 = r->Parameter(4);
   
   //  mytext->Draw();
   

 Char_t message[200];
 
 sprintf(message,"#chi^{2}/NDF = %.1f / %.d",int_fit_pt->GetChisquare(),int_fit_pt->GetNDF());
 
 TPaveText *mytext = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
 mytext->SetTextSize(0.035);
 mytext->SetFillColor(0);
 mytext->SetTextAlign(12);
 mytext->AddText(message);
 mytext->Draw();
 
 
   TLatex l3;
   l3.SetTextSize(0.06);
   l3.SetTextAlign(13);
   l3.SetTextColor(4);
   
   char text4[100];
   sprintf(text4,"p + Au");
   
   l3.SetTextAlign(12);
   l3.DrawLatexNDC(0.75, 0.93, text4); //4.4,150
   

   // write pT integrated array for Run15pAu
 if(write_minbias && !south_arm)
    {
      std::fstream fit_array_file("pt_int/Run15pAu_reco_pt_int_N.C", std::ofstream::out);
      fit_array_file << "double Run15pAu_reco_pt_int_N[100] = { " ;
      
      for(int bin = 0; bin < 100; bin++)
	{
	  if(bin < 99)
	    fit_array_file <<  Run15pAu_minbias[1][bin]  <<  ", ";
	  if(bin == 99)
	    fit_array_file <<  Run15pAu_minbias[1][bin] << "};" << endl;
	  
	}
	  
      fit_array_file << endl;
      
    } // write_minbias north
  
 if(write_minbias && south_arm)
    {
      std::fstream fit_array_file("pt_int/Run15pAu_reco_pt_int_S.C", std::ofstream::out);
      fit_array_file << "double Run15pAu_reco_pt_int_S[100] = { " ;
      
      for(int bin = 0; bin < 100; bin++)
	{
	  if(bin < 99)
	    fit_array_file <<  Run15pAu_minbias[0][bin]  <<  ", ";
	  if(bin == 99)
	    fit_array_file <<  Run15pAu_minbias[0][bin] << "};" << endl;
	}

      fit_array_file << endl;

    } // write_minbias south
  

 if(print_screen == true)
    {
      cout << "Run15pAu_pt_int_array[0][" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	  cout << pt_int_array[0][i] << ", ";
      	}
      cout << pt_int_array[0][bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
     
    } // print out
  
  
 
  
}// end void 


