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

void yuehang_minbias_fx_macro_pp_xuan_ratio()
{

  bool pt_direct_int = false;
  bool pt_mass_ratio = true;
 
  double x_array[100] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,5.05,5.15,5.25,5.35,5.45,5.55,5.65,5.75,5.85,5.95,6.05,6.15,6.25,6.35,6.45,6.55,6.65,6.75,6.85,6.95,7.05,7.15,7.25,7.35,7.45,7.55,7.65,7.75,7.85,7.95,8.05,8.15,8.25,8.35,8.45,8.55,8.65,8.75,8.85,8.95,9.05,9.15,9.25,9.35,9.45,9.55,9.65,9.75,9.85,9.95};

  // initial fits evaluated at each mass bin
#include "yuehang_pt_int/yuehang_direct_pt_int_500mev_N.C" 
#include "yuehang_pt_int/yuehang_direct_pt_int_500mev_S.C" 
#include "yuehang_pt_int/yuehang_direct_pt_int_1gev_N.C" 
#include "yuehang_pt_int/yuehang_direct_pt_int_1gev_S.C" 

// refits evaluated at each mass bin
#include "yuehang_pt_int/fit_array_200mev_N.C" 
#include "yuehang_pt_int/fit_array_200mev_S.C" 
#include "yuehang_pt_int/fit_array_500mev_N.C" 
#include "yuehang_pt_int/fit_array_500mev_S.C" 
#include "yuehang_pt_int/fit_array_1gev_N.C" 
#include "yuehang_pt_int/fit_array_1gev_S.C" 
  
  //bool south_arm = false;
  bool south_arm = true;

  bool print_screen = true;
  
  TF1 *int_fit_pt;
  TFitResultPtr r;

  double pt_slices[3] = {5,10,7};

  // read in area under curve for 200, 500 and 1,000 mev binwidths
  char basename2[6][500] = {"yuehang_pt_int/weighted_fits_200mev_S_","yuehang_pt_int/weighted_fits_500mev_S_","yuehang_pt_int/weighted_fits_1gev_S_",
			    "yuehang_pt_int/weighted_fits_200mev_N_","yuehang_pt_int/weighted_fits_500mev_N_","yuehang_pt_int/weighted_fits_1gev_N_"};

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
      counter2 = 12;
     
      for(int series = 0; series < 3; series++) 
	{
	  for(int pt = 0; pt < 14; pt++)  
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
		  if( (series == 1) && (pt < 10) )  // pt is place in array file that is read in
		    {
		      if( (pt > 1) && (pt < 10) )  // counter is place in array file being created
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
  
  int counter = 0;

  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < 14; k++)
	{
	  cout << "arm " << arm << ", bin: " << k+1 << endl;
	  cout << area_corrbg[arm][k] << endl;
	}
    }
 
  double x_errors[2][100] = {0};  
  double y_errors[2][100] = {0};
  
  double num[2][100] = {0};
  double fit_array[14][2][100] = {0};  
  double Run15pp_minbias[2][100] = {0};

  int bins = 100;
  
    for(int arm = 0; arm < 2; arm++)
   {
     for(int k = 0; k < 14; k++)
       {
	 sum_area[arm]+= area_corrbg[arm][k];
       }
   }
    
 cout << "sum_area S: " << sum_area[0] << endl;
 cout << "sum_area N: " << sum_area[1] << endl;

 double scale[2] = {0.01305,0.01393};
 
 for(int k = 0; k < 14; k++)
   {
     for(int bin = 0; bin < 100; bin++)
       {
	 if(k < 4) // elements 0,1,2,3 == elements 0,1,2,3 (since 400-600 was skipped)
	   {
	     fit_array[k][0][bin] = fit_array_200mev_S[k][bin];
	     fit_array[k][1][bin] = fit_array_200mev_N[k][bin];
	   }
	 if((k < 12) && (k > 3))  // elements 4,5,6,7,8,9,10,11 == elements 2,3,4,5,6,7,8,9   --> k-2
	   {
	     fit_array[k][0][bin] = fit_array_500mev_S[k-2][bin];
	     fit_array[k][1][bin] = fit_array_500mev_N[k-2][bin];
	   }
	 if(k > 11)  // elements 12,13 == elements 5,6  --> k-7
	   {
	     fit_array[k][0][bin] = fit_array_1gev_S[k-7][bin];
	     fit_array[k][1][bin] = fit_array_1gev_N[k-7][bin];
	
	     //  cout << "fit_array_1gev_N : " << fit_array_1gev_N[k-7][bin] << ", for k-7: " << k-7 << ", and bin: " << bin << endl;
	      cout << "fit_array_1gev_S : " << fit_array_1gev_S[k-7][bin] << ", for k-7: " << k-7 << ", and bin: " << bin << endl;
	   }
	 
	 num[0][bin] += fit_array[k][0][bin] * area_corrbg[0][k];
	 num[1][bin] += fit_array[k][1][bin] * area_corrbg[1][k];

	   } // for loop bin
       } // for loop k
 
 for(int arm  = 0; arm < 2; arm++)
   {
     for(int bin = 0; bin < 100; bin++)
       {
	  Run15pp_minbias[arm][bin] =   num[arm][bin] /sum_area[arm];
	  y_errors[arm][bin] = sqrt(scale[arm]*Run15pp_minbias[arm][bin]);
       }
   }


 //PLOT MINBIAS AND FIT
 
   TCanvas *c20 = new TCanvas("c20","Corr BG pT Integrated",200,10,700,500);
   gPad->SetGrid();
     
   TGraphErrors *corrbg_int = new TGraphErrors(bins,x_array,Run15pp_minbias[arm_low],x_errors[arm_low],y_errors[arm_low]);
   
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
       f2 = 10;
       f3 = 10;
       f4 = 10;
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

   int_fit_pt = new TF1("int_fit_pt"," [2]/ pow( ( (exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",1,5);

   int_fit_pt->SetParameter(0, f0);
   int_fit_pt->SetParameter(1, f1);
   int_fit_pt->SetParameter(2, f2);
   int_fit_pt->SetParameter(3, f3);
   int_fit_pt->SetParameter(4, f4);
   
   int_fit_pt->SetLineColor(kRed);
   int_fit_pt->SetLineStyle(5);
   int_fit_pt->SetLineWidth(2);

 
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);  
 
   r = corrbg_int->Fit(int_fit_pt,"S R");
   
//////////////////////////////////////// begin xuan ratio plot ////////////////////////////////////////
  
   double Run15pp_ratio[2][100] ={0};
   double num_ratio[2][100] ={0};
   double denom[2][100] = {0};
   double ratio_errors[2][100] = {0};
   double denom_errors[2][100] = {0};
   double denom_200[2][100] = {0};
   double denom_500[2][100] = {0};
   double denom_1000[2][100] = {0};

   for(int arm = 0; arm < 2; arm++)
     {
       for(int bin = 0; bin < 100; bin++)
	 {
	   num_ratio[arm][bin] = Run15pp_minbias[arm][bin];
	   pt_int_array[arm][bin] = int_fit_pt->Eval(x_array[bin]);
	  
	 }
     }

   cout << "South num: " << Run15pp_minbias[0][26] << endl;
   cout << "North num: " << Run15pp_minbias[1][26] << endl;
   
   for(int arm = 0; arm < 2; arm++)
     {
       for(int pt = 0; pt < 14; pt++)
	 {
	   for(int bin = 0; bin < 100; bin++)
	     {
	       if(pt < 4 && arm == 0)
		 denom_200[arm][bin] +=  fit_array_200mev_S[pt][bin];
	       if(pt < 4 && arm == 1)
		 denom_200[arm][bin] += fit_array_200mev_N[pt][bin];
	       if(pt > 3 && pt < 12 && arm == 0)
		 denom_500[arm][bin] += yuehang_direct_pt_int_500mev_S[pt][bin];
	       if(pt > 3 && pt < 12 && arm == 1)
		 denom_500[arm][bin] += yuehang_direct_pt_int_500mev_N[pt][bin];
	       if(pt > 11 && arm == 0)
		 denom_1000[arm][bin] += yuehang_direct_pt_int_1gev_S[pt][bin];
	       if(pt > 11 && arm == 1)
		 denom_1000[arm][bin] += yuehang_direct_pt_int_1gev_N[pt][bin];
	     }

	 }
     }


   // printed out values for both distributions at bin 26 (right around peak)
     
   double scale_peak[2] = {1855.33/(317.361),1837.03/(358.544)};
  
   // read in Yue hang pT integrated from sum of pt slices  

   for(int arm = 0; arm < 2; arm++)
     { 
       for(int bin = 0; bin < 100; bin++)
	 {
	   denom[arm][bin] = (denom_200[arm][bin] + denom_500[arm][bin] + denom_1000[arm][bin])/scale_peak[arm];
	 }
      }

   cout << "North denom: " << denom[1][26] << endl;
   cout << "South denom: " << denom[0][26] << endl;
   
   for(int arm = 0; arm < 2; arm++)
     {   
       for(int bin = 0; bin < 100; bin++)
	 {
	   if(denom[arm][bin] != 0)
	     {
	       Run15pp_ratio[arm][bin] = num_ratio[arm][bin]/denom[arm][bin];
	       ratio_errors[arm][bin] = 0.02*Run15pp_ratio[arm][bin];
	     }
	 }
     }
   
   //create ratio plot
   // comment out for pt_direct_int
   // TCanvas *c2 = new TCanvas("c2","Corr BG pT Integrated",200,10,700,500);
   //  gPad->SetGrid();
   
   TGraphErrors *corrbg_direct_pt_int = new TGraphErrors(bins,x_array,denom[arm_low],x_errors[arm_low],denom_errors[arm_low]);
   TGraphErrors *corrbg_int_ratio = new TGraphErrors(bins,x_array,Run15pp_ratio[arm_low],x_errors[arm_low],ratio_errors[arm_low]);
   
   
   if(pt_direct_int)// if plotting only pt direct integrated
     {
       if(south_arm)
	 corrbg_direct_pt_int->SetTitle("Corr BG pT Integrated, South");  
       else
	 corrbg_direct_pt_int->SetTitle("Corr BG pT Integrated, North");  
       corrbg_direct_pt_int->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
       corrbg_direct_pt_int->GetXaxis()->SetLabelSize(0.04);
       corrbg_direct_pt_int->GetXaxis()->SetTitleSize(0.04);
       corrbg_direct_pt_int->GetXaxis()->SetTitleOffset(0.9);
       corrbg_direct_pt_int->GetXaxis()->SetLimits(2,5);
       corrbg_direct_pt_int->GetYaxis()->SetLabelSize(0.04);
       corrbg_direct_pt_int->GetYaxis()->SetTitleSize(0.1); 
       corrbg_direct_pt_int->GetYaxis()->SetTitleOffset(0.52);
       //corrbg_direct_pt_int->SetMaximum(0.3);
       // corrbg_direct_pt_int->SetMinimum(pow());
       
         if(south_arm)
	 corrbg_direct_pt_int->SetMarkerColor(kViolet);  
       else
	 corrbg_direct_pt_int->SetMarkerColor(kAzure);  
       corrbg_direct_pt_int->SetMarkerSize(0.75);
       corrbg_direct_pt_int->SetMarkerStyle(20);
       corrbg_direct_pt_int->Draw("P");
     }
 
   if(pt_mass_ratio)
     {
       if(south_arm)
	 corrbg_int_ratio->SetTitle("Corr BG p_{T} Integrated Mass Ratio (refit/raw), South");  
       else
	 corrbg_int_ratio->SetTitle("Corr BG p_{T} Integrated Mass Ratio (refit/raw), North");  
       corrbg_int_ratio->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
       corrbg_int_ratio->GetXaxis()->SetLabelSize(0.04);
       corrbg_int_ratio->GetXaxis()->SetTitleSize(0.04);
       corrbg_int_ratio->GetXaxis()->SetTitleOffset(0.9);
       corrbg_int_ratio->GetXaxis()->SetLimits(2,5);
       corrbg_int_ratio->GetYaxis()->SetLabelSize(0.04);
       corrbg_int_ratio->GetYaxis()->SetTitleSize(0.1); 
       corrbg_int_ratio->GetYaxis()->SetTitleOffset(0.52);
       corrbg_int_ratio->SetMaximum(1.4);
       corrbg_int_ratio->SetMinimum(0.6);

       
       if(south_arm)
	 corrbg_int_ratio->SetMarkerColor(kViolet);  
       else
	 corrbg_int_ratio->SetMarkerColor(kAzure);  
       corrbg_int_ratio->SetMarkerSize(0.75);
       corrbg_int_ratio->SetMarkerStyle(20);
       corrbg_int_ratio->Draw("AP");
     }

   TLegend *leg3 = new TLegend(0.51, 0.6, 0.8, 0.9);   
   
   if(pt_direct_int)
     {
       leg3 = new TLegend(0.51, 0.6, 0.8, 0.9);  
       leg3->SetFillColor(0); 
       leg3->SetTextSize(0.035);
       
       leg3->AddEntry(corrbg_int,"Reconstructed p_{T} integrated (f(x))", "p");
       leg3->AddEntry(int_fit_pt,"Fit to Reconstructed p_{T} integrated", "l");
       leg3->AddEntry(corrbg_direct_pt_int,"Raw p_{T} integrated", "p");
       
       leg3->Draw(); 
     }

 TF1 *fit;
 TF1 *fit_S;

 fit = new TF1("fit","[0]",2,5);
 fit->SetParameter(0,1);
 fit->SetLineStyle(10);
 
 corrbg_int_ratio->Fit(fit,"R"); 
 
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
   sprintf(text4,"p + p");
   
   l3.SetTextAlign(12);
   l3.DrawLatexNDC(0.75, 0.93, text4); //4.4,150
   
   Char_t message3[200];
   
   sprintf(message3,"y = %.3f +/- %.3f",fit->GetParameter(0),fit->GetParError(0));
   

 TPaveText *mytext3 = new TPaveText(0.7,0.8,0.9,0.7,"NDC"); // x0,y0,x1,y1
  mytext3->SetTextSize(0.035);
  mytext3->SetFillColor(0);
  mytext3->SetTextAlign(12);
  mytext3->AddText(message3);
  mytext3->Draw();
 
 


  // bins = 100;

 if(print_screen == true)
    {
      cout << "Run15pp_pt_int_array[0][" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	  cout << pt_int_array[0][i] << ", ";
      	}
      cout << pt_int_array[0][bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
     
    } // print out
  








  /*

   cout << "Chi Square/NDF: " << int_fit_pt->GetChisquare() << "/" << int_fit_pt->GetNDF() << "" << endl;

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
   sprintf(text4,"p + p");
   
   l3.SetTextAlign(12);
   l3.DrawLatexNDC(0.75, 0.93, text4); //4.4,150
   
   
   for(int pt = 0; pt < 7; pt++)
     {
       for(int bin = 0; bin < 100; bin++)
	 {
	   //cout << "area_corrbg[arm][pt]: " << area_corrbg[arm][pt] << endl;
	   //cout << "fit coefficients:  " << fit_array[pt][arm][bin] << endl;
	   
	   //cout << "fit array 500 mev S: " << fit_array_500mev_S[pt][bin] << ", for pt: " << pt << ", and bin: " << bin << endl;
	   //	 cout << "fit array 1 gev S: " << fit_array_1gev_S[pt][bin] << ", for pt: " << pt << ", and bin: " << bin << endl;
	   
	   //  cout << "fit array 500 mev N: " << fit_array_500mev_N[pt][bin] << ", for pt: " << pt << ", and bin: " << bin << endl;
	   //cout << "fit array 1 gev N: " << fit_array_1gev_N[pt][bin] << ", for pt: " << pt << ", and bin: " << bin << endl;
	   
	 }
     }
   
   
 
 cout <<"length of N coefficient array 200 mev should be 400 and is: " << sizeof(fit_array_200mev_N)/sizeof(fit_array_200mev_N[0][0]) << endl;
cout <<"length of N coefficient array 500 mev should be 1000 and is: " << sizeof(fit_array_500mev_N)/sizeof(fit_array_500mev_N[0][0]) << endl;
cout <<"length of N coefficient array 1 gev should be 700 and is: " << sizeof(fit_array_1gev_N)/sizeof(fit_array_1gev_N[0][0]) << endl;
cout <<"length of S coefficient array 200 mev should be 400 and is: " << sizeof(fit_array_200mev_S)/sizeof(fit_array_200mev_S[0][0]) << endl;
cout <<"length of S coefficient array 200 mev should be 1000 and is: " << sizeof(fit_array_500mev_S)/sizeof(fit_array_500mev_S[0][0]) << endl;
cout <<"length of S coefficient array 200 mev should be 700 and is: " << sizeof(fit_array_1gev_S)/sizeof(fit_array_1gev_S[0][0]) << endl;
 */
  
}// end void 


