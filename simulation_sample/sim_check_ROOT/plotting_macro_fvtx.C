
#include <TF1.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <TVirtualFitter.h>
#include <TMatrixDSym.h>
#include <TFile.h>
#include <TLine.h>
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
# include <fstream>
# include <iostream>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TNtuple.h>
#include <TPaveText.h>


using namespace std;
  
void plotting_macro_fvtx()
{
 
   int arm = 0;

   // cout << "Enter the arm ('0' for South, '1' for North" << endl;
   // cin >> arm;

   bool north_arm = true;

  // if(arm == 0)
  // bool  north_arm = false;

  int plot = 1;
  //cout << "Enter the plot number " << endl;
  //cin >> plot;

  int color = 1;
  // cout << "Enter marker color (1-7)" << endl;
  // cin >> color;

  char name100[500];
  char name150[500];

  const int narm = 2;
  const int nvars = 6;
  const int nsteps = 100; 

  double min_S = 0;
  double temp1_S = 2;
  double temp1_N = 2;
  double temp2_S = 0;
  double temp2_N = 0;
 
  double min_N = 0;
  double array[100] = {0};

  double var_scan[narm][nvars][nsteps] = {0};

  int counter = 0;
  int k = 0;
  
 for(int arm = 0; arm < 2; arm++)
    {
      if(arm == 0)
      	{
	  //sprintf(name100,"fvtx_south_loop/cat_all_S.txt");
	 // sprintf(name100,"sngtrk_south_loop/cat_all_S.txt");
	  //sprintf(name100,"fvtx_south_sim_211/cat_all_S.txt");

	  // sprintf(name100,"new_MC_S/cat_all_S.txt");
	  sprintf(name100,"newlib_MC/cat_all_S.txt");

	  // sprintf(name100,"with_bg_south/cat_all_S.txt");
	 // sprintf(name100,"fvtxsngtrk_south_loop/cat_all_S.txt");
	 // sprintf(name100,"/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/Jpsis/Run15pp/fit_NOFVTX/krista_fits/CB_parameter_scan_vs_sigma_2nd/test_fvtx_south_loop/cat_all_S.txt");
      	}
      if(arm == 1)
	{
	  // sprintf(name100,"fvtx_north_loop/cat_all_N.txt");
	  //sprintf(name100,"sngtrk_north_loop/cat_all_N.txt");

	  // sprintf(name100,"new_MC_N/cat_all_N.txt");
	  sprintf(name100,"newlib_MC/cat_all_N.txt");

	  // sprintf(name100,"with_bg_north/cat_all_N.txt");
	  //sprintf(name100,"fvtxsngtrk_north_loop/cat_all_N.txt");
	}

      counter = 0;
    
  ifstream Run15pp(name100);
      if(Run15pp)
	{
	  do 
	    {
	      double Njpsi, err_Njpsi, Npsi2s, err_Npsi2s, ratio_int, chi2, chi2ndf, fixed, CB_mean, N_2nd, X_mean, x_mean, alpha, n, sigma, N_CB;
	      Run15pp >> Njpsi >> err_Njpsi >> fixed >> Npsi2s >> err_Npsi2s; // >>
	       
	      var_scan[arm][0][counter] = Njpsi;
	      var_scan[arm][1][counter] = err_Njpsi;
	      var_scan[arm][2][counter] = counter;
	      var_scan[arm][3][counter] = Npsi2s;
	      var_scan[arm][4][counter] = err_Npsi2s;
	      var_scan[arm][5][counter] = 0;
	    
	      cout <<   var_scan[arm][2][counter]  << endl;

	      Njpsi = 0;  // [0]
	      err_Njpsi= 0; // [1]
	      Npsi2s = 0; // [2]
	      err_Npsi2s = 0; // [3]
	      ratio_int = 0;  // [4]
	      chi2 = 0; // [5] 
	      chi2ndf = 0;  // [6]
	      fixed = 0;  // [7]
	      CB_mean = 0; // [10]
	      alpha = 0;
	      n = 0;
	      sigma = 0;
	      N_CB = 0;
	      N_2nd = 0;

	      counter++;
	      
	    }while(counter < 100);  
	} // if
      else
	cout << "file does not exist for arm " << arm << endl;   
    } // for arm
     
     

 double chi2ndf_N = 0;
 double chi2ndf_S = 0;
 double N_ratio = 0;
 double fixed_N = 0;
 double fixed_S = 0;
 double sigma_S = 0;
 double sigma_N = 0;
 

  for(int i = 0; i < 100; i++)
    {
      for(int arm = 0; arm < 2; arm++)
  	{
	  
  	  if(arm == 0)
  	    {
  	      chi2ndf_S = var_scan[arm][1][i];
	      // N_ratio = var_scan[arm][9][i];
	      fixed_S = var_scan[arm][2][i];
	      sigma_S = var_scan[arm][3][i];
	    
	      // if(chi2ndf_S < 2.84159) // fvtx psi2s
	      // if(chi2ndf_S < 2.80988) // fvtx jpsi
		   //	 cout << "South arm minimum: " << chi2ndf_S << ", at " << fixed_S << endl;

  	    }
  	  if(arm == 1)
  	    {
  	      chi2ndf_N = var_scan[arm][1][i];
	      // N_ratio = var_scan[arm][9][i];
  	      fixed_N = var_scan[arm][2][i];
	      sigma_N = var_scan[arm][3][i];

	      // if(chi2ndf_N < 2.20480)  // fvtx psi2s
	      // if(chi2ndf_N < 1.84555)  // fvtx jpsi
	      //	cout << "North arm minimum: " << chi2ndf_N << ", at " << fixed_N << endl;
  	    	      
	    }

	  if( chi2ndf_N < 1.52 && chi2ndf_N > 1.51)
	    cout  << fixed << endl;
	  if( chi2ndf_S < 1.52 && chi2ndf_S > 1.51)
	    cout  << fixed << endl;
  	}
    }


  TGraphErrors *gr1;
  TGraphErrors *gr2;
 


  
 
  if(north_arm == false)
    {
      if(plot == 1)
	{
	  gr1 = new TGraphErrors(100,var_scan[0][2],var_scan[0][1],var_scan[0][0],var_scan[0][1]);
	   gr2 = new TGraphErrors(100,var_scan[1][2],var_scan[1][1],var_scan[0][0],var_scan[0][1]);
	}
      // if(plot == 2)
      // 	gr1 = new TGraph(100,var_scan[0][7],var_scan[0][0]);   
      // if(plot == 3)
      // 	gr1 = new TGraph(100,var_scan[0][7],var_scan[0][4]);   
      // if(plot == 4)
      // 	gr1 = new TGraph(100,var_scan[0][8],var_scan[0][0]);  
      // if(plot == 5) 
      // 	gr1 = new TGraph(100,var_scan[0][8],var_scan[0][0]); 
      // if(plot == 6) 
      // 	gr1 = new TGraph(100,var_scan[0][7],var_scan[0][4]);  
      // if(plot == 7) 
      // 	gr1 = new TGraph(100,var_scan[0][4],var_scan[0][1]);   // ratio N vs. fixed
    }
  else
    {
      if(plot == 1)
	{
	  // gr1 = new TGraphErrors(100,var_scan[1][2],var_scan[1][3],var_scan[1][5],var_scan[1][4]);   // psi2s
	  //gr2 = new TGraphErrors(100,var_scan[0][2],var_scan[0][3],var_scan[0][5],var_scan[0][4]);

	 gr1 = new TGraphErrors(100,var_scan[1][2],var_scan[1][0],var_scan[1][5],var_scan[1][1]);  // jpsi
	 gr2 = new TGraphErrors(100,var_scan[0][2],var_scan[0][0],var_scan[0][5],var_scan[0][1]);
	}
      // if(plot == 2)
      // 	gr1 = new TGraph(100,var_scan[1][7],var_scan[1][9]);   
      // if(plot == 3)
      // 	gr1 = new TGraph(100,var_scan[1][7],var_scan[1][6]);   
      // if(plot == 4)
      // 	gr1 = new TGraph(100,var_scan[1][8],var_scan[1][9]);  
      // if(plot == 5) 
      // 	gr1 = new TGraph(100,var_scan[1][8],var_scan[1][6]); 
      // if(plot == 6) 
      // 	gr1 = new TGraph(100,var_scan[1][7],var_scan[1][4]);  
      // if(plot == 7) 
      // 	gr1 = new TGraph(100,var_scan[1][4],var_scan[1][6]);   // chi2 vs n
     
    } 
    


  TCanvas *c1 = new TCanvas("c1","c1",0,0,690,545);
  
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->cd();
  c1->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.04);
  gPad->SetBottomMargin(0.13);
 
  if(color == 1)
    gr1->SetMarkerColor(kViolet+3);
  if(color == 2)
    gr1->SetMarkerColor(kBlue+3);
  if(color == 3)
    gr1->SetMarkerColor(kMagenta+3);
  if(color == 4)
    gr1->SetMarkerColor(kCyan+3);
  if(color == 5)
    gr1->SetMarkerColor(kSpring+3);
  if(color == 6)
    gr1->SetMarkerColor(kRed+3);

  gr1->SetTitle("");
  gr1->SetLineWidth(1);
  gr1->SetLineColor(kViolet+3);
  gr1->SetMarkerStyle(6);
  
  gr1->GetXaxis()->SetTitleOffset(0.99);
  gr1->GetXaxis()->SetTitleFont(102);
  gr1->GetXaxis()->SetLabelFont(102);

  if(plot == 1 || plot == 2 || plot == 4 || plot == 5 || plot == 6)
    {
      // gr1->GetXaxis()->SetTitle("#sigma_{2nd}");
      gr1->GetXaxis()->SetTitle("Throws");
      //gr1->GetXaxis()->SetTitle("#sigma_{CB}");
      // gr1->GetXaxis()->SetTitle("#bar{x}_{CB}");
      //gr1->GetYaxis()->SetTitle("N_{CB}/N_{Jpsi}");
      // gr1->GetYaxis()->SetTitle("#chi^{2}/NDF, J/#psi sim");
      gr1->GetYaxis()->SetTitle("Counts");
      // gr1->GetYaxis()->SetTitle("Int Ratio");//, J/#psi sim");
      gr1->GetXaxis()->SetRangeUser(0,100);
       gr1->GetXaxis()->SetNdivisions(8);
    }
  if(plot == 3)
    {
      gr1->GetYaxis()->SetTitle("CB n");
      gr1->GetXaxis()->SetTitle("#sigma_{2nd}");
      gr1->GetXaxis()->SetLabelSize(0.035);
      gr1->GetXaxis()->SetTitleSize(0.065);
    if(north_arm)
      {
	gr1->GetXaxis()->SetNdivisions(7);
	//gr1->GetXaxis()->SetRangeUser(0.0,0.04);
	gr1->GetYaxis()->SetRangeUser(5,20);
      }
    else
      {
	gr1->GetXaxis()->SetRangeUser(0.0,0.04);
	gr1->GetXaxis()->SetNdivisions(7);
	//	gr1->GetYaxis()->SetRangeUser(0.0,40);
      }
    //gr1->GetYaxis()->SetRangeUser(1.4,2);
    }
 
  gr1->GetXaxis()->SetTitleSize(0.055);
  gr1->GetXaxis()->SetTitleOffset(0.85);
    gr1->GetXaxis()->SetLabelSize(0.045);
 
  gr1->GetYaxis()->SetTitleOffset(1.25);
  gr1->GetYaxis()->SetTitleFont(102);
  gr1->GetYaxis()->SetLabelFont(102);
  gr1->GetYaxis()->SetTitleSize(0.055);
  gr1->GetYaxis()->SetLabelSize(0.045);
  if(plot == 1)
    {
      //gr1->GetYaxis()->SetTitle("#chi^{2}/NDF");
      //gr1->GetYaxis()->SetTitle("N_{2nd}");///N_{CB}");
      // gr1->GetYaxis()->SetRangeUser(1,6);
      // gr1->GetYaxis()->SetRangeUser(1,20);

       gr1->GetYaxis()->SetRangeUser(4700,5200);
      //gr1->GetYaxis()->SetRangeUser(100,250);

      //gr1->GetXaxis()->SetRangeUser(400,1200);
      //gr1->GetYaxis()->SetRangeUser(0,4);
    }
  if(plot == 2)
    {
      gr1->GetYaxis()->SetTitle("CB alpha");
      gr1->GetYaxis()->SetRangeUser(0.02,1.3);
      gr1->GetXaxis()->SetLabelSize(0.035);
      gr1->GetYaxis()->SetLabelSize(0.035);
      gr1->GetXaxis()->SetTitleSize(0.065);
     
    }
  if(plot == 4)
    {
      gr1->GetYaxis()->SetTitle("CB alpha");
      gr1->GetXaxis()->SetTitle("#bar{X} 2nd");
      gr1->GetYaxis()->SetRangeUser(0.02,1.3);
      gr1->GetXaxis()->SetRangeUser(3,3.3);
      gr1->GetXaxis()->SetNdivisions(4);

    }
  if(plot == 5)
    {
      gr1->GetXaxis()->SetTitle("#bar{X} 2nd");
      gr1->GetYaxis()->SetTitle("CB n");
      gr1->GetYaxis()->SetRangeUser(5,20);
      gr1->GetXaxis()->SetRangeUser(3,3.3);
      gr1->GetXaxis()->SetNdivisions(4);

    }
  if(plot == 6)
    {
      gr1->GetYaxis()->SetTitle("ratio");
      gr1->GetXaxis()->SetTitle("#sigma_{2nd}");
      gr1->GetXaxis()->SetLabelSize(0.035);
      gr1->GetXaxis()->SetTitleSize(0.065);
      gr1->GetYaxis()->SetLabelSize(0.035);
      gr1->GetYaxis()->SetTitleSize(0.065);
      //gr1->GetYaxis()->SetRangeUser(3800,4200);
      
    }
  if(plot == 7)
   {
      gr1->GetYaxis()->SetTitle("N_CB/N_2nd");
      gr1->GetXaxis()->SetTitle("#sigma_{2nd}");
      gr1->GetXaxis()->SetLabelSize(0.035);
      gr1->GetXaxis()->SetTitleSize(0.065);
      gr1->GetYaxis()->SetLabelSize(0.035);
      gr1->GetYaxis()->SetTitleSize(0.065);
      gr1->GetYaxis()->SetRangeUser(0.05,0.2);
      
    }

  //gr1->GetYaxis()->SetNdivisions(7);
  gr1->Draw("AP");
  if(plot == 1)
    {
      gr1->GetXaxis()->SetTitleSize(0.065);
      gr1->SetMarkerStyle(6);
      gr2->SetMarkerStyle(6);
      gr2->SetLineWidth(1);
      gr1->SetMarkerColor(kBlue+2);
      gr2->SetMarkerColor(kRed+2);
      gr2->SetLineColor(kRed+2);
      gr1->SetLineColor(kBlue+2);
      //gr1->GetXaxis()->SetNdivisions(8);
      gr2->Draw("P");

    }

   // TLine *line = new TLine(0.257,1,0.257,6);
  // TLine *line = new TLine(0.220,1,0.220,6);
  TLine *line = new TLine(0,5000,100,5000);
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw("L");

  TLine *line2 = new TLine(0,150,100,150);
  line2->SetLineWidth(2);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(2);
   line2->Draw("L");

  //TLine *line3 = new TLine(0.23,1,0.23,6);
  //TLine *line3 = new TLine(0.284,1,0.284,6);
  TLine *line3 = new TLine(3.1313,1,3.1313,20);
  line3->SetLineWidth(2);
  line3->SetLineColor(kBlack);
  line3->SetLineStyle(2);
  line3->Draw("L");

 TF1 *y1 = new TF1("y1","[0]",0,100);
  y1->SetParameter(0,0);
   y1->SetLineColor(kBlue+2);
  y1->SetLineStyle(1);
  y1->SetLineWidth(2);

  gr1->Fit(y1, "R");

TF1 *y2 = new TF1("y2","[0]",0,100);
  y2->SetParameter(0,0);
   y2->SetLineColor(kRed+2);
  y2->SetLineStyle(1);
  y2->SetLineWidth(2);

  gr2->Fit(y2, "R");

  TLatex l2;
  l2.SetTextSize(0.05);
  l2.SetTextAlign(12);
  l2.SetTextColor(kBlack);
  
  char text2[100];
  char text3[100];
  char text4[100];
  char text5[100];
  char text6[100];

  if(plot !=1)
    {
      if(north_arm)
  	sprintf(text2," 1.2 < y < 2.2");
      else
  	sprintf(text2,"-2.2 < y < -1.2");
    }
  cout << "here" << endl;

  sprintf(text3, " Sng+Dbl Tracks");
  sprintf(text4, "#sigma_{2nd} = 0.219");
   sprintf(text6, "Run15pp");
  sprintf(text5, "MC pairs "); 
  // sprintf(text5, "#bar{x}_{CB}^{min} = 3.1265");
  l2.SetTextFont(102);  
  if(plot !=6)
    {
      l2.DrawLatexNDC(0.565, 0.735, text3); // x0, y0
      // l2.DrawLatexNDC(0.55, 0.825, text2); // x0, y0
      l2.DrawLatexNDC(0.18, 0.885, text4); // x0, y0
      //  l2.DrawLatexNDC(0.18, 0.825, text5); // x0, y0
        l2.DrawLatexNDC(0.18, 0.825, text6); // x0, y0
    }
  else
    {
      l2.DrawLatexNDC(0.28, 0.745, text3); // x0, y0
      // l2.DrawLatexNDC(0.28, 0.825, text2); // x0, y0
    }
  TLegend *leg2 = new TLegend(0.50, 0.77, 0.75, 0.92);  //(start x, start y, end x, end y)
  leg2->SetFillColor(19); 
  leg2->SetFillStyle(3003); 
  leg2->SetLineWidth(1);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.05); 
  leg2->SetTextFont(102);

  leg2->AddEntry(gr1," 1.2 < y < 2.2", "l");
  //leg2->AddEntry(gr1," South, J/#psi Sim", "l");

  if(plot == 1)
    {
      //leg2->AddEntry(gr2," South, dimu02", "l");
      leg2->AddEntry(gr2," -2.2 < y < -1.2", "l");
      leg2->Draw();
    }
    
   
  Char_t message2[200];
  Char_t message3[200];
 
  sprintf(message2,"North J/#psi fit = %.3f  +/- %.3f",y1->GetParameter(0),y1->GetParError(0));
  sprintf(message3,"South J/#psi fit = %.3f +/- %.3f",y2->GetParameter(0), y2->GetParError(0));
  // sprintf(message2,"North #psi(2S) fit = %.3f  +/- %.3f",y1->GetParameter(0),y1->GetParError(0));
  // sprintf(message3,"South #psi(2S) fit = %.3f +/- %.3f",y2->GetParameter(0), y2->GetParError(0));
  
  TPaveText *mytext2 = new TPaveText(0.18,0.16,0.64,0.27,"NDC"); // x0,y0,x1,y1
  mytext2->SetTextSize(0.035);
  mytext2->SetFillColor(0);
  mytext2->SetTextAlign(12);
  mytext2->AddText(message2);
  mytext2->AddText(message3);
  mytext2->Draw();
 
 
 

} // void end macro




