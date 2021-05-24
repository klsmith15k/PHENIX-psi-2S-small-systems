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

//#define CC
//#define BB
#define DY
//#define CH_UL

using namespace std;

void yuehang_bg_macro_comp()
{

  //int bg_comp = 0;  // for cc component
  //int bg_comp = 1; // for bb component
  int bg_comp = 2; // for dy
  // int bg_comp  = 3; // for UL corr hadrons
  // int bg_comp = 4; // for LS corr hadrons


  bool all_fits = true;
  bool write = false;
  bool south_arm = true;
  //bool south_arm = false;

  int pt_l = 2;
  int pt_h = 3;

  double weight[2][3] = {1,1.03,1,
			 1,0.985,1};
		
  double scale[2] = {2.415*pow(10,11),4.5*pow(10,11)}; // for LS background 

  int pt_binwidth;
  cout << "What pt binning are you running?  Enter '0' for 100 MeV binwidth, enter '1' for 200 MeV binwidth, enter '2' for 500 MeV binwidth, '3' for 1 GeV binwidth, '4' for 2 GeV binwidth" << endl;
  cin >> pt_binwidth;
  
  int pt_slices[5] = {10,6,8,10,5}; 
    
#ifdef CC
  std::string filename_comp[5][2] = {"yuehang_Run15pAu_S_rebinned_500mev_cc.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_cc.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_cc.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root",

				     "yuehang_Run15pAu_S_rebinned_500mev_cc.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_cc.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_cc.root"};
  
#endif
  
#ifdef BB
  std::string filename_comp[5][2] = {"yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  				     "yuehang_Run15pAu_N_rebinned_500mev_bb.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_bb.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_bb.root",

				     "yuehang_Run15pAu_S_rebinned_500mev_bb.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_bb.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  				   "yuehang_Run15pAu_N_rebinned_500mev_bb.root",
				     
				     "yuehang_Run15pAu_S_rebinned_500mev_bb.root",
				     "yuehang_Run15pAu_N_rebinned_500mev_bb.root"};
#endif
  
  
#ifdef DY
  std::string filename_comp[5][2] = {"yuehang_Run15pp_S_rebinned_500_mev.root",
				     "yuehang_Run15pp_N_rebinned_500_mev.root",
				     
				     "yuehang_Run15pp_S_rebinned_500_mev.root",
				     "yuehang_Run15pp_N_rebinned_500_mev.root",
				     
				     "yuehang_Run15pp_S_rebinned_500_mev.root",
				     "yuehang_Run15pp_N_rebinned_500_mev.root",
				     
				     "yuehang_Run15pp_S_rebinned_500_mev.root",
				     "yuehang_Run15pp_N_rebinned_500_mev.root",
				     
				     "yuehang_Run15pp_S_rebinned_500_mev.root",
				     "yuehang_Run15pp_N_rebinned_500_mev.root"};
#endif
  
   
  double pt_width[5][38] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
			    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			    1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  double pt_center[5][38] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,
			     0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     
			     0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			     0.5,2.0,4.0,6.0,8.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


  TFile *file_comp;
  TH1D *bg[3][2][38];
  
 
/////////////

  // add other filenames once they are prepared...
  // std::string filename_bb_LS[2] = {"yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  //std::string filename_bb_UL[2] = {"yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  // std::string filename_corr_had_UL[2] = {"yuehang_Run15pAu_S_rebinned_500mev_bb.root",
  //std::string filename_corr_had_LS[2] = {"yuehang_Run15pAu_S_rebinned_500mev_bb.root",

  //////////////
 

  std::string obj_filename[3][2][38];  //[number of files][arm][pt slices]
  
  char unique1[800] = {0};
  char unique2[800] ={0};
  char unique3[800] ={0};
   
  int bins = 100;
  
  double pt_low;
  double pt_high;
  
  double bin_low;
  double bin_high;
  
  double delta_pt = 0.001;
  double a_array[38];
  double b_array[38];
  
  for(int arm = 0; arm < 2; arm++)
    {
      for(int k = 0; k < pt_slices[pt_binwidth]; k++)
	{ 
	  pt_low = pt_center[pt_binwidth][k] - pt_width[pt_binwidth][k]/2 + delta_pt; 
	  pt_high = pt_center[pt_binwidth][k] + pt_width[pt_binwidth][k]/2 - delta_pt;
  
	  double a = (pt_low - delta_pt)*1000; 
	  double b = (pt_high + delta_pt)*1000;

	  if(bg_comp == 0)
	    {
	      sprintf(unique1,"cc_pp%.0f_%.0f",a,b);                  // cc unmodified is [0]
	      sprintf(unique2,"cc_powheg_%.0f_%.0f",a,b);        // pp unmodified is [1]
	      sprintf(unique3,"cc_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) is [2]  
	  //   --->   cc modifed = (pp modified/pp unmodfied)*cc unmodified = ( [2] / [1] ) * [0]
	    }
	  if(bg_comp == 1)
	    {
	      sprintf(unique1,"bb_pp%.0f_%.0f",a,b);                  // bb unmodified is [0]
	      sprintf(unique2,"bb_powheg_%.0f_%.0f",a,b);        // pp unmodified is [1]
	      sprintf(unique3,"bb_suppression_%.0f_%.0f",a,b);   // pp modified (RdAu) is [2]  
	  //   --->   bb modifed = (pp modified/pp unmodfied)*bb unmodified = ( [2] / [1] ) * [0]
	    }
	  if(bg_comp == 2)
	     sprintf(unique1,"dy_%.0f_%.0f",a,b);                 

	  a_array[k] = a/1000;
	  b_array[k] = b/1000;

	  obj_filename[0][arm][k] = unique1;
	  obj_filename[1][arm][k] = unique2;
	  obj_filename[2][arm][k] = unique3;
	}
    }
  
  for(int i_histo = 0; i_histo < 3; i_histo++)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	    {
	      cout << obj_filename[i_histo][arm][pt].c_str() << " " << endl;
	    }
	}
    }
  
 
  for(int i_histo = 0; i_histo < 3; i_histo++)
    { 
      for(int arm = 0; arm < 2; arm++)
	{
	  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++) 
	    {
	      file_comp = TFile::Open(filename_comp[pt_binwidth][arm].c_str()); 
	      file_comp->GetObject(obj_filename[i_histo][arm][pt].c_str(),bg[i_histo][arm][pt]);
	    }
	}	  
    }

  TAxis *xaxis;
  TAxis *yaxis;

  double bin_two;
  double bin_five;

  double x_array[100];
  double x_array_mass_ratio[100];
  double x_low[100];
  double x_high[100];

  Double_t binCenterX;
  Double_t binLowX;
  Double_t binWidthX;
  Int_t numBinsX;
 
  Double_t binCenterY;
  Double_t binLowY;
  Double_t binWidthY;
  Int_t numBinsY;

 
  for(int arm = 0; arm < 2; arm++)
    {
      if(bg_comp == 2)
	xaxis = bg[0][arm][0]->GetXaxis(); 
      else
	xaxis = bg[1][arm][0]->GetXaxis();  // powheg location

      binCenterX = xaxis->GetBinCenter(1);
      binLowX = xaxis->GetBinLowEdge(1);
      binWidthX = xaxis->GetBinWidth(1);
      numBinsX = xaxis->GetNbins();
	     
      bin_two = xaxis->FindBin(2.001);
      bin_five = xaxis->FindBin(4.999);
      cout << " ______________________________ " << endl;
      cout << "powheg mass bin 2.0: " << bin_two << ", and bin 5.0: " << bin_five << endl;
      cout << "x axis starts at: " << binLowX <<  endl;
      cout << "x binwidth: " << binWidthX <<  endl;
      cout << "x axis bin 1 center: " << binCenterX << endl;
      cout << "x axis total bins: " << numBinsX << endl;
      cout << " ______________________________ " << endl;

      if(bg_comp != 2)
	{
	  yaxis = bg[2][arm][0]->GetYaxis(); // dAu suppression location
	  binCenterY = yaxis->GetBinCenter(1);
	  binLowY = yaxis->GetBinLowEdge(1);
	  binWidthY = yaxis->GetBinWidth(1);
	  numBinsY = yaxis->GetNbins();
	  
	  cout << " ______________________________ " << endl;
	  cout << "dAu suppression mass bin 2.0: " << bin_two << ", and bin 5.0: " << bin_five << endl;
	  cout << "y axis starts at: " << binLowY <<  endl;
	  cout << "y binwidth: " << binWidthY <<  endl;
	  cout << "y axis bin 1 center: " << binCenterY << endl;
	  cout << "y axis total bins: " << numBinsY << endl;
	  cout << " ______________________________ " << endl;
	}
    }
  
  if(bg_comp != 2)      
    xaxis = bg[1][0][0]->GetXaxis();  
  else
    xaxis = bg[0][0][0]->GetXaxis(); 
 
  for(int i = 0; i < 100; i++)
    {
      x_array[i]  = xaxis->GetBinCenter(i+1);
      // cout << "x array: " << x_array[i] << endl;
    }
  
  int arm_low;
  int arm_high;
  
  if(south_arm == true)
    {
      arm_low = 0;
      arm_high = 1;
    }
  else
    {
      arm_low = 1;
      arm_high = 2;
    }
  
  // scale cc/bb/dy components from pp histograms (these were invariant yields, not the counts)
  for(int arm = 0; arm < 2; arm++)
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
	{
	  bg[0][arm][pt]->Scale(scale[arm]);
	}
    }

  char uniqueh[38][500];
  char uniqued[38][500];
  int i = 0;
  
  for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
    {
      i++;
      sprintf(uniqueh[pt],"h_comp_%d",i); 
      sprintf(uniqued[pt],"h_RdAu_%d",i); 
    }
  
  TH1D *h1 = new TH1D("h1", "cc distribution", 100, 0, 10);
  TH1D *h2 = new TH1D("h2", "RdAu distribution", 100, 0, 10);

  TH1D *h_comp[2][38];
  TH1D *h_RdAu[2][38];
   
  double comp_1[2][38][100];
  double comp_2[2][38][100];
  double comp_3[2][38][100];
  double comp_pAu[2][38][100];
  double R_pAu[2][38][100];

  double comp_err_1[2][38][100];
  double comp_err_2[2][38][100];
  double comp_err_3[2][38][100];
  
  double quad[2][38][100];
  double quad_norm[2][38][100];
  double temp[2][100];
  double x_errors[2][100];
 
  for(int arm = 0; arm < 2; arm++)  
    {
      for(int pt = 0; pt < pt_slices[pt_binwidth]; pt++)
 	{
	  for(int i = 0; i < bins; i++)
	    {
	      if(bg_comp == 2)
		h1->SetBinContent(i+1,( bg[0][arm][pt]->GetBinContent(i+1) ) );
	      else
		{
		  // Multiply dAu Suppression by unmodified distribution to get modified distribution for pAu
		  if(bg[1][arm][pt]->GetBinContent(i+1) != 0)
		    h1->SetBinContent(i+1,( bg[0][arm][pt]->GetBinContent(i+1)*( bg[2][arm][pt]->GetBinContent(i+1)/bg[1][arm][pt]->GetBinContent(i+1) ) )*weight[arm][bg_comp] );
		  else
		    h1->SetBinContent(i+1, 0);
		  
		  if(bg[1][arm][pt]->GetBinContent(i+1) != 0)
		    h2->SetBinContent(i+1, bg[2][arm][pt]->GetBinContent(i+1)/bg[1][arm][pt]->GetBinContent(i+1) ); 
		  else
		    h2->SetBinContent(i+1, 0);
		  
		  // component arrays for print out, not plotting
		  comp_1[arm][pt][i] = bg[0][arm][pt]->GetBinContent(i+1);
		  comp_2[arm][pt][i] = bg[1][arm][pt]->GetBinContent(i+1);
		  comp_3[arm][pt][i] = bg[2][arm][pt]->GetBinContent(i+1);
		  
		  comp_pAu[arm][pt][i] = comp_1[arm][pt][i] *( comp_3[arm][pt][i] / comp_2[arm][pt][i] )*weight[arm][bg_comp];
		  R_pAu[arm][pt][i] = (comp_3[arm][pt][i] / comp_2[arm][pt][i]) * weight[arm][bg_comp];
		  
		  comp_err_1[arm][pt][i] = (bg[0][arm][pt]->GetBinError(i+1));  //  ~10% 
		  comp_err_2[arm][pt][i] = (bg[1][arm][pt]->GetBinError(i+1));  //  ~10% 
		  comp_err_3[arm][pt][i] = (bg[2][arm][pt]->GetBinError(i+1));  //  ~10% 
		  
		  //quad[arm][pt][i] = sqrt( pow(comp_err_1[arm][pt][i],2) +  pow(comp_err_2[arm][pt][i],2) +  pow(comp_err_3[arm][pt][i],2) );
		  quad[arm][pt][i] = ( pow(comp_err_1[arm][pt][i]/comp_1[arm][pt][i],2) +  pow(comp_err_2[arm][pt][i]/comp_2[arm][pt][i],2) +  pow(comp_err_3[arm][pt][i]/comp_3[arm][pt][i],2) )  * pow( (comp_1[arm][pt][i]*comp_3[arm][pt][i])/comp_2[arm][pt][i] ,2) ;
		  
		  if(h1->GetBinContent(i+1) > 0)
		    {
		      h1->SetBinError(i+1,sqrt (quad[arm][pt][i] ) );
		      // h1->SetBinError(i+1,(h1->GetBinContent(i+1))*per_err[arm][pt]); 
		    }
		  if(h2->GetBinContent(i+1) > 0)
		    {
		      h2->SetBinError(i+1,(h2->GetBinContent(i+1)*0.02)); // this is just B/C, not A*B/C   
		    }
		}
		  
	      x_errors[arm][i] = 0.0;
	      
	      h_comp[arm][pt] = (TH1D *) h1->Clone(uniqueh[pt]); 
	      h_comp[arm][pt]->SetBinContent(i+1,h1->GetBinContent(i+1));

	      if(bg_comp != 2)
		{
		  h_RdAu[arm][pt] = (TH1D *) h2->Clone(uniqued[pt]); 
		  h_RdAu[arm][pt]->SetBinContent(i+1,h2->GetBinContent(i+1));
		}
	    	      
	    }
	}
    }
  
  if(bg_comp !=2)
    {
      for(int arm = 0; arm < 2; arm++)
	{
	  cout << " _____________________________ " << endl;
	  cout << " For mass bin 30 and pt slice 0.5 - 1.0 GeV/c in arm : " << arm << endl;
	  if(bg_comp == 0)
	    cout << " cc value in Run15pp: " << comp_1[arm][1][30] << endl;
	  if(bg_comp == 1)
	    cout << " bb value in Run15pp: " << comp_1[arm][1][30] << endl;
	  cout << " powheg (unmodified) value: " <<  comp_2[arm][1][30] <<  endl;
	  cout << " RpAu value: " <<  R_pAu[arm][1][30] <<  endl;
	  if(bg_comp == 0)
	    cout << " cc modified value: " <<  comp_pAu[arm][1][30] << endl;
	  if(bg_comp == 1)
	    cout << " bb modified value: " <<  comp_pAu[arm][1][30] << endl;
	}
    }

/////////////// CUSTOMIZED PLOTS ///////////////////////
  
// Histogram comp
  TCanvas *c6 = new TCanvas("c6","Corr BG comp ",200,10,700,500);
  gPad->SetLeftMargin(0.15); 
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gPad->SetGrid();



  if(all_fits)
    {
      if(bg_comp == 0)
	{
	  h_comp[0][0]->SetTitle("Run15pAu Corr BG ( cc(UL )), South");
	  h_comp[1][0]->SetTitle("Run15pAu Corr BG ( cc(UL )), North");
	}
      if(bg_comp == 1)
	{
	  h_comp[0][0]->SetTitle("Run15pAu Corr BG ( bb(UL )), South");
	  h_comp[1][0]->SetTitle("Run15pAu Corr BG ( bb(UL )), North");
	}
      if(bg_comp == 2)
	{
	  h_comp[0][0]->SetTitle("Run15pAu Drell Yan, South");
	  h_comp[1][0]->SetTitle("Run15pAu Drell Yan, North");
	}
      h_comp[arm_low][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_comp[arm_low][0]->GetXaxis()->SetLabelSize(0.04);
      h_comp[arm_low][0]->GetXaxis()->SetTitleSize(0.04);
      h_comp[arm_low][0]->GetXaxis()->SetTitleOffset(0.9);
      h_comp[arm_low][0]->GetYaxis()->SetLabelSize(0.04);
      h_comp[arm_low][0]->GetYaxis()->SetTitleSize(0.1); 
      h_comp[arm_low][0]->GetYaxis()->SetTitleOffset(0.52);
      h_comp[arm_low][0]->SetMinimum(0);
      h_comp[arm_low][0]->SetMaximum(1000);
      h_comp[arm_low][0]->SetAxisRange(0,5);
    }
  else
    {
      if(bg_comp == 0)
	{
	  h_comp[0][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) ), South");
	  h_comp[1][pt_l]->SetTitle("Run15pAu Corr BG (cc(UL)), North");
	}
      if(bg_comp == 1)
	{
	  h_comp[0][pt_l]->SetTitle("Run15pAu Corr BG ( bb(UL) ), South");
	  h_comp[1][pt_l]->SetTitle("Run15pAu Corr BG (bb(UL)), North");
	}
      if(bg_comp == 2)
	{
	  h_comp[0][0]->SetTitle("Run15pAu Drell Yan, South");
	  h_comp[1][0]->SetTitle("Run15pAu Drell Yan, North");
	}
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
      h_comp[arm_low][pt_l]->GetXaxis()->SetLabelSize(0.04);
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitleSize(0.04);
      h_comp[arm_low][pt_l]->GetXaxis()->SetTitleOffset(0.9);
      h_comp[arm_low][pt_l]->GetYaxis()->SetLabelSize(0.04);
      h_comp[arm_low][pt_l]->GetYaxis()->SetTitleSize(0.1); 
      h_comp[arm_low][pt_l]->GetYaxis()->SetTitleOffset(0.52);
      h_comp[arm_low][pt_l]->SetMinimum(0);
      h_comp[arm_low][pt_l]->SetMaximum(1200);
      h_comp[arm_low][pt_l]->SetAxisRange(0,5);
    }

  for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
    {
      h_comp[arm_low][pt]->SetMarkerColor(pt+1);  
      if(pt+1==10)
	h_comp[arm_low][pt]->SetMarkerColor(30); 
      if(pt+1==2)
	h_comp[arm_low][pt]->SetMarkerColor(28); 
      if(pt+1 == 5)
	h_comp[arm_low][pt]->SetMarkerColor(40); 
      h_comp[arm_low][pt]->SetMarkerSize(0.7);
      h_comp[arm_low][pt]->SetMarkerStyle(20);

      if(all_fits)
	{
	  if(pt == 0)
	    h_comp[arm_low][pt]->Draw();
	  else
	    h_comp[arm_low][pt]->Draw("SAME");
	}
      else
	{
	  if(pt == pt_l)
	    h_comp[arm_low][pt_l]->Draw();
	}
    }

 /////////////////////////////////////////////////////////////////////////////////
 /// Make legend for h_comp histogram

  char unique_l[800];
  TLegend *leg_sum[2];
 
  double for_begin; 
  double for_end; 

  if(all_fits)
    {
      for_begin = 0;
      for_end = pt_slices[pt_binwidth];
    }
  else
    {
      for_begin = pt_l;
      for_end = pt_h;
    }
 
  for(int arm = 0; arm < 2; arm++)
    {
      leg_sum[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
      leg_sum[arm]->SetFillColor(0); 
      leg_sum[arm]->SetTextSize(0.035);
	  	 
      for(int pt = for_begin; pt < for_end; pt++)
	{
	  if(bg_comp == 0)
	    sprintf(unique_l,"cc corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  if(bg_comp == 1)
	    sprintf(unique_l,"bb corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  if(bg_comp == 2)
	    sprintf(unique_l,"dy corr bg %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	  leg_sum[arm]->AddEntry(h_comp[arm][pt],unique_l, "p");
	}
    }

  leg_sum[arm_low]->Draw();



  if(bg_comp !=2)
    {
      /////////////////////////////////////////////////////////////////  Histogram RdAu //////////////////////////////////////////////////////////////////////////////
      TCanvas *c1 = new TCanvas("c1","RdAu Distribution",200,10,700,500);
      gPad->SetLeftMargin(0.15); 
      gPad->SetBottomMargin(0.15);
      gStyle->SetOptStat(0);
      // gPad->SetLogy();
      gPad->SetGrid();
  

      if(all_fits)
	{
	  if(bg_comp == 0)
	    {
	      h_RdAu[0][0]->SetTitle("RdAu Distribution ( cc(UL) ), South");
	      h_RdAu[1][0]->SetTitle("RdAu Distribution ( cc(UL) ), North");
	    }
	  if(bg_comp == 1)
	    {
	      h_RdAu[0][0]->SetTitle("RdAu Distribution ( bb(UL) ), South");
	      h_RdAu[1][0]->SetTitle("RdAu Distribution ( bb(UL) ), North");
	    }
	  if(bg_comp == 2)
	    {
	      h_RdAu[0][0]->SetTitle("RdAu Distribution ( drell yan(UL) ), South");
	      h_RdAu[1][0]->SetTitle("RdAu Distribution ( drell yan(UL) ), North");
	    }

	  h_RdAu[arm_low][0]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
	  h_RdAu[arm_low][0]->GetXaxis()->SetLabelSize(0.04);
	  h_RdAu[arm_low][0]->GetXaxis()->SetTitleSize(0.04);
	  h_RdAu[arm_low][0]->GetXaxis()->SetTitleOffset(0.9);
	  h_RdAu[arm_low][0]->GetYaxis()->SetLabelSize(0.04);
	  h_RdAu[arm_low][0]->GetYaxis()->SetTitleSize(0.1); 
	  h_RdAu[arm_low][0]->GetYaxis()->SetTitleOffset(0.52);
	  h_RdAu[arm_low][0]->SetMinimum(0);
	  h_RdAu[arm_low][0]->SetMaximum(5);
	  h_RdAu[arm_low][0]->SetAxisRange(2,5);
	}
      else
	{
	  if(bg_comp == 0)
	    {
	      h_RdAu[0][pt_l]->SetTitle("Run15pAu Corr BG ( cc(UL) ), South");
	      h_RdAu[1][pt_l]->SetTitle("Run15pAu Corr BG (cc(UL)), North");
	    }
	  if(bg_comp == 1)
	    {
	      h_RdAu[0][pt_l]->SetTitle("Run15pAu Corr BG ( bb(UL) ), South");
	      h_RdAu[1][pt_l]->SetTitle("Run15pAu Corr BG (bb(UL)), North");
	    }
	  if(bg_comp == 2)
	    {
	      h_RdAu[0][pt_l]->SetTitle("Run15pAu Corr BG ( dy(UL) ), South");
	      h_RdAu[1][pt_l]->SetTitle("Run15pAu Corr BG ( dy(UL) ), North");
	    }
	  h_RdAu[arm_low][pt_l]->GetXaxis()->SetTitle("Inv. Mass (GeV/c^{2})");
	  h_RdAu[arm_low][pt_l]->GetXaxis()->SetLabelSize(0.04);
	  h_RdAu[arm_low][pt_l]->GetXaxis()->SetTitleSize(0.04);
	  h_RdAu[arm_low][pt_l]->GetXaxis()->SetTitleOffset(0.9);
	  h_RdAu[arm_low][pt_l]->GetYaxis()->SetLabelSize(0.04);
	  h_RdAu[arm_low][pt_l]->GetYaxis()->SetTitleSize(0.1); 
	  h_RdAu[arm_low][pt_l]->GetYaxis()->SetTitleOffset(0.52);
	  h_RdAu[arm_low][0]->SetMinimum(0);
	  h_RdAu[arm_low][0]->SetMaximum(5);
	  h_RdAu[arm_low][0]->SetAxisRange(2,5);
	}

      for(int pt = 0; pt < pt_slices[pt_binwidth];pt++)
	{
	  h_RdAu[arm_low][pt]->SetMarkerColor(pt+1);  
	  if(pt+1==10)
	    h_RdAu[arm_low][pt]->SetMarkerColor(30); 
	  if(pt+1==2)
	    h_RdAu[arm_low][pt]->SetMarkerColor(28); 
	  if(pt+1 == 5)
	    h_RdAu[arm_low][pt]->SetMarkerColor(40); 
	  h_RdAu[arm_low][pt]->SetMarkerSize(0.7);
	  h_RdAu[arm_low][pt]->SetMarkerStyle(20);

	  if(all_fits)
	    {
	      if(pt == 0)
		h_RdAu[arm_low][pt]->Draw();
	      else
		h_RdAu[arm_low][pt]->Draw("SAME");
	    }
	  else
	    {
	      if(pt == pt_l)
		h_RdAu[arm_low][pt_l]->Draw();
	    }
	}

      /////////////////////////////////////////////////////////////////////////////////
      /// Make legend for RdAu histogram

      char unique_r[800];
      TLegend *leg_h[2];
 
      if(all_fits)
	{
	  for_begin = 0;
	  for_end = pt_slices[pt_binwidth];
	}
      else
	{
	  for_begin = pt_l;
	  for_end = pt_h;
	}
 
      for(int arm = 0; arm < 2; arm++)
	{
	  leg_h[arm] = new TLegend(0.51, 0.5, 0.8, 0.9);  
	  leg_h[arm]->SetFillColor(0); 
	  leg_h[arm]->SetTextSize(0.035);
	  	 
	  for(int pt = for_begin; pt < for_end; pt++)
	    {
	      sprintf(unique_r,"RdAu distribution %.1f - %.1f GeV",a_array[pt],b_array[pt]); 
	      leg_h[arm]->AddEntry(h_RdAu[arm][pt],unique_r, "p");
	    }
	}

      leg_h[arm_low]->Draw();
    }

 
}// void end macro



      
  
