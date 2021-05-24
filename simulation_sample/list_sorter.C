
#include <TF1.h>
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
#include <TPaveText.h>
#include<bits/stdc++.h>

using namespace std;


void list_sorter()
{
  double size = 0;

  char const *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup_Run14HeAu/sng_tracks/sanghoon_simulations/psi2s_dimu_pythia_embedA_301.txt";
  //  char const *infile = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/sanghoon_simulations/psi2s_dimu_02_embedC_301.txt";

 int total = 0;

  ifstream runlist(infile);
  if(runlist)
    {
      do 
	{

	  char filepath [1000];
	  runlist.getline(filepath,1000);
	 
	  TFile *f = new TFile(filepath, "read");

	  if(f)
	    {
	      size = f->GetNbytesInfo();
	      
	      //  cout <<  filepath << endl;
	      
	      if(size > 10000) // unpopulated HEPMC containers have ~7,000 bytes.  Files with HEPMC have ~15,000 bytes
		{	      
		  if(!f)
		    cout << filepath << endl;
		  else
		    {
		      std::ofstream  outfile;
		      outfile.open ("sanghoon_simulations/psi2s_dimu_pythia_embedA_301_all.txt",std::ios_base::app);
		  
		      outfile << filepath << std::endl; 

		     total += size;
		    }
		}
	      
	      // if(f->IsZombie())
		//	cout << "yes: "  << filepath << endl;
	    }
	 	 
	  f->Close();

	}while(runlist.good());  
    } 
  else
    cout << "file does not exist" << endl;   

  runlist.close();

  cout << "total bytes: " << total << endl;

} 
  


