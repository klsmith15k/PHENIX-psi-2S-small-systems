// sanghoon's macro to make dimuon container without using analysis train

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

void make_picodstobj_dimu(const int nevt,
		     const char* infile,
		     const char* outfile
){
  gSystem->Load("libfvtx_subsysreco.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libMWGOO");
  gSystem->Load("libmutrg");
  gSystem->Load("librpc_subsysreco");
  gSystem->Load("librpc_muotrackreco");
  gSystem->Load("librecal");
  //gSystem->Load("libpicodst_object.so");

  gSystem->Load("/direct/phenix+u/klsmith/install/lib/libpicodst_object.so");

  // gSystem->Load("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/09.re_association/picoDST_object/install/lib/libpicodst_object.so");

/////////////////////////////////////////////////////////////////
  //  Server...
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  FvtxReadbackDST *fvtxoo = new FvtxReadbackDST();
  fvtxoo->Verbosity(0);
  se->registerSubsystem(fvtxoo);

	se->registerSubsystem(new MuonReadbackDST());

	//DB access control

  // counter
  MuonCounter* counter = new MuonCounter();
  counter->set_event_dump( 10000 );
  se->registerSubsystem( counter );

  // new
 // MasterRecalibrator *recal = new MasterRecalibrator();
 //  recal->Unlock();
 //  recal->UseOnly("Run15pAu200GeVCentralityReco");
 //  se->registerSubsystem(recal);

	/*
	mFvtxPostProductionAlignment* align = new mFvtxPostProductionAlignment();
	align->set_align_what(mFvtxPostProductionAlignment::DCA);
	//align->set_align_what(mFvtxPostProductionAlignment::NONE);
	//align->set_sim(true);
	se->registerSubsystem(align);
	*/

  mFillDiMuonContainer* mdi = new mFillDiMuonContainer(true);
	//mdi->set_bbcz_cut(15); //use fvtxz for mixed event
	mdi->set_bbcz_cut(30); //use bbcz for mixed event
	mdi->set_mass_cut(1.0); // 
	//mdi->set_last_gap_cut(3);        
	//	mdi->set_nidhits_cut(14);   
	mdi->set_is_sim(false);
	//mdi->set_is_sim(true);

	mdi->set_is_pp(false);
	mdi->set_make_dimuon_vertex_origin_widerange(false);
	se->registerSubsystem(mdi);

  /////////////////////////////////////////////////////////////////
  //  Input Managers...
  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("input","DST");
  se->registerInputManager(in);

  Fun4AllOutputManager *outpico = new Fun4AllDstOutputManager("Outee",outfile);
  outpico->AddNode("SingleMuonContainer"); 
  // outpico->AddNode("MutrRefitSingleMuonContainer");   // commented out  Feb 14th 2020 to use T->Project and not TTree Reader
  outpico->AddNode("DiMuonContainer");
  outpico->AddNode("TrigLvl1");
	outpico->AddNode("EventHeader");
	outpico->AddNode("RunHeader");
	//	outpico->AddNode("MCHepMCParticleContainer");
  outpico->AddNode("Sync");
  outpico->AddNode("VtxOut");
  outpico->AddNode("PHGlobal");
  outpico->AddEventSelector("mFillDiMuonContainer");
  se->registerOutputManager(outpico);

	if ( strstr(infile,"root") ){
		in->AddFile(infile);
	}else{
		in->AddListFile(infile);
	}

	gSystem->ListLibraries();

  se->run(nevt);

  se->End();

}
