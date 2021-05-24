void Run_fill_ana_tree(
		       const char *outFile = "sanghoon_simulations/sngdbl_psi2s_dimu_02_MC_embedC_301.root",
		       
		       // const char *filein = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/my_test_simulations/txt_files/jpsi_9398_S_noembedC_all.txt",
		       //const char *filein = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/sanghoon_simulations/jpsi_dimu_02_embedC_201_all.txt",
		       const char *filein = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/sng_dbl_tracks/sanghoon_simulations/psi2s_dimu_02_embedC_301_all.txt",
		       // const char *filein = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/sanghoon_simulations/psi2s_dimu_02_embedC_301_all.txt",
		       //const char *filein = "/gpfs/mnt/gpfs02/phenix/hhj/hhj1/klsmith/simulation/sanghoon_setup/fvtx_tracks/my_test_simulations/psi2S_embedc_969798_all.txt",
		       // bool is_sim = false
		      
){
  //gSystem->Load("/direct/phenix+u/klsmith/install/lib/libfill_ana_tree.so");
  // source /gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/phenix_new/setup.csh;
  gSystem->Load("libfill_ana_tree.so");
	gSystem->Load("libfun4all");

  // gROOT->LoadMacro("g4libs.C");
  // g4libs();

  // gSystem->Load("libfun4all.so");	// framework + reco modules
  // gSystem->Load("libfvtx_subsysreco.so");
  // gSystem->Load("libsimreco.so");
  // gSystem->Load("libPHHepMCNode.so");
  // gSystem->Load("libPHG3toG4.so");

	//	gSystem->Load("libpicodst_obj.so");

  Fun4AllServer *se = Fun4AllServer::instance();

  //int run = 431846,

    //recoConsts *rc = recoConsts::instance();
    //rc->set_IntFlag("RUNNUMBER", run);

  mFillAnaTree* mtree = new mFillAnaTree("mFillAnaTree",outFile);
 	mtree->set_dataset("Run15pp200");
	se->registerSubsystem(mtree);

	Fun4AllInputManager *inMan = new Fun4AllDstInputManager("Background","DST","TOP");  // embed?
	//Fun4AllInputManager *inMan = new Fun4AllDstInputManager("Signal","DST","TOP");  // no embed?
	se->registerInputManager(inMan);
	inMan->AddListFile(filein);

	se->run(0);
	se->End();

	delete se;
}

