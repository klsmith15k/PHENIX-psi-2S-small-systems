void Run_fill_ana_tree(
		       const char *outFile = "checking_info.root",
		       const char *filein = "test.txt",
		       bool is_sim = true
){
	gSystem->Load("/direct/phenix+u/klsmith/install/lib/libfill_ana_tree.so");
	gSystem->Load("libfun4all");

  Fun4AllServer *se = Fun4AllServer::instance();

  mFillAnaTree* mtree = new mFillAnaTree("mFillAnaTree",outFile, is_sim);
	mtree->set_dataset("Run15pp200");
	se->registerSubsystem(mtree);

	Fun4AllInputManager *inMan = new Fun4AllDstInputManager("Background","DST","TOP");
	//Fun4AllInputManager *inMan = new Fun4AllDstInputManager("Signal","DST","TOP");
	se->registerInputManager(inMan);
	inMan->AddListFile(filein);

	se->run(0);
	se->End();

	delete se;
}

