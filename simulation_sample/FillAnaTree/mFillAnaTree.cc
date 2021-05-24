#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <MWGConsts.h>
#include <Tools.h>
#include <Fun4AllServer.h>
#include <SyncObject.h>
#include <RunHeader.h>
#include <VtxOut.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <MCHepMCParticleContainer.h>
#include <MCHepMCParticle.h>
#include <SingleMuonContainer.h>
#include <MCSingleMuonFvtxContainer.h>
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <PHGlobal.h>
#include <PHTFileServer.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <PHPoint.h>
#include <PHGeometry.h>
#include <TFvtxCompactTrkMap.h>
#include <TFvtxCompactTrk.h>
#include <SvxSegmentList.h>
#include <SvxSegment.h>
#include <TrigLvl1.h>
#include <fstream>
#include <iostream>
#include <TrigRunLvl1.h>
#include <PHPythiaContainer.h>
#include <TMCParticle.h>
#include <MutrRefitSingleMuonContainer.h>
#include <TriggerEmulatorContainer.h>
#include "mFillAnaTree.h"

using namespace std;

//___________________________________________________________________
mFillAnaTree::mFillAnaTree(const char *name, const char *rootfile):
  SubsysReco(name)
{

	_filename = rootfile;
	_dataset = "Run15pp200";
	//	_is_sim = is_sim;
	//	_is_sim = is_sim;
	_is_upsilon = false;
	_is_phi = false;
	_target_pid = 100443;
	_nevent = 0;
	_nevent_skip = 0;

	return;
}

//___________________________________________________________________
mFillAnaTree::~mFillAnaTree()
{

	return;
}

//___________________________________________________________________
int mFillAnaTree::Init(PHCompositeNode *top_node)
{  

	MUTOO::PRINT(std::cout,"mFillAnaTree::Init");

	first = true;

	file_out = new TFile(_filename.data(),"RECREATE");
	std::cout << "opening TFile " << _filename.data() << std::endl;

	//
	hpT_eta_zvtx = new TH3F("hpT_eta_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx = new TH3F("hpT_y_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_z = new TH3F("hpT_y_zvtx_w_z","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_dAu = new TH3F("hpT_y_zvtx_w_dAu","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_dAu_y = new TH3F("hpT_y_zvtx_w_dAu_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_dAu_pT = new TH3F("hpT_y_zvtx_w_dAu_pT","",48,0,12,88,-2.2,2.2,60,-30,30);

	hpT_y_zvtx_w_pp = new TH3F("hpT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_pp_y = new TH3F("hpT_y_zvtx_w_pp_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_w_pp_pT = new TH3F("hpT_y_zvtx_w_pp_pT","",48,0,12,88,-2.2,2.2,60,-30,30);

	/*
	hpT_eta_zvtx_S = new TH3F("hpT_eta_zvtx_S","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S = new TH3F("hpT_y_zvtx_S","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_z = new TH3F("hpT_y_zvtx_S_w_z","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_dAu = new TH3F("hpT_y_zvtx_S_w_dAu","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_dAu_y = new TH3F("hpT_y_zvtx_S_w_dAu_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_dAu_pT = new TH3F("hpT_y_zvtx_S_w_dAu_pT","",48,0,12,88,-2.2,2.2,60,-30,30);

	hpT_y_zvtx_S_w_pp = new TH3F("hpT_y_zvtx_S_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_pp_y = new TH3F("hpT_y_zvtx_S_w_pp_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_S_w_pp_pT = new TH3F("hpT_y_zvtx_S_w_pp_pT","",48,0,12,88,-2.2,2.2,60,-30,30);

	hpT_eta_zvtx_N = new TH3F("hpT_eta_zvtx_N","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N = new TH3F("hpT_y_zvtx_N","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_z = new TH3F("hpT_y_zvtx_N_w_z","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_dAu = new TH3F("hpT_y_zvtx_N_w_dAu","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_dAu_y = new TH3F("hpT_y_zvtx_N_w_dAu_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_dAu_pT = new TH3F("hpT_y_zvtx_N_w_dAu_pT","",48,0,12,88,-2.2,2.2,60,-30,30);

	hpT_y_zvtx_N_w_pp = new TH3F("hpT_y_zvtx_N_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_pp_y = new TH3F("hpT_y_zvtx_N_w_pp_y","",48,0,12,88,-2.2,2.2,60,-30,30);
	hpT_y_zvtx_N_w_pp_pT = new TH3F("hpT_y_zvtx_N_w_pp_pT","",48,0,12,88,-2.2,2.2,60,-30,30);
	*/

	hp_pT_y_zvtx = new TH3F("hp_pT_y_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);  // this is one we want
	hp_pT_eta_zvtx = new TH3F("hp_pT_eta_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);
	//	hd_pT_y_zvtx = new TH3F("hd_pT_y_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);  // for MB

	hd_pT_y_zvtx_N = new TH3F("hd_pT_y_zvtx_N","",48,0,12,88,-2.2,2.2,60,-30,30);  // 
	hd_pT_y_zvtx_S = new TH3F("hd_pT_y_zvtx_S","",48,0,12,88,-2.2,2.2,60,-30,30);  // 

	hmu_pT_y_zvtx = new TH3F("hmu_pT_y_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);
	hmu_pT_y_zvtx_w_dAu = new TH3F("hmu_pT_y_zvtx_w_dAu","",48,0,12,88,-2.2,2.2,60,-30,30);
	//hmu_pT_y_zvtx_w_pp = new TH3F("hmu_pT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
        hmu_pT_y_zvtx_w_pp = new TH3F("hmu_pT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
	hmu_y0_y1 = new TH2F("hmu_y0_y1","",88,-2.2,2.2,88,-2.2,2.2);

        hdimu_pT_y_zvtx = new TH3F("hdimu_pT_y_zvtx","",48,0,12,88,-2.2,2.2,60,-30,30);
	hdimu_pT_y_zvtx_cuts = new TH3F("hdimu_pT_y_zvtx_cuts","",48,0,12,88,-2.2,2.2,60,-30,30);
	hdimu_ul_pT_y_zvtx_w_pp = new TH3F("hdimu_ul_pT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);

	hdimu_pp_pT_y_zvtx_w_pp = new TH3F("hdimu_pp_pT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);
	hdimu_mm_pT_y_zvtx_w_pp = new TH3F("hdimu_mm_pT_y_zvtx_w_pp","",48,0,12,88,-2.2,2.2,60,-30,30);

	//////////////////////////

	hmass_N_fvtx = new TH3F("hmass_N_fvtx","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp = new TH3F("hmass_N_fvtx_pp","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm = new TH3F("hmass_N_fvtx_mm","",300,0,15,48,0,12,60,-30,30);

	// 0-10
	hmass_N_fvtx_010 = new TH3F("hmass_N_fvtx_010","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_010 = new TH3F("hmass_N_fvtx_pp_010","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_010 = new TH3F("hmass_N_fvtx_mm_010","",300,0,15,48,0,12,60,-30,30);

	// 10-30
	hmass_N_fvtx_1030 = new TH3F("hmass_N_fvtx_1030","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_1030 = new TH3F("hmass_N_fvtx_pp_1030","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_1030 = new TH3F("hmass_N_fvtx_mm_1030","",300,0,15,48,0,12,60,-30,30);

	//30-50
	hmass_N_fvtx_3050 = new TH3F("hmass_N_fvtx_3050","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_3050 = new TH3F("hmass_N_fvtx_pp_3050","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_3050 = new TH3F("hmass_N_fvtx_mm_3050","",300,0,15,48,0,12,60,-30,30);

	// 40-60
	hmass_N_fvtx_4060 = new TH3F("hmass_N_fvtx_4060","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_4060 = new TH3F("hmass_N_fvtx_pp_4060","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_4060 = new TH3F("hmass_N_fvtx_mm_4060","",300,0,15,48,0,12,60,-30,30);

	// 50-84
	hmass_N_fvtx_5084 = new TH3F("hmass_N_fvtx_5084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_5084 = new TH3F("hmass_N_fvtx_pp_5084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_5084 = new TH3F("hmass_N_fvtx_mm_5084","",300,0,15,48,0,12,60,-30,30);

	// 0-20
	hmass_N_fvtx_020 = new TH3F("hmass_N_fvtx_020","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_020 = new TH3F("hmass_N_fvtx_pp_020","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_020 = new TH3F("hmass_N_fvtx_mm_020","",300,0,15,48,0,12,60,-30,30);

	// 0-5
	hmass_N_fvtx_2040 = new TH3F("hmass_N_fvtx_2040","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_2040 = new TH3F("hmass_N_fvtx_pp_2040","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_2040 = new TH3F("hmass_N_fvtx_mm_2040","",300,0,15,48,0,12,60,-30,30);

	// 5-10
	hmass_N_fvtx_4084 = new TH3F("hmass_N_fvtx_4084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_4084 = new TH3F("hmass_N_fvtx_pp_4084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_4084 = new TH3F("hmass_N_fvtx_mm_4084","",300,0,15,48,0,12,60,-30,30);

	// 40-84
	hmass_N_fvtx_6084 = new TH3F("hmass_N_fvtx_6084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_pp_6084 = new TH3F("hmass_N_fvtx_pp_6084","",300,0,15,48,0,12,60,-30,30);
	hmass_N_fvtx_mm_6084 = new TH3F("hmass_N_fvtx_mm_6084","",300,0,15,48,0,12,60,-30,30);

	//////////////////

	hmass_S_fvtx = new TH3F("hmass_S_fvtx","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp = new TH3F("hmass_S_fvtx_pp","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm = new TH3F("hmass_S_fvtx_mm","",300,0,15,48,0,12,60,-30,30);

	// 0-10
	hmass_S_fvtx_010 = new TH3F("hmass_S_fvtx_010","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_010 = new TH3F("hmass_S_fvtx_pp_010","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_010 = new TH3F("hmass_S_fvtx_mm_010","",300,0,15,48,0,12,60,-30,30);

	// 10-30
	hmass_S_fvtx_1030 = new TH3F("hmass_S_fvtx_1030","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_1030 = new TH3F("hmass_S_fvtx_pp_1030","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_1030 = new TH3F("hmass_S_fvtx_mm_1030","",300,0,15,48,0,12,60,-30,30);

	//30-50
	hmass_S_fvtx_3050 = new TH3F("hmass_S_fvtx_3050","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_3050 = new TH3F("hmass_S_fvtx_pp_3050","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_3050 = new TH3F("hmass_S_fvtx_mm_3050","",300,0,15,48,0,12,60,-30,30);

	// 40-60
	hmass_S_fvtx_4060 = new TH3F("hmass_S_fvtx_4060","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_4060 = new TH3F("hmass_S_fvtx_pp_4060","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_4060 = new TH3F("hmass_S_fvtx_mm_4060","",300,0,15,48,0,12,60,-30,30);

	// 50-84
	hmass_S_fvtx_5084 = new TH3F("hmass_S_fvtx_5084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_5084 = new TH3F("hmass_S_fvtx_pp_5084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_5084 = new TH3F("hmass_S_fvtx_mm_5084","",300,0,15,48,0,12,60,-30,30);

	// 0-20
	hmass_S_fvtx_020 = new TH3F("hmass_S_fvtx_020","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_020 = new TH3F("hmass_S_fvtx_pp_020","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_020 = new TH3F("hmass_S_fvtx_mm_020","",300,0,15,48,0,12,60,-30,30);

	// 0-5
	hmass_S_fvtx_2040 = new TH3F("hmass_S_fvtx_2040","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_2040 = new TH3F("hmass_S_fvtx_pp_2040","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_2040 = new TH3F("hmass_S_fvtx_mm_2040","",300,0,15,48,0,12,60,-30,30);

	// 5-10
	hmass_S_fvtx_4084 = new TH3F("hmass_S_fvtx_4084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_4084 = new TH3F("hmass_S_fvtx_pp_4084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_4084 = new TH3F("hmass_S_fvtx_mm_4084","",300,0,15,48,0,12,60,-30,30);

	// 40-84
	hmass_S_fvtx_6084 = new TH3F("hmass_S_fvtx_6084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_pp_6084 = new TH3F("hmass_S_fvtx_pp_6084","",300,0,15,48,0,12,60,-30,30);
	hmass_S_fvtx_mm_6084 = new TH3F("hmass_S_fvtx_mm_6084","",300,0,15,48,0,12,60,-30,30);

	/////////////////

	hmass_fvtx_below_3pT = new TH3F("hmass_fvtx_below_3pT","",300,0,15,48,0,12,88,-2.2,2.2);
	hmass_fvtx_above_3pT = new TH3F("hmass_fvtx_above_3pT","",300,0,15,48,0,12,88,-2.2,2.2);

	hmass_mismatch = new TH3F("hmass_mismatch","",300,0,15,2000,-1000,1000,88,-2.2,2.2);

	file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","READ");
	//file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive_pythia6.root","READ");
	fw_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
	fw_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
	fw_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");

	fw_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
	fw_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");

	if ( !fw_pT_dAu_fwd || !fw_pT_dAu_bwd || !fw_y_dAu || !fw_y_pp || !fw_pT_pp ){
		cout << "CAN NOT FIND WEIGHT-FUNCTION!!" << endl;
		exit(1);
	}

	if ( _dataset=="Run15pp200" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pp200.root","READ");
	}else if ( _dataset=="Run14HeAu200" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run14HeAu200.root","READ");
	}else if ( _dataset=="Run15pAu200" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pAu200.root","READ");
	}else if ( _dataset=="Run15pAl200" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pAl200.root","READ");
	}else if ( _dataset=="Run9pp200" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/acceff/outfile_jpsi_zweight_Run9pp200.root","READ");
	}else if ( _dataset=="Run9pp500" ){
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/acceff/outfile_jpsi_zweight_Run9pp500.root","READ");
	}else{
		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pp200.root","READ");
	}
	fw_z = (TF1*)file_weight1->Get("fwt_z");
	//fw_z_S = (TF1*)file_weight1->Get("fwt_z_S");
	//fw_z_N = (TF1*)file_weight1->Get("fwt_z_N");

	if ( !fw_z ){
		cout << "CAN NOT FIND WEIGHT-FUNCTION!!" << endl;
		exit(1);
	}else{
		cout << "OPEN: " << file_weight1->GetName() << endl;
		fw_z->Print();
	}

	cout << "OPEN weight file: " << file_weight0->GetName() << endl;
	cout << "OPEN weight file: " << file_weight1->GetName() << endl;


	cout << "SET SIM: " << _is_sim << endl;
	cout << "IS UPSILON: " << _is_upsilon << endl;
	cout << "IS PHI: " << _is_phi << endl;
	cout << "TARGET PID: " << _target_pid << endl;
	MUTOO::PRINT(std::cout,"**");
  return 0;
}

//______________________________________________________
int mFillAnaTree::InitRun(PHCompositeNode *top_node)
{
	MUTOO::PRINT(std::cout,"mFillAnaTree::InitRun");

	if ( first ){
		//MCHepMCParticleContainer
		hepmc_con = findNode::getClass<MCHepMCParticleContainer>(top_node,"MCHepMCParticleContainer");
		if (!hepmc_con){
			std::cout << PHWHERE << "mFillAnaTree:: MCHepMCParticleContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//SingleMuonContainer
		sngmucon = findNode::getClass<SingleMuonContainer>(top_node,"SingleMuonContainer");
		if (!sngmucon){
			std::cout << PHWHERE << "mFillAnaTree:: SingleMuonContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//FVTX MCSingleMuonContainer
		fvtxmc_sngmucon = findNode::getClass<MCSingleMuonFvtxContainer>(top_node,"MCSingleMuonFvtxContainer");
		if (!fvtxmc_sngmucon){
			std::cout << PHWHERE << "mFillAnaTree:: MCSingleMuonFvtxContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//DiMuonContainer
		dimucon = findNode::getClass<DiMuonContainer>(top_node,"DiMuonContainer");
		if (!dimucon){
			std::cout << PHWHERE << "mFillAnaTree:: DiMuonContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		// TrigRunLvl1
		trigrun= findNode::getClass<TrigRunLvl1>(top_node,"TrigRunLvl1");
		if (!trigrun){
			std::cout << PHWHERE << "mFillAnaTree:: TrigRunLvl1 not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//VtxOut
		vtx = findNode::getClass<VtxOut>(top_node,"VtxOut");
		if (!vtx){
			std::cout << PHWHERE << "mFillAnaTree:: VtxOut not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//TrigLvl1
		triglvl1 = findNode::getClass<TrigLvl1>(top_node,"TrigLvl1");
		if (!triglvl1){
			std::cout << PHWHERE << "mFillAnaTree:: TrigLvl1 not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		//FvtxCompactTrkMap
		//fvtxtrk_map = findNode::getClass<TFvtxCompactTrkMap>(top_node,"TFvtxCompactTrkMap");
	      	if (!fvtxtrk_map){
		        std::cout << PHWHERE << "mFillAnaTree:: TFvtxCompactTrkMap not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		// Trigger Emulator
		trigemcon= findNode::getClass<TriggerEmulatorContainer>(top_node,"TriggerEmulatorContainer");
		if (!trigemcon){
			std::cout << PHWHERE << "mFillAnaTree:: TriggerEmulatorContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}


		//SvxSegmentList
		//	svxseg_list = findNode::getClass<SvxSegmentList>(top_node,"SvxSegmentList");
	      	if (!svxseg_list){
	                std::cout << PHWHERE << "mFillAnaTree:: SvxSegmentList not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		phpythia_con = findNode::getClass<PHPythiaContainer>(top_node,"PHPythia");
		if (!phpythia_con){
			std::cout << PHWHERE << "mFillAnaTree:: PHPythiaContainer not in Node Tree" << std::endl;
			//return ABORTRUN;
		}

		// refit single muons
		sngrefitmucon = findNode::getClass<MutrRefitSingleMuonContainer>(top_node,"MutrRefitSingleMuonContainer");
		if(!sngrefitmucon){
		  std::cout << PHWHERE << "mFillAnaTree:  RefitSingleMuonContainer not in Node tree" << std::endl;
		}

		MUTOO::PRINT(std::cout,"mFillAnaTree::GetNodes");
		first = false;
	}


	MUTOO::PRINT(std::cout,"**");
	return 0;
}

//______________________________________________________
int mFillAnaTree::process_event(PHCompositeNode *top_node)
{

	_nevent++;

	//float mult_fvtx_map = fvtxtrk_map->count();

	float mult_fvtx_prim_cut = 0.0;
	float mult_svx_prim_cut = 0.0;
	float mult_fvtx_prim_cut_eta25 = 0.0;

	//float mult_fvtxN = sngmucon->get_Evt_Mult_FVTXN();
	//float mult_fvtxS = sngmucon->get_Evt_Mult_FVTXS();
	float mult_svx	= sngmucon->get_Evt_Mult_SVX();

	//scan TFvtxCompactTrkMap
	if ( fvtxtrk_map ){
		TFvtxCompactTrkMap::iterator iter( fvtxtrk_map->range() );
		while( TFvtxCompactTrkMap::const_pointer fvtx_ptr = iter.next() ) {
			float fvtx_chi2 = fvtx_ptr->get()->get_chi2_ndf();
			if ( fvtx_chi2>4.0 || isnan(fvtx_chi2) ) continue;
			if ( fvtx_ptr->get()->get_nhits()<2.5 ) continue;

			mult_fvtx_prim_cut++;

			float rapidity = fvtx_ptr->get()->get_fvtx_eta();
			if ( fabs(rapidity)<2.5 ) mult_fvtx_prim_cut_eta25++;
		} 
	}

	//scan SvxSegmentList
	if ( svxseg_list ){
		for(int isvx=0; isvx<svxseg_list->get_nSegments(); ++isvx) {
			SvxSegment *segment = svxseg_list->get_segment(isvx);

			if (isnan(segment->get3Momentum(0))) continue;
			if (isnan(segment->get3Momentum(1))) continue;
			if (isnan(segment->get3Momentum(2))) continue;
			if (isnan(segment->getMomentum())) continue;

			if (isnan(segment->getInnerMostProjectedPosition(0))) continue;
			if (isnan(segment->getInnerMostProjectedPosition(1))) continue;
			if (isnan(segment->getInnerMostProjectedPosition(2))) continue;

			mult_svx_prim_cut++;
		}    
	}else{
		mult_svx_prim_cut = mult_svx;
	}

	/*
	bool trig_S = false;
	bool trig_N = false;

	if ( triglvl1 ){
		unsigned int trig_scaled = triglvl1->get_lvl1_trigscaled();

		if ( (trig_scaled&0x00200000) || (trig_scaled&0x00800000) || (trig_scaled&0x02000000) || (trig_scaled&0x08000000) ) trig_S = true; 
		if ( (trig_scaled&0x00100000) || (trig_scaled&0x00400000) || (trig_scaled&0x01000000) || (trig_scaled&0x04000000) ) trig_N = true; 
	}
	*/


	//cout << "Fvtx Mult. pdst : " << mult_fvtxN+mult_fvtxS << ", prim. cut : " << mult_fvtx_prim_cut << endl;

	//get vertex information
	PHPoint sim_vertex = vtx->get_Vertex("SIM");
	PHPoint fvtx_vertex = vtx->get_Vertex("FVTX"); 
	//PHPoint fvtx_vertex_err = vtx->get_VertexError("FVTX");
		float simZ = sim_vertex.getZ();
	//float simX = sim_vertex.getX();
	//float simY = sim_vertex.getY();

	//float fvtxZ = sngmucon->get_Evt_fvtxZ();
	//float fvtxX = sngmucon->get_Evt_fvtxX();
	//float fvtxY = sngmucon->get_Evt_fvtxY();
	//float fvtxZ_err = sngmucon->get_Evt_fvtxZ_Err();
	//float fvtxX_err = sngmucon->get_Evt_fvtxX_Err();
	//float fvtxY_err = sngmucon->get_Evt_fvtxY_Err();

		//	float fvtxZ = dimucon->get_Evt_fvtxZ();
	// float fvtxX = dimucon->get_Evt_fvtxX();
	// float fvtxY = dimucon->get_Evt_fvtxY();
	// float fvtxZ_err = dimucon->get_Evt_fvtxZ_Err();
	// float fvtxX_err = dimucon->get_Evt_fvtxX_Err();
	// float fvtxY_err = dimucon->get_Evt_fvtxY_Err();

	vector<int> mu_vertex_begin; 
	mu_vertex_begin.clear();

	vector<float> mu_pT;
	vector<float> mu_rap;
	vector<float> mu_eta;
	vector<float> jpsi_pT;
	vector<float> jpsi_rap;

	int count = 0;

	if ( hepmc_con ){
		int n_hepmc_p = hepmc_con->get_nMCHepMCParticles();
		//Check double-Jpsi event
		for (int ip=0; ip<n_hepmc_p; ip++){
			MCHepMCParticle *hepmc = hepmc_con->get_MCHepMCParticle(ip);
			int mc_pid 		= hepmc->get_mc_pid();
			int mc_status = hepmc->get_mc_status();

			if ( mc_status==1 && abs(mc_pid)==13 ){
				mu_vertex_begin.push_back(hepmc->get_mc_vertex_begin());

				float mc_px 	= hepmc->get_mc_px();
				float mc_py 	= hepmc->get_mc_py();
				float mc_pz 	= hepmc->get_mc_pz();
				float mc_e		= hepmc->get_mc_e();

				float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
				float mc_ptot	= sqrt(mc_pT*mc_pT + mc_pz*mc_pz); 
				float mc_eta	= 0.5*log((mc_ptot+mc_pz)/(mc_ptot-mc_pz));
				float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

				mu_pT.push_back(mc_pT);
				mu_rap.push_back(mc_y);
				mu_eta.push_back(mc_eta);
			}    

			if ( mc_status==2 && mc_pid==_target_pid ){
				float mc_px 	= hepmc->get_mc_px();
				float mc_py 	= hepmc->get_mc_py();
				float mc_pz 	= hepmc->get_mc_pz();
				float mc_e		= hepmc->get_mc_e();

				float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
				float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

				//cout << "mc_pT: " << mc_pT << endl;

				jpsi_pT.push_back(mc_pT);
				jpsi_rap.push_back(mc_y);

				count++;
			}
		}//n_hepmc_p
	}else if ( phpythia_con ){

		for (unsigned int ip=0; ip<phpythia_con->size(); ip++){
			TMCParticle *part = phpythia_con->getParticle(ip);

			/*
			if ( part->GetKF()==443 || abs(part->GetKF())==13 ){
				cout
					<< part->GetKF() << " "
					<< part->GetKS() << " "
					<< part->GetPx() << " "
					<< part->GetPy() << " "
					<< part->GetPz() << " "
					<< endl;
			}
			*/

			if ( abs(part->GetKF())==13 && part->GetKS()==1 ){
				float mc_px = part->GetPx();
				float mc_py = part->GetPy();
				float mc_pz = part->GetPz();
				float mc_e = part->GetEnergy();

				float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
				float mc_p		= sqrt(mc_pT*mc_pT + mc_pz*mc_pz);
				float mc_eta	= 0.5*log((mc_p+mc_pz)/(mc_p-mc_pz));
				float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

				mu_pT.push_back(mc_pT);
				mu_rap.push_back(mc_y);
				mu_eta.push_back(mc_eta);
			}


			if ( part->GetKF()==_target_pid && part->GetKS()==11 ){
				float mc_px = part->GetPx();
				float mc_py = part->GetPy();
				float mc_pz = part->GetPz();
				float mc_e = part->GetEnergy();

				float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
				float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

				jpsi_pT.push_back(mc_pT);
				jpsi_rap.push_back(mc_y);

				count++;
			}
		}//ip

	}

	if ( count>1 ){
		cout << "Skip event containing multiple " << _target_pid << ":" << count << ", muon: " << int(mu_pT.size()) << endl;
		cout << "pT0 : " << jpsi_pT[0] << ", y0: " << jpsi_rap[0] << ", pT1: " << jpsi_pT[1] << ", y1: " << jpsi_rap[1] << endl;
		_nevent_skip++;
		return ABORTEVENT;
	}

	if ( count>0 && (_is_sim==3 || _is_sim==4) ){
		cout << "Skip event containing J/psi in cc-bar/bb-bar sim: " << count << endl;
		return ABORTEVENT;
	}  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(dimucon)
	  {

	    int ndim = dimucon->get_nDiMuons();

	    for (int id=0; id<ndim; id++)
	      {

		DiMuon *dimu = dimucon->get_DiMuon(id);
	
		int run_num = trigrun->get_run_number();
		
		//	float dimu_pT = dimu->get_pT();
	
		float hits_rap = dimu->get_rapidity();
		//float centrality = dimucon->get_Evt_Cent();

		///////////////
		//	float mass = dimu->get_mass();
		//	float mass_sngtrk = dimu->get_mass_sngtrk();
		float mass_fvtxmutr = dimu->get_mass_fvtxmutr();
		///////////////
		
		float Tr0_xst1 = dimu->get_Tr0_xst1();
		float Tr0_yst1 = dimu->get_Tr0_yst1();
		float Tr1_xst1 = dimu->get_Tr1_xst1();
		float Tr1_yst1 = dimu->get_Tr1_yst1();

		float Tr0_xst2 = dimu->get_Tr0_xst2();
		float Tr0_yst2 = dimu->get_Tr0_yst2();
		float Tr1_xst2 = dimu->get_Tr1_xst2();
		float Tr1_yst2 = dimu->get_Tr1_yst2();

		float Tr0_xst3 = dimu->get_Tr0_xst3();
		float Tr0_yst3 = dimu->get_Tr0_yst3();
		float Tr1_xst3 = dimu->get_Tr1_xst3();
		float Tr1_yst3 = dimu->get_Tr1_yst3()
		  ;

		float Tr0_dr_fvtx = dimu->get_Tr0_dr_fvtx();
		float Tr1_dr_fvtx = dimu->get_Tr1_dr_fvtx();
		float Tr0_dtheta_fvtx = dimu->get_Tr0_dtheta_fvtx();
		float Tr1_dtheta_fvtx = dimu->get_Tr1_dtheta_fvtx();
		float Tr0_dphi_fvtx = dimu->get_Tr0_dphi_fvtx();
		float Tr1_dphi_fvtx = dimu->get_Tr1_dphi_fvtx();
		
		float Tr0_DG0 = dimu->get_Tr0_DG0();
		float Tr1_DG0 = dimu->get_Tr1_DG0();
		float Tr0_DDG0 = dimu->get_Tr0_DDG0();
		float Tr1_DDG0 = dimu->get_Tr1_DDG0();

		float nhits_fvtx_Tr0 = dimu->get_Tr0_nhits_fvtx();
		float nhits_fvtx_Tr1 = dimu->get_Tr1_nhits_fvtx();
		float Tr0_rap = dimu->get_Tr0_rapidity();
		float Tr1_rap = dimu->get_Tr1_rapidity();

		float Charge = dimu->get_charge();
		float idhits_Tr0 = dimu->get_Tr0_idhits();
		float idhits_Tr1 = dimu->get_Tr1_idhits();
		float chi2_Tr0 = dimu->get_Tr0_trchi2();
		float chi2_Tr1 = dimu->get_Tr1_trchi2();
		short last_gap_0 = dimu->get_Tr0_lastgap();
		short last_gap_1 = dimu->get_Tr1_lastgap();

		float px_Tr0 = dimu->get_Tr0_px_fvtxmutr();
		float px_Tr1 = dimu->get_Tr1_px_fvtxmutr();
		float py_Tr0 = dimu->get_Tr0_py_fvtxmutr();
		float py_Tr1 = dimu->get_Tr1_py_fvtxmutr();
		float pz_Tr0 = dimu->get_Tr0_pz_fvtxmutr();
		float pz_Tr1 = dimu->get_Tr1_pz_fvtxmutr();

		// float mutr_px_Tr0 = dimu->get_Tr0_px();
		// float mutr_px_Tr1 = dimu->get_Tr1_px();
		// float mutr_py_Tr0 = dimu->get_Tr0_py();
		// float mutr_py_Tr1 = dimu->get_Tr1_py();
		// float mutr_pz_Tr0 = dimu->get_Tr0_pz();
		// float mutr_pz_Tr1 = dimu->get_Tr1_pz();

		//	unsigned int trig_scale = triglvl1->get_lvl1_trigscaled();
		bool same_evt = dimu->get_same_event();

		//	float chi2_fvtxmutr_Tr0 = dimu->get_Tr0_chi2_fvtxmutr();
		//	float chi2_fvtxmutr_Tr1 = dimu->get_Tr1_chi2_fvtxmutr();
		
		//	float chi2_fvtx_Tr0 = dimu->get_Tr0_chi2_fvtx();
		//	float chi2_fvtx_Tr1 = dimu->get_Tr1_chi2_fvtx();

		// float weight_pT_dAu = 1.0;
		// float weight_y_dAu = 1.0;
		float weight_pT_pp = 1.0;
		float weight_y_pp = 1.0;
		float weight_z = 1.0;

		//	float fvtx_prob1 = 0;
		//	float fvtx_prob2 = 0;

		double dr_cut_N = 0;
		double dr_cut_S = 0;
		double dtheta_cut_N = 0;
		double dtheta_cut_S = 0;
		double dphi_cut_N = 0;
		double dphi_cut_S = 0;

		// double track0_x = 0;
		// double track0_y = 0;
		// double track1_x = 0;
		// double track1_y = 0;
		// double dx = 0;
		// double dy = 0;
  
		double DG0_N = 0;
		double DG0_S = 0;
		double DDG0_N = 0;
		double DDG0_S = 0;

		float dimu_pT = sqrt( (px_Tr0+px_Tr1)*(px_Tr0+px_Tr1) + (py_Tr0+py_Tr1)*(py_Tr0+py_Tr1) );

		double octant = abs(  ( ((int((atan2(Tr0_xst1,Tr0_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) - ( ((int((atan2(Tr1_xst1,Tr1_yst1)+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8 ) ); 		
	
		// December 18th cuts using root5 - cuts from all runs
		// DG0_N = Tr0_DG0<(8.49898+100.736/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DG0<(8.49898+100.736/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// DG0_S = Tr0_DG0<(8.79546+376.293/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DG0<(8.79546+376.293/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// DDG0_N = Tr0_DDG0<(4.51376+137.931/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DDG0<(4.51376+137.931/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// DDG0_S = Tr0_DDG0<(4.71542+119.735/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DDG0<(4.71542+119.735/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dr_cut_N = Tr0_dr_fvtx<(1.69833+24.7377/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dr_fvtx<(1.69833+24.7377/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dr_cut_S = Tr0_dr_fvtx<(1.54454+27.5314/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dr_fvtx<(1.54454+27.5314/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dtheta_cut_N = Tr0_dtheta_fvtx<(0.0482243+0.732457/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dtheta_fvtx<(0.0482243+0.732457/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dtheta_cut_S = Tr0_dtheta_fvtx<(0.0559102+0.580277/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dtheta_fvtx<(0.0559102+0.580277/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dphi_cut_N = Tr0_dphi_fvtx<(0.117266+1.53911/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dphi_fvtx<(0.117266+1.53911/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		// dphi_cut_S = Tr0_dphi_fvtx<(0.113169+1.6207/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dphi_fvtx<(0.113169+1.6207/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));

		// cuts December 30th - using only good runs
		dr_cut_N = Tr0_dr_fvtx<(1.76048+23.9016/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dr_fvtx<(1.76048+23.9016/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		dr_cut_S = Tr0_dr_fvtx<(1.39017+29.279/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dr_fvtx<(1.39017+29.279/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		dtheta_cut_N = Tr0_dtheta_fvtx<(0.0499702+0.694368/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dtheta_fvtx<(0.0499702+0.694368/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		dtheta_cut_S = Tr0_dtheta_fvtx<(0.0555553+0.557225/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dtheta_fvtx<(0.0555553+0.557225/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		dphi_cut_N = Tr0_dphi_fvtx<(0.11936+1.44372/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dphi_fvtx<(0.11936+1.44372/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		dphi_cut_S = Tr0_dphi_fvtx<(0.116105+1.51064/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_dphi_fvtx<(0.116105+1.51064/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		DG0_N = Tr0_DG0<(8.14581+106.858/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DG0<(8.14581+106.858/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		DG0_S = Tr0_DG0<(8.93709+366.856/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DG0<(8.93709+366.856/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		DDG0_N = Tr0_DDG0<(4.3602+142.651/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DDG0<(4.3602+142.651/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));
		DDG0_S = Tr0_DDG0<(4.66658+120.627/(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0)) && Tr1_DDG0<(4.66658+120.627/(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1));

		double asym = abs( (sqrt(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0 ) - sqrt(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1 ))/ (sqrt(px_Tr0 * px_Tr0 + py_Tr0*py_Tr0 + pz_Tr0*pz_Tr0 ) + sqrt(px_Tr1 * px_Tr1 + py_Tr1*py_Tr1 + pz_Tr1*pz_Tr1 )) ) ;
		
		// PPG188 fid cut N
		double fiducial_N = !(  atan2(Tr0_yst3,Tr0_xst3)>0.35 && atan2(Tr0_yst3,Tr0_xst3)<0.77  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>0.77 && atan2(Tr0_yst3,Tr0_xst3)<1.15 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>160 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst3,Tr0_xst3)>-2.75 && atan2(Tr0_yst3,Tr0_xst3)<-1.92 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)>237 && sqrt(Tr0_yst3*Tr0_yst3 + Tr0_xst3*Tr0_xst3)<275  ) && !(  atan2(Tr0_yst2,Tr0_xst2)>-2.75 && atan2(Tr0_yst2,Tr0_xst2)<-1.92 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)>140 && sqrt(Tr0_yst2*Tr0_yst2 + Tr0_xst2*Tr0_xst2)<163  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>-2.75 && atan2(Tr0_yst1,Tr0_xst1)<-2 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>76 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<80  )  && !(  atan2(Tr0_yst1,Tr0_xst1)>0.52 && atan2(Tr0_yst1,Tr0_xst1)<1.11 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)>60 && sqrt(Tr0_yst1*Tr0_yst1 + Tr0_xst1*Tr0_xst1)<75 )  && !(  atan2(Tr1_yst3,Tr1_xst3)>0.35 && atan2(Tr1_yst3,Tr1_xst3)<0.77  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>0.77 && atan2(Tr1_yst3,Tr1_xst3)<1.15 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>160 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst3,Tr1_xst3)>-2.75 && atan2(Tr1_yst3,Tr1_xst3)<-1.92 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)>237 && sqrt(Tr1_yst3*Tr1_yst3 + Tr1_xst3*Tr1_xst3)<275  ) && !(  atan2(Tr1_yst2,Tr1_xst2)>-2.75 && atan2(Tr1_yst2,Tr1_xst2)<-1.92 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)>140 && sqrt(Tr1_yst2*Tr1_yst2 + Tr1_xst2*Tr1_xst2)<163  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>-2.75 && atan2(Tr1_yst1,Tr1_xst1)<-2 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>76 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<80  )  && !(  atan2(Tr1_yst1,Tr1_xst1)>0.52 && atan2(Tr1_yst1,Tr1_xst1)<1.11 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)>60 && sqrt(Tr1_yst1*Tr1_yst1 + Tr1_xst1*Tr1_xst1)<75 );

		// ppg228 cut (no ppg188 fid cut S exists)
		//	double fiducial_S = !( atan2(Tr0_yst1,Tr0_xst1)>0 && atan2(Tr0_yst1,Tr0_xst1)<0.4 ) && !( atan2(Tr0_yst2,Tr0_xst2)>0 && atan2(Tr0_yst2,Tr0_xst2)<0.4 ) && !( atan2(Tr0_yst3,Tr0_xst3)>0 && atan2(Tr0_yst3,Tr0_xst3)<0.4 ) && !( atan2(Tr0_yst1,Tr0_xst1)>1.55 && atan2(Tr0_yst1,Tr0_xst1)<1.9 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>1.55 && atan2(Tr0_yst2,Tr0_xst2)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>1.55 && atan2(Tr0_yst3,Tr0_xst3)<1.9 ) && !( atan2(Tr0_yst3,Tr0_xst3)>-0.60 && atan2(Tr0_yst3,Tr0_xst3)<-0.40 ) && !( atan2(Tr1_yst1,Tr1_xst1)>0 && atan2(Tr1_yst1,Tr1_xst1)<0.4 ) && !( atan2(Tr1_yst2,Tr1_xst2)>0 && atan2(Tr1_yst2,Tr1_xst2)<0.4 ) && !( atan2(Tr1_yst3,Tr1_xst3)>0 && atan2(Tr1_yst3,Tr1_xst3)<0.4 ) && !( atan2(Tr1_yst1,Tr1_xst1)>1.55 && atan2(Tr1_yst1,Tr1_xst1)<1.9 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>1.55 && atan2(Tr1_yst2,Tr1_xst2)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>1.55 && atan2(Tr1_yst3,Tr1_xst3)<1.9 ) && !( atan2(Tr1_yst3,Tr1_xst3)>-0.60 && atan2(Tr1_yst3,Tr1_xst3)<-0.40 ) && !( atan2(Tr0_yst1,Tr0_xst1)>-2.35 && atan2(Tr0_yst1,Tr0_xst1)<-1.95 ) &&  !( atan2(Tr0_yst2,Tr0_xst2)>-2.35 && atan2(Tr0_yst2,Tr0_xst2)<-1.95 ) &&  !( atan2(Tr0_yst3,Tr0_xst3)>-2.35 && atan2(Tr0_yst3,Tr0_xst3)<-1.95 ) && !( atan2(Tr1_yst1,Tr1_xst1)>-2.35 && atan2(Tr1_yst1,Tr1_xst1)<-1.95 ) &&  !( atan2(Tr1_yst2,Tr1_xst2)>-2.35 && atan2(Tr1_yst2,Tr1_xst2)<-1.95 ) &&  !( atan2(Tr1_yst3,Tr1_xst3)>-2.35 && atan2(Tr1_yst3,Tr1_xst3)<-1.95 ); 
	
				
		weight_z = fw_z->Eval(simZ);
			
		
		if ( fabs(hits_rap)>0.6 && fabs(hits_rap)<3.0 ){

		  // if ( dimu_rap<0 ){
		  //   //	weight_y_dAu = fw_y_dAu->Eval(dimu_rap);
		  // 	weight_pT_dAu = (dimu_pT<8.0) ? fw_pT_dAu_bwd->Eval(dimu_pT) : fw_pT_dAu_bwd->Eval(8.0);
		  // }else{
		  // 	weight_y_dAu = fw_y_dAu->Eval(dimu_rap);
		  // 	weight_pT_dAu = (dimu_pT<8.0) ? fw_pT_dAu_fwd->Eval(dimu_pT) : fw_pT_dAu_fwd->Eval(8.0);
		  // }

		  weight_y_pp = fw_y_pp->Eval(hits_rap);
		  weight_pT_pp = (dimu_pT<8.0) ? fw_pT_pp->Eval(dimu_pT) : fw_pT_pp->Eval(8.0);

		}


		// if (nhits_fvtx_Tr0>0){
		//   int ndf = 2*nhits_fvtx_Tr0 - 1; 
		//   float chi2 = chi2_fvtx_Tr0;
		//   fvtx_prob1 = TMath::Prob(chi2*ndf, ndf);
		// }

		// if (nhits_fvtx_Tr1>0){
		//   int ndf = 2*nhits_fvtx_Tr1 - 1; 
		//   float chi2 = chi2_fvtx_Tr1;
		//   fvtx_prob2 = TMath::Prob(chi2*ndf, ndf);
		// }

		int badrun= 0;


		//cout << "run num: " << run_num << endl;

		std::string  filename = "badruns_N.txt";
	       

		ifstream run_files(filename.c_str() );
		if(run_files)
		  {
		    do 
		      {
			badrun = 0;
			run_files >> badrun;
			//run_array[counter] = run;
			//	counter++;
			if(badrun == run_num)
			  continue;

		      }while(run_files.good());

		  } // end if

		if(badrun == run_num)
		  continue;


		// // 	///////////////////////////////////////////////////////////////////////
		int nsng = sngmucon->get_nSingleMuons();
		
		for (int is=0; is<nsng; is++)
		  {
		    SingleMuon *sngmu1 = sngmucon->get_SingleMuon(is);
		 
		    // bool sng1_MUID1D = sngmu1->get_MUID1D();
		    float sng1_px = sngmu1->get_px_fvtxmutr();
		    float sng1_py = sngmu1->get_py_fvtxmutr();
		    float sng1_pz = sngmu1->get_pz_fvtxmutr();
		    // float mutr_sng1_px= sngmu1->get_px();
		    // float mutr_sng1_py = sngmu1->get_py();
		    // float mutr_sng1_pz = sngmu1->get_pz();
		   
		    float cent1 = sngmucon->get_Evt_Cent();
  		      
		    MCSingleMuonFvtx *fvtxmc_sngmu1 = fvtxmc_sngmucon->get_MCSingleMuonFvtx(is); //"with the same index"
		    float mc_pid1 = fvtxmc_sngmu1->get_mc_g_pid();

		    // float nhits_fvtx1 = fvtxmc_sngmu1->get_mc_hits_fvtx();
		    // float nhits_svx1 = fvtxmc_sngmu1->get_mc_hits_svx();
		    // float nhits_g_fvtx1 = fvtxmc_sngmu1->get_mc_hits_fvtx_true();
		    // float nhits_g_svx1 = fvtxmc_sngmu1->get_mc_hits_svx_true();

		    // float mismatch_num1 = nhits_g_fvtx1 + nhits_g_svx1;
		    // float mismatch_denom1 = nhits_fvtx1 + nhits_svx1;
		    // float mismatch_ratio1 = mismatch_num1 / mismatch_denom1;

		    for (int it=0; it<nsng; it++)
		      {
			SingleMuon *sngmu2 = sngmucon->get_SingleMuon(it);
			
			//	bool sng2_MUID1D = sngmu2->get_MUID1D();
			float sng2_px = sngmu2->get_px_fvtxmutr();
			float sng2_py = sngmu2->get_py_fvtxmutr();
			float sng2_pz = sngmu2->get_pz_fvtxmutr();

			// float mutr_sng2_px= sngmu2->get_px();
		    	// float mutr_sng2_py = sngmu2->get_py();
		    	// float mutr_sng2_pz = sngmu2->get_pz();
		
			//	float cent2 = sngmucon->get_Evt_Cent();

			MCSingleMuonFvtx *fvtxmc_sngmu2 = fvtxmc_sngmucon->get_MCSingleMuonFvtx(it); // "with the same index"
			float mc_pid2 = fvtxmc_sngmu2->get_mc_g_pid();
		
		

			float centrality = 0;
			
			//	if(cent1 == cent2 )
			centrality = cent1;

			
			//	int trig1D_N =  trigemcon->GetLL1Mu1DNorth();
			//	int trig1D_S = trigemcon->GetLL1Mu1DSouth();
 		      
			//	if(fvtx_prob1<0.05 && fvtx_prob2<0.05)  
			//	continue;
			//////////////////////////////////////////////

			//	unsigned int trig_scale = triglvl1->get_lvl1_trigscaled();

			//	if( ( refit_sng1_px == mutr_sng1_px ) && ( refit_sng1_py == mutr_sng1_py ) && ( refit_sng1_pz == mutr_sng1_pz ) ) 
			{
			  // if( (refit_sng2_px == mutr_sng2_px ) && ( refit_sng2_py == mutr_sng2_py ) && ( refit_sng2_pz == mutr_sng2_pz ))
			  {		    
			    if( (px_Tr0 == sng1_px) &&  (py_Tr0 == sng1_py)  &&  (pz_Tr0 == sng1_pz))// && (sng1_MUID1D) ) 
			      {
				if( abs(mc_pid1) == 13)	
				  {
				    if( (px_Tr1 == sng2_px) &&  (py_Tr1 == sng2_py)  &&  (pz_Tr1 == sng2_pz))// && (sng2_MUID1D) ) 
				      {
					if( abs(mc_pid2) == 13)	
					  {
					    //////////////////////////////////////////////////////////////////////////////////////////////////
					    ///////////// analysis cuts currently for PPG228 or PPG188
					    //	if (fvtx_prob1<0.05 && fvtx_prob2<0.05)  
					    // continue;

					    //cout << "inside loop: " << mc_pid1 << " , " << mc_pid2 << endl;

					    //  if( abs(mc_pid1) != 13)
					    //   hmass_mismatch->Fill(mass_fvtxmutr,mc_pid1,hits_rap);

					    // if( abs(mc_pid2) != 13)
					    //   hmass_mismatch->Fill(mass_fvtxmutr,mc_pid2,hits_rap);

					    {
					      if( (nhits_fvtx_Tr0 > 0) && (nhits_fvtx_Tr1 > 0) )  
						{
						  //	  if(chi2_fvtx_Tr0<10 && chi2_fvtx_Tr1<10 && chi2_fvtx_Tr0>0 && chi2_fvtx_Tr1>0)
						  {
						    //   if(chi2_fvtxmutr_Tr0<10 && chi2_fvtxmutr_Tr0>0 && chi2_fvtxmutr_Tr1<10 && chi2_fvtxmutr_Tr1>0)
						    {
						      if(mass_fvtxmutr > 1.5)  // if(mass > 1.5)
							{
							  if(same_evt == 1)
							    {
							      if( (idhits_Tr0 > 14) &&  (idhits_Tr1 > 14) )
								{
								  if(abs(simZ) < 30)   //	if(abs(fvtxZ) < 40) 
								    {
								      if( asym < 1 )
									{
									  if( (chi2_Tr0 < 23) && (chi2_Tr1 < 23) )
									    {
									      if(octant > 0)
										{
										  if( (last_gap_0 >= 3) &&  (last_gap_1 >= 3) )
										    {
										      //////////////////////////////// 
										
										      if( (hits_rap >=1.2) && (hits_rap <=2.2))// && (trig_scale &0x00100000) && (trig1D_N == 1) )//  && Tr0_rap >=1.2 && Tr0_rap <=2.4 && Tr1_rap >=1.2 && Tr1_rap <=2.4 )   // North
											{
											  if(dr_cut_N && dtheta_cut_N && dphi_cut_N)
											    {
											      if(DG0_N && DDG0_N && ( abs(pz_Tr0) > 2) && (abs(pz_Tr1) > 2) && (pz_Tr0 > 0) && (pz_Tr1 > 0) && fiducial_N) 
												{
												  ///////////////////////////////////////////////////////////////////////////
												  //  if( (trig_scale &0x00000010) && (abs(simZ < 10) ) )
												    {
												      if(Charge == 0)
													{
													  //  if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_ul_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);  
													      hdimu_pT_y_zvtx_cuts->Fill(dimu_pT, hits_rap, simZ);
													      hmass_N_fvtx->Fill(mass_fvtxmutr,dimu_pT,simZ);
													    }

													  if( (centrality>20) && (centrality<=40))
													    hmass_N_fvtx_2040->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>60) && (centrality<=84))
													    hmass_N_fvtx_6084->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>10) && (centrality<=30))
													    hmass_N_fvtx_1030->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>0) && (centrality<=10))
													    hmass_N_fvtx_010->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>0) && (centrality<=20))
													    hmass_N_fvtx_020->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>30) && (centrality<=50))
													    hmass_N_fvtx_3050->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>40) && (centrality<=60))
													    hmass_N_fvtx_4060->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
													  if( (centrality>50) && (centrality<=84))
													    hmass_N_fvtx_5084->Fill(mass_fvtxmutr,dimu_pT,simZ);

													  if( (centrality>40) && (centrality<=84))
													    hmass_N_fvtx_4084->Fill(mass_fvtxmutr,dimu_pT,simZ);

													  if(dimu_pT < 3)
													    hmass_fvtx_below_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);  // to check if 2nd gaussian is pT dep.
													  else
													    hmass_fvtx_above_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);

													  // if( (mismatch_ratio1 < 0.5) || (mismatch_ratio2 < 0.5) )
													  //	hmass_mismatch->Fill(mass_fvtxmutr, dimu_pT, hits_rap);

													}
												      if(Charge > 0)
													{

													  // if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_pp_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);  
													      hmass_N_fvtx_pp->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													    }

													  if( (centrality>20) && (centrality<=40))
													    hmass_N_fvtx_pp_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>60) && (centrality<=84))
													    hmass_N_fvtx_pp_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>10) && (centrality<=30))
													    hmass_N_fvtx_pp_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>0) && (centrality<=10))
													    hmass_N_fvtx_pp_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>0) && (centrality<=20))
													    hmass_N_fvtx_pp_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>30) && (centrality<=50))
													    hmass_N_fvtx_pp_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>40) && (centrality<=60))
													    hmass_N_fvtx_pp_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>50) && (centrality<=84))
													    hmass_N_fvtx_pp_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

													  if( (centrality>40) && (centrality<=84))
													    hmass_N_fvtx_pp_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

													}
												      if(Charge < 0)
													{

													  // if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);   
													      hmass_N_fvtx_mm->Fill(mass_fvtxmutr,dimu_pT,simZ);
													    }

													  if( (centrality>20) && (centrality<=40))
													    hmass_N_fvtx_mm_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>60) && (centrality<=84))
													    hmass_N_fvtx_mm_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>10) && (centrality<=30))
													    hmass_N_fvtx_mm_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>0) && (centrality<=10))
													    hmass_N_fvtx_mm_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>0) && (centrality<=20))
													    hmass_N_fvtx_mm_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>30) && (centrality<=50))
													    hmass_N_fvtx_mm_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>40) && (centrality<=60))
													    hmass_N_fvtx_mm_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
													  if( (centrality>50) && (centrality<=84))
													    hmass_N_fvtx_mm_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

													  if( (centrality>40) && (centrality<=84))
													    hmass_N_fvtx_mm_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													}
												    }
												  ///////////////////////////////////////////////////////////////////////////
												  // if( ( (trig_scale &0x00000001) ||  (trig_scale &0x00000002) ) && (abs(simZ < 30) ) )
												  //   {
												  //     if(Charge == 0)
												  // 	{
												  // 	  // if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_ul_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);  
												  // 	      hdimu_pT_y_zvtx_cuts->Fill(dimu_pT, hits_rap, simZ);
												  // 	      hmass_N_fvtx->Fill(mass_fvtxmutr,dimu_pT,simZ);
												  // 	    }

												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_N_fvtx_2040->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_N_fvtx_6084->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_N_fvtx_1030->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_N_fvtx_010->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_N_fvtx_020->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_N_fvtx_3050->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_N_fvtx_4060->Fill(mass_fvtxmutr,dimu_pT,simZ);
													
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_N_fvtx_5084->Fill(mass_fvtxmutr,dimu_pT,simZ);

												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_N_fvtx_4084->Fill(mass_fvtxmutr,dimu_pT,simZ);

												  // 	  if(dimu_pT < 3)
												  // 	    hmass_fvtx_below_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);  // to check if 2nd gaussian is pT dep.
												  // 	  else
												  // 	    hmass_fvtx_above_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);

												  // 	  // if( (mismatch_ratio1 < 0.5) || (mismatch_ratio2 < 0.5) )
												  // 	  //	hmass_mismatch->Fill(mass_fvtxmutr, dimu_pT, hits_rap);

												  // 	}
												  //     if(Charge > 0)
												  // 	{

												  // 	  //	  if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_pp_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);  
												  // 	      hmass_N_fvtx_pp->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												  // 	    }

												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_N_fvtx_pp_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_N_fvtx_pp_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_N_fvtx_pp_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_N_fvtx_pp_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_N_fvtx_pp_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_N_fvtx_pp_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_N_fvtx_pp_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_N_fvtx_pp_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_N_fvtx_pp_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

												  // 	}
												  //     if(Charge < 0)
												  // 	{

												  // 	  //	  if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);   
												  // 	      hmass_N_fvtx_mm->Fill(mass_fvtxmutr,dimu_pT,simZ);
												  // 	    }

												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_N_fvtx_mm_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_N_fvtx_mm_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_N_fvtx_mm_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_N_fvtx_mm_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_N_fvtx_mm_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_N_fvtx_mm_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_N_fvtx_mm_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
													
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_N_fvtx_mm_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 

												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_N_fvtx_mm_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												  // 	}
												  //   }
												  //////////////////////////////////
												}
											    }
											}
										      
										      if( (hits_rap <=-1.2) && (hits_rap >=-2.2)) // && (trig_scale &0x00200000)  )  //Tr0_rap <=-1.2 && Tr0_rap >=-2.2  && Tr1_rap <=-1.2 && Tr1_rap >=-2.2)  //South
											{
											  if(DG0_S && DDG0_S && ( abs(pz_Tr0) > 2) && (abs(pz_Tr1) > 2) && (pz_Tr0 < 0) && (pz_Tr1 < 0))// && fiducial_S)
											    {
											      if(dr_cut_S && dtheta_cut_S && dphi_cut_S)
												{
												  ////////////////////////////////////////////////////
												  //  if( (trig_scale &0x00000010) && (abs(simZ < 10) ) )
												    {
												      if(Charge == 0)
													{
												      
													  // if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z); 
													      hmass_S_fvtx->Fill(mass_fvtxmutr,dimu_pT,simZ);  
													    }

													  if( (centrality>20) && (centrality<=40))
													    hmass_S_fvtx_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>60) && (centrality<=84))
													    hmass_S_fvtx_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>10) && (centrality<=30))
													    hmass_S_fvtx_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>0) && (centrality<=10))
													    hmass_S_fvtx_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>0) && (centrality<=20))
													    hmass_S_fvtx_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>30) && (centrality<=50))
													    hmass_S_fvtx_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>40) && (centrality<=60))
													    hmass_S_fvtx_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>50) && (centrality<=84))
													    hmass_S_fvtx_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if( (centrality>40) && (centrality<=84))
													    hmass_S_fvtx_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
													  if(dimu_pT < 3)
													    hmass_fvtx_below_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);
													  else
													    hmass_fvtx_above_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);
												      
													  //   if( (mismatch_ratio1 < 0.5) || (mismatch_ratio2 < 0.5) )
													  //	hmass_mismatch->Fill(mass_fvtxmutr, dimu_pT, hits_rap);
												      
													}
												      if(Charge > 0)
													{
												      
													  // if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_pp_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);   
													      hmass_S_fvtx_pp->Fill(mass_fvtxmutr,dimu_pT,simZ);
													    }
												      
													  if( (centrality>20) && (centrality<=40))
													    hmass_S_fvtx_pp_2040->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>60) && (centrality<=84))
													    hmass_S_fvtx_pp_6084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>10) && (centrality<=30))
													    hmass_S_fvtx_pp_1030->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>0) && (centrality<=10))
													    hmass_S_fvtx_pp_010->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>0) && (centrality<=20))
													    hmass_S_fvtx_pp_020->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>30) && (centrality<=50))
													    hmass_S_fvtx_pp_3050->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>40) && (centrality<=60))
													    hmass_S_fvtx_pp_4060->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>50) && (centrality<=84))
													    hmass_S_fvtx_pp_5084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>40) && (centrality<=84))
													    hmass_S_fvtx_pp_4084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													}
												      if(Charge < 0)
													{
												      
													  // if( (centrality>0) && (centrality<=100))
													    {
													      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z); 
													      hmass_S_fvtx_mm->Fill(mass_fvtxmutr,dimu_pT,simZ);  
													    }
												      
													  if( (centrality>20) && (centrality<=40))
													    hmass_S_fvtx_mm_2040->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>60) && (centrality<=84))
													    hmass_S_fvtx_mm_6084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>10) && (centrality<=30))
													    hmass_S_fvtx_mm_1030->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>0) && (centrality<=10))
													    hmass_S_fvtx_mm_010->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>0) && (centrality<=20))
													    hmass_S_fvtx_mm_020->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>30) && (centrality<=50))
													    hmass_S_fvtx_mm_3050->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>40) && (centrality<=60))
													    hmass_S_fvtx_mm_4060->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>50) && (centrality<=84))
													    hmass_S_fvtx_mm_5084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													  if( (centrality>40) && (centrality<=84))
													    hmass_S_fvtx_mm_4084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
													} // charge
												    } // simX south
												      ///////////////////////////////////////////////////
												  // if( ( (trig_scale &0x00000001) ||  (trig_scale &0x00000002) ) && (abs(simZ < 30) ) )
												  //   {
												  //     if(Charge == 0)
												  // 	{
												      
												  // 	  //if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z); 
												  // 	      hmass_S_fvtx->Fill(mass_fvtxmutr,dimu_pT,simZ);  
												  // 	    }

												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_S_fvtx_2040->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_S_fvtx_6084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_S_fvtx_1030->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_S_fvtx_010->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_S_fvtx_020->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_S_fvtx_3050->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_S_fvtx_4060->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_S_fvtx_5084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_S_fvtx_4084->Fill(mass_fvtxmutr,dimu_pT,simZ); 
												      
												  // 	  if(dimu_pT < 3)
												  // 	    hmass_fvtx_below_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);
												  // 	  else
												  // 	    hmass_fvtx_above_3pT->Fill(mass_fvtxmutr,dimu_pT,hits_rap);
												      
												  // 	  //   if( (mismatch_ratio1 < 0.5) || (mismatch_ratio2 < 0.5) )
												  // 	  //	hmass_mismatch->Fill(mass_fvtxmutr, dimu_pT, hits_rap);
												      
												  // 	}
												  //     if(Charge > 0)
												  // 	{
												      
												  // 	  // if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_pp_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);   
												  // 	      hmass_S_fvtx_pp->Fill(mass_fvtxmutr,dimu_pT,simZ);
												  // 	    }
												      
												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_S_fvtx_pp_2040->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_S_fvtx_pp_6084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_S_fvtx_pp_1030->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_S_fvtx_pp_010->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_S_fvtx_pp_020->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_S_fvtx_pp_3050->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_S_fvtx_pp_4060->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_S_fvtx_pp_5084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_S_fvtx_pp_4084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	}
												  //     if(Charge < 0)
												  // 	{
												      
												  // 	  // if( (centrality>0) && (centrality<=100))
												  // 	    {
												  // 	      hdimu_mm_pT_y_zvtx_w_pp->Fill(dimu_pT, hits_rap, simZ, weight_pT_pp*weight_y_pp*weight_z); 
												  // 	      hmass_S_fvtx_mm->Fill(mass_fvtxmutr,dimu_pT,simZ);  
												  // 	    }
												      
												  // 	  if( (centrality>20) && (centrality<=40))
												  // 	    hmass_S_fvtx_mm_2040->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>60) && (centrality<=84))
												  // 	    hmass_S_fvtx_mm_6084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>10) && (centrality<=30))
												  // 	    hmass_S_fvtx_mm_1030->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>0) && (centrality<=10))
												  // 	    hmass_S_fvtx_mm_010->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>0) && (centrality<=20))
												  // 	    hmass_S_fvtx_mm_020->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>30) && (centrality<=50))
												  // 	    hmass_S_fvtx_mm_3050->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>40) && (centrality<=60))
												  // 	    hmass_S_fvtx_mm_4060->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>50) && (centrality<=84))
												  // 	    hmass_S_fvtx_mm_5084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	  if( (centrality>40) && (centrality<=84))
												  // 	    hmass_S_fvtx_mm_4084->Fill(mass_fvtxmutr, dimu_pT, simZ);
												      
												  // 	} // charge
												  //   } // simX south
												      ///////////////////////////////////////////////
												} // dr
											    } // DG0
											} // hits_rap
											  /////////////////////////////////////
										    }
										}
									    }
									} // bbc < 30
								    } // octant
								} // // last gap
							    } // chi2
							}
						    }
						  }
						}
					    }
					  }
				      }
				  }
			      }
			  }
			}
		

			
			if( (hits_rap >=1.2 && hits_rap <=2.2 ) && Tr0_rap >=1.2 && Tr0_rap <=2.4 && Tr1_rap >=1.2 && Tr1_rap <=2.4 )   //  rapidity cut North
			  hdimu_pT_y_zvtx->Fill(dimu_pT, hits_rap, simZ);
			if( (hits_rap <=-1.2 && hits_rap >=-2.2 )&& Tr0_rap <=-1.2 && Tr0_rap >=-2.2  && Tr1_rap <=-1.2 && Tr1_rap >=-2.2)  // rapidity cut South
			  hdimu_pT_y_zvtx->Fill(dimu_pT, hits_rap, simZ);
			
			
			
		      }// int sngmu1 loop
		  }// int sngmu2 loop
	      } //int dimuon loop
	    
	    
	    
	    
	    
	    // hdimu_y_zvtx_w_dAu->Fill(dimu_pT, dimu_rap, simZ, weight_pT_dAu*weight_y_dAu*weight_z);
	    // hdimu_y_zvtx_w_dAu_pT->Fill(dimu_pT, dimu_rap, simZ, weight_pT_dAu*weight_z);
	    // hdimu_y_zvtx_w_dAu_y->Fill(dimu_pT, dimu_rap, simZ, weight_y_dAu*weight_z);
	    
	    //  if ( fabs(dimu_rap)>1.2 && fabs(dimu_rap)<2.2 ){
	    //	hdimu_pT_y_zvtx_w_pp->Fill(dimu_pT, dimu_rap, simZ, weight_pT_pp*weight_y_pp*weight_z);   // forward and backward rpaiidty no cut
	    //	hdimu_y_zvtx_w_pp_pT->Fill(dimu_pT, dimu_rap, simZ, weight_pT_pp*weight_z);
	    //	hdimu_y_zvtx_w_pp_y->Fill(dimu_pT, dimu_rap, simZ, weight_y_pp*weight_z);
	    //  }
	    
	  } // dimucon
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	count = 0;
	    
	//float eta_Jpsi = -999;
	float y_Jpsi = -999;
	//int Jpsi_vertex_begin = 0;

	if ( hepmc_con ){
		int n_hepmc_p = hepmc_con->get_nMCHepMCParticles();

		//Fill generated Jpsi information
		for (int ip=0; ip<n_hepmc_p; ip++){

			MCHepMCParticle *hepmc = hepmc_con->get_MCHepMCParticle(ip);

			int mc_pid 		= hepmc->get_mc_pid();
			int mc_status = hepmc->get_mc_status();

			if ( !(mc_status==2 && mc_pid==_target_pid) ) continue;

			float mc_px 	= hepmc->get_mc_px();
			float mc_py 	= hepmc->get_mc_py();
			float mc_pz 	= hepmc->get_mc_pz();
			float mc_e		= hepmc->get_mc_e();

			float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
			float mc_ptot	= sqrt(mc_pT*mc_pT + mc_pz*mc_pz); 
			float mc_eta	= 0.5*log((mc_ptot+mc_pz)/(mc_ptot-mc_pz));
			float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

			//Jpsi_vertex_begin = hepmc->get_mc_vertex_begin();

			//eta_Jpsi = mc_eta;
			y_Jpsi = mc_y;

			count++;

			hpT_eta_zvtx->Fill(mc_pT, mc_eta, simZ);
			hpT_y_zvtx->Fill(mc_pT, mc_y, simZ);

			if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 ){
			  if( (y_Jpsi > 1.2) && (y_Jpsi < 2.2))
			    hd_pT_y_zvtx_N->Fill(mc_pT, mc_y, simZ);
			  if( (y_Jpsi > -2.2) && (y_Jpsi < -1.2))
			    hd_pT_y_zvtx_S->Fill(mc_pT, mc_y, simZ);
			}

			float weight_pT_dAu = 1.0;
			float weight_y_dAu = 1.0;
			float weight_pT_pp = 1.0;
			float weight_y_pp = 1.0;
			float weight_z = 1.0;

			weight_z = fw_z->Eval(simZ);

			//	if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 )
			if ( fabs(y_Jpsi)>0.6 && fabs(y_Jpsi)<3.0 )
			  {
			    
			    if ( y_Jpsi<0 )
			      {
				weight_y_dAu = fw_y_dAu->Eval(y_Jpsi);
				weight_pT_dAu = (mc_pT<8.0) ? fw_pT_dAu_bwd->Eval(mc_pT) : fw_pT_dAu_bwd->Eval(8.0);
			      }
			    else
			      {
				weight_y_dAu = fw_y_dAu->Eval(y_Jpsi);
				weight_pT_dAu = (mc_pT<8.0) ? fw_pT_dAu_fwd->Eval(mc_pT) : fw_pT_dAu_fwd->Eval(8.0);
			      }
			    
			    weight_y_pp = fw_y_pp->Eval(y_Jpsi);
			    weight_pT_pp = (mc_pT<8.0) ? fw_pT_pp->Eval(mc_pT) : fw_pT_pp->Eval(8.0);
			    
			    if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 )
			      {
			      for (unsigned int ii=0; ii<mu_pT.size(); ii++)
				{
				  hmu_pT_y_zvtx->Fill(mu_pT[ii],mu_rap[ii],simZ);
				  hmu_pT_y_zvtx_w_dAu->Fill(mu_pT[ii],mu_rap[ii],simZ,weight_pT_dAu*weight_y_dAu*weight_z);
				  hmu_pT_y_zvtx_w_pp->Fill(mu_pT[ii],mu_rap[ii],simZ,weight_pT_pp*weight_y_pp*weight_z);
			      }
			      
			      if ( mu_pT.size()==2 )
				{
				  hmu_y0_y1->Fill(mu_rap[0],mu_rap[1]);
				}
			    }
			  }
			
			hpT_y_zvtx_w_z->Fill(mc_pT, mc_y, simZ, weight_z);

			hpT_y_zvtx_w_dAu->Fill(mc_pT, mc_y, simZ, weight_pT_dAu*weight_y_dAu*weight_z);
			hpT_y_zvtx_w_dAu_pT->Fill(mc_pT, mc_y, simZ, weight_pT_dAu*weight_z);
			hpT_y_zvtx_w_dAu_y->Fill(mc_pT, mc_y, simZ, weight_y_dAu*weight_z);

			hpT_y_zvtx_w_pp->Fill(mc_pT, mc_y, simZ, weight_pT_pp*weight_y_pp*weight_z);
			hpT_y_zvtx_w_pp_pT->Fill(mc_pT, mc_y, simZ, weight_pT_pp*weight_z);
			hpT_y_zvtx_w_pp_y->Fill(mc_pT, mc_y, simZ, weight_y_pp*weight_z);
		}//n_hepmc_p
	}else if ( phpythia_con ){

		for (unsigned int ip=0; ip<phpythia_con->size(); ip++){
			TMCParticle *part = phpythia_con->getParticle(ip);

			if ( part->GetKF()!=_target_pid || part->GetKS()!=11 ) continue;

			float mc_px = part->GetPx();
			float mc_py = part->GetPy();
			float mc_pz = part->GetPz();
			float mc_e = part->GetEnergy();

			float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
			float mc_ptot	= sqrt(mc_pT*mc_pT + mc_pz*mc_pz); 
			float mc_eta	= 0.5*log((mc_ptot+mc_pz)/(mc_ptot-mc_pz));
			float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));

			y_Jpsi = mc_y;

			count++;

			hpT_eta_zvtx->Fill(mc_pT, mc_eta, simZ);
			hpT_y_zvtx->Fill(mc_pT, mc_y, simZ);

			if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 ){
			  if( (y_Jpsi > 1.2) && (y_Jpsi < 2.2))
			    hd_pT_y_zvtx_N->Fill(mc_pT, mc_y, simZ);
			  if( (y_Jpsi > -2.2) && (y_Jpsi < -1.2))
			    hd_pT_y_zvtx_S->Fill(mc_pT, mc_y, simZ);
			}

			float weight_pT_dAu = 1.0;
			float weight_y_dAu = 1.0;
			float weight_pT_pp = 1.0;
			float weight_y_pp = 1.0;
			float weight_z = 1.0;

			weight_z = fw_z->Eval(simZ);

			//if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 ){
			if ( fabs(y_Jpsi)>0.6 && fabs(y_Jpsi)<3.0 ){

				if ( y_Jpsi<0 ){
					weight_y_dAu = fw_y_dAu->Eval(y_Jpsi);
					weight_pT_dAu = (mc_pT<8.0) ? fw_pT_dAu_bwd->Eval(mc_pT) : fw_pT_dAu_bwd->Eval(8.0);
				}else{
					weight_y_dAu = fw_y_dAu->Eval(y_Jpsi);
					weight_pT_dAu = (mc_pT<8.0) ? fw_pT_dAu_fwd->Eval(mc_pT) : fw_pT_dAu_fwd->Eval(8.0);
				}

				weight_y_pp = fw_y_pp->Eval(y_Jpsi);
				weight_pT_pp = (mc_pT<8.0) ? fw_pT_pp->Eval(mc_pT) : fw_pT_pp->Eval(8.0);

				if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 ){
					for (unsigned int ii=0; ii<mu_pT.size(); ii++){
						hmu_pT_y_zvtx->Fill(mu_pT[ii],mu_rap[ii],simZ);
						hmu_pT_y_zvtx_w_dAu->Fill(mu_pT[ii],mu_rap[ii],simZ,weight_pT_dAu*weight_y_dAu*weight_z);
						hmu_pT_y_zvtx_w_pp->Fill(mu_pT[ii],mu_rap[ii],simZ,weight_pT_pp*weight_y_pp*weight_z);
					}

					if ( mu_pT.size()==2 ){
						hmu_y0_y1->Fill(mu_rap[0],mu_rap[1]);
					}
				}
			}

			hpT_y_zvtx_w_z->Fill(mc_pT, mc_y, simZ, weight_z);

			hpT_y_zvtx_w_dAu->Fill(mc_pT, mc_y, simZ, weight_pT_dAu*weight_y_dAu*weight_z);
			hpT_y_zvtx_w_dAu_pT->Fill(mc_pT, mc_y, simZ, weight_pT_dAu*weight_z);
			hpT_y_zvtx_w_dAu_y->Fill(mc_pT, mc_y, simZ, weight_y_dAu*weight_z);

			hpT_y_zvtx_w_pp->Fill(mc_pT, mc_y, simZ, weight_pT_pp*weight_y_pp*weight_z);
			hpT_y_zvtx_w_pp_pT->Fill(mc_pT, mc_y, simZ, weight_pT_pp*weight_z);
			hpT_y_zvtx_w_pp_y->Fill(mc_pT, mc_y, simZ, weight_y_pp*weight_z);

		}//ip

	}

	if ( count>1 ){
		cout << "##################################################" << endl;
		cout << "FIND event containing multiple Jpsi after check: " << count << endl;
		cout << "##################################################" << endl;
	}

	/*
	if ( 0 ){
		//Find mother
		for (int ip=0; ip<n_hepmc_p; ip++){
			MCHepMCParticle *hepmc = hepmc_con->get_MCHepMCParticle(ip);

			int mc_vertex_end = hepmc->get_mc_vertex_end();

			if ( mc_vertex_end!= Jpsi_vertex_begin ) continue;
			int mc_pid 		= hepmc->get_mc_pid();

			if ( (fabs(mc_pid)>500 && fabs(mc_pid)<600) || (fabs(mc_pid)>5000 && fabs(mc_pid)<6000) ){
				float mc_px 	= hepmc->get_mc_px();
				float mc_py 	= hepmc->get_mc_py();
				float mc_pz 	= hepmc->get_mc_pz();

				float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
				//float mc_ptot	= sqrt(mc_pT*mc_pT + mc_pz*mc_pz); 
				float mc_e		= hepmc->get_mc_e();
				float mc_y		= 0.5*log((mc_e+mc_pz)/(mc_e-mc_pz));
				//float mc_eta	= 0.5*log((mc_ptot+mc_pz)/(mc_ptot-mc_pz));

				if ( fabs(y_Jpsi)>1.2 && fabs(y_Jpsi)<2.2 ){
					hp_pT_y_zvtx->Fill(mc_pT, mc_y, simZ);
				}
			}

			//cout << "MC_PID: " << mc_pid << endl;
		}
	}
	*/

	//if ( fabs(eta_Jpsi)<1.2 || fabs(eta_Jpsi)>2.2 ) return EVENT_OK;

	//rapidity distribution
	/*
	for (int ip=0; ip<n_hepmc_p; ip++){

		MCHepMCParticle *hepmc = hepmc_con->get_MCHepMCParticle(ip);

		int mc_pid 		= hepmc->get_mc_pid();
		int mc_status = hepmc->get_mc_status();

		if ( !(mc_status==1) ) continue;
		if ( mc_pid==22 ) continue;

		float mc_px 	= hepmc->get_mc_px();
		float mc_py 	= hepmc->get_mc_py();
		float mc_pz 	= hepmc->get_mc_pz();
		
		float mc_pT		= sqrt(mc_px*mc_px + mc_py*mc_py);
		float mc_ptot	= sqrt(mc_pT*mc_pT + mc_pz*mc_pz); 
		float mc_eta	= 0.5*log((mc_ptot+mc_pz)/(mc_ptot-mc_pz));

		//
		heta_simz->Fill(mc_eta, simZ);

		//if ( (mult_fvtxN+mult_fvtxS+mult_svx)>=3 && (mult_fvtxN+mult_fvtxS+mult_svx)<=5 ){
		if ( mult_svx<0.5 ){
			heta_simz_svx0->Fill(mc_eta, simZ);
			if ( !isnan(fvtxZ) ){
				heta_simz_svx0_zreco->Fill(mc_eta, simZ);
			}
		}//total mult 3-5
	}//n_hepmc_p
	*/

	/*
	if ( mult_svx<0.5 ){
		hsimz_svx0->Fill(simZ);
		if ( !isnan(fvtxZ) ){
			hsimz_svx0_zreco->Fill(simZ);
		}
	}
	*/

	/*
	if ( mult_svx<0.5 && mult_fvtx_prim_cut>=2 && mult_fvtx_prim_cut<=3 ){
		hsimz_svx0_fvtx2345->Fill(simZ);
		if ( !isnan(fvtxZ) ){
			hsimz_svx0_fvtx2345_zreco->Fill(simZ);
		}

		TFvtxCompactTrkMap::iterator iter( fvtxtrk_map->range() );
		while( TFvtxCompactTrkMap::const_pointer fvtx_ptr = iter.next() ) {
			float fvtx_chi2 = fvtx_ptr->get()->get_chi2_ndf();
			if ( fvtx_chi2>4.0 || isnan(fvtx_chi2) ) continue;
			if ( fvtx_ptr->get()->get_nhits()<2.5 ) continue;

			float rapidity = fvtx_ptr->get()->get_fvtx_eta();

			heta_simz_svx0_fvtx2345->Fill(rapidity, simZ);

			if ( !isnan(fvtxZ) ){
				heta_simz_svx0_fvtx2345_zreco->Fill(rapidity, simZ);
			}
		}//FvtxCompactTrk 
	}//multiplicity condition
	*/

	/*
	if ( fabs(simZ)<10.0 ){
			//h2mult_fvtx_svx->Fill(mult_fvtxN+mult_fvtxS, mult_svx);
			//h2mult_fvtx_svx->Fill(mult_fvtx_prim_cut, mult_svx);
			h2mult_fvtx_svx->Fill(mult_fvtx_prim_cut_eta25, mult_svx);
		if ( !isnan(fvtxZ) ){
			//h2mult_fvtx_svx_zreco->Fill(mult_fvtxN+mult_fvtxS, mult_svx);
			h2mult_fvtx_svx_zreco->Fill(mult_fvtx_prim_cut_eta25, mult_svx);
		}
	}
	*/

	//Fill event information
	//hmult_simz_fvtx->Fill(mult_fvtxN+mult_fvtxS, simZ);
	//hmult_simz_all->Fill(mult_fvtxN+mult_fvtxS+mult_svx, simZ);
	//hmult_simz_svx->Fill(mult_svx, simZ);

	//hmult_cut_simz_fvtx->Fill(mult_fvtx_prim_cut_eta25, simZ);
	//hmult_cut_simz_all->Fill(mult_fvtx_prim_cut_eta25+mult_svx_prim_cut, simZ);
	//hmult_cut_simz_svx->Fill(mult_svx_prim_cut, simZ);

	/*
	if ( isnan(fvtxZ) ){

	}else{

		if ( fabs(simZ)<10 ){
			hfvtx_zerr_res_svx_mult->Fill(fvtxZ_err, fvtxZ-simZ, mult_svx);
			hfvtx_xerr_res_svx_mult->Fill(fvtxX_err, fvtxX-simX, mult_svx);
			hfvtx_yerr_res_svx_mult->Fill(fvtxY_err, fvtxY-simY, mult_svx);

			hfvtx_zerr_res_all_mult->Fill(fvtxZ_err, fvtxZ-simZ, mult_svx+mult_fvtx_prim_cut_eta25);
			hfvtx_xerr_res_all_mult->Fill(fvtxX_err, fvtxX-simX, mult_svx+mult_fvtx_prim_cut_eta25);
			hfvtx_yerr_res_all_mult->Fill(fvtxY_err, fvtxY-simY, mult_svx+mult_fvtx_prim_cut_eta25);

			if ( mult_svx>0 ){
				hfvtx_zerr_res->Fill(fabs(fvtxZ-simZ), fvtxZ_err);
			}
		}

		hmult_simz_zreco_fvtx->Fill(mult_fvtxN+mult_fvtxS, simZ);
		hmult_simz_zreco_all->Fill(mult_fvtxN+mult_fvtxS+mult_svx, simZ);
		hmult_simz_zreco_svx->Fill(mult_svx, simZ);

		hmult_cut_simz_zreco_fvtx->Fill(mult_fvtx_prim_cut_eta25, simZ);
		hmult_cut_simz_zreco_all->Fill(mult_fvtx_prim_cut_eta25+mult_svx_prim_cut, simZ);
		hmult_cut_simz_zreco_svx->Fill(mult_svx_prim_cut, simZ);

	}//isnan(fvtxZ)
	*/


	return EVENT_OK;
}

//______________________________________________________
int mFillAnaTree::End(PHCompositeNode *top_node)
{
  MUTOO::PRINT(std::cout,"mFillAnaTree::End");

	cout << "Number of processed events: " << _nevent << endl;
	cout << "here I am" << endl;
	cout << "Number of skipped events: " << _nevent_skip << endl;
	cout << "Fraction of skipped events: " << _nevent_skip*1.0/_nevent << endl;

	file_out->cd();
	file_out->Write();

	delete file_out;

	MUTOO::PRINT(std::cout,"**");
  return 0;
}

