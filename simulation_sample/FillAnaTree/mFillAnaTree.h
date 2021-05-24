#include <SubsysReco.h>

class VtxOut;
class SingleMuonContainer;
class MCSingleMuonFvtxContainer;
class TTree;
class TH3F;
class TH2F;
class TH1F;
class TFile;
class TF1;
class MCHepMCParticleContainer;
class TFvtxCompactTrkMap;
class SvxSegmentList;
class TrigLvl1;
class PHPythiaContainer;
class DiMuonContainer;
class TrigRunLvl1;
class TriggerEmulatorContainer;
class MutrRefitSingleMuonContainer;
class TriggerEmulatorContainer;

class mFillAnaTree : public SubsysReco
{ 
 public:

  //! default constructor
  mFillAnaTree(const char *name, const char *filename);
  
  //! destructor
  virtual 
  ~mFillAnaTree();

  //! global initialization
  int Init(PHCompositeNode *topNode);
  
  //! Run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event method
  int process_event(PHCompositeNode *topNode);
  
  //! global termination
  int End(PHCompositeNode *topNode);

	void set_sim(int sim_num){ _is_sim = sim_num; }

	void set_upsilon(bool upsilon){ _is_upsilon = upsilon; if (upsilon) _target_pid = 553; }

	void set_phi(bool phi){ _is_phi = phi; if (phi) _target_pid = 333; }

	void set_dataset(const char *dataset){ _dataset = dataset; }

 private:

	VtxOut *vtx;
	SingleMuonContainer *sngmucon;
	MCSingleMuonFvtxContainer *fvtxmc_sngmucon;
	DiMuonContainer *dimucon;
	MCHepMCParticleContainer *hepmc_con;
	TFvtxCompactTrkMap *fvtxtrk_map;
	SvxSegmentList *svxseg_list;
	TrigLvl1 *triglvl1;
	PHPythiaContainer *phpythia_con;
	TrigRunLvl1 *trigrun;
	MutrRefitSingleMuonContainer *sngrefitmucon;
	TriggerEmulatorContainer *trigemcon;

	bool first;
	int _is_sim;
	bool _is_upsilon;
	bool _is_phi;
	int _target_pid;

	int _nevent;
	int _nevent_skip;

	TFile *file_weight0;
	TFile *file_weight1;
	TF1 *fw_pT_dAu_fwd;
	TF1 *fw_pT_dAu_bwd;
	TF1 *fw_y_dAu;
	TF1 *fw_pT_pp;
	TF1 *fw_y_pp;

	TF1 *fw_z;
	TF1 *fw_z_S;
	TF1 *fw_z_N;

	TFile *file_out;
	TTree *ana_tree;

	TH3F *hpT_eta_zvtx;
	TH3F *hpT_y_zvtx;
	TH3F *hpT_y_zvtx_w_z;

	TH3F *hpT_y_zvtx_w_dAu;
	TH3F *hpT_y_zvtx_w_dAu_pT;
	TH3F *hpT_y_zvtx_w_dAu_y;

	TH3F *hpT_y_zvtx_w_pp;
	TH3F *hpT_y_zvtx_w_pp_pT;
	TH3F *hpT_y_zvtx_w_pp_y;

	TH3F *hpT_eta_zvtx_S;
	TH3F *hpT_y_zvtx_S;
	TH3F *hpT_y_zvtx_S_w_z;

	TH3F *hpT_y_zvtx_S_w_dAu;
	TH3F *hpT_y_zvtx_S_w_dAu_pT;
	TH3F *hpT_y_zvtx_S_w_dAu_y;

	TH3F *hpT_y_zvtx_S_w_pp;
	TH3F *hpT_y_zvtx_S_w_pp_pT;
	TH3F *hpT_y_zvtx_S_w_pp_y;

	TH3F *hpT_eta_zvtx_N;
	TH3F *hpT_y_zvtx_N;
	TH3F *hpT_y_zvtx_N_w_z;

	TH3F *hpT_y_zvtx_N_w_dAu;
	TH3F *hpT_y_zvtx_N_w_dAu_pT;
	TH3F *hpT_y_zvtx_N_w_dAu_y;

	TH3F *hpT_y_zvtx_N_w_pp;
	TH3F *hpT_y_zvtx_N_w_pp_pT;
	TH3F *hpT_y_zvtx_N_w_pp_y;

	TH3F *hmu_pT_y_zvtx;
	TH3F *hmu_pT_y_zvtx_w_dAu;
	TH3F *hmu_pT_y_zvtx_w_pp;
	TH2F *hmu_y0_y1;

	TH3F *hdimu_pT_y_zvtx;
	TH3F *hdimu_pT_y_zvtx_cuts;
	TH3F *hdimu_ul_pT_y_zvtx_w_pp;

	TH3F *hdimu_pp_pT_y_zvtx_w_pp;
	TH3F *hdimu_mm_pT_y_zvtx_w_pp;

	TH3F *hp_pT_y_zvtx;
	TH3F *hp_pT_eta_zvtx;
	TH3F *hd_pT_y_zvtx;

	TH3F *hfvtx_zerr_res_svx_mult;
	TH3F *hfvtx_xerr_res_svx_mult;
	TH3F *hfvtx_yerr_res_svx_mult;

	TH3F *hfvtx_zerr_res_all_mult;
	TH3F *hfvtx_xerr_res_all_mult;
	TH3F *hfvtx_yerr_res_all_mult;

	TH2F *hfvtx_zerr_res;

	TH2F *h2mult_zreco;

	TH2F *hmult_simz_fvtx;
	TH2F *hmult_simz_svx;
	TH2F *hmult_simz_all;

	TH2F *hmult_simz_zreco_fvtx;
	TH2F *hmult_simz_zreco_svx;
	TH2F *hmult_simz_zreco_all;

	//cuts in prim. vtx finder
	TH2F *hmult_cut_simz_fvtx;
	TH2F *hmult_cut_simz_svx;
	TH2F *hmult_cut_simz_all;

	TH2F *hmult_cut_simz_zreco_fvtx;
	TH2F *hmult_cut_simz_zreco_svx;
	TH2F *hmult_cut_simz_zreco_all;

	TH2F *heta_simz;

	TH2F *h2mult_fvtx_svx;
	TH2F *h2mult_fvtx_svx_zreco;

	TH2F *heta_simz_svx0;
	TH2F *heta_simz_svx0_zreco;
	TH1F *hsimz_svx0;
	TH1F *hsimz_svx0_zreco;

	TH2F *heta_simz_svx0_fvtx2345;
	TH2F *heta_simz_svx0_fvtx2345_zreco;
	TH1F *hsimz_svx0_fvtx2345;
	TH1F *hsimz_svx0_fvtx2345_zreco;

	///////////////////////////////


	TH3F *hmass_fvtx_below_3pT;
	TH3F *hmass_fvtx_above_3pT;
	TH3F *hmass_mismatch;

	TH3F *hmass_N_fvtx;
	TH3F *hmass_N_fvtx_pp;
	TH3F *hmass_N_fvtx_mm;

	TH3F * hmass_N_fvtx_4084;
	TH3F *hmass_N_fvtx_pp_4084;
	TH3F *hmass_N_fvtx_mm_4084;

	TH3F * hmass_N_fvtx_2040;
	TH3F *hmass_N_fvtx_pp_2040;
	TH3F *hmass_N_fvtx_mm_2040;

	TH3F * hmass_N_fvtx_010;
	TH3F *hmass_N_fvtx_pp_010;
	TH3F *hmass_N_fvtx_mm_010;

	TH3F * hmass_N_fvtx_020;
	TH3F *hmass_N_fvtx_pp_020;
	TH3F *hmass_N_fvtx_mm_020;

	TH3F * hmass_N_fvtx_3050;
	TH3F *hmass_N_fvtx_pp_3050;
	TH3F *hmass_N_fvtx_mm_3050;

	TH3F * hmass_N_fvtx_4060;
	TH3F *hmass_N_fvtx_pp_4060;
	TH3F *hmass_N_fvtx_mm_4060;

	TH3F * hmass_N_fvtx_1030;
	TH3F *hmass_N_fvtx_pp_1030;
	TH3F *hmass_N_fvtx_mm_1030;

	TH3F * hmass_N_fvtx_5084;
	TH3F *hmass_N_fvtx_pp_5084;
	TH3F *hmass_N_fvtx_mm_5084;

	TH3F * hmass_N_fvtx_6084;
	TH3F *hmass_N_fvtx_pp_6084;
	TH3F *hmass_N_fvtx_mm_6084;

	/////////////////////////////////

	TH3F *hmass_S_fvtx;
	TH3F *hmass_S_fvtx_pp;
	TH3F *hmass_S_fvtx_mm;

	TH3F * hmass_S_fvtx_4084;
	TH3F *hmass_S_fvtx_pp_4084;
	TH3F *hmass_S_fvtx_mm_4084;

	TH3F * hmass_S_fvtx_2040;
	TH3F *hmass_S_fvtx_pp_2040;
	TH3F *hmass_S_fvtx_mm_2040;

	TH3F * hmass_S_fvtx_010;
	TH3F *hmass_S_fvtx_pp_010;
	TH3F *hmass_S_fvtx_mm_010;

	TH3F * hmass_S_fvtx_020;
	TH3F *hmass_S_fvtx_pp_020;
	TH3F *hmass_S_fvtx_mm_020;

	TH3F * hmass_S_fvtx_3050;
	TH3F *hmass_S_fvtx_pp_3050;
	TH3F *hmass_S_fvtx_mm_3050;

	TH3F * hmass_S_fvtx_4060;
	TH3F *hmass_S_fvtx_pp_4060;
	TH3F *hmass_S_fvtx_mm_4060;

	TH3F * hmass_S_fvtx_1030;
	TH3F *hmass_S_fvtx_pp_1030;
	TH3F *hmass_S_fvtx_mm_1030;

	TH3F * hmass_S_fvtx_5084;
	TH3F *hmass_S_fvtx_pp_5084;
	TH3F *hmass_S_fvtx_mm_5084;

	TH3F * hmass_S_fvtx_6084;
	TH3F *hmass_S_fvtx_pp_6084;
	TH3F *hmass_S_fvtx_mm_6084;

	TH3F *hd_pT_y_zvtx_N;
	TH3F *hd_pT_y_zvtx_S;

	std::string _filename;
	std::string _dataset;

};
