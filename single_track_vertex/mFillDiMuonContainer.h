
#include <SubsysReco.h>
#include <set>
#include <vector>

#ifndef __CINT__
#include <boost/smart_ptr.hpp>
#include <MWGVertex.h>
#endif
 
class SingleMuon;
class SingleMuonContainer;
class SingleMuonContainer_v16;
class DiMuon;
class DiMuon_v12;
class PHPoint;
class PHCompositeNode;

class mFillDiMuonContainer : public SubsysReco
{
 public:

  //! default constructor
 mFillDiMuonContainer(const bool _make_event_mix=true, const float _z_vertex_error=0.5, const bool _only_same_arm=true)
   :SubsysReco("mFillDiMuonContainer")
    {
      make_event_mix = _make_event_mix;
      only_same_arm = _only_same_arm;
      z_vertex_error = _z_vertex_error;
      mass_cut = 0.0;
      is_pp = false;
      is_sim = false;
      make_dimuon_vertex_origin_widerange = false;
      dimuon_vtxr_cut = 100.0;
      dimuon_vtxchi2_cut = 999.0;
      bbczcut = 30.0;
      vtxerror = PHPoint(1.0, 1.0, z_vertex_error);
      nevents = 0;
      naccepted_N = 0;
      naccepted_S = 0;
      naccepted_N_mix = 0;
      naccepted_S_mix = 0;
      for (int i=0; i<NZ; i++)
	for (int j=0; j<NCENT; j++)
	  for (int k=0; k<NRP; k++)
	    {
	      nbuff[i][j][k] = 0;
	      for (int ii=0; ii<NBUF; ii++)
		for (int jj=0; jj<NBUF; jj++)
		  {
		    used_comb[i][j][k][ii][jj] = false;
		  }
	    }
      _vtx = NULL;
    
			_max_fvtx_vertexes = 10;
			reset_fvtx_vertex_names();
			_vtx_top_node = "";
    }
  
  //! destructor
  virtual ~mFillDiMuonContainer() {_fvtx_vertex_names.clear();}
  
  //! global initialization
  int Init(PHCompositeNode *topNode);
  
  //! Run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event method
  int process_event(PHCompositeNode *topNode);
  
  //! global termination
  int End(PHCompositeNode *topNode);
  
  void set_mass_cut(float a) {mass_cut = a;}
  void set_is_sim(bool a) {is_sim = a;}
  void set_is_pp(bool a) {is_pp = a;}
  void set_make_dimuon_vertex_origin_widerange(bool a) {make_dimuon_vertex_origin_widerange = a;}
  bool get_make_dimuon_vertex_origin_widerange() {return make_dimuon_vertex_origin_widerange;}
  void set_bbcz_cut(float a) {bbczcut = a;}
  void set_dimuon_vtxr_cut(float a) {dimuon_vtxr_cut = a;}
  void set_dimuon_vtxchi2_cut(float a) {dimuon_vtxchi2_cut = a;}
  void set_only_same_arm(bool a) {only_same_arm = a;}

  void set_vtx_top_node(std::string s){_vtx_top_node = s;}
 
	//! follow the naming rule in the FVTXPrimVertex module
	void reset_fvtx_vertex_names()
	{
		_fvtx_vertex_names.resize(_max_fvtx_vertexes);
		for (int i = 0;i<_max_fvtx_vertexes;i++)
		{
			if(i==0) _fvtx_vertex_names[i] = "FVTX";
			else if(i==1) _fvtx_vertex_names[i] = "FVTX_SECOND";
			else _fvtx_vertex_names[i] = Form("FVTX_%i",i+1);
		}
	}

	//!
	void set_fvtx_vertex_names(int i = 0, std::string name = "FVTX")
	{
		if(i<_max_fvtx_vertexes){
			_fvtx_vertex_names[i] = name;
		}else{
			return;
		}
	}

	//!
	const std::string get_fvtx_vertex_names(int i = 0) const
	{
		if(i < _max_fvtx_vertexes){
			return _fvtx_vertex_names[i];
		}else{
			std::string dummy("DUMMY");
			return dummy;
		}
	}

 protected:

  DiMuon_v12 make_dimuon(int imuon1, int imuon2, SingleMuonContainer* muons1, SingleMuonContainer* muons2);
  void fillFirstMuon(DiMuon_v12* dimuon, SingleMuon *muon);
  void fillSecondMuon(DiMuon_v12* dimuon, SingleMuon *muon);
  void calculate_from_fvtx(SingleMuon* muon1, SingleMuon* muon2, DiMuon_v12* dimuon);
  void calculate_from_sngtrk12(SingleMuon* muon1, SingleMuon* muon2, DiMuon_v12* dimuon);
  void calculate_from_sngtrk21(SingleMuon* muon2, SingleMuon* muon1, DiMuon_v12* dimuon);
  void get_dimuon_vertex(SingleMuon* muon1, SingleMuon* muon2, PHPoint& vtx, PHPoint& vtx_error);

 private:

	enum {NZ = 40};
	enum {NCENT = 20};
	enum {NRP = 20};
	enum {NBUF = 4};  // 4

#ifndef __CINT__
  MWGVertex get_vertex(unsigned short imu1, unsigned short imu2,
		       SingleMuonContainer* muon1, SingleMuonContainer* muon2);
  MWGVertex get_vertex_from_fvtx(unsigned short imu1, unsigned short imu2,
				 SingleMuonContainer* muon1, SingleMuonContainer* muon2);
  MWGVertex get_vertex_from_sngtrk12(unsigned short imu1, unsigned short imu2,
                                 SingleMuonContainer* muon1, SingleMuonContainer* muon2);
 MWGVertex get_vertex_from_sngtrk21(unsigned short imu2, unsigned short imu1,
                                 SingleMuonContainer* muon2, SingleMuonContainer* muon1);

  void fill_mom_vtx(MWGVertex &vertex, DiMuon_v12 &dimuon);
  void fill_mom_vtx_fvtx(MWGVertex &vertex_fvtx, DiMuon_v12 &dimuon);
  void fill_mom_vtx_sngtrk(MWGVertex &vertex_fvtx, DiMuon_v12 &dimuon);

#endif

  float mass_cut;
  float z_vertex_error;
  bool only_same_arm;
  
	//!
	int _max_fvtx_vertexes;

	//! naming rule of the FVTX vertices
	std::vector<std::string> _fvtx_vertex_names;

  //! counters
  int nevents;
  int naccepted_N;
  int naccepted_S;
  int naccepted_N_mix;
  int naccepted_S_mix;

#ifndef __CINT__
  typedef boost::shared_ptr<SingleMuonContainer_v16> SingleMuonContainer_ptr;
  SingleMuonContainer_ptr single_queue[NZ][NCENT][NRP][NBUF];
#endif
  
  int nbuff[NZ][NCENT][NRP];
  
  bool used_comb[NZ][NCENT][NRP][NBUF][NBUF];

  float bbczcut;
  unsigned int trig_scale;
  bool is_sim;
  bool is_pp;
  bool make_dimuon_vertex_origin_widerange;
  float dimuon_vtxr_cut;
  float dimuon_vtxchi2_cut;
  bool make_event_mix;
  PHPoint vtxerror;
  bool fvtx_vtx;
  bool sngtrk12_vtx;

  VtxOut* _vtx;
  std::string _vtx_top_node;
  PHCompositeNode *vtxTopNode;
};
