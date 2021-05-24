
#include <Fun4AllServer.h>
#include <exception>
#include <PHGlobal.h>
#include <ReactionPlaneObject.h>
#include <PHMuoTracksOut.h>
#include <PHGeometry.h>
#include <RunHeader.h>
#include <MUTOO.h>
#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <TrigLvl1.h>
#include <VtxOut.h>
#include <MWG.h>
#include <MWGConsts.h>
#include <MWGVersion.h>
#include <Tools.h>
#include <SingleMuon.h>
#include <SingleMuonContainer.h>
#include <SingleMuonContainer_v16.h>
#include <DiMuon.h>
#include <DiMuon_v12.h>
#include <DiMuonContainer.h>
#include <DiMuonContainer_v12.h>
#include <mFillSingleMuonContainer.h>

#include <mFillDiMuonContainer.h>

#define E_CMS 200

using namespace std;
typedef PHIODataNode<PHObject> PHObjectNode_t;

//___________________________________________________________________
int mFillDiMuonContainer::Init(PHCompositeNode *top_node)
{
//create output node in the node tree
PHNodeIterator iter(top_node);
PHCompositeNode *dstNode
= static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode","DST"));

if(!dstNode) {
cout << "dstNode not found"<<endl;
} else {
cout << "dstNode is found"<<endl;
}

DiMuonContainer_v12 *dimuons = new DiMuonContainer_v12();
if(dimuons) {
PHObjectNode_t *muonNode =
  new PHIODataNode<PHObject>(dimuons,"DiMuonContainer","PHObject");
dstNode->addNode(muonNode);
cout << "DiMuonContainer is added" <<endl;
} else {
cout << ThisName << " Init() failed to create output object"<<endl;
return ABORTRUN;
}

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
for (int l=0; l<NBUF; l++)
  single_queue[i][j][k][l] = boost::make_shared<SingleMuonContainer_v16>();
}
  
for (int i=0; i<NZ; i++)
  for (int j=0; j<NCENT; j++)
    for (int k=0; k<NRP; k++)
      for (int ii=0; ii<NBUF; ii++)
	for (int jj=0; jj<NBUF; jj++)
	  {
used_comb[i][j][k][ii][jj] = false;
}

if(!_vtx_top_node.empty())
  {
Fun4AllServer *se = Fun4AllServer::instance();
vtxTopNode = se->topNode(_vtx_top_node.c_str());
cout << "mFillSingleMuonContainer::Init() - Setting vtx top_node to " << _vtx_top_node << endl;
if(!vtxTopNode) return ABORTRUN;
}

return 0;
}

//______________________________________________________
int mFillDiMuonContainer::InitRun(PHCompositeNode *top_node)
{

// Make sure global db parameters are set:
TFvtxGlobalParCntrl::init_run();

MUTOO::PRINT( cout, "mFillDiMuonContainer::InitRun" );

return 0;
}

//______________________________________________________
  int mFillDiMuonContainer::process_event(PHCompositeNode *top_node)
  {

    // ///////////////////////////////////////////////////////////////////////////////////////////// added 1/15/2021
     int Run_Number = 0;
    
    RunHeader* runh = findNode::getClass<RunHeader>(top_node,"RunHeader");
    if (runh)
      Run_Number = runh->get_RunNumber();
    ///////////////////////////////////////////////////////////////////////////////////////////// 
    

    ///////////////////////////////////////////////////////////////////////////////////////////// added 11/15/2020
    TrigLvl1* trig = findNode::getClass<TrigLvl1>(top_node,"TrigLvl1");
    unsigned int trig_scale = trig->get_lvl1_trigscaled();
    //////unsigned int trig_live = trig->get_lvl1_triglive();
    ///////////////////////////////////////////////////////////////////////////////////////////// 

    // for Run15pp (runnumber range: 421815-432008
     if( (Run_Number > 421814) && (Run_Number < 432009) )
      {
    	///////////////////////////////////////////////////////////////////////////////////////////// added 11/15/2020
    	if( ! ((trig_scale&0x00400000) ||  (trig_scale&0x00800000) ||  (trig_scale&0x00100000) ||  (trig_scale&0x00200000)) )
    	  return DISCARDEVENT;
    	///////////////////////////////////////////////////////////////////////////////////////////// 
      }

    // for Run15pAu (runnumber range: 432639-436647)
     if( (Run_Number > 432638) && (Run_Number < 436648) )
      {
    	///////////////////////////////////////////////////////////////////////////////////////////// added 11/15/2020
    	if( ! ((trig_scale&0x00400000) ||  (trig_scale&0x00800000) ||  (trig_scale&0x00100000) ||  (trig_scale&0x00200000)) )
    	  return DISCARDEVENT;
    	///////////////////////////////////////////////////////////////////////////////////////////// 
      }
    // // for Run14HeAu (runnumber range: 415751-416842)
    //  if( (Run_Number > 415750) && (Run_Number < 416843) )
    //   {
    //  	///////////////////////////////////////////////////////////////////////////////////////////// added 1/15/2021
    //  	if( ! ((trig_scale&0x00100000) ||  (trig_scale&0x00200000)  ||  (trig_scale&0x00040000) ||  (trig_scale&0x00080000)) )
    //  	  return DISCARDEVENT;
    //  	/////////////////////////////////////////////////////////////////////////////////////////////
    //   }

  // // for Run15pAl (runnumber range: 436759 - 438422)
     // if( (Run_Number > 436758) && (Run_Number < 438423) )
     //  {
     // 	///////////////////////////////////////////////////////////////////////////////////////////// added 2/18/2021
     // 	if( ! ((trig_scale&0x00400000) ||  (trig_scale&0x00800000) ||  (trig_scale&0x00100000) ||  (trig_scale&0x00200000)) )
     // 	  return DISCARDEVENT;
     // 	/////////////////////////////////////////////////////////////////////////////////////////////
     //  }

    // cout << "scaled trig1D N: " << (trig_scale&0x00400000) << ", trig1D S: " << (trig_scale&0x00800000)  << ", trig 2D N : " << (trig_scale&0x00100000) << ", trig 2D S: " <<   (trig_scale&0x00200000) << endl;
    //cout << "lve trig1D N: " << (trig_live&0x00400000) << ", trig1D S: " << (trig_live&0x00800000)  << ", trig 2D N : " << (trig_live&0x00100000) << ", trig 2D S: " <<   (trig_live&0x00200000) << endl;

    if(nevents%100000==0){
      cout << "Event: " << nevents << endl;
    }
    nevents++;
    //cout << "nevents: " << nevents << endl;   ///// added Nov. 16th

    vtxerror = PHPoint(1.0, 1.0, z_vertex_error);
    // VTX
    if(_vtx_top_node.empty()) _vtx = findNode::getClass<VtxOut>(top_node, "VtxOut" );
    else _vtx = findNode::getClass<VtxOut>(vtxTopNode, "VtxOut" );
    if(!_vtx){ cout << "mFillDiMuonContainer - Could not find VtxOut on " << vtxTopNode->getName() << endl; return ABORTRUN;}

    if (_vtx)
      {
	if (fabs(_vtx->get_Vertex(_fvtx_vertex_names[0].data()).getZ())<200)
	  vtxerror = _vtx->get_VertexError(_fvtx_vertex_names[0].data());
	else if (fabs(_vtx->get_Vertex("SVX_PRECISE").getZ())<200)
	  vtxerror = _vtx->get_VertexError("SVX_PRECISE");
	else 
	  vtxerror = _vtx->get_VertexError("BBC");
      }
    else
      {
	if (nevents<2)
	  cout << PHWHERE << "mFillSingleMuonContainer:: VtxOut not in Node Tree" << endl;
      }
  
    // new framework MWG tracks
    SingleMuonContainer* muons = findNode::getClass<SingleMuonContainer>(top_node,"SingleMuonContainer");
    if (!muons)
      muons = findNode::getClass<SingleMuonContainer>(top_node,"SingleMuonContainer_filtered");
    if (!muons)
      {
	cout << "mFillDiMuonContainer:: SingleMuonContainer not in Node Tree" << endl;
	return ABORTRUN;
      }

    if (verbosity)
      cout << muons->get_Evt_fvtxX() << " " << muons->get_Evt_fvtxY() << " " << muons->get_Evt_fvtxZ() << endl;

    int npart = muons->get_nSingleMuons();

      
    //  if ( fabs(muons->get_Evt_vtxZ()) > bbczcut && fabs(muons->get_Evt_bbcZ())>bbczcut) return DISCARDEVENT;
    //  if ( fabs(muons->get_Evt_vtxZ()) > bbczcut) return DISCARDEVENT;
    if ( fabs(muons->get_Evt_bbcZ()) > bbczcut) return DISCARDEVENT;

    if ( !is_sim && !is_pp )
      {
	if (muons->get_Evt_Cent() < 0 || muons->get_Evt_Cent() > 98) return DISCARDEVENT;
	//      if (muons->get_rx_NS() < -pi/2 || muons->get_rx_NS() > pi/2) return DISCARDEVENT;
      }

    float rx_NS = -1000.;
    ReactionPlaneObject* rp = findNode::getClass<ReactionPlaneObject>(top_node,"ReactionPlaneObject");
    if (rp)
      {
	// Check whether rp is implemented to avoid invalid error messages
	const TString v = rp->ClassName();
	if ((v == TString("ReactionPlaneObject"))
	    || (v == TString("ReactionPlaneObjectv4")))
	  {
	    static bool once = true;

	    if (once)
	      {

		cout
		  << "mFillDiMuonContainer::process_event - WARNING - Incompatible version of ReactionPlaneObject ("
		  << v << ") received from the TOP node. Ignore it." << endl;

		once = false;
	      }

	  }
	else
	  {
	    rx_NS = rp->getRXNrp18();
	  }
      }

    DiMuonContainer_v12 *dimuons = findNode::getClass<DiMuonContainer_v12>(top_node, "DiMuonContainer" );
    if (!dimuons)
      {
	cout << "mFillDiMuonContainer:: DiMuonContainer not in Node Tree" << endl;
	return ABORTRUN;
      }
    dimuons->Reset();
    if (!is_pp && !is_sim)
      dimuons->set_Evt_Cent(-999.9);

    //==============================================================

    // Condition block for mixed events
    if (make_event_mix)
      {
	//__________________________________________________
	// fill single muon pool for mixed event
	int i_dz = (int)(NZ*(muons->get_Evt_bbcZ() + bbczcut)/(2*bbczcut)); // calcualtes bin number for BBCZ
	int i_cent = (int)(NCENT*muons->get_Evt_Cent()/100); // finds correct bin number
	int i_rp = (int)(NRP*(rx_NS + M_PI/2)/M_PI);  // finds correct bin number
      
	if ( is_sim || is_pp)
	  {
	    i_cent = NCENT/2;
	    i_rp = NRP/2;
	  }  

	if (i_dz<0 || i_dz>=NZ)
	  {
	    cout << "i_dz=" << i_dz << " for BBCZ=" << muons->get_Evt_bbcZ() << endl;
	    return DISCARDEVENT;
	  }
	if (i_cent<0 || i_cent>=NCENT)
	  {
	    cout << "i_cent=" << i_cent << " for Centrality=" << muons->get_Evt_Cent() << endl;
	    return DISCARDEVENT;
	  }
	if (i_rp<0 || i_rp>=NRP)
	  {
	    //      cout << "i_rp=" << i_rp << " for RP=" << rx_NS << endl;
	    //      return DISCARDEVENT;
	    i_rp = 0;
	  }

	if (nbuff[i_dz][i_cent][i_rp]>=NBUF)
	  nbuff[i_dz][i_cent][i_rp]=0;
	int nbuff_ = nbuff[i_dz][i_cent][i_rp];
	single_queue[i_dz][i_cent][i_rp][nbuff_]->Reset();
	for (int ii=0; ii<NBUF; ii++)
	  for (int jj=0; jj<NBUF; jj++)
	    if (ii==nbuff_ || jj==nbuff_)
	      used_comb[i_dz][i_cent][i_rp][ii][jj] = false;
  
	single_queue[i_dz][i_cent][i_rp][nbuff_]->copyfrom(muons);
	nbuff[i_dz][i_cent][i_rp]++;

	//__________________________________________________
	// Loop over single dimuons to fill container with
	// mixed events
	for (int ievent=0; ievent<NBUF; ievent++)
	  for (int jevent=ievent+1; jevent<NBUF; jevent++)
	    {
	      if (used_comb[i_dz][i_cent][i_rp][ievent][jevent])
		continue;
	      else
		used_comb[i_dz][i_cent][i_rp][ievent][jevent] = true;

		
	      SingleMuonContainer_ptr muons1 = single_queue[i_dz][i_cent][i_rp][ievent];
	      SingleMuonContainer_ptr muons2 = single_queue[i_dz][i_cent][i_rp][jevent];

	      int nsingles1 = muons1->get_nSingleMuons();
	      int nsingles2 = muons2->get_nSingleMuons();
	   
	      if (nsingles1==0 || nsingles2==0) continue;
	    
	
	      dimuons->set_Evt_vtxX(muons1->get_Evt_vtxX());
	      dimuons->set_Evt_vtxY(muons1->get_Evt_vtxY());
	      dimuons->set_Evt_vtxZ(muons1->get_Evt_vtxZ());
	      //VTX vertex error.
	      dimuons->set_Evt_vtxX_Err(muons1->get_Evt_vtxX_Err());
	      dimuons->set_Evt_vtxY_Err(muons1->get_Evt_vtxY_Err());
	      dimuons->set_Evt_vtxZ_Err(muons1->get_Evt_vtxZ_Err());
	    
	      dimuons->set_Evt_fvtxX(muons1->get_Evt_fvtxX());
	      dimuons->set_Evt_fvtxY(muons1->get_Evt_fvtxY());
	      dimuons->set_Evt_fvtxZ(muons1->get_Evt_fvtxZ());

	      //FVTX vertex error.                                
	      dimuons->set_Evt_fvtxX_Err(muons1->get_Evt_fvtxX_Err());
	      dimuons->set_Evt_fvtxY_Err(muons1->get_Evt_fvtxY_Err());
	      dimuons->set_Evt_fvtxZ_Err(muons1->get_Evt_fvtxZ_Err());

	      // Fill x,y with beam-average position with slope correction if it was used in reconstruction:
	      if (TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy")){
		dimuons->set_Evt_fvtxX(TFvtxGlobalParCntrl::get_float_par("beam_x_seed") + dimuons->get_Evt_fvtxZ() * TFvtxGlobalParCntrl::get_float_par("beam_dxdz"));
		dimuons->set_Evt_fvtxY(TFvtxGlobalParCntrl::get_float_par("beam_y_seed") + dimuons->get_Evt_fvtxZ() * TFvtxGlobalParCntrl::get_float_par("beam_dydz"));
	      }
	    
	      dimuons->set_Evt_fvtxX2(muons1->get_Evt_fvtxX2());
	      dimuons->set_Evt_fvtxY2(muons1->get_Evt_fvtxY2());
	      dimuons->set_Evt_fvtxZ2(muons1->get_Evt_fvtxZ2());
	    
	      dimuons->set_Evt_bbcZ(muons1->get_Evt_bbcZ());
	      dimuons->set_Evt_Cent(muons1->get_Evt_Cent());

	      dimuons->set_Evt_Mult_FVTXN(muons1->get_Evt_Mult_FVTXN());
	      dimuons->set_Evt_Mult_FVTXS(muons1->get_Evt_Mult_FVTXS());
	      dimuons->set_Evt_Mult_sFVTXN(muons1->get_Evt_Mult_sFVTXN());
	      dimuons->set_Evt_Mult_sFVTXS(muons1->get_Evt_Mult_sFVTXS());

	      if ( is_sim || is_pp)
		dimuons->set_Evt_Cent(50.0);
	    
	      for (int iimu=0; iimu<nsingles1; iimu++)
		for (int jjmu=0; jjmu<nsingles2; jjmu++)  
		  {
		      
		    {
		      DiMuon_v12 dimuon = make_dimuon(iimu, jjmu, muons1.get(), muons2.get());
		      if (dimuon.mass < mass_cut) continue;
		      if (dimuon.Evt_vtxchi2 > dimuon_vtxchi2_cut ) continue;
		      float dimuon_vtxr = sqrt(dimuon.X0*dimuon.X0 + dimuon.Y0*dimuon.Y0);
		      if (dimuon_vtxr > dimuon_vtxr_cut ) continue;
			
		      dimuon.same_event = false;
		      dimuons->AddDiMuon(dimuon);
		      if (dimuon.Pz > 0) naccepted_N_mix++;
		      if (dimuon.Pz < 0) naccepted_S_mix++;
			
		   	
		    }
		  }
	    }
      } // End of making mixed events
      
    //==============================================================
      
    if (npart<2) return DISCARDEVENT;
 
    // same event
  dimuons->set_Evt_vtxX(muons->get_Evt_vtxX());
  dimuons->set_Evt_vtxY(muons->get_Evt_vtxY());
  dimuons->set_Evt_vtxZ(muons->get_Evt_vtxZ());

  //VTX vertex error.                                                 
  dimuons->set_Evt_vtxX_Err(muons->get_Evt_vtxX_Err());
  dimuons->set_Evt_vtxY_Err(muons->get_Evt_vtxY_Err());
  dimuons->set_Evt_vtxZ_Err(muons->get_Evt_vtxZ_Err());

  dimuons->set_Evt_fvtxX(muons->get_Evt_fvtxX());
  dimuons->set_Evt_fvtxY(muons->get_Evt_fvtxY());
  dimuons->set_Evt_fvtxZ(muons->get_Evt_fvtxZ());

  //FVTX vertex error.                                                 
  dimuons->set_Evt_fvtxX_Err(muons->get_Evt_fvtxX_Err());
  dimuons->set_Evt_fvtxY_Err(muons->get_Evt_fvtxY_Err());
  dimuons->set_Evt_fvtxZ_Err(muons->get_Evt_fvtxZ_Err());

  // Fill x,y with beam-average position with slope correction if it was used in reconstruction:
  if (TFvtxGlobalParCntrl::get_bool_par("beam_use_average_xy")){
    dimuons->set_Evt_fvtxX(TFvtxGlobalParCntrl::get_float_par("beam_x_seed") + dimuons->get_Evt_fvtxZ() * TFvtxGlobalParCntrl::get_float_par("beam_dxdz"));
    dimuons->set_Evt_fvtxY(TFvtxGlobalParCntrl::get_float_par("beam_y_seed") + dimuons->get_Evt_fvtxZ() * TFvtxGlobalParCntrl::get_float_par("beam_dydz"));
  }

  dimuons->set_Evt_fvtxX2(muons->get_Evt_fvtxX2());
  dimuons->set_Evt_fvtxY2(muons->get_Evt_fvtxY2());
  dimuons->set_Evt_fvtxZ2(muons->get_Evt_fvtxZ2());

  dimuons->set_Evt_bbcZ(muons->get_Evt_bbcZ());
  dimuons->set_Evt_Cent(muons->get_Evt_Cent());
  if ( is_sim || is_pp)
    dimuons->set_Evt_Cent(50.0);

  for (int imu1=0; imu1<npart; imu1++)
    for (int imu2=imu1+1; imu2<npart; imu2++)
      {
	DiMuon_v12 dimuon = make_dimuon(imu1, imu2, muons, muons);
	float dimuon_vtxr = sqrt(dimuon.X0*dimuon.X0 + dimuon.Y0*dimuon.Y0);
	if (dimuon.get_mass() < mass_cut) continue;
	if (dimuon.Evt_vtxchi2 > dimuon_vtxchi2_cut ) continue;
	if (dimuon_vtxr > dimuon_vtxr_cut ) continue;
	dimuon.same_event = true;
	dimuons->AddDiMuon(dimuon);
	if (dimuon.Pz > 0) naccepted_N++;
	if (dimuon.Pz < 0) naccepted_S++;
      }

  if (dimuons->get_nDiMuons() > 0) return EVENT_OK;

  return DISCARDEVENT;
 }

/////////////////////////////////////// this is part that needs to change
DiMuon_v12 mFillDiMuonContainer::make_dimuon(int imuon1, int imuon2,
					     SingleMuonContainer* muons1,
					     SingleMuonContainer* muons2)
{
  SingleMuon* muon1 = muons1->get_SingleMuon(imuon1);
  SingleMuon* muon2 = muons2->get_SingleMuon(imuon2);
  
  DiMuon_v12 dimuon;
  dimuon.Reset();

  if (only_same_arm && (muon1->get_pz()/fabs(muon1->get_pz()) !=  muon2->get_pz()/fabs(muon2->get_pz())))
    return dimuon;

  MWGVertex vertex;
  try
    {
      vertex = MWGVertex(get_vertex(imuon1, imuon2, muons1, muons2));
    }
  catch(exception& e)
    {
      cout << PHWHERE << e.what() << endl;
      vertex.print();
      exit(0);
    }

  if (vertex.get_mass() < mass_cut)
    return dimuon;
  dimuon.mass = vertex.get_mass();
  fill_mom_vtx(vertex,dimuon);
  int charge1 = (muon1->get_charge()>0) ? 1 : -1;
  int charge2 = (muon2->get_charge()>0) ? 1 : -1;
  dimuon.charge = charge1 + charge2;


  float fvtx_prob1 = 0;
  float fvtx_prob2 = 0;

  if (muon1->get_nhits_fvtx()>0){
    int ndf = 2*muon1->get_nhits_fvtx() - 1;
    float chi2 = muon1->get_chi2_fvtx();
    fvtx_prob1 = TMath::Prob(chi2*ndf, ndf);
  }

  if (muon2->get_nhits_fvtx()>0){
    int ndf = 2*muon2->get_nhits_fvtx() - 1;
    float chi2 = muon2->get_chi2_fvtx();
    fvtx_prob2 = TMath::Prob(chi2*ndf, ndf);
  }
  
 
  if(  (fabs(muon1->get_dphi_fvtx()) < 10 ) && (fabs(muon2->get_dphi_fvtx()) > 10) )  // muon1 is fvtx and muon2 is not fvtx
    {
      //if (muon1->get_nhits_fvtx()>0 && muon2->get_nhits_fvtx()==0)
      if( (fvtx_prob1>0.05) && (fvtx_prob2<0.05) )  
      {
	calculate_from_sngtrk12(muon1, muon2, &dimuon);

	//if (dimuon.get_mass_sngtrk() > mass_cut ) 
	if (dimuon.get_mass_fvtx() > mass_cut ) 
	  {
	    MWGVertex vertex_sngtrk;
	    try
	      {
		vertex_sngtrk = MWGVertex(get_vertex_from_sngtrk12(imuon1, imuon2, muons1, muons2));
		//dimuon.mass_sngtrk = vertex_sngtrk.get_mass();   
		dimuon.mass_fvtxmutr = vertex_sngtrk.get_mass();   
		fill_mom_vtx_sngtrk(vertex_sngtrk,dimuon);
		sngtrk12_vtx = true;
	      }
	    catch(exception& e)
	      {
		cout << PHWHERE << e.what() << endl;
		//              exit(0);
	      }
	  }
      }
    }
     
  /////////////////////////////////////////////
  if(  (fabs(muon2->get_dphi_fvtx()) < 10 ) && (fabs(muon1->get_dphi_fvtx()) > 10) )  // muon2 is fvtx and muon1 is not fvtx
    {
      //  if (muon1->get_nhits_fvtx()==0 && muon2->get_nhits_fvtx()>0)
      if( (fvtx_prob2>0.05) && (fvtx_prob1<0.05) )  
      {
	calculate_from_sngtrk21(muon2, muon1, &dimuon);

	//if (dimuon.get_mass_sngtrk() > mass_cut ) 
	if (dimuon.get_mass_fvtx() > mass_cut ) 
	  {
	    MWGVertex vertex_sngtrk;
	    try
	      {
		vertex_sngtrk = MWGVertex(get_vertex_from_sngtrk21(imuon2, imuon1, muons2, muons1));
		//dimuon.mass_sngtrk = vertex_sngtrk.get_mass();
		dimuon.mass_fvtxmutr = vertex_sngtrk.get_mass();
		fill_mom_vtx_sngtrk(vertex_sngtrk,dimuon);
	      }
	    catch(exception& e)
	      {
		cout << PHWHERE << e.what() << endl;
		//              exit(0);
	      }
	  }
      }
    }
  ///////////////////////////////////////////////

  if (fabs(muon1->get_dphi_fvtx())<10 && fabs(muon2->get_dphi_fvtx())<10)  // fvtx+fvtx
    //if (muon1->get_nhits_fvtx()>0 && muon2->get_nhits_fvtx()>0)
    {
      if( (fvtx_prob1>0.05) && (fvtx_prob2>0.05) )
      {
	calculate_from_fvtx(muon1, muon2, &dimuon);

	if (dimuon.get_mass_fvtx() > mass_cut ) 
	  {
	    MWGVertex vertex_fvtx;
	    try
	      {
		vertex_fvtx = MWGVertex(get_vertex_from_fvtx(imuon1, imuon2, muons1, muons2)); 
		dimuon.mass_fvtxmutr = vertex_fvtx.get_mass();  
		fill_mom_vtx_fvtx(vertex_fvtx,dimuon);
		fvtx_vtx = true;
	      }
	    catch(exception& e)
	      {
		cout << PHWHERE << e.what() << endl;
		//	      exit(0);
	      }
	  }
      }
    }
     

  PHPoint vtx, vtx_error;

  if ( _vtx )
    get_dimuon_vertex(muon1, muon2, vtx, vtx_error);
  else if ( fvtx_vtx || sngtrk12_vtx )
    {
      vtx = PHPoint(muons1->get_Evt_fvtxX(), muons1->get_Evt_fvtxY(), muons1->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  else
    {
      vtx = PHPoint(muons2->get_Evt_fvtxX(), muons2->get_Evt_fvtxY(), muons2->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }

  dimuon.Evt_vtxoor = sqrt(pow(muon1->get_x0()-muon2->get_x0(),2)+ 
			   pow(muon1->get_y0()-muon2->get_y0(),2));

  float P0[3] = {0,0,0};
  float P1[3] = {0,0,0};
  float beam[3] = {0,0,E_CMS/2};
  
  if (muon1->get_charge()>0 || dimuon.charge<0)
    {
      fillFirstMuon(&dimuon, muon1);
      fillSecondMuon(&dimuon, muon2);
      P0[0] = vertex.get_px(0);
      P0[1] = vertex.get_py(0);
      P0[2] = vertex.get_pz(0);
      P1[0] = vertex.get_px(1);
      P1[1] = vertex.get_py(1);
      P1[2] = vertex.get_pz(1);
    }
  else
    {
      fillFirstMuon(&dimuon, muon2);
      fillSecondMuon(&dimuon, muon1);	  
      P0[0] = vertex.get_px(1);
      P0[1] = vertex.get_py(1);
      P0[2] = vertex.get_pz(1);
      P1[0] = vertex.get_px(0);
      P1[1] = vertex.get_py(0);
      P1[2] = vertex.get_pz(0);
    }
  
  dimuon.costhCS = Tools::costhetaCS(MU_MASS, &P0[0], &P1[0] );
  dimuon.costhGJ = Tools::costhetaGJ(MU_MASS, &P0[0], &P1[0], &beam[0]);
  dimuon.costhHX = Tools::costhetaHX(MU_MASS, &P0[0], &P1[0] );
  dimuon.phiCS = Tools::phiCS(MU_MASS, &P0[0], &P1[0], &beam[0] );
  dimuon.phiGJ = Tools::phiGJ(MU_MASS, &P0[0], &P1[0], &beam[0] );
  dimuon.phiHX = Tools::phiHX(MU_MASS, &P0[0], &P1[0], &beam[0] );
  
  return dimuon;
}


//______________________________________________________
int mFillDiMuonContainer::End(PHCompositeNode *top_node)
{

  MUTOO::PRINT( cout, "mFillDiMuonContainer::End" );

  cout << "events: " << nevents << endl;
  cout << "North Arm accepted dimuons:  " << naccepted_N << endl;
  cout << "South Arm accepted dimuons:  " << naccepted_S << endl;
  cout << "North Arm accepted dimuons in mixed events:  " << naccepted_N_mix << endl;
  cout << "South Arm accepted dimuons in mixed events:  " << naccepted_S_mix << endl;

  return 0 ;
}


void mFillDiMuonContainer::fillFirstMuon(DiMuon_v12* dimuon, SingleMuon *muon)
{
  dimuon->Tr0_x0 = muon->get_x0();
  dimuon->Tr0_y0 = muon->get_y0();
  
  dimuon->Tr0_px = muon->get_px();
  dimuon->Tr0_py = muon->get_py();
  dimuon->Tr0_pz = muon->get_pz();

  dimuon->Tr0_rapidity = muon->get_rapidity();
  dimuon->Tr0_trhits   = muon->get_trhits();
  dimuon->Tr0_idhits   = muon->get_idhits();
  dimuon->Tr0_idpanel  = muon->get_idpanel();
  dimuon->Tr0_DG0      = muon->get_DG0();
  dimuon->Tr0_DDG0     = muon->get_DDG0();
  dimuon->Tr0_DS3      = muon->get_DS3();
  dimuon->Tr0_trchi2   = muon->get_trchi2();
  dimuon->Tr0_idchi2   = muon->get_idchi2();
  dimuon->Tr0_ntrhits  = muon->get_ntrhits();
  dimuon->Tr0_nidhits  = muon->get_nidhits();
  dimuon->Tr0_lastgap  = muon->get_lastgap();
  dimuon->Tr0_xst1     = muon->get_xst1();
  dimuon->Tr0_yst1     = muon->get_yst1();
  dimuon->Tr0_xst2     = muon->get_xst2();
  dimuon->Tr0_yst2     = muon->get_yst2();
  dimuon->Tr0_xst3     = muon->get_xst3();
  dimuon->Tr0_yst3     = muon->get_yst3();
  dimuon->Tr0_idx      = muon->get_idx();
  dimuon->Tr0_idy      = muon->get_idy();
  dimuon->Tr0_st1px    = muon->get_st1px();
  dimuon->Tr0_st1py    = muon->get_st1py();
  dimuon->Tr0_st1pz    = muon->get_st1pz();
  dimuon->Tr0_MUID1D   = muon->get_MUID1D();
  dimuon->Tr0_MUID1S   = muon->get_MUID1S();
  dimuon->Tr0_dca_z    = muon->get_dca_z();
  dimuon->Tr0_dca_r    = muon->get_dca_r();
  dimuon->Tr0_dca_phi    = muon->get_dca_phi();

  dimuon->Tr0_px_fvtx = muon->get_px_fvtx();
  dimuon->Tr0_py_fvtx = muon->get_py_fvtx();
  dimuon->Tr0_pz_fvtx = muon->get_pz_fvtx();
  dimuon->Tr0_x0_fvtx = muon->get_x0_fvtx();
  dimuon->Tr0_y0_fvtx = muon->get_y0_fvtx();
  dimuon->Tr0_z0_fvtx = muon->get_z0_fvtx();
  dimuon->Tr0_px_fvtxmutr = muon->get_px_fvtxmutr();
  dimuon->Tr0_py_fvtxmutr = muon->get_py_fvtxmutr();
  dimuon->Tr0_pz_fvtxmutr = muon->get_pz_fvtxmutr();
  dimuon->Tr0_x0_fvtxmutr = muon->get_x0_fvtxmutr();
  dimuon->Tr0_y0_fvtxmutr = muon->get_y0_fvtxmutr();
  dimuon->Tr0_z0_fvtxmutr = muon->get_z0_fvtxmutr();
  dimuon->Tr0_nhits_fvtx = muon->get_nhits_fvtx();
  dimuon->Tr0_dphi_fvtx = muon->get_dphi_fvtx();
  dimuon->Tr0_dtheta_fvtx = muon->get_dtheta_fvtx();
  dimuon->Tr0_dr_fvtx = muon->get_dr_fvtx();
  dimuon->Tr0_chi2_fvtx = muon->get_chi2_fvtx();
  dimuon->Tr0_chi2_fvtxmutr = muon->get_chi2_fvtxmutr();
  for (int i=0; i<4; i++)
    {
      dimuon->Tr0_fvtx_strip[i] = muon->get_fvtx_strip(i);
    }
  dimuon->Tr0_nfvtx_tracklets_cone = muon->get_nfvtx_tracklets_cone();
  dimuon->Tr0_nfvtx_clusters_cone = muon->get_nfvtx_clusters_cone();
  dimuon->Tr0_clusters_size1 = muon->get_clusters_size1();

  dimuon->Tr0_n_matchings = muon->get_n_matchings();
  dimuon->Tr0_vtx_index = muon->get_vtx_index();
  dimuon->Tr0_hit_pattern = muon->get_hit_pattern();
  dimuon->Tr0_maxres_sigma = muon->get_maxres_sigma();
  dimuon->Tr0_track_id = muon->get_track_id();
  dimuon->Tr0_is_sfvtx = muon->get_is_sfvtx();
  dimuon->Tr0_best_fvtxmutr_match = muon->get_best_fvtxmutr_match();
  dimuon->Tr0_x_fvtxproj = muon->get_x_fvtxproj(0);
  dimuon->Tr0_y_fvtxproj = muon->get_y_fvtxproj(0);
}

void mFillDiMuonContainer::fillSecondMuon(DiMuon_v12* dimuon, SingleMuon *muon)
{
  dimuon->Tr1_x0 = muon->get_x0();
  dimuon->Tr1_y0 = muon->get_y0();

  dimuon->Tr1_px = muon->get_px();
  dimuon->Tr1_py = muon->get_py();
  dimuon->Tr1_pz = muon->get_pz();
  
  dimuon->Tr1_rapidity = muon->get_rapidity();
  dimuon->Tr1_trhits   = muon->get_trhits();
  dimuon->Tr1_idhits   = muon->get_idhits();
  dimuon->Tr1_idpanel  = muon->get_idpanel();  
  dimuon->Tr1_DG0      = muon->get_DG0();
  dimuon->Tr1_DDG0     = muon->get_DDG0();
  dimuon->Tr1_DS3      = muon->get_DS3();
  dimuon->Tr1_trchi2   = muon->get_trchi2();
  dimuon->Tr1_idchi2   = muon->get_idchi2();
  dimuon->Tr1_ntrhits  = muon->get_ntrhits();
  dimuon->Tr1_nidhits  = muon->get_nidhits();
  dimuon->Tr1_lastgap  = muon->get_lastgap();
  dimuon->Tr1_xst1     = muon->get_xst1();
  dimuon->Tr1_yst1     = muon->get_yst1();
  dimuon->Tr1_xst2     = muon->get_xst2();
  dimuon->Tr1_yst2     = muon->get_yst2();
  dimuon->Tr1_xst3     = muon->get_xst3();
  dimuon->Tr1_yst3     = muon->get_yst3();
  dimuon->Tr1_idx      = muon->get_idx();
  dimuon->Tr1_idy      = muon->get_idy();
  dimuon->Tr1_st1px    = muon->get_st1px();
  dimuon->Tr1_st1py    = muon->get_st1py();
  dimuon->Tr1_st1pz    = muon->get_st1pz();
  dimuon->Tr1_MUID1D   = muon->get_MUID1D();
  dimuon->Tr1_MUID1S   = muon->get_MUID1S();
  dimuon->Tr1_dca_z    = muon->get_dca_z();
  dimuon->Tr1_dca_r    = muon->get_dca_r();
  dimuon->Tr1_dca_phi  = muon->get_dca_phi();

  dimuon->Tr1_px_fvtx = muon->get_px_fvtx();
  dimuon->Tr1_py_fvtx = muon->get_py_fvtx();
  dimuon->Tr1_pz_fvtx = muon->get_pz_fvtx();
  dimuon->Tr1_x0_fvtx = muon->get_x0_fvtx();
  dimuon->Tr1_y0_fvtx = muon->get_y0_fvtx();
  dimuon->Tr1_z0_fvtx = muon->get_z0_fvtx();
  dimuon->Tr1_px_fvtxmutr = muon->get_px_fvtxmutr();
  dimuon->Tr1_py_fvtxmutr = muon->get_py_fvtxmutr();
  dimuon->Tr1_pz_fvtxmutr = muon->get_pz_fvtxmutr();
  dimuon->Tr1_x0_fvtxmutr = muon->get_x0_fvtxmutr();
  dimuon->Tr1_y0_fvtxmutr = muon->get_y0_fvtxmutr();
  dimuon->Tr1_z0_fvtxmutr = muon->get_z0_fvtxmutr();
  dimuon->Tr1_nhits_fvtx = muon->get_nhits_fvtx();
  dimuon->Tr1_dphi_fvtx = muon->get_dphi_fvtx();
  dimuon->Tr1_dtheta_fvtx = muon->get_dtheta_fvtx();
  dimuon->Tr1_dr_fvtx = muon->get_dr_fvtx();
  dimuon->Tr1_chi2_fvtx = muon->get_chi2_fvtx();
  dimuon->Tr1_chi2_fvtxmutr = muon->get_chi2_fvtxmutr();
  for (int i=0; i<4; i++)
    {
      dimuon->Tr1_fvtx_strip[i] = muon->get_fvtx_strip(i);
    }
  dimuon->Tr1_nfvtx_tracklets_cone = muon->get_nfvtx_tracklets_cone();
  dimuon->Tr1_nfvtx_clusters_cone = muon->get_nfvtx_clusters_cone();
  dimuon->Tr1_clusters_size1 = muon->get_clusters_size1();

  dimuon->Tr1_n_matchings = muon->get_n_matchings();
  dimuon->Tr1_vtx_index = muon->get_vtx_index();
  dimuon->Tr1_hit_pattern = muon->get_hit_pattern();
  dimuon->Tr1_maxres_sigma = muon->get_maxres_sigma();
  dimuon->Tr1_track_id = muon->get_track_id();
  dimuon->Tr1_is_sfvtx = muon->get_is_sfvtx();
  dimuon->Tr1_best_fvtxmutr_match = muon->get_best_fvtxmutr_match();
  dimuon->Tr1_x_fvtxproj = muon->get_x_fvtxproj(0);
  dimuon->Tr1_y_fvtxproj = muon->get_y_fvtxproj(0);
}


MWGVertex mFillDiMuonContainer::get_vertex(unsigned short imu1, unsigned short imu2,
					   SingleMuonContainer* muons1, SingleMuonContainer* muons2)
{
  SingleMuon* muon1 = muons1->get_SingleMuon(imu1);
  SingleMuon* muon2 = muons2->get_SingleMuon(imu2);

  static PHMuoTracksOut *particle = MWG::newPHMuoTracksFvtx();
  particle->Reset();
  particle->set_npart(2);
  particle->set_TClonesArraySize(2);

  particle->AddPHParticle(0);
  particle->set_charge(0,(muon1->get_charge()>0) ? 1 : -1);
  particle->set_px(0,0,muon1->get_px());
  particle->set_py(0,0,muon1->get_py());
  particle->set_pz(0,0,muon1->get_pz());
  particle->set_xpos(0,0,muon1->get_x0());
  particle->set_ypos(0,0,muon1->get_y0());
  particle->set_zpos(0,0,muon1->get_z0());
  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,0,muon1->get_cov(i,j));

  particle->AddPHParticle(1);
  particle->set_charge(1, (muon2->get_charge()>0) ? 1 : -1);
  particle->set_px(0,1,muon2->get_px());
  particle->set_py(0,1,muon2->get_py());
  particle->set_pz(0,1,muon2->get_pz());
  particle->set_xpos(0,1,muon2->get_x0());
  particle->set_ypos(0,1,muon2->get_y0());
  particle->set_zpos(0,1,muon2->get_z0());
  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,1,muon2->get_cov(i,j));

  MWGVertex vertex;

  vertex.add_track(0,(PHMuoTracksOut*)particle);
  vertex.add_track(1,(PHMuoTracksOut*)particle);

  // Retrieve vertex assigned to each track from mFillSingleMuonContainer:
  PHPoint vtx, vtx_error;
  if (make_dimuon_vertex_origin_widerange)
    {
      vtx.setX(0.);
      vtx.setY(0.);
      vtx.setZ(0.);

      vtx_error.setX(1.);
      vtx_error.setY(1.);
      vtx_error.setZ(100. );
    }
  else if ( _vtx )
    get_dimuon_vertex(muon1, muon2, vtx, vtx_error);
  else if ( fvtx_vtx || sngtrk12_vtx )  
    {
      vtx = PHPoint(muons1->get_Evt_fvtxX(), muons1->get_Evt_fvtxY(), muons1->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  else  // sngtrk21
    {
      vtx = PHPoint(muons2->get_Evt_fvtxX(), muons2->get_Evt_fvtxY(), muons2->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  vertex.add_vertex(vtx, vtx_error);
    
  vertex.fit();
  return vertex;
}

MWGVertex mFillDiMuonContainer::get_vertex_from_fvtx(unsigned short imu1, unsigned short imu2,
						     SingleMuonContainer* muons1, SingleMuonContainer* muons2)
{
  SingleMuon* muon1 = muons1->get_SingleMuon(imu1);
  SingleMuon* muon2 = muons2->get_SingleMuon(imu2);
  
  static PHMuoTracksOut *particle = MWG::newPHMuoTracksFvtx();
  particle->Reset();
  particle->set_npart(2);
  particle->set_TClonesArraySize(2);

  particle->AddPHParticle(0);
  particle->set_charge(0, (muon1->get_charge()>0) ? 1 : -1);
  particle->set_px(0,0,muon1->get_px_fvtxmutr());
  particle->set_py(0,0,muon1->get_py_fvtxmutr());
  particle->set_pz(0,0,muon1->get_pz_fvtxmutr());
  particle->set_xpos(0,0,muon1->get_x0_fvtxmutr());
  particle->set_ypos(0,0,muon1->get_y0_fvtxmutr());
  particle->set_zpos(0,0,muon1->get_z0_fvtxmutr());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,0,muon1->get_cov_fvtxmutr(i,j));

  particle->AddPHParticle(1);
  particle->set_charge(1, (muon2->get_charge()>0) ? 1 : -1);
  particle->set_px(0,1,muon2->get_px_fvtxmutr());
  particle->set_py(0,1,muon2->get_py_fvtxmutr());
  particle->set_pz(0,1,muon2->get_pz_fvtxmutr());
  particle->set_xpos(0,1,muon2->get_x0_fvtxmutr());
  particle->set_ypos(0,1,muon2->get_y0_fvtxmutr());
  particle->set_zpos(0,1,muon2->get_z0_fvtxmutr());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,1,muon2->get_cov_fvtxmutr(i,j));

  MWGVertex vertex;
  vertex.add_track(0,(PHMuoTracksOut*)particle);
  vertex.add_track(1,(PHMuoTracksOut*)particle);

  // Retrieve vertex assigned to each track from mFillSingleMuonContainer:
  PHPoint vtx, vtx_error;
  if ( _vtx )
    get_dimuon_vertex(muon1, muon2, vtx, vtx_error);
  else if ( fvtx_vtx || sngtrk12_vtx )  
    {
      vtx = PHPoint(muons1->get_Evt_fvtxX(), muons1->get_Evt_fvtxY(), muons1->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  else 
    {
      vtx = PHPoint(muons2->get_Evt_fvtxX(), muons2->get_Evt_fvtxY(), muons2->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  vertex.add_vertex(vtx, vtx_error);

  vertex.fit();

  return vertex;
}

MWGVertex mFillDiMuonContainer::get_vertex_from_sngtrk12(unsigned short imu1, unsigned short imu2,
							 SingleMuonContainer* muons1, SingleMuonContainer* muons2) 
{
  SingleMuon* muon1 = muons1->get_SingleMuon(imu1);
  SingleMuon* muon2 = muons2->get_SingleMuon(imu2);
  
  static PHMuoTracksOut *particle = MWG::newPHMuoTracksFvtx();
  particle->Reset();
  particle->set_npart(2);
  particle->set_TClonesArraySize(2);

  particle->AddPHParticle(0);
  particle->set_charge(0, (muon1->get_charge()>0) ? 1 : -1);
  particle->set_px(0,0,muon1->get_px_fvtxmutr());
  particle->set_py(0,0,muon1->get_py_fvtxmutr());
  particle->set_pz(0,0,muon1->get_pz_fvtxmutr());
  particle->set_xpos(0,0,muon1->get_x0_fvtxmutr());
  particle->set_ypos(0,0,muon1->get_y0_fvtxmutr());
  particle->set_zpos(0,0,muon1->get_z0_fvtxmutr());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,0,muon1->get_cov_fvtxmutr(i,j));

  particle->AddPHParticle(1);
  particle->set_charge(1, (muon2->get_charge()>0) ? 1 : -1);
  particle->set_px(0,1,muon2->get_px());
  particle->set_py(0,1,muon2->get_py());
  particle->set_pz(0,1,muon2->get_pz());
  particle->set_xpos(0,1,muon2->get_x0());
  particle->set_ypos(0,1,muon2->get_y0());
  particle->set_zpos(0,1,muon2->get_z0());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,1,muon2->get_cov(i,j));

  MWGVertex vertex;
  vertex.add_track(0,(PHMuoTracksOut*)particle);
  vertex.add_track(1,(PHMuoTracksOut*)particle);

  // Retrieve vertex assigned to each track from mFillSingleMuonContainer:
  PHPoint vtx, vtx_error;
  if ( _vtx )
    get_dimuon_vertex(muon1, muon2, vtx, vtx_error);
  else if ( fvtx_vtx || sngtrk12_vtx )  
    {
      vtx = PHPoint(muons1->get_Evt_fvtxX(), muons1->get_Evt_fvtxY(), muons1->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  else 
    {
      vtx = PHPoint(muons2->get_Evt_fvtxX(), muons2->get_Evt_fvtxY(), muons2->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  vertex.add_vertex(vtx, vtx_error);

  vertex.fit();

  return vertex;
}

MWGVertex mFillDiMuonContainer::get_vertex_from_sngtrk21(unsigned short imu2, unsigned short imu1,
							 SingleMuonContainer* muons2, SingleMuonContainer* muons1) 
{
  SingleMuon* muon2 = muons2->get_SingleMuon(imu2);
  SingleMuon* muon1 = muons1->get_SingleMuon(imu1);
  
  static PHMuoTracksOut *particle = MWG::newPHMuoTracksFvtx();
  particle->Reset();
  particle->set_npart(2);
  particle->set_TClonesArraySize(2);

  particle->AddPHParticle(0);
  particle->set_charge(0, (muon2->get_charge()>0) ? 1 : -1);
  particle->set_px(0,0,muon2->get_px_fvtxmutr());
  particle->set_py(0,0,muon2->get_py_fvtxmutr());
  particle->set_pz(0,0,muon2->get_pz_fvtxmutr());
  particle->set_xpos(0,0,muon2->get_x0_fvtxmutr());
  particle->set_ypos(0,0,muon2->get_y0_fvtxmutr());
  particle->set_zpos(0,0,muon2->get_z0_fvtxmutr());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,0,muon2->get_cov_fvtxmutr(i,j));

  particle->AddPHParticle(1);
  particle->set_charge(1, (muon1->get_charge()>0) ? 1 : -1);
  particle->set_px(0,1,muon1->get_px());
  particle->set_py(0,1,muon1->get_py());
  particle->set_pz(0,1,muon1->get_pz());
  particle->set_xpos(0,1,muon1->get_x0());
  particle->set_ypos(0,1,muon1->get_y0());
  particle->set_zpos(0,1,muon1->get_z0());

  for (unsigned short i=0; i<5; i++)
    for (unsigned short j=0; j<5; j++)
      particle->set_cov(i,j,1,muon1->get_cov(i,j));

  MWGVertex vertex;
  vertex.add_track(0,(PHMuoTracksOut*)particle);
  vertex.add_track(1,(PHMuoTracksOut*)particle);

  // Retrieve vertex assigned to each track from mFillSingleMuonContainer:
  PHPoint vtx, vtx_error;
  if ( _vtx )
    get_dimuon_vertex(muon1, muon2, vtx, vtx_error);
  else if ( fvtx_vtx || sngtrk12_vtx )  
    {
      vtx = PHPoint(muons1->get_Evt_fvtxX(), muons1->get_Evt_fvtxY(), muons1->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
      fvtx_vtx = false; 
      sngtrk12_vtx = false; 
    }
  else 
    {
      vtx = PHPoint(muons2->get_Evt_fvtxX(), muons2->get_Evt_fvtxY(), muons2->get_Evt_fvtxZ());
      vtxerror = PHPoint(0.01, 0.01, 0.05);
    }
  vertex.add_vertex(vtx, vtx_error);

  vertex.fit();

  return vertex;
}

void mFillDiMuonContainer::get_dimuon_vertex(SingleMuon* muon1, SingleMuon* muon2, PHPoint& vtx, PHPoint& vtx_error)
{

  PHPoint vtx1, vtx2;
  PHPoint vtx1_error, vtx2_error;

  if ( is_sim ){
    vtx1 = _vtx->get_Vertex("SIM");
    vtx1_error = _vtx->get_VertexError("SIM");
    if ( vtx1_error.getZ() == 0 )
      vtx1_error.setZ(0.005);
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::FVTX){
    vtx1 = _vtx->get_Vertex(_fvtx_vertex_names[0].data());
    vtx1_error = _vtx->get_VertexError(_fvtx_vertex_names[0].data());
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::FVTX_SECOND){
    vtx1 = _vtx->get_Vertex(_fvtx_vertex_names[1].data());
    vtx1_error = _vtx->get_VertexError(_fvtx_vertex_names[1].data());
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::FVTX_3){
    vtx1 = _vtx->get_Vertex(_fvtx_vertex_names[2].data());
    vtx1_error = _vtx->get_VertexError(_fvtx_vertex_names[2].data());
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::FVTX_4){
    vtx1 = _vtx->get_Vertex(_fvtx_vertex_names[3].data());
    vtx1_error = _vtx->get_VertexError(_fvtx_vertex_names[3].data());
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::SVX_PRECISE){
    vtx1 = _vtx->get_Vertex("SVX_PRECISE");
    vtx1_error = _vtx->get_VertexError("SVX_PRECISE");
  }
  else if (muon1->get_vtx_index() == mFillSingleMuonContainer::BBC){
    vtx1 = _vtx->get_Vertex("BBC");
    vtx1_error = _vtx->get_VertexError("BBC");
  }

  if ( is_sim ){
    vtx2 = _vtx->get_Vertex("SIM");
    vtx2_error = _vtx->get_VertexError("SIM");
    if ( vtx2_error.getZ() == 0 )
      vtx2_error.setZ(0.5);
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::FVTX){
    vtx2 = _vtx->get_Vertex(_fvtx_vertex_names[0].data());
    vtx2_error = _vtx->get_VertexError(_fvtx_vertex_names[0].data());
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::FVTX_SECOND){
    vtx2 = _vtx->get_Vertex(_fvtx_vertex_names[1].data());
    vtx2_error = _vtx->get_VertexError(_fvtx_vertex_names[1].data());
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::FVTX_3){
    vtx2 = _vtx->get_Vertex(_fvtx_vertex_names[2].data());
    vtx2_error = _vtx->get_VertexError(_fvtx_vertex_names[2].data());
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::FVTX_4){
    vtx2 = _vtx->get_Vertex(_fvtx_vertex_names[3].data());
    vtx2_error = _vtx->get_VertexError(_fvtx_vertex_names[3].data());
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::SVX_PRECISE){
    vtx2 = _vtx->get_Vertex("SVX_PRECISE");
    vtx2_error = _vtx->get_VertexError("SVX_PRECISE");
  }
  else if (muon2->get_vtx_index() == mFillSingleMuonContainer::BBC){
    vtx2 = _vtx->get_Vertex("BBC");
    vtx2_error = _vtx->get_VertexError("BBC");
  }

  vtx.setX((vtx1.getX()+vtx2.getX())/2);
  vtx.setY((vtx1.getY()+vtx2.getY())/2);
  vtx.setZ((vtx1.getZ()+vtx2.getZ())/2);

  vtx_error.setX( sqrt((vtx1_error.getX()*vtx1_error.getX() + vtx2_error.getX()*vtx2_error.getX())/2.0) );
  vtx_error.setY( sqrt((vtx1_error.getY()*vtx1_error.getY() + vtx2_error.getY()*vtx2_error.getY())/2.0) );
  vtx_error.setZ( sqrt((vtx1_error.getZ()*vtx1_error.getZ() + vtx2_error.getZ()*vtx2_error.getZ())/2.0) );
}

void mFillDiMuonContainer::calculate_from_fvtx(SingleMuon* muon1, SingleMuon* muon2, DiMuon_v12* dimuon)
{
  double p1 = muon1->get_p();
  PHVector v1( muon1->get_px_fvtxmutr()/muon1->get_pz_fvtxmutr(),
	       muon1->get_py_fvtxmutr()/muon1->get_pz_fvtxmutr(), 1 );
  v1.normalize();
  if( muon1->get_pz()<0 ) v1 = v1*(-1);
  double phi1 = atan2(muon1->get_py_fvtxmutr(), muon1->get_px_fvtxmutr());
  double theta1 = atan(sqrt(v1.getX()*v1.getX() + v1.getY()*v1.getY()));
  double px1 = p1*cos(phi1)*sin(theta1);
  double py1 = p1*sin(phi1)*sin(theta1);
  
  double p2 = muon2->get_p();
  
  PHVector v2( muon2->get_px_fvtxmutr()/muon2->get_pz_fvtxmutr(),
	       muon2->get_py_fvtxmutr()/muon2->get_pz_fvtxmutr(), 1 );
  v2.normalize();
  if( muon2->get_pz()<0 ) v2 = v2*(-1);
  double phi2 = atan2(muon2->get_py_fvtxmutr(), muon2->get_px_fvtxmutr());
  double theta2 = atan(sqrt(v2.getX()*v2.getX() + v2.getY()*v2.getY()));
  double px2 = p2*cos(phi2)*sin(theta2);
  double py2 = p2*sin(phi2)*sin(theta2);
  
  double costh = cos(v1.angle(v2));
  double E1 = sqrt(p1*p1 + MU_MASS*MU_MASS);
  double E2 = sqrt(p2*p2 + MU_MASS*MU_MASS);
  double dot = p1*p2*costh;
  dimuon->mass_fvtx = sqrt(2*(E1*E2 + MU_MASS*MU_MASS - dot));
  dimuon->Px = px1+px2;
  dimuon->Py = py1+py2;
  
  if (costh!=1.00)
    {
      double r1 = sqrt(pow(muon1->get_x0_fvtxmutr(),2) + pow(muon1->get_y0_fvtxmutr(),2));
      PHCylPoint p1cyl(r1, phi1, muon1->get_z0_fvtxmutr());
      PHPoint pnt1 = PHPoint();
      PHGeometry::cylindricalToCartesian(p1cyl, pnt1);
      PHLine l1(pnt1, v1);
      
      double r2 = sqrt(pow(muon2->get_x0_fvtxmutr(),2) + pow(muon2->get_y0_fvtxmutr(),2));
      PHCylPoint p2cyl(r2, phi2, muon2->get_z0_fvtxmutr());
      PHPoint pnt2 = PHPoint();
      PHGeometry::cylindricalToCartesian(p2cyl, pnt2);
      PHLine l2(pnt2, v2);
      
      PHPoint CA = PHGeometry::closestApproachLineLine(l1, l2);      

      dimuon->X0_fvtx = CA.getX();
      dimuon->Y0_fvtx = CA.getY();
      dimuon->Z0_fvtx = CA.getZ();
      
      PHPoint CA1 = PHGeometry::closestApproachLinePoint(l1,CA);
      PHPoint CA2 = PHGeometry::closestApproachLinePoint(l2,CA);

      dimuon->dca_z = CA1.getZ() - CA2.getZ();
      
      PHPoint pz1 = PHPoint();
      PHPoint pz2 = PHPoint();

      PHVector vzref = PHVector(0,0,1);
      PHPoint pzref = PHPoint(0,0,CA.getZ());
      PHPlane plzref = PHPlane(pzref,vzref);

      if (PHGeometry::intersectionLinePlane(l1,plzref,pz1) && PHGeometry::intersectionLinePlane(l2,plzref,pz2))
	{
	  dimuon->dca_r = PHGeometry::distancePointToPoint(pz1,pz2);
	  phi1 = atan2(pz1.getY()-CA.getY(), pz1.getX()-CA.getX());
	  phi2 = atan2(pz2.getY()-CA.getY(), pz2.getX()-CA.getX());
	  dimuon->phi_dca = phi1 - phi2;
	  dimuon->phi_dca = 0.5*atan2(sin(2*dimuon->phi_dca),cos(2*dimuon->phi_dca));
	}
      else
	dimuon->dca_r = 1000;
    }
}


void mFillDiMuonContainer::calculate_from_sngtrk12(SingleMuon* muon1, SingleMuon* muon2, DiMuon_v12* dimuon) 
{
  double p1 = muon1->get_p(); 
  PHVector v1( muon1->get_px_fvtxmutr()/muon1->get_pz_fvtxmutr(),
	       muon1->get_py_fvtxmutr()/muon1->get_pz_fvtxmutr(), 1 );
  v1.normalize();
  if( muon1->get_pz()<0 ) 
    v1 = v1*(-1);
  double phi1 = atan2(muon1->get_py_fvtxmutr(), muon1->get_px_fvtxmutr());
  double theta1 = atan(sqrt(v1.getX()*v1.getX() + v1.getY()*v1.getY()));
  double px1 = p1*cos(phi1)*sin(theta1);
  double py1 = p1*sin(phi1)*sin(theta1);
  
  double p2 = muon2->get_p();
  
  PHVector v2( muon2->get_px()/muon2->get_pz(),
               muon2->get_py()/muon2->get_pz(), 1 );
  v2.normalize();
  if( muon2->get_pz()<0 ) v2 = v2*(-1);
  double phi2 = atan2(muon2->get_py(), muon2->get_px());
  double theta2 = atan(sqrt(v2.getX()*v2.getX() + v2.getY()*v2.getY()));
  double px2 = p2*cos(phi2)*sin(theta2);
  double py2 = p2*sin(phi2)*sin(theta2);
  
  double costh = cos(v1.angle(v2));
  double E1 = sqrt(p1*p1 + MU_MASS*MU_MASS);
  double E2 = sqrt(p2*p2 + MU_MASS*MU_MASS);
  double dot = p1*p2*costh;
  //dimuon->mass_sngtrk = sqrt(2*(E1*E2 + MU_MASS*MU_MASS - dot));
  dimuon->mass_fvtx = sqrt(2*(E1*E2 + MU_MASS*MU_MASS - dot));
  dimuon->Px = px1+px2;
  dimuon->Py = py1+py2;
  
  if (costh!=1.00)
    {
      double r1 = sqrt(pow(muon1->get_x0_fvtxmutr(),2) + pow(muon1->get_y0_fvtxmutr(),2));
      PHCylPoint p1cyl(r1, phi1, muon1->get_z0_fvtxmutr());
      PHPoint pnt1 = PHPoint();
      PHGeometry::cylindricalToCartesian(p1cyl, pnt1);
      PHLine l1(pnt1, v1);
      
      double r2 = sqrt(pow(muon2->get_x0(),2) + pow(muon2->get_y0(),2));
      PHCylPoint p2cyl(r2, phi2, muon2->get_z0());
      PHPoint pnt2 = PHPoint();
      PHGeometry::cylindricalToCartesian(p2cyl, pnt2);
      PHLine l2(pnt2, v2);
      
      PHPoint CA = PHGeometry::closestApproachLineLine(l1, l2);      

      //dimuon->X0_sngtrk = CA.getX();  
      //dimuon->Y0_sngtrk = CA.getY();
      //dimuon->Z0_sngtrk = CA.getZ();
      dimuon->X0_fvtx = CA.getX();  
      dimuon->Y0_fvtx = CA.getY();
      dimuon->Z0_fvtx = CA.getZ();
      
      PHPoint CA1 = PHGeometry::closestApproachLinePoint(l1,CA);
      PHPoint CA2 = PHGeometry::closestApproachLinePoint(l2,CA);

      dimuon->dca_z = CA1.getZ() - CA2.getZ();
      
      PHPoint pz1 = PHPoint();
      PHPoint pz2 = PHPoint();

      PHVector vzref = PHVector(0,0,1);
      PHPoint pzref = PHPoint(0,0,CA.getZ());
      PHPlane plzref = PHPlane(pzref,vzref);

      if (PHGeometry::intersectionLinePlane(l1,plzref,pz1) && PHGeometry::intersectionLinePlane(l2,plzref,pz2))
        {
          dimuon->dca_r = PHGeometry::distancePointToPoint(pz1,pz2);
          phi1 = atan2(pz1.getY()-CA.getY(), pz1.getX()-CA.getX());
          phi2 = atan2(pz2.getY()-CA.getY(), pz2.getX()-CA.getX());
          dimuon->phi_dca = phi1 - phi2;
          dimuon->phi_dca = 0.5*atan2(sin(2*dimuon->phi_dca),cos(2*dimuon->phi_dca));
        }
      else
        dimuon->dca_r = 1000;
    }
}

void mFillDiMuonContainer::calculate_from_sngtrk21(SingleMuon* muon2, SingleMuon* muon1, DiMuon_v12* dimuon) 
{
  double p2 = muon2->get_p();  // muon2 with fvtx 
  PHVector v2( muon2->get_px_fvtxmutr()/muon2->get_pz_fvtxmutr(),
	       muon2->get_py_fvtxmutr()/muon2->get_pz_fvtxmutr(), 1 );
  v2.normalize();
  if( muon2->get_pz()<0 ) 
    v2 = v2*(-1);
  double phi2 = atan2(muon2->get_py_fvtxmutr(), muon2->get_px_fvtxmutr());
  double theta2 = atan(sqrt(v2.getX()*v2.getX() + v2.getY()*v2.getY()));
  double px2 = p2*cos(phi2)*sin(theta2);
  double py2 = p2*sin(phi2)*sin(theta2);
  
  double p1 = muon1->get_p();
  
  PHVector v1( muon1->get_px()/muon1->get_pz(),
               muon1->get_py()/muon1->get_pz(), 1 );
  v1.normalize();
  if( muon1->get_pz()<0 ) v1 = v1*(-1);
  double phi1 = atan2(muon1->get_py(), muon1->get_px());
  double theta1 = atan(sqrt(v1.getX()*v1.getX() + v1.getY()*v1.getY()));
  double px1 = p1*cos(phi1)*sin(theta1);
  double py1 = p1*sin(phi1)*sin(theta1);
  
  double costh = cos(v2.angle(v1));
  double E2 = sqrt(p2*p2 + MU_MASS*MU_MASS);
  double E1 = sqrt(p1*p1 + MU_MASS*MU_MASS);
  double dot = p2*p1*costh;
  //dimuon->mass_sngtrk = sqrt(2*(E2*E1 + MU_MASS*MU_MASS - dot));
  dimuon->mass_fvtx = sqrt(2*(E2*E1 + MU_MASS*MU_MASS - dot));
  dimuon->Px = px2+px1;
  dimuon->Py = py2+py1;
  
  if (costh!=1.00)
    {
      double r2 = sqrt(pow(muon2->get_x0_fvtxmutr(),2) + pow(muon2->get_y0_fvtxmutr(),2));
      PHCylPoint p2cyl(r2, phi2, muon2->get_z0_fvtxmutr());
      PHPoint pnt2 = PHPoint();
      PHGeometry::cylindricalToCartesian(p2cyl, pnt2);
      PHLine l2(pnt2, v2);
      
      double r1 = sqrt(pow(muon1->get_x0(),2) + pow(muon1->get_y0(),2));
      PHCylPoint p1cyl(r1, phi1, muon1->get_z0());
      PHPoint pnt1 = PHPoint();
      PHGeometry::cylindricalToCartesian(p1cyl, pnt1);
      PHLine l1(pnt1, v1);
      
      PHPoint CA = PHGeometry::closestApproachLineLine(l2, l1);      

      //dimuon->X0_sngtrk = CA.getX();
      //dimuon->Y0_sngtrk = CA.getY();
      //dimuon->Z0_sngtrk = CA.getZ();
      dimuon->X0_fvtx = CA.getX();
      dimuon->Y0_fvtx = CA.getY();
      dimuon->Z0_fvtx = CA.getZ();
      
      PHPoint CA2 = PHGeometry::closestApproachLinePoint(l2,CA);
      PHPoint CA1 = PHGeometry::closestApproachLinePoint(l1,CA);

      dimuon->dca_z = CA2.getZ() - CA1.getZ();
      
      PHPoint pz2 = PHPoint();
      PHPoint pz1 = PHPoint();

      PHVector vzref = PHVector(0,0,1);
      PHPoint pzref = PHPoint(0,0,CA.getZ());
      PHPlane plzref = PHPlane(pzref,vzref);

      if (PHGeometry::intersectionLinePlane(l2,plzref,pz2) && PHGeometry::intersectionLinePlane(l1,plzref,pz1))
        {
          dimuon->dca_r = PHGeometry::distancePointToPoint(pz2,pz1);
          phi2 = atan2(pz2.getY()-CA.getY(), pz2.getX()-CA.getX());
          phi1 = atan2(pz1.getY()-CA.getY(), pz1.getX()-CA.getX());
          dimuon->phi_dca = phi2 - phi1;
          dimuon->phi_dca = 0.5*atan2(sin(2*dimuon->phi_dca),cos(2*dimuon->phi_dca));
        }
      else
        dimuon->dca_r = 1000;
    }
}
        
void
mFillDiMuonContainer::fill_mom_vtx(MWGVertex &vertex, DiMuon_v12 &dimuon)
{
  dimuon.Px = vertex.get_px(0) + vertex.get_px(1);
  dimuon.Py = vertex.get_py(0) + vertex.get_py(1);
  dimuon.Pz = vertex.get_pz(0) + vertex.get_pz(1);
  
  dimuon.rapidity = vertex.get_rapidity();
  dimuon.Evt_vtxchi2 = vertex.get_chisquare() / vertex.get_ndf();
  
  dimuon.X0 = vertex.get_vtx_x();
  dimuon.Y0 = vertex.get_vtx_y();
  dimuon.Z0 = vertex.get_vtx_z();
  
}

void
mFillDiMuonContainer::fill_mom_vtx_fvtx(MWGVertex &vertex_fvtx, DiMuon_v12 &dimuon)
{
  dimuon.Px_fvtxmutr = vertex_fvtx.get_px(0) + vertex_fvtx.get_px(1);
  dimuon.Py_fvtxmutr = vertex_fvtx.get_py(0) + vertex_fvtx.get_py(1);
  dimuon.Pz_fvtxmutr = vertex_fvtx.get_pz(0) + vertex_fvtx.get_pz(1);

  dimuon.rapidity_fvtxmutr = vertex_fvtx.get_rapidity();
  dimuon.Evt_vtxchi2_fvtxmutr = vertex_fvtx.get_chisquare() / vertex_fvtx.get_ndf();

  dimuon.X0_fvtxmutr = vertex_fvtx.get_vtx_x();
  dimuon.Y0_fvtxmutr = vertex_fvtx.get_vtx_y();
  dimuon.Z0_fvtxmutr = vertex_fvtx.get_vtx_z();

}

void mFillDiMuonContainer::fill_mom_vtx_sngtrk(MWGVertex &vertex_sngtrk, DiMuon_v12 &dimuon) 
{
  /*
    dimuon.Px_sngtrk = vertex_sngtrk.get_px(0) + vertex_sngtrk.get_px(1);  
    dimuon.Py_sngtrk = vertex_sngtrk.get_py(0) + vertex_sngtrk.get_py(1);   
    dimuon.Pz_sngtrk = vertex_sngtrk.get_pz(0) + vertex_sngtrk.get_pz(1);   

    dimuon.rapidity_sngtrk = vertex_sngtrk.get_rapidity();
    dimuon.Evt_vtxchi2_sngtrk = vertex_sngtrk.get_chisquare() / vertex_sngtrk.get_ndf();   

    dimuon.X0_sngtrk = vertex_sngtrk.get_vtx_x();   
    dimuon.Y0_sngtrk = vertex_sngtrk.get_vtx_y();   
    dimuon.Z0_sngtrk = vertex_sngtrk.get_vtx_z();   
  */

  dimuon.Px_fvtxmutr = vertex_sngtrk.get_px(0) + vertex_sngtrk.get_px(1);  
  dimuon.Py_fvtxmutr = vertex_sngtrk.get_py(0) + vertex_sngtrk.get_py(1);   
  dimuon.Pz_fvtxmutr = vertex_sngtrk.get_pz(0) + vertex_sngtrk.get_pz(1);   

  dimuon.rapidity_fvtxmutr = vertex_sngtrk.get_rapidity();
  dimuon.Evt_vtxchi2_fvtxmutr = vertex_sngtrk.get_chisquare() / vertex_sngtrk.get_ndf();   

  dimuon.X0_fvtxmutr = vertex_sngtrk.get_vtx_x();   
  dimuon.Y0_fvtxmutr = vertex_sngtrk.get_vtx_y();   
  dimuon.Z0_fvtxmutr = vertex_sngtrk.get_vtx_z();   

}    


