/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/* Simple class for working with truth DAODs */
/* Goal is to validate truth DAOD contents:
    MET_Truth
    Truth particles with aux: TruthElectrons, TruthMuons, TruthPhotons, TruthTaus, TruthNeutrinos, TruthBSM, TruthBottom, TruthTop, TruthBoson
    Full collections with aux: TruthWbosonWithDecayParticles, TruthWbosonWithDecayVertices
    AntiKt4TruthDressedWZJets.GhostCHadronsFinalCount.GhostBHadronsFinalCount.pt.HadronConeExclTruthLabelID.ConeTruthLabelID.PartonTruthLabelID.TrueFlavor // with Aux
    AntiKt10TruthTrimmedPtFrac5SmallR20Jets.pt.Tau1_wta.Tau2_wta.Tau3_wta.D2 // with Aux
    TruthEvents.Q.XF1.XF2.PDGID1.PDGID2.PDFID1.PDFID2.X1.X2.crossSection
    xAOD::TruthMetaDataContainer#TruthMetaData // with Aux
   Plus extra features
*/

// Setup for reading ATLAS data
#ifdef XAOD_STANDALONE
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#else
#include "POOLRootAccess/TEvent.h"
#include "GaudiKernel/SmartIF.h"
#endif

// Boost for argument parsing
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Truth weight tool
#include "AsgTools/MessageCheck.h"
#include "AsgTools/AnaToolHandle.h"
#include "PMGAnalysisInterfaces/IPMGTruthWeightTool.h"
//#include "PMGTools/PMGTruthWeightTool.h"

// ATLAS data dependencies
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/TruthMetaDataContainer.h"

// ROOT dependencies
#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>

// STL includes
#include <iostream>
#include <cmath>

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//	 		Helper Functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Helper macro for retrieving containers from the event
#define CHECK_RETRIEVE( container , name ) { \
    if (!event.retrieve( container , name ).isSuccess()){ \
      Error( APP_NAME , "Could not load event %s from the file!" , name ); \
      throw; \
    } \
  }

// Helper for counting children
int countChildren( const xAOD::TruthParticle * p ){
  if (!p) return 0;
  int children = 0;
  for (size_t n=0;n<p->nChildren();++n){
    children += countChildren( p->child(n) );
  }
  return children;
}

// Helper for counting parents
int countParents( const xAOD::TruthParticle * p ){
  if (!p) return 0;
  int parents = 0;
  for (size_t n=0;n<p->nParents();++n){
    parents += countParents( p->parent(n) );
  }
  return parents;
}

// Helper for dPhi
float getdPhi(float phi1, float phi2){
  float dPhi = fabs(phi1 - phi2);
   if(dPhi > 3.14)  dPhi -= 2*3.14;
   if(dPhi < -3.14) dPhi += 2*3.14;
   return dPhi;
}

// Helper for mT
float GetMt(TLorentzVector v1, TLorentzVector v2, float met_met, float met_phi){
 	TLorentzVector jj = v1 + v2;
	TLorentzVector met_v;
	met_v.SetPtEtaPhiM(met_met,0,met_phi,0.0);
 	float dijetEt = sqrt(pow(jj.M(),2) + pow(jj.Pt(),2));
	float mT2 = pow(dijetEt + met_v.Et(),2) - pow((jj+met_v).Pt(),2);	
 	return sqrt(mT2);
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//			Main Routine
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Main routine... here we go!
int main(int argc, char **argv) {

  // The application's name:
  const char* APP_NAME = argv[ 0 ];

  // Make sure things know we are not in StatusCode land
  using namespace asg::msgUserCode;
  ANA_CHECK_SET_TYPE (int);

  // Setup for reading -- if this fails, we have major problems
#ifdef XAOD_STANDALONE
  if ( ! xAOD::Init().isSuccess() ) {
    throw std::runtime_error("Cannot initialise xAOD access !");
  }
  ANA_MSG_INFO("Using xAOD access");
#else
  SmartIF<IAppMgrUI> app = POOL::Init();
  ANA_MSG_INFO("Using POOL access");
#endif

  // Get and parse the command-line arguments
  po::options_description desc("Running a simple truth analysis");
  std::string outputName,inputName;
  desc.add_options()
    ("help,h", "print usage and exit")
    ("output,o",    po::value<std::string>(&outputName), "Output file name")
    ("input,i",     po::value<std::string>(&inputName), "Input file name")
    ("nevents", po::value<int>()->default_value(-1), "number of events to run on (set to -1 to ignore it")
    ;

  po::variables_map vm;
  po::store( po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  long long nevents = -1;
  if (vm.count("nevents")) {
      nevents = vm["nevents"].as<int>();
  }

  Info( "Reading from %s and writing to %s" , inputName.c_str() , outputName.c_str() );

  // Input chain
  std::unique_ptr< TFile > ifile( TFile::Open( inputName.c_str() , "READ" ) );
  ANA_CHECK( ifile.get() );
#ifdef XAOD_STANDALONE
  xAOD::TEvent event(xAOD::TEvent::kClassAccess);
#else
  POOL::TEvent event(POOL::TEvent::kClassAccess);
#endif
  ANA_CHECK( event.readFrom( ifile.get() ) );

  // Load metadata
  event.getEntries();

  // Make some temporary variables that we'll get out during the event loop
  const size_t nParticleContainers = 9;
  std::string particleKeyList[nParticleContainers] = { "TruthElectrons", "TruthMuons", "TruthPhotons", "TruthTaus", "TruthNeutrinos", "TruthBSM", "TruthBottom", "TruthTop", "TruthBoson" };
  const xAOD::TruthEventContainer* xTruthEventContainer(nullptr);
  const xAOD::EventInfo* xEventInfo(nullptr);
  const xAOD::JetContainer* smallRJets(nullptr);
  const xAOD::JetContainer* largeRJets(nullptr);
  const xAOD::TruthMetaDataContainer* truthMeta(nullptr);
  const xAOD::MissingETContainer* truthMET(nullptr);
  //const xAOD::TruthParticleContainer * truthParticles[nParticleContainers];
  const xAOD::TruthParticleContainer * truthParticles(nullptr);
  //const xAOD::TruthParticleContainer * truthBSM(nullptr);
  //for (size_t n=0;n<nParticleContainers;++n) truthParticles[n] = nullptr;

/*
    Missing full collections with aux: TruthBosonsWithDecayParticles, TruthBosonsWithDecayVertices
*/

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //			Define Histograms
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // MET histograms
  TH1D* h_metNonInt = new TH1D("MET_NonInt","MET; E^{T}_{miss}",50,0,800.);
  TH1D* h_metInt = new TH1D("MET_Int","",50,0,200.);
  TH1D* h_metIntOut = new TH1D("MET_IntOut","",50,0,100.);
  TH1D* h_metIntMuons = new TH1D("MET_IntMuons","",50,0,300.);
  // Particle collections: pT and connections for all
  //TH1D* h_partPt[nParticleContainers];
  //for (size_t n=0;n<nParticleContainers;++n) h_partPt[n] = new TH1D((particleKeyList[n]+"_pT").c_str(),"",50,0.,500.);
  //TH1D* h_partConn[nParticleContainers];
  //for (size_t n=0;n<nParticleContainers;++n) h_partConn[n] = new TH1D((particleKeyList[n]+"_Connections").c_str(),"",35,-10,25);
  //Dark Quarks
  TH1D* h_xdPt = new TH1D("xdPt", "Dark Quark Pt; p_{T} [GeV]", 70, 0, 4000.);
  TH1D* h_xdEta = new TH1D("xdEta", "Dark Quark Eta; #eta", 50, -5., 5.);
  TH1D* h_xdPhi = new TH1D("xdPhi", "Dark Quark Phi; #phi", 50, -3.14, 3.14);
  TH1D* h_xdM = new TH1D("xdM", "Dark Quark Mass; M [GeV]", 20, 5, 15);
  TH1D* h_xdDPhi = new TH1D("xdDPhi", "#Delta#phi Dark Quarks; #Delta#phi(xd,xd)", 20, 0,3.3);

  TH1D* h_zpPt = new TH1D("zpPt", "Z' Pt; p_{T} [GeV]", 70, 0, 4000.);
  TH1D* h_zpEta = new TH1D("zpEta", "Z' Eta; #eta", 50, -5., 5.);
  TH1D* h_zpPhi = new TH1D("zpPhi", "Z' Phi; #phi", 50, -3.14, 3.14);
  TH1D* h_zpM = new TH1D("zpM", "Z' Mass; M [GeV]", 80, 0, 8000);

  TH1D* h_r_inv = new TH1D("r_inv", "r_{inv} fraction; r_{inv}", 20, 0, 1);
  TH1D* h_children = new TH1D("children", "children of all dark hadrons; pdgId",80,0,80);
  TH1D* h_children_111 = new TH1D("children_111", "children of 4900111; PDGID",80,0,80);
  TH1D* h_children_113 = new TH1D("children_113", "children of 4900113; PDGID",80,0,80);
  TH1D* h_children_211 = new TH1D("children_211", "children of 4900211; PDGID",80,0,80);
  TH1D* h_children_213 = new TH1D("children_213", "children of 4900213; PDGID",80,0,80);
  TH1D* h_hadrons = new TH1D("hadrons", "Hadron Decays; ID", 16,0,16);
  h_hadrons->GetXaxis()->SetBinLabel(1,"4900111");
  h_hadrons->GetXaxis()->SetBinLabel(4,"No decay");
  h_hadrons->GetXaxis()->SetBinLabel(3,"Invisible");
  h_hadrons->GetXaxis()->SetBinLabel(2,"Visible");
  h_hadrons->GetXaxis()->SetBinLabel(5,"4900113");
  h_hadrons->GetXaxis()->SetBinLabel(8,"No decay");
  h_hadrons->GetXaxis()->SetBinLabel(7,"Invisible");
  h_hadrons->GetXaxis()->SetBinLabel(6,"Visible");
  h_hadrons->GetXaxis()->SetBinLabel(9,"4900211");
  h_hadrons->GetXaxis()->SetBinLabel(12,"No decay");
  h_hadrons->GetXaxis()->SetBinLabel(11,"Invisible");
  h_hadrons->GetXaxis()->SetBinLabel(10,"Visible");
  h_hadrons->GetXaxis()->SetBinLabel(13,"4900213");
  h_hadrons->GetXaxis()->SetBinLabel(16,"No decay");
  h_hadrons->GetXaxis()->SetBinLabel(15,"Invisible");
  h_hadrons->GetXaxis()->SetBinLabel(14,"Visible");

  TH1D* h_xdxdM = new TH1D("xdxdM", "Dark Quark Invariant Mass; M_{xd,xd} [GeV]",80, 0, 8000);
  TH1D* h_nJetsMatched = new TH1D("nJetsMatched", "N jets satisfying dR(j,xd) < 0.4; nJets", 10, 0, 10);
  TH1D* h_dRxdj = new TH1D("dRxdj", "DeltaR(quark, closest jet); #Delta R", 20,0,5);
  TH1D* h_dR_MET = new TH1D("dR_MET", "DeltaR(quark, MET aligned jet); #Delta R", 25,0,5);
  TH1D* h_dR_aMET = new TH1D("dR_aMET", "DeltaR(quark, MET anti-aligned jet); #Delta R", 25,0,5);
  TH1D* h_dRxdj1 = new TH1D("dRxdj1", "DeltaR(quark, leading jet); #Delta R", 25,0,5);
  TH1D* h_dRxdj2 = new TH1D("dRxdj2", "DeltaR(quark, subleading jet); #Delta R", 25,0,5);
  TH1D* h_mjj = new TH1D("mjj", "Invariant Mass 2 Closest Jets; M_{jj} [GeV]", 50, 0, 4000);
  TH1D* h_mT_12 = new TH1D("mT_12", "mT Sum (2 Leading + MET); m_{T} [GeV]", 80, 0, 5000);
  TH1D* h_mT_jj = new TH1D("mT_jj", "mT Sum (2 Matched Jets + MET); m_{T} [GeV]", 50, 0, 2500);
  TH1D* h_xdj_match_idx = new TH1D("xdj_match_idx", "Index of (matched) Jet Closest To Quark", 10, 0, 10);  
  TH1D* h_xdj_idx = new TH1D("xdj_idx", "Index of (any) Jet Closest To Quark", 15, 0, 15);  
  TH1D* h_dPhi_j_MET = new TH1D("dPhi_j_MET", "#Delta#phi MET and Closest Jet; #Delta#phi(MET,j_{closest})", 20, 0, 3.3);
  TH1D* h_dPhi_xdj_MET = new TH1D("dPhi_xdj_MET", "#Delta#phi MET and Closest dR Matched Jet; #Delta#phi(MET,xdj_{closest})", 20, 0, 3.3);
  TH1D* h_dPhi_xd_MET = new TH1D("dPhi_xd_MET", "#Delta#phi MET and Closest Dark Quark; #Delta#phi(MET,xd_{closest})", 20, 0, 3.3);

  // Small-R jets
  TH1D* h_nSmallR = new TH1D("nSmallR", "nJets (Small-R); N_{jets}", 15, 0, 15);
  TH1D* h_jetPt = new TH1D("JetPt","Leading Jet Pt (Small-R); p_{T} [GeV]",80,0,1000.);
  TH1D* h_jetEta = new TH1D("JetEta","Leading Jet Eta (Small-R); #eta",50,-5.,5.);
  TH1D* h_jetPhi = new TH1D("JetPhi","Leading Jet Phi (Small-R); #phi",50,-3.14,3.14);
  TH1D* h_jjdPhi = new TH1D("jjDPhi","#Delta#phi dR Matched Jets (Small-R); #Delta#phi(j,j)", 20, 0, 3.3);

  // Large R : pT, tau2, D2
  TH1D* h_nLargeR = new TH1D("nLargeR", "nJets (Large-R); N_{jets}", 10, 0, 10);
  TH1D* h_jetLRPt = new TH1D("JetLRPt","Jet Pt (Large-R); p_{T} [GeV]",80,0,3000.);
  TH1D* h_jetLREta = new TH1D("JetLREta","Jet Eta (Large-R); #eta",50,-5.,5.);
  TH1D* h_jetLRPhi = new TH1D("JetLRPhi","Jet Phi (Large-R); #phi",50,-3.14,3.14);
  // Truth events: PDF info , cross section , filter outcome
  TH1D* h_x1 = new TH1D("X1","",50,0,1.);
  TH1D* h_x2 = new TH1D("X2","",50,0,1.);
  TH1D* h_xsec = new TH1D("XSec","",50,0,1.);
  TH1D* h_HTFilt = new TH1D("HTFilt","",50,0,1000.);
  // Lazy for truth weights
  asg::AnaToolHandle< PMGTools::IPMGTruthWeightTool > weightTool("PMGTools::PMGTruthWeightTool/PMGTruthWeightTool");
  ANA_CHECK(weightTool.retrieve());

  // For the weights in the file
  //std::vector<TH1D*> h_weights;
  // For testing the truth metadata
  uint32_t channelNumber = 0;
  std::vector<std::string> weightNames;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //			Event Loop
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Into the event loop, with an optional cap
  long long nEvents = nevents > 0 ? std::min(static_cast<long long>(event.getEntries()), nevents) : event.getEntries();
  Info( APP_NAME , "Beginning event loop over %llu events." , nEvents );
  for (long evt=0;evt<nEvents;++evt){
    // Grab the data
    event.getEntry(evt);
    if (evt%1000==0) Info( APP_NAME , "Working on entry %lu" , evt );

    // Get the containers out
    CHECK_RETRIEVE( xTruthEventContainer , "TruthEvents" )
    CHECK_RETRIEVE( xEventInfo , "EventInfo" )
    CHECK_RETRIEVE( smallRJets , "AntiKt4TruthDressedWZJets" )
    CHECK_RETRIEVE( largeRJets , "AntiKt10TruthTrimmedPtFrac5SmallR20Jets" )
    CHECK_RETRIEVE( truthMET , "MET_Truth" )
    //CHECK_RETRIEVE( truthBSM, "TruthBSM" )
    CHECK_RETRIEVE( truthParticles, "TruthParticles" )

    //for (size_t n=0;n<nParticleContainers;++n){
    //  CHECK_RETRIEVE( truthParticles[n] , particleKeyList[n].c_str() )
    //}

    if (!event.retrieveMetaInput( truthMeta , "TruthMetaData" ).isSuccess()){
      Error( APP_NAME , "Could not load event TruthMetaData from the file!" );
      throw;
    }

    // Metadata check
    if (truthMeta->size()>1){
      Error( APP_NAME , "Truth metadata size: %lu . No one will look past item 0!" , truthMeta->size() );
      throw;
    }
    if (channelNumber==0){
      channelNumber = (*truthMeta)[0]->mcChannelNumber();
      weightNames = std::vector<std::string>((*truthMeta)[0]->weightNames());
      Info( APP_NAME , "Got channel number %u and %lu weight names" , channelNumber , weightNames.size() );
    }
    if (channelNumber != (*truthMeta)[0]->mcChannelNumber()){
      Error( APP_NAME , "Channel number changed mid-file! Was: %u now: %u " , channelNumber , (*truthMeta)[0]->mcChannelNumber() );
      throw;
    }
    if (weightNames != (*truthMeta)[0]->weightNames() ){
      Error( APP_NAME , "Weights have changed!" );
      for (size_t n=0;n<std::max( weightNames.size() , (*truthMeta)[0]->weightNames().size() );++n){
        std::cerr << "   " << n << " ";
        if (n<weightNames.size()) std::cerr << weightNames[n] << " ";
        else                      std::cerr << "- ";
        if (n<(*truthMeta)[0]->weightNames().size()) std::cerr << (*truthMeta)[0]->weightNames()[n] << " ";
        else                                         std::cerr << "- ";
        std::cerr << std::endl;
      }
      throw;
    }

    // Truth particles
    float x1=0.,x2=0.;
    (*xTruthEventContainer)[0]->pdfInfoParameter( x1 , xAOD::TruthEvent::X1 );
    (*xTruthEventContainer)[0]->pdfInfoParameter( x2 , xAOD::TruthEvent::X2 );
    h_x1->Fill( x1 );
    h_x2->Fill( x2 );
    h_xsec->Fill( (*xTruthEventContainer)[0]->crossSection() );
    h_HTFilt->Fill( xEventInfo->auxdata<float>("GenFiltHT") );

    // For MET: NonInt, Int, IntOut, IntMuons
    h_metNonInt->Fill( (*truthMET)["NonInt"]->met()*0.001 );
    //h_metNonInt->Fill( (*truthMET)["Int"]->met()*0.001 );
    //h_metNonInt->Fill( (*truthMET)["IntOut"]->met()*0.001 );
    //h_metNonInt->Fill( (*truthMET)["IntMuons"]->met()*0.001 );
    // Truth particles
    //for (size_t n=0;n<nParticleContainers;++n){ // PT and connections for all
    //  for (const auto * p : *truthParticles[n]){
    //    h_partPt[n]->Fill( p->pt() );
    //    h_partConn[n]->Fill( countChildren( p ) );
    //    h_partConn[n]->Fill( -1-countParents( p ) );
    //  }
    //}
    std::vector<TLorentzVector> quarks;
    int n_vis = 0;
    int n_invs = 0;
    float r_inv;
    for (const auto * bsm: *truthParticles){
      if ( fabs(bsm->pdgId()) == 4900101 && bsm->status() == 23){
      /*if ( fabs(bsm->pdgId()) == 4900101){
	std::cout << "pdgId = " << bsm->pdgId() << ", nParents = " << bsm->nParents() << ", status = " << bsm->status() << " | ";
	for (size_t n=0; n<bsm->nParents();n++){
		const xAOD::TruthParticle *p = bsm->parent(n);
		std::cout << "parent pdgId = " << p->pdgId() << ", parent status = " << p->status();
        }
        std::cout<<std::endl;*/
        TLorentzVector v_xd(0,0,0,0);
        v_xd.SetPxPyPzE(bsm->px()*0.001, bsm->py()*0.001, bsm->pz()*0.001, bsm->e()*.001);
        h_xdPt->Fill( v_xd.Pt() );
        h_xdM->Fill( v_xd.M() );
        h_xdEta->Fill( v_xd.Eta() );
        h_xdPhi->Fill( v_xd.Phi() );
        quarks.push_back(v_xd);
      }
      if ( bsm->pdgId() == 5000001 ){
        TLorentzVector v_zp(0,0,0,0);
        v_zp.SetPxPyPzE(bsm->px()*0.001, bsm->py()*0.001, bsm->pz()*0.001, bsm->e()*.001);
	h_zpPt->Fill( v_zp.Pt() );
	h_zpEta->Fill( v_zp.Eta() );
	h_zpPhi->Fill( v_zp.Phi() );
	h_zpM->Fill( v_zp.M() );
      }
      if (fabs(bsm->pdgId()) <4900214 && fabs(bsm->pdgId()) > 4900110){
      //if (fabs(bsm->pdgId()) == 4900211){
	// calculate r_inv
	bool visible = false;
        for (size_t n=0; n<bsm->nChildren();n++){
           const xAOD::TruthParticle *c = bsm->child(n); 
           //if (fabs(c->pdgId())< 23){visible = true; /*n_vis++;*/}
           //else if (fabs(c->pdgId()) != 51 && fabs(c->pdgId()) != 53) std::cout << "Unknown child ID " << c->pdgId() << " for hadron " << bsm->pdgId() << std::endl;
	   if (fabs(c->pdgId()) != 51 && fabs(c->pdgId()) != 53) visible = true;
           h_children->Fill(fabs(c->pdgId()));
	   if (fabs(bsm->pdgId()) == 4900111) h_children_111->Fill(fabs(c->pdgId()));
	   if (fabs(bsm->pdgId()) == 4900113) h_children_113->Fill(fabs(c->pdgId()));
	   if (fabs(bsm->pdgId()) == 4900211) h_children_211->Fill(fabs(c->pdgId()));
	   if (fabs(bsm->pdgId()) == 4900213) h_children_213->Fill(fabs(c->pdgId()));
	}
        if (visible) n_vis++;
	//else if (!visible && bsm->nChildren() > 0) n_invs++;
	else n_invs++;

	int factor=-16;
        if (fabs(bsm->pdgId()) == 4900111) factor = 0;
        if (fabs(bsm->pdgId()) == 4900113) factor = 4;
        if (fabs(bsm->pdgId()) == 4900211) factor = 8;
        if (fabs(bsm->pdgId()) == 4900213) factor = 12;
	h_hadrons->Fill(factor);
	if(bsm->nChildren() == 0) h_hadrons->Fill(factor+3);
	if(!visible && bsm->nChildren() > 0) h_hadrons->Fill(factor+2);
	if(visible) h_hadrons->Fill(factor+1);
      }//dark hadron if
    }//bsm loop
    r_inv = (float)n_invs/((float)n_invs+(float)n_vis);
    h_r_inv->Fill(r_inv);

    //if(quarks.size() != 2){std::cout << "No quarks found" << std::endl; continue;}
    h_xdxdM->Fill((quarks[0]+quarks[1]).M());
    h_xdDPhi->Fill(fabs(quarks[0].DeltaPhi(quarks[1])));

    //h_nSmallR->Fill(smallRJets->size());
    int nSmallR = 0;
    for (const auto * j : *smallRJets){
      if (j->pt() > 100000 && fabs(j->eta()) < 2.5) nSmallR++;
    }
    h_nSmallR->Fill(nSmallR);

    // MET Vec
    TLorentzVector v_met(0,0,0,0);
    v_met.SetPtEtaPhiM((*truthMET)["NonInt"]->met()*0.001,0,(*truthMET)["NonInt"]->phi(),0.0);

    // Small R Jets
    std::vector<TLorentzVector> jets;
    int j_xd1_idx = -1;
    int j_xd2_idx = -1;
    float dR_xd1_min = 10;
    float dR_xd2_min = 10;
    float dPhi_any_j_min = 10;
    float dPhi_any_j_max = 0;
    float dPhi_matched_j_min = 10;
    TLorentzVector j_xd1, j_xd2, j_MET, j_aMET, jxd_MET, j1, j2;
    int j_idx = -1;
    bool j1pass = false;
    bool j2pass = false;
    int atleast_2_jets = 0;
    for (const auto * j : *smallRJets){ // Small-R jets
    //for (const auto * j : *largeRJets){ // Small-R jets
      // small_R jet selections
      j_idx++;
      if (j->pt()*0.001 < 30) continue;
      if (fabs(j->eta()) > 2.8) continue;
      // large_R jet selections
      //if (j->pt()*0.001 < 100) continue;
      //if (fabs(j->eta()) > 2.8) continue;
      //if (j->m()*0.001 < 30 && j->pt()*0.001 < 1000) continue;
      atleast_2_jets++;
      if (j_idx==0) j1pass = true;
      if (j_idx==1) j2pass = true;
      if(j_idx == 0){h_jetPt->Fill( j->pt()*0.001 );
        h_jetEta->Fill( j->eta());
        h_jetPhi->Fill( j->phi());}
      TLorentzVector v_j(0,0,0,0);
      v_j.SetPtEtaPhiE(j->pt()*0.001, j->eta(), j->phi(), j->e()*0.001);

      // Save leading and subleading
      if (j_idx == 0) j1 = v_j;
      if (j_idx == 1) j2 = v_j;

      // Find jet closest to each dark quark
      if (v_j.DeltaR(quarks[0]) < dR_xd1_min){ dR_xd1_min = v_j.DeltaR(quarks[0]); j_xd1_idx = j_idx; j_xd1 = v_j;} 
      if (v_j.DeltaR(quarks[1]) < dR_xd2_min){ dR_xd2_min = v_j.DeltaR(quarks[1]); j_xd2_idx = j_idx; j_xd2 = v_j;} 

      // Find jet (of all jets) closest to MET
      if (fabs(v_j.DeltaPhi(v_met)) < dPhi_any_j_min){ dPhi_any_j_min = fabs(v_j.DeltaPhi(v_met)); j_MET = v_j;}

      // Find jet (of all jets) furthest from MET
      if (fabs(v_j.DeltaPhi(v_met)) > dPhi_any_j_max){ dPhi_any_j_max = fabs(v_j.DeltaPhi(v_met)); j_aMET = v_j;}
 
      // Find dR matched jet closest to MET
      if (fabs(v_j.DeltaPhi(v_met)) < dPhi_matched_j_min && (v_j.DeltaR(quarks[0]) < 0.4 || v_j.DeltaR(quarks[1]) < 0.4)){ dPhi_matched_j_min = fabs(v_j.DeltaPhi(v_met)); jxd_MET = v_j;}

      // Count dR matched jets
      if (v_j.DeltaR(quarks[0]) < 0.4 || v_j.DeltaR(quarks[1]) < 0.4) jets.push_back(v_j);
    }
    if (atleast_2_jets < 2) continue;
    h_nJetsMatched->Fill(jets.size());
    h_dPhi_j_MET->Fill(dPhi_any_j_min);
    if (jets.size() > 0) h_dPhi_xdj_MET->Fill(dPhi_matched_j_min);

    // Calculate mT with leading, subleading
    if (j1pass && j2pass) h_mT_12->Fill( GetMt(j1,j2, (*truthMET)["NonInt"]->met()*0.001, (*truthMET)["NonInt"]->phi()) );

    // Check dR leading, subleading
    h_dRxdj1->Fill( std::min(j1.DeltaR(quarks[0]), j1.DeltaR(quarks[1])) );
    h_dRxdj2->Fill( std::min(j2.DeltaR(quarks[0]), j2.DeltaR(quarks[1])) );

    // Check dR MET aligned, MET anti-aligned
    h_dR_MET->Fill( std::min(j_MET.DeltaR(quarks[0]), j_aMET.DeltaR(quarks[1])) );
    h_dR_aMET->Fill( std::min(j_aMET.DeltaR(quarks[0]), j_aMET.DeltaR(quarks[1])) );

    // Find angle between dark quark and MET;
    if(fabs(v_met.DeltaPhi(quarks[0])) < fabs(v_met.DeltaPhi(quarks[1]))) h_dPhi_xd_MET->Fill(fabs(v_met.DeltaPhi(quarks[0])));
    else h_dPhi_xd_MET->Fill(fabs(v_met.DeltaPhi(quarks[1])));

    if (smallRJets->size() > 1 && j_xd1_idx != j_xd2_idx){
    //if (largeRJets->size() > 1 && j_xd1_idx != j_xd2_idx){
	h_dRxdj->Fill(quarks[0].DeltaR(j_xd1));
	h_dRxdj->Fill(quarks[1].DeltaR(j_xd2));
        h_xdj_idx->Fill(j_xd1_idx);
        h_xdj_idx->Fill(j_xd2_idx);
        if (dR_xd1_min < 0.4) h_xdj_match_idx->Fill(j_xd1_idx);
        if (dR_xd2_min < 0.4) h_xdj_match_idx->Fill(j_xd2_idx);
        if (dR_xd1_min < 0.4 && dR_xd2_min < 0.4){ 
          h_mjj->Fill((j_xd1+j_xd2).M());
          h_mT_jj->Fill( GetMt(j_xd1,j_xd2, (*truthMET)["NonInt"]->met()*0.001, (*truthMET)["NonInt"]->phi()) );
          h_jjdPhi->Fill(fabs(j_xd1.DeltaPhi(j_xd2)));
        }
    }


    //h_nLargeR->Fill(largeRJets->size());
    int nLargeR = 0;
    for (const auto * j : *largeRJets){
      if (j->pt() > 150000 && fabs(j->eta()) < 2.5) nLargeR++;
    }
    h_nLargeR->Fill(nLargeR);

    for (const auto * j : *largeRJets){ // Large-R jets
      h_jetLRPt->Fill( j->pt()*0.001 );
      h_jetLREta->Fill( j->eta());
      h_jetLRPhi->Fill( j->phi());
    }
    


  } // End of event loop
  Info( APP_NAME , "Done with event loop" );

  // Output file  
  TFile * oRoot = new TFile(outputName.c_str(),"RECREATE");
  oRoot->cd();

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //			Write Histogrmas
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // MET histograms
  h_metNonInt->Write();
  h_metInt->Write();
  h_metIntOut->Write();
  h_metIntMuons->Write();
  // Truth particle histograms
  //for (size_t n=0;n<nParticleContainers;++n) h_partPt[n]->Write();
  //for (size_t n=0;n<nParticleContainers;++n) h_partConn[n]->Write();
  h_xdPt->Write();
  h_xdM->Write();
  h_xdEta->Write();
  h_xdPhi->Write();

  h_zpPt->Write();
  h_zpM->Write();
  h_zpEta->Write();
  h_zpPhi->Write();

  h_r_inv->Write();
  h_children->Write();
  h_children_111->Write();
  h_children_113->Write();
  h_children_211->Write();
  h_children_213->Write();
  h_hadrons->Write();

  h_xdxdM->Write();
  h_xdDPhi->Write();
  h_dRxdj->Write();
  h_xdj_idx->Write();
  h_xdj_match_idx->Write();
  h_dPhi_j_MET->Write();
  h_dPhi_xdj_MET->Write();
  h_dPhi_xd_MET->Write();  
  h_nJetsMatched->Write();
  h_dRxdj1->Write();
  h_dRxdj2->Write();
  h_dR_MET->Write();
  h_dR_aMET->Write();
  h_mjj->Write();
  h_mT_jj->Write();
  h_mT_12->Write();

  // Truth jet histograms
  h_jetPt->Write();
  h_jetEta->Write();
  h_jetPhi->Write();
  h_jjdPhi->Write();

  h_jetLRPt->Write();
  h_jetLREta->Write();
  h_jetLRPhi->Write();
  h_nSmallR->Write();
  h_nLargeR->Write();
  // Event histograms
  h_x1->Write();
  h_x2->Write();
  h_xsec->Write();
  h_HTFilt->Write();
  //for (auto * h : h_weights) h->Write();

  // Close up -- all done!
  oRoot->Close();

  // trigger finalization of all services and tools created by the Gaudi Application
#ifndef XAOD_STANDALONE
  app->finalize();
#endif

  Info( APP_NAME , "All done -- goodbye." );

  return 0;
}
