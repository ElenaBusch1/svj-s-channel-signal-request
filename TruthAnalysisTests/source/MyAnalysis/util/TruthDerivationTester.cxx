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
  const xAOD::TruthParticleContainer * truthBSM(nullptr);
  //for (size_t n=0;n<nParticleContainers;++n) truthParticles[n] = nullptr;

/*
    Missing full collections with aux: TruthBosonsWithDecayParticles, TruthBosonsWithDecayVertices
*/
  // Make histograms
  // MET histograms
  TH1D* h_metNonInt = new TH1D("MET_NonInt","MET; E^{T}_{miss}",50,0,2000.);
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

  TH1D* h_xdxdM = new TH1D("xdxdM", "Dark Quark Invariant Mass; M_{xd,xd} [GeV]",80, 0, 8000);
  TH1D* h_nJetsMatched = new TH1D("nJetsMatched", "N jets satisfying dR(j,xd) < 0.4; nJets", 10, 0, 10);
  TH1D* h_dRxdj = new TH1D("dRxdj", "DeltaR(quark, closest jet); #Delta R", 20,0,5);
  TH1D* h_mjj = new TH1D("mjj", "Invariant Mass 2 Closest Jets; M_{jj} [GeV]", 50, 0, 4000);
  TH1D* h_mT = new TH1D("mT", "System Transverse Mass; m_{T} [GeV]", 50, 0, 4600);
  TH1D* h_xdj_idx = new TH1D("xdj_idx", "Index of Jet Closest To Quark", 15, 0, 15);  
  TH1D* h_dPhi_j_MET = new TH1D("dPhi_j_MET", "#Delta#phi MET and Closest Jet; #Delta#phi(MET,j_{closest})", 20, 0, 3.3);
  TH1D* h_dPhi_xdj_MET = new TH1D("dPhi_xdj_MET", "#Delta#phi MET and Closest dR Matched Jet; #Delta#phi(MET,xdj_{closest})", 20, 0, 3.3);
  TH1D* h_dPhi_xd_MET = new TH1D("dPhi_xd_MET", "#Delta#phi MET and Closest Dark Quark; #Delta#phi(MET,xd_{closest})", 20, 0, 3.3);

  // Small-R jets
  TH1D* h_nSmallR = new TH1D("nSmallR", "nJets (Small-R); N_{jets}", 15, 0, 15);
  TH1D* h_jetPt = new TH1D("JetPt","Jet Pt (Small-R); p_{T} [GeV]",80,0,2000.);
  TH1D* h_jetEta = new TH1D("JetEta","Jet Eta (Small-R); #eta",50,-5.,5.);
  TH1D* h_jetPhi = new TH1D("JetPhi","Jet Phi (Small-R); #phi",50,-3.14,3.14);
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
    CHECK_RETRIEVE( truthBSM, "TruthBSM" )
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
    // Event weight handling
    /*
    if (h_weights.size()==0){
      // First event, init!
      for (auto & name : weightNames ){
        h_weights.push_back( new TH1D(("h_W_"+name).c_str(),"",100,-10.,10.) );
      }
      h_weights.push_back( new TH1D("h_W_nominalTest","",100,-10.,10.) );
    }
    for (size_t n=0;n<weightNames.size();++n) h_weights[n]->Fill( weightTool->getWeight(weightNames[n]) );
    */
    // Eventually this should be the nominal weight without needing to give an explicit name
    //h_weights[weightNames.size()]->Fill( weightTool->getWeight(" nominal ") );
    // Event info
    float x1=0.,x2=0.;
    (*xTruthEventContainer)[0]->pdfInfoParameter( x1 , xAOD::TruthEvent::X1 );
    (*xTruthEventContainer)[0]->pdfInfoParameter( x2 , xAOD::TruthEvent::X2 );
    h_x1->Fill( x1 );
    h_x2->Fill( x2 );
    h_xsec->Fill( (*xTruthEventContainer)[0]->crossSection() );
    h_HTFilt->Fill( xEventInfo->auxdata<float>("GenFiltHT") );
    // For MET: NonInt, Int, IntOut, IntMuons
    h_metNonInt->Fill( (*truthMET)["NonInt"]->met()*0.001 );
    h_metNonInt->Fill( (*truthMET)["Int"]->met()*0.001 );
    h_metNonInt->Fill( (*truthMET)["IntOut"]->met()*0.001 );
    h_metNonInt->Fill( (*truthMET)["IntMuons"]->met()*0.001 );
    // Truth particles
    //for (size_t n=0;n<nParticleContainers;++n){ // PT and connections for all
    //  for (const auto * p : *truthParticles[n]){
    //    h_partPt[n]->Fill( p->pt() );
    //    h_partConn[n]->Fill( countChildren( p ) );
    //    h_partConn[n]->Fill( -1-countParents( p ) );
    //  }
    //}
    std::vector<TLorentzVector> quarks;
    for (const auto * bsm: *truthBSM){
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
    }
    //if(quarks.size() != 2){std::cout << "No quarks found" << std::endl; continue;}
    h_xdxdM->Fill((quarks[0]+quarks[1]).M());
    h_xdDPhi->Fill(fabs(quarks[0].DeltaPhi(quarks[1])));
    h_nSmallR->Fill(smallRJets->size());

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
    float dPhi_matched_j_min = 10;
    TLorentzVector j_xd1, j_xd2, j_MET, jxd_MET;
    int j_idx = -1;
    for (const auto * j : *smallRJets){ // Small-R jets
    //for (const auto * j : *largeRJets){ // Small-R jets
      j_idx++;
      h_jetPt->Fill( j->pt()*0.001 );
      h_jetEta->Fill( j->eta());
      h_jetPhi->Fill( j->phi());
      TLorentzVector v_j(0,0,0,0);
      v_j.SetPtEtaPhiE(j->pt()*0.001, j->eta(), j->phi(), j->e()*0.001);
      // Find jet closest to each dark quark
      if (v_j.DeltaR(quarks[0]) < dR_xd1_min){ dR_xd1_min = v_j.DeltaR(quarks[0]); j_xd1_idx = j_idx; j_xd1 = v_j;} 
      if (v_j.DeltaR(quarks[1]) < dR_xd2_min){ dR_xd2_min = v_j.DeltaR(quarks[1]); j_xd2_idx = j_idx; j_xd2 = v_j;} 

      // Find jet (of all jets) closest to MET
      if (fabs(v_j.DeltaPhi(v_met)) < dPhi_any_j_min){ dPhi_any_j_min = fabs(v_j.DeltaPhi(v_met)); j_MET = v_j;}
      
      // Find dR matched jet closest to MET
      if (fabs(v_j.DeltaPhi(v_met)) < dPhi_matched_j_min && (v_j.DeltaR(quarks[0]) < 0.4 || v_j.DeltaR(quarks[1]) < 0.4)){ dPhi_matched_j_min = fabs(v_j.DeltaPhi(v_met)); jxd_MET = v_j;}

      // Count dR matched jets
      if (v_j.DeltaR(quarks[0]) < 0.4 || v_j.DeltaR(quarks[1]) < 0.4) jets.push_back(v_j);
    }
    h_nJetsMatched->Fill(jets.size());
    h_dPhi_j_MET->Fill(dPhi_any_j_min);
    h_dPhi_xdj_MET->Fill(dPhi_matched_j_min);

    // Find angle between dark quark and MET;
    if(fabs(v_met.DeltaPhi(quarks[0])) < fabs(v_met.DeltaPhi(quarks[1]))) h_dPhi_xd_MET->Fill(fabs(v_met.DeltaPhi(quarks[0])));
    else h_dPhi_xd_MET->Fill(fabs(v_met.DeltaPhi(quarks[1])));


    if (smallRJets->size() > 1 && j_xd1_idx != j_xd2_idx){
    //if (largeRJets->size() > 1 && j_xd1_idx != j_xd2_idx){
	h_dRxdj->Fill(quarks[0].DeltaR(j_xd1));
	h_dRxdj->Fill(quarks[1].DeltaR(j_xd2));
        h_xdj_idx->Fill(j_xd1_idx);
        h_xdj_idx->Fill(j_xd2_idx);
        if (dR_xd1_min < 0.4 && dR_xd2_min < 0.4){ 
          h_mjj->Fill((j_xd1+j_xd2).M());
          h_mT->Fill( GetMt(j_xd1,j_xd2, (*truthMET)["NonInt"]->met()*0.001, (*truthMET)["NonInt"]->phi()) );
          h_jjdPhi->Fill(fabs(j_xd1.DeltaPhi(j_xd2)));
        }
    }


    h_nLargeR->Fill(largeRJets->size());
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

  // Write histograms
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

  h_xdxdM->Write();
  h_xdDPhi->Write();
  h_dRxdj->Write();
  h_xdj_idx->Write();
  h_dPhi_j_MET->Write();
  h_dPhi_xdj_MET->Write();
  h_dPhi_xd_MET->Write();  
  h_nJetsMatched->Write();

  h_mjj->Write();
  h_mT->Write();

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
