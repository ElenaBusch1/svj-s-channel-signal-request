#include "string.h"
#include "TPave.h"
#include <TArrayD.h>
#include <math.h>

vector<Color_t> mycolors    = { kBlue+1, kAzure+7, kGreen+3, kSpring, kRed+1, kOrange-3, kPink+10, kPink+1, kYellow, kYellow-3 };
vector<Color_t> mycolors2   = { kBlue+1, kOrange+7, kGreen+3, kRed+1, kAzure+7, kPink+10, kYellow+1};
vector<Color_t> rainbow    = { kRed, kOrange-3, kYellow+1, kSpring-1, kAzure, kBlue+2, kViolet};
string save_dir = "/eos/user/e/ebusch/SVJ/Plots";


map<string,TH1D*> GetHistograms(vector<string> variables, string input_file){
  map<string,TH1D*> hist_map;
  TFile *file_data = new TFile(Form("%s",input_file.c_str()), "READ"); 
  for(auto &var: variables){
	hist_map[var] = (TH1D*) file_data->Get(Form("%s",var.c_str()));
  }
  return hist_map;
}

void Plot_Histograms(string var, vector<TH1D*> hists, vector<string> legend_names){
  TCanvas *canv = new TCanvas("canv","canv",1600,1200);
  TLegend *leg;
  if(var.find("DPhi") != string::npos) leg = new TLegend(0.05,0.6,0.35,0.94);
  else leg = new TLegend(0.65,0.6,0.96,0.94);
  float max = 0;
  gStyle->SetOptStat(0);

  for(int i=0; i<hists.size();i++){
    //hists[i]->SetLineColor(rainbow[i]);
    hists[i]->SetLineColor(mycolors2[i]);
    hists[i]->SetLineWidth(2);
    hists[i]->Scale(1./hists[i]->Integral());
    hists[i]->GetYaxis()->SetTitle("A.U.");
    if (hists[i]->GetMaximum() > max) max = hists[i]->GetMaximum();
    //leg->AddEntry(hists[i], Form("%s (%i) ",legend_names[i].c_str(), int(hists[i]->GetEntries())));
    leg->AddEntry(hists[i], Form("%s  ",legend_names[i].c_str()));
  }
  for(int i=0; i<hists.size();i++){
    hists[i]->SetMaximum(max*1.2);
    hists[i]->Draw("same hist");
  }
  if (var != "nSmallR" && var != "nLargeR"  && var != "nJetsMatched") canv->SetLogy();
  leg->Draw();
  canv->SaveAs(Form("%s/%s_xJ_pt100eta25_750_2.png", save_dir.c_str(), var.c_str()));
  delete canv;
} 

void myPlotter(){


  vector<string> my_vars = {
		//"JetPt", "JetEta", "JetPhi",
		//"JetLRPt", "JetLREta", "JetLRPhi", 
		//"MET_NonInt", 
		"nSmallR", "nLargeR"
		//"xdPt", "xdM", "xdPhi", "xdEta",
		//"xdxdM",
		//"xdDPhi", "jjDPhi"
		//"zpPt", "zpM", "zpPhi", "zpEta"
		//"mjj", "mT", 
		//"dPhi_xdj_MET", "dPhi_xd_MET", "dPhi_j_MET",
		//"dRxdj", "xdj_idx", "nJetsMatched"
		};


  vector<string> input_files = {"hists_output_750_2_ej_pt25eta25.root", "hists_output_750_2_nej_pt25eta25.root"};
  /*
  vector<string> input_files = {
		"hists_output_750_2_LR71.root",
		"hists_output_750_8_LR71.root",
		//"hists_output_1500_2_zp.root",
		//"hists_output_2500_2_zp.root",
		"hists_output_3500_2_LR71.root",
		"hists_output_3500_8_LR71.root",
		//"hists_output_4000_2.root",
		//"hists_output_4500_2_zp.root",
		//"hists_output_5000_2.root",
		//"hists_output_5500_2_zp.root",
		"hists_output_6000_2_LR71.root",
		"hists_output_6000_8_LR71.root"
		//"hists_output_6500_2_zp.root"
		//"hists_output_7000_2.root"
		//"hists_output_7000_8.root"
		};
  */

  vector<map<string,TH1D*>> all_hists;
  for (auto & f: input_files){
    map<string,TH1D*> hist_map = GetHistograms(my_vars, Form("%s", f.c_str()));
    all_hists.push_back(hist_map);
  }
  //map<string,TH1D*> extra_jets = GetHistograms(my_vars, "test_with_extra_jets.root");
  //map<string,TH1D*> normal_jets = GetHistograms(my_vars, "test_no_extra_jets.root");

  vector<string> legend_names = {"Extra Jets", "No Extra Jets"};
  /*
  vector<string> legend_names = {
		"M_{Z'}=750 GeV | r_{inv}=0.2",
		"M_{Z'}=750 GeV | r_{inv}=0.8",
		//"M_{Z'}=1500 GeV | r_{inv}=0.2",
		//"M_{Z'}=1500 GeV | r_{inv}=0.8",
		//"M_{Z'}=2500 GeV | r_{inv}=0.2",
		//"M_{Z'}=2500 GeV | r_{inv}=0.8",
		"M_{Z'}=3500 GeV | r_{inv}=0.2",
		"M_{Z'}=3500 GeV | r_{inv}=0.8",
		//"M_{Z'}=4000 GeV | r_{inv}=0.2",
		//"M_{Z'}=4500 GeV | r_{inv}=0.2",
		//"M_{Z'}=5000 GeV | r_{inv}=0.2",
		//"M_{Z'}=5500 GeV | r_{inv}=0.2",
		"M_{Z'}=6000 GeV | r_{inv}=0.2",
		"M_{Z'}=6000 GeV | r_{inv}=0.8"
		//"M_{Z'}=6500 GeV | r_{inv}=0.2"
		//"M_{Z'}=7000 GeV | r_{inv}=0.2"
		//"M_{Z'}=7000 GeV | r_{inv}=0.8"
		};
  */

  for(auto &var: my_vars){
    vector<TH1D*> plot_hists;
    for(auto &mp: all_hists){
      plot_hists.push_back(mp[var]);
    }
    //plot_hists.push_back(extra_jets[var]);
    //plot_hists.push_back(normal_jets[var]);
    Plot_Histograms(var, plot_hists, legend_names);
  }

  /*
  TFile *file_data_j  = new TFile("test_with_extra_jets.root", "READ");

  TH1D* pt_j = (TH1D*) file_data_j->Get("JetPt");
  TH1D* eta_j = (TH1D*) file_data_j->Get("JetEta");

  TFile *file_data  = new TFile("test_no_extra_jets.root", "READ");
  TH1D* pt = (TH1D*) file_data->Get("JetPt");
  TH1D* eta = (TH1D*) file_data->Get("JetEta");

  
  TCanvas *c[10];
  c[0] = new TCanvas("c0","c0",1600,1200);
  pt->Draw();
  pt_j->SetLineColor(kRed);
  pt_j->Draw("same");
  c[0]->SaveAs(Form("%s/pT.png", save_dir.c_str()));

  //MET
  c[0] = new TCanvas("c0","c0",1600,1200);
  TH1F* h_met = new TH1F("h_met", "Truth MET; MET [GeV]; nEvents", 50, 0, 2000);
  tree0->Draw("MET_TruthAuxDyn.sumet/1000.>>h_met");
  c[0]->SetLogy();
  c[0]->SaveAs(Form("%s/MET.png", save_dir.c_str()));

  //Z'
  c[1] = new TCanvas("c1","c1",1600,1200);
  TH1F* h_mz = new TH1F("h_mz", "Mass Z' (5000001); M_{Z'} [GeV]; nEvents", 50, 600, 2200);
  tree0->Draw("TruthBSMWithDecayParticlesAuxDyn.m/1000.>>h_mz", "TruthBSMWithDecayParticlesAuxDyn.pdgId==5000001");
  c[1]->SaveAs(Form("%s/m_zp.png", save_dir.c_str()));

  //xd
  c[2] = new TCanvas("c2","c2",1600,1200);
  TH1F* h_xd1_m = new TH1F("h_xd1_m", "Dark Quark Mass (4900101); M_{xd} [GeV]; nEvents", 50, 0, 1000);
  tree0->Draw("TruthBSMWithDecayParticlesAuxDyn.m/1000.>>h_xd1_m", "TruthBSMWithDecayParticlesAuxDyn.pdgId==51 || TruthBSMWithDecayParticlesAuxDyn.pdgId==53");
  c[2]->SaveAs(Form("%s/xd1_m.png", save_dir.c_str()));

  //jets
  c[3] = new TCanvas("c3","c3",1600,1200);
  TH1F* h_pt_j1 = new TH1F("h_pt_j1", "Jet pT; p_{T,j} [GeV]; nEvents", 50, 0, 1700);
  tree0->Draw("AntiKt10TruthTrimmedPtFrac5SmallR20JetsAux.pt/1000.>>h_pt_j1");
  c[3]->SaveAs(Form("%s/pt_j.png", save_dir.c_str()));  

  // jet multiplicity
  c[4] = new TCanvas("c4","c4",1600,1200);
  TH1F* h_njet = new TH1F("h_njet", "Jet Multiplicity; n_{jets}; nEvents", 20, 0, 20); 
  const xAOD::JetContainer* jets = nullptr; 
  */
 
}
