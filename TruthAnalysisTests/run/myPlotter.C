#include "string.h"
#include "TPave.h"
#include <TArrayD.h>
#include <math.h>

vector<Color_t> mycolors    = { kBlue+1, kAzure+7, /*kGreen+3, kSpring,*/ kRed+1, kOrange-3, kPink+10, kPink+1, kYellow, kYellow-3 };
vector<Color_t> mycolors2   = { kRed, kBlue, kGreen+3, kRed+1, kAzure+7, kPink+10, kYellow+1};
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
  TCanvas *canv = new TCanvas("canv","canv",1400,1200);
  TLegend *leg;
  //if(var.find("dPhi") != string::npos) leg = new TLegend(0.1,0.1,0.4,0.4);
  //else leg = new TLegend(0.65,0.7,0.95,0.95);
  leg = new TLegend(0.65,0.7,0.95,0.95);
  float max = 0;
  gStyle->SetOptStat(0);

  for(int i=0; i<hists.size();i++){
    //hists[i]->SetLineColor(rainbow[i]);
    hists[i]->SetLineColor(mycolors[i]);
    hists[i]->SetLineWidth(3);
    //hists[i]->Scale(1./hists[i]->Integral());
    //hists[i]->GetYaxis()->SetTitle("A.U.");
    if (hists[i]->GetMaximum() > max) max = hists[i]->GetMaximum();
    leg->AddEntry(hists[i], Form("%s",legend_names[i].c_str()));
    //leg->AddEntry(hists[i], Form("%s (%i) ",legend_names[i].c_str(), int(hists[i]->GetEntries())));
    //leg->AddEntry(hists[i], Form("%s; mean = %.2f",legend_names[i].c_str(), hists[i]->GetMean()));
    if (var.find("dR") != string::npos) cout << var << " " << legend_names[i] << " percent dR matched: " << hists[i]->Integral(1,3) << endl;
  }
  for(int i=0; i<hists.size();i++){
    hists[i]->SetMaximum(max*1.2);
    hists[i]->Draw("same hist");
  }
  if (var != "nSmallR" && var != "nLargeR"  && var != "nJetsMatched") canv->SetLogy();
  leg->Draw();
  canv->SaveAs(Form("%s/%s_validation.png", save_dir.c_str(), var.c_str()));
  delete canv;
} 

void myPlotter(){


  vector<string> my_vars = {
		"JetPt", //"JetEta", "JetPhi",
		//"JetLRPt", "JetLREta", "JetLRPhi", 
		"MET_NonInt", 
		"nSmallR", //"nLargeR"
		"xdPt", //"xdM", "xdPhi", "xdEta",
		"xdxdM",
		"xdDPhi", "jjDPhi",
		//"zpPt", "zpM", "zpPhi", "zpEta"
		"mT_jj", "mT_12", "mjj",
                //"dR_MET", "dR_aMET",
		//"dRxdj1", "dRxdj2",
                "nJetsMatched",
                "r_inv", "hadrons",
		"dPhi_xdj_MET", "dPhi_xd_MET", "dPhi_j_MET",
		"xdj_idx", "xdj_match_idx" //"dRxdj", "nJetsMatched"

		};


  //vector<string> input_files = {"hists_output_750_2_ej_pt25eta25.root", "hists_output_750_2_nej_pt25eta25.root"};
  
  vector<string> input_files = {
		//"hists_output_750_2_LR.root",
		//"hists_output_750_8_LR.root",
		//"hists_output_1500_2_zp.root",
		//"hists_output_2500_2_zp.root",
		//"hists_output_3500_2_LR.root",
		//"hists_output_3500_8_LR.root",
		//"hists_output_4000_2.root",
		//"hists_output_4500_2_zp.root",
		//"hists_output_5000_2.root",
		//"hists_output_5500_2_zp.root",
		//"hists_output_6000_2_LR.root",
		//"hists_output_6000_8_LR.root"
		//"hists_output_6500_2_zp.root"
		//"hists_output_7000_2.root"
		//"hists_output_7000_8.root"
		"hists_snowmass58_750_2.root",
		//"hists_cms_750_2.root",
		"hists_snowmass58_750_8.root",
		//"hists_cms_750_8.root",
		"hists_snowmass58_4000_2.root",
		//"hists_cms_4000_2.root"
		"hists_snowmass58_4000_8.root",
		};
  

  vector<map<string,TH1D*>> all_hists;
  for (auto & f: input_files){
    map<string,TH1D*> hist_map = GetHistograms(my_vars, Form("%s", f.c_str()));
    all_hists.push_back(hist_map);
  }
  //map<string,TH1D*> extra_jets = GetHistograms(my_vars, "test_with_extra_jets.root");
  //map<string,TH1D*> normal_jets = GetHistograms(my_vars, "test_no_extra_jets.root");

  //vector<string> legend_names = {"Extra Jets", "No Extra Jets"};
  
  vector<string> legend_names = {
		"M_{Z'}=750 GeV | r_{inv}=0.2",
		"M_{Z'}=750 GeV | r_{inv}=0.8",
		//"M_{Z'}=1500 GeV | r_{inv}=0.2",
		//"M_{Z'}=1500 GeV | r_{inv}=0.8",
		//"M_{Z'}=2500 GeV | r_{inv}=0.2",
		//"M_{Z'}=2500 GeV | r_{inv}=0.8",
		//"M_{Z'}=3500 GeV | r_{inv}=0.2",
		//"M_{Z'}=3500 GeV | r_{inv}=0.8",
		"M_{Z'}=4000 GeV | r_{inv}=0.2",
		"M_{Z'}=4000 GeV | r_{inv}=0.8",
		//"M_{Z'}=4500 GeV | r_{inv}=0.2",
		//"M_{Z'}=5000 GeV | r_{inv}=0.2",
		//"M_{Z'}=5500 GeV | r_{inv}=0.2",
		//"M_{Z'}=6000 GeV | r_{inv}=0.2",
		//"M_{Z'}=6000 GeV | r_{inv}=0.8"
		//"M_{Z'}=6500 GeV | r_{inv}=0.2"
		//"M_{Z'}=7000 GeV | r_{inv}=0.2"
		//"M_{Z'}=7000 GeV | r_{inv}=0.8"
		};
  

  for(auto &var: my_vars){
    vector<TH1D*> plot_hists;
    for(auto &mp: all_hists){
      plot_hists.push_back(mp[var]);
    }
    //plot_hists.push_back(extra_jets[var]);
    //plot_hists.push_back(normal_jets[var]);
    Plot_Histograms(var, plot_hists, legend_names);
  }


}
