#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/RVec.hxx"

TString era = "2017";
TFile*f=TFile::Open("data/TriggerSF_"+era+"UL.root");
TH2D*h1_ee=(TH2D*)f->Get("h2D_SF_ee_lep1pteta");
TH2D*h2_ee=(TH2D*)f->Get("h2D_SF_ee_lep2pteta");
TH2D*h1_mm=(TH2D*)f->Get("h2D_SF_mumu_lep1pteta");
TH2D*h2_mm=(TH2D*)f->Get("h2D_SF_mumu_lep2pteta");
TH2D*h1_em=(TH2D*)f->Get("h2D_SF_emu_lep1pteta");
TH2D*h2_em=(TH2D*)f->Get("h2D_SF_emu_lep2pteta");

TFile*f_cf=TFile::Open("data/ChargeFlipSF_" + era + "_MLE.root");
TH2F*h_OS=(TH2F*)f_cf->Get("OS_ChargeFlip_SF");
TH2F*h_SS = (TH2F*)f_cf->Get("SS_ChargeFlip_SF_sys"); 

TFile*f_cfregion=TFile::Open("data/ChargeFlipProbability_" + era + "_MLE.root");
TH2F*h_data = (TH2F*) f_cfregion->Get("data_CFRate");

TFile*fele=TFile::Open("data/EleIDSF_" + era + ".root");
TH2D*h_eleSF=(TH2D*)fele->Get("EleIDDataEff");

TFile*fctag=TFile::Open("data/DeepJet_ctagSF_Summer20UL17_interp.root");
TH2D*h_flavc=(TH2D*)fctag->Get("SFc_hist");
TH2D*h_flavb=(TH2D*)fctag->Get("SFb_hist");
TH2D*h_flavl=(TH2D*)fctag->Get("SFl_hist");

float eleID_sf_ee(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>500) l1_pt=499.;
	if(l2_pt>500) l2_pt=499.;
	float sf_l1=h_eleSF->GetBinContent(h_eleSF->FindBin(l1_pt,fabs(l1_eta)));
	float sf_l2=h_eleSF->GetBinContent(h_eleSF->FindBin(l2_pt,fabs(l2_eta)));
	return sf_l1*sf_l2;
}

float trigger_sf_ee(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_ee->GetBinContent(h1_ee->FindBin(l1_pt,fabs(l1_eta)));
	float sf_l2=h2_ee->GetBinContent(h2_ee->FindBin(l2_pt,fabs(l2_eta)));
	return sf_l1;
}

float trigger_sf_mm(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_mm->GetBinContent(h1_mm->FindBin(l1_pt,fabs(l1_eta)));
	float sf_l2=h2_mm->GetBinContent(h2_mm->FindBin(l2_pt,fabs(l2_eta)));
	return sf_l1;
}

float trigger_sf_em(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_em->GetBinContent(h1_em->FindBin(l1_pt,fabs(l1_eta)));
	float sf_l2=h2_em->GetBinContent(h2_em->FindBin(l2_pt,fabs(l2_eta)));
	return sf_l1;
}

int kinematic(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
  std::vector<Float_t> pt_region  = {20., 40., 60., 100., 100000000000.};
  std::vector<Float_t> eta_region = {0.,  0.8, 1.479, 2.4};
  int pt_bins  = pt_region.size() - 1;
  int eta_bins = eta_region.size() - 1;
  int l1_pt_index = -1;
  int l2_pt_index = -1;
  int l1_eta_index = -1;
  int l2_eta_index = -1;
  for(int i = 0; i < pt_bins; i++){
    if(pt_region[i] <= l1_pt && l1_pt < pt_region[i+1]) l1_pt_index = i;
    if(pt_region[i] <= l2_pt && l2_pt < pt_region[i+1]) l2_pt_index = i;
  }
  for(int i = 0; i < eta_bins; i++){
    if(eta_region[i] <= fabs(l1_eta) && fabs(l1_eta) < eta_region[i+1]) l1_eta_index = i;
    if(eta_region[i] <= fabs(l2_eta) && fabs(l2_eta) < eta_region[i+1]) l2_eta_index = i;
  }
  if(l1_pt_index <0 || l2_pt_index <0 || l1_eta_index <0 || l2_eta_index<0) return -1;
  return l2_eta_index + eta_bins*(l2_pt_index) + eta_bins*pt_bins*(l1_eta_index) + eta_bins*pt_bins*eta_bins*(l1_pt_index);
}
     

float chargeflip_sf(int kinematic_region, int charge1, int charge2, float sigma){
    int index = kinematic_region+1;
    float sf = 1.0;
    if((charge1*charge2)<0){
      float OS_SF = h_OS->GetBinContent(index);
      float OS_sigma = h_OS->GetBinError(index);
      sf = OS_SF - sigma*OS_sigma;
    }
    else if((charge1*charge2)>0){
      float SS_SF = h_SS->GetBinContent(index);
      float SS_sigma = h_SS->GetBinError(index);
      sf = SS_SF + sigma*SS_sigma;
    }
    if(sf<0.) sf=0.;
    return sf;
}
//n_tight_jet,tightJets_id_in24,Jet_puId,Jet_hadronFlavour,Jet_btagDeepFlavCvL,Jet_btagDeepFlavCvB
float CtagSF(int n_tight_jet, ROOT::VecOps::RVec<Int_t> tightJets_id_in24,ROOT::VecOps::RVec<Int_t> JetpuId, ROOT::VecOps::RVec<Float_t> Jet_pt, ROOT::VecOps::RVec<Int_t> Jet_hadflav, ROOT::VecOps::RVec<Float_t> CvsL, ROOT::VecOps::RVec<Float_t> CvsB){
  float sf = 1.0;
  for(int i = 0; i < n_tight_jet; i++){
    int index = tightJets_id_in24[i];
    if (JetpuId[index] == 0 && Jet_pt[index] < 50) continue;
    if(Jet_hadflav[index]==4){
      sf *= h_flavc->GetBinContent(h_flavc->FindBin(CvsL[index],CvsB[index]));
    }
    else if(Jet_hadflav[index]==5){
      sf *= h_flavb->GetBinContent(h_flavb->FindBin(CvsL[index],CvsB[index]));
    }
    else{
      sf *= h_flavl->GetBinContent(h_flavl->FindBin(CvsL[index],CvsB[index]));
    }
  }
  return sf;

}
