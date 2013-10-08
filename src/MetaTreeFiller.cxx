#include "include/MetaTreeFiller.h"
#include "include/EventCalc.h"
#include "include/TopFitCalc.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include "TH3F.h"



using namespace std;


float IsoMuonTree( Muon mu, float deltaR ){
  float chargedHadronIso=0;
  float neutralHadronIso=0;
  float photonIso=0;
  float puiso=0;

  EventCalc* calc = EventCalc::Instance();

  unsigned int npfp=calc->GetIsoPFParticles()->size();
  for(unsigned int j=0; j<npfp; ++j){
    PFParticle pfp = calc->GetIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.0 && pfp.deltaR(mu)>0.0001 ) chargedHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eH0 && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) neutralHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eGamma && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) photonIso += pfp.pt();
    }
  }
  
  unsigned int npfppu=calc->GetPUIsoPFParticles()->size();
  for(unsigned int j=0; j<npfppu; ++j){
    PFParticle pfp = calc->GetPUIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01 ) puiso += pfp.pt();
    }
  }
  
  return (chargedHadronIso + std::max( 0.0, neutralHadronIso + photonIso - 0.5*puiso));
};

float relIsoParticleTree( Particle mu, float deltaR ){
  float chargedHadronIso=0;
  float neutralHadronIso=0;
  float photonIso=0;
  float puiso=0;

  EventCalc* calc = EventCalc::Instance();

  unsigned int npfp=calc->GetIsoPFParticles()->size();
  for(unsigned int j=0; j<npfp; ++j){
    PFParticle pfp = calc->GetIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.0 && pfp.deltaR(mu)>0.0001 ) chargedHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eH0 && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) neutralHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eGamma && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) photonIso += pfp.pt();
    }
  }
  
  unsigned int npfppu=calc->GetPUIsoPFParticles()->size();
  for(unsigned int j=0; j<npfppu; ++j){
    PFParticle pfp = calc->GetPUIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01 ) puiso += pfp.pt();
    }
  }
  
  return (chargedHadronIso + std::max( 0.0, neutralHadronIso + photonIso - 0.5*puiso))/mu.pt();
};








MetaTreeFiller::MetaTreeFiller(const char* name,  HypothesisDiscriminator *discr, TTree* MetaOutputTree) : BaseHists(name), m_tree(MetaOutputTree)
{
  // named default constructor
  m_discr = discr;
}

MetaTreeFiller::~MetaTreeFiller()
{
  // default destructor, does nothing
}

void MetaTreeFiller::Init()
{
  
  //book output tree branches
  
  //m_tree->Branch("data");
  
  m_tree->Branch("weight"    , &m_weight     , "weight/D");
  m_tree->Branch("pT_mu"     , &m_pT_mu      , "pT_mu/D");
  m_tree->Branch("HT"        , &m_HT         , "HT/D");
  m_tree->Branch("2D"        , &m_2D         , "2D/D");
  m_tree->Branch("mtt"       , &m_mtt        , "mtt/D");
  m_tree->Branch("met_pt"    , &m_met_pt     , "met_pt/D");
  m_tree->Branch("Jet_pt_max", &m_Jet_pt_max , "Jet_pt_max/D");
  m_tree->Branch("btag"      , &m_btag       , "btag/D");


  m_tree->Branch("Iso02"   , &m_Iso02   , "Iso02/D");
  m_tree->Branch("Iso018"  , &m_Iso018  , "Iso018/D");
  m_tree->Branch("Iso016"  , &m_Iso016  , "Iso016/D");
  m_tree->Branch("Iso014"  , &m_Iso014  , "Iso014/D");
  m_tree->Branch("Iso012"  , &m_Iso012  , "Iso012/D");
  m_tree->Branch("Iso01"   , &m_Iso01   , "Iso01/D");
  m_tree->Branch("Iso008"  , &m_Iso008  , "Iso008/D");
  m_tree->Branch("Iso006"  , &m_Iso006  , "Iso006/D");
  m_tree->Branch("Iso004"  , &m_Iso004  , "Iso004/D");
  m_tree->Branch("Iso002"  , &m_Iso002  , "Iso002/D");

  m_tree->Branch("Iso04"   , &m_Iso04   , "Iso04/D");
  m_tree->Branch("Iso05"   , &m_Iso05   , "Iso05/D");
}

void MetaTreeFiller::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  TopFitCalc* fitcalc = TopFitCalc::Instance();


  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  MET* met = calc->GetMET();

  LuminosityHandler* lumih = calc->GetLumiHandler();

  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  // important: get the event weight
  double weight = calc->GetWeight();

  std::vector<Muon>* muons = calc->GetMuons();
  
  std::vector<TopJet>* cajets = calc->GetCAJets();
  std::vector<Jet>* antikjets = calc->GetJets();

  /*
  cout<<1<<endl;
  for(int i=0; i<antikjets->size();++i)
    cout<<antikjets->at(i).pt()<<endl;
  */

  Particle neutrino; 
  neutrino.set_v4(hyp->neutrino_v4());

  Particle wlep;
  wlep.set_v4(hyp->wlep_v4());

  LorentzVector toplep = hyp->toplep_v4();
  LorentzVector tophad = hyp->tophad_v4();

  double HT = calc->GetHT();
  Muon muon = muons->at(0);


  std::vector< PFParticle >* IsoPFParticle = calc->GetIsoPFParticles();

  int Nmuons = muons->size();

  m_weight= weight;
  m_pT_mu = muon.pt(); 
  m_HT    = HT;

  if(deltaRmin(&muon, antikjets)<0.5 && pTrel(&muon, antikjets)<25)
    m_2D = 0.;
  else 
    m_2D = 1.;

  m_mtt =  (toplep+tophad).isTimelike()?(toplep+tophad).M(): -sqrt((toplep+tophad).M2()) ;
 
  m_met_pt = met->pt();
  m_Jet_pt_max = antikjets->at(0).pt();

  m_Iso05  = IsoMuonTree(muon,0.5); 
  m_Iso04  = IsoMuonTree(muon,0.4);  
 
  m_Iso02  = IsoMuonTree(muon,0.2); 
  m_Iso018 = IsoMuonTree(muon,0.18); 
  m_Iso016 = IsoMuonTree(muon,0.16); 
  m_Iso014 = IsoMuonTree(muon,0.14); 
  m_Iso012 = IsoMuonTree(muon,0.12);  
  m_Iso01  = IsoMuonTree(muon,0.1);  
  m_Iso008 = IsoMuonTree(muon,0.08); 
  m_Iso006 = IsoMuonTree(muon,0.06); 
  m_Iso004 = IsoMuonTree(muon,0.04); 
  m_Iso002 = IsoMuonTree(muon,0.02);  
 

  m_btag = 0; 

  for(int i= 0; i<antikjets->size(); ++i)
    if(IsTagged(antikjets->at(i),e_CSVT)) m_btag += 1; 

  m_tree->Fill();
  

 
}

void MetaTreeFiller::Finish()
{
  // final calculations, like division and addition of certain histograms

}

