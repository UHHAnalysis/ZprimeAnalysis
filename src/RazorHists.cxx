#include "include/RazorHists.h"
#include "include/EventCalc.h"
#include <iostream>
#include "include/SelectionModules.h"
#include "include/TopFitCalc.h"
#include "SFrameTools/include/SubJetTagger.h"
#include "TH3F.h"

using namespace std;

RazorHists::RazorHists(const char* name, HypothesisDiscriminator *discr) : BaseHists(name)
{
  // named default constructor
  m_discr = discr;
   
}

RazorHists::~RazorHists()
{
  // default destructor, does nothing
}

void RazorHists::Init()
{
  // book all histograms here

  Book( TH1F( "mrt_w_tophad", "m_{r,T} with W and Top_{had} ", 100,0,5000) );
  Book( TH1F( "mrt_toplep_tophad", "m_{r,T} with Top_{lep} and Top_{had} ", 100,0,5000) );
  Book( TH1F( "mr_w_tophad", "m_{r} with W and Top_{had} ", 100,0,5000) );
  Book( TH1F( "mr_toplep_tophad", "m_{r} with Top_{lep} and Top_{had} ", 100,0,5000) );
  Book( TH1F( "beta", "#beta^{R}_{L}", 100,0,1));
  Book( TH1F( "beta_ly", "#beta^{R}_{L}", 100,0,1));

  Book( TH1F( "boost_phi_neutrino_blep", "W-Boson CM #Delta #phi (#nu ,b-jet)", 100,0,4));
  Book( TH1F( "boost_delR_neutrino_blep", "W-Boson CM #Delta R(#nu ,b-jet)", 100,0,7));
  Book( TH1F( "boost_delR_muon_nu", "W-Boson CM #Delta #phi(#nu ,muon)", 100,0,7));
  Book( TH1F( "boost_phi_muon_nu", "W-Boson CM #Delta R(#nu ,muon)", 100,0,4));

  Book( TH1F( "topboost_phi_wlep_blep", "W-Boson CM #Delta #phi (#nu ,b-jet)", 100,0,4));
  Book( TH1F( "topboost_delR_wlep_blep", "W-Boson CM #Delta R(#nu ,b-jet)", 100,0,7));
  Book( TH1F( "topboost_delR_muon_nu", "W-Boson CM #Delta #phi(#nu ,muon)", 100,0,7));
  Book( TH1F( "topboost_phi_muon_nu", "W-Boson CM #Delta R(#nu ,muon)", 100,0,4));
  Book( TH1F( "topboost_phi_neutrino_blep", "W-Boson CM #Delta #phi (#nu ,b-jet)", 100,0,4));
  Book( TH1F( "topboost_delR_neutrino_blep", "W-Boson CM #Delta R(#nu ,b-jet)", 100,0,7));

  Book( TH2F("mr_w_tophad_beta" ,"mr_w_tophad beta",100,0,5000,100,0,1));
  Book( TH2F("mr_toplep_tophad_beta" ,"mr_toplep_tophad beta",100,0,5000,100,0,1));

  Book( TH2F("mrt_w_tophad_beta" ,"mr_w_tophad beta",100,0,5000,100,0,1));
  Book( TH2F("mrt_toplep_tophad_beta" ,"",100,0,5000,100,0,1));

  Book( TH2F("mr_w_tophad_mttbar" ,"mr_w_tophad m_{tt#bar}",100,0,5000,100,0,5000));
  Book( TH2F("mr_toplep_tophad_mttbar" ,"mr_toplep_tophad mttbar",100,0,5000,100,0,5000));

  Book( TH2F("mrt_w_tophad_mttbar" ,"mr_w_tophad mttbar",100,0,5000,100,0,5000));
  Book( TH2F("mrt_toplep_tophad_mttbar" ,"mr_toplep_tophad mttbar",100,0,5000,100,0,5000));

  Book( TH2F("beta_mttbar" ,"mr_toplep_tophad mttbar",100,0,1,100,0,5000));

  Book( TH2F("mrt_w_tophad_pTrel" ,"mr_w_tophad beta",100,0,5000,100,0,1));

  Book( TH2F("dalitz_lep" ,"",100,0,400000,100,0,400000));
  //Book( TH2F("dalitz_had" ,"",100,0,1000,100,0,1000));
  Book( TH2F("dalitz_tag" ,"",100,2000,30000,100,0,400000));

}

void RazorHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  TopFitCalc* fitcalc = TopFitCalc::Instance();
  double weight = calc -> GetWeight();
  bool IsRealData = calc->IsRealData();

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  std::vector<TopJet>* cajets = calc->GetCAJets();
  std::vector<Jet>* antikjets = calc->GetJets();
  
  MET* met = calc->GetMET();

  Muon muon  = bcc->muons->at(0);

  Particle neutrino;
  neutrino.set_v4(hyp->neutrino_v4());
  Particle wlep;
  wlep.set_v4(hyp->wlep_v4());
  Particle bjet;
  bjet.set_v4(hyp->blep_v4());

  LorentzVector toplep = hyp->toplep_v4();
  LorentzVector tophad = hyp->tophad_v4();


  TVector3 plepT;
  plepT.SetPtEtaPhi(toplep.pt(),0,toplep.phi());
  TVector3 phadT;
  phadT.SetPtEtaPhi(tophad.pt(),0,tophad.phi());
  TVector3 pwlepT; 
  pwlepT.SetPtEtaPhi(wlep.pt(),0,wlep.phi());
  TVector3 met3Vec; 
  met3Vec.SetPtEtaPhi(met->pt(),0,met->phi());
  TVector3 blepT;
  blepT.SetPtEtaPhi(bjet.pt(),0,bjet.phi());
  TVector3 muonT;
  muonT.SetPtEtaPhi(muon.pt(),0,muon.phi());




  double mttbar_mass = (toplep+tophad).isTimelike() ? (toplep+tophad).M() : -sqrt((toplep+tophad).M2());


  double beta = abs(toplep.Pz()+tophad.Pz())/(toplep.E()+tophad.E());

  ((TH2F*)Hist("beta_mttbar"))->Fill(beta,mttbar_mass,weight);
  Hist("beta_ly")->Fill(beta,weight);
  Hist("beta")->Fill(beta,weight);

  double mrazorT = sqrt(met->pt()*(toplep.pt()+tophad.pt())- met3Vec*(plepT+phadT));
  double mrazor = sqrt(pow(toplep.E()+tophad.E(),2)-pow(toplep.Pz()+tophad.Pz(),2));

  Hist("mrt_toplep_tophad")->Fill(mrazorT,weight);
  Hist("mr_toplep_tophad")->Fill(mrazor,weight);
  
  ((TH2F*)Hist("mr_toplep_tophad_beta"))->Fill(mrazor,beta,weight);
  ((TH2F*)Hist("mrt_toplep_tophad_beta"))->Fill(mrazorT,beta,weight);
  ((TH2F*)Hist("mr_toplep_tophad_mttbar"))->Fill(mrazor,mttbar_mass,weight);
  ((TH2F*)Hist("mrt_toplep_tophad_mttbar"))->Fill(mrazorT,mttbar_mass,weight);
  
  mrazorT = sqrt(met->pt()*(wlep.pt()+tophad.pt())- met3Vec*(pwlepT+phadT));
  mrazor = sqrt(pow(wlep.v4().E()+tophad.E(),2)-pow(wlep.v4().Pz()+tophad.Pz(),2));


  Hist("mrt_w_tophad")->Fill(mrazorT,weight);
  Hist("mr_w_tophad")->Fill(mrazor,weight);

  ((TH2F*)Hist("mr_w_tophad_beta"))->Fill(mrazor,beta,weight);
  ((TH2F*)Hist("mrt_w_tophad_beta"))->Fill(mrazorT,beta,weight);
  ((TH2F*)Hist("mr_w_tophad_mttbar"))->Fill(mrazor,mttbar_mass,weight);
  ((TH2F*)Hist("mrt_w_tophad_mttbar"))->Fill(mrazorT,mttbar_mass,weight);
  ((TH2F*)Hist("mrt_w_tophad_pTrel"))->Fill(mrazor,pTrel(&muon,antikjets),weight);


  //((TH2F*)Hist(""))->Fill(mrazor,,weight);



  ((TH2F*)Hist("dalitz_lep"))->Fill((hyp->wlep_v4()).M2(),(muon.v4()+hyp->blep_v4()).M2(),weight);

  //cout<<(hyp->wlep_v4()).M2()<<" "<<(muon.v4()+hyp->blep_v4()).M2() <<endl;


  CMSTopTagger toptag;
  toptag.SetTau32Cut();

  for(unsigned int i = 0; i<cajets->size(); ++i){
    if(toptag.Tag(cajets->at(i))){
      std::vector<Particle> topsub = cajets->at(i).subjets();

      //cout<<(topsub.at(0).v4()+topsub.at(1).v4()).M2()<<" "<<(topsub.at(1).v4()+topsub.at(2).v4()).M2()<<endl;

      ((TH2F*)Hist("dalitz_tag"))->Fill((topsub.at(0).v4()+topsub.at(1).v4()).M2(),(topsub.at(1).v4()+topsub.at(2).v4()).M2(),weight);
    }
  }


  TLorentzVector neutrino_t;
  neutrino_t.SetPtEtaPhiE(hyp->neutrino_v4().pt(),hyp->neutrino_v4().eta(),hyp->neutrino_v4().phi(),hyp->neutrino_v4().E());
  TLorentzVector wlep_t;
  wlep_t.SetPtEtaPhiE(wlep.pt(),wlep.eta(),wlep.phi(),wlep.v4().E());
  TLorentzVector bjet_t;
  bjet_t.SetPtEtaPhiE(bjet.pt(),bjet.eta(),bjet.phi(),bjet.v4().E());
  TLorentzVector muon_t;
  muon_t.SetPtEtaPhiE(muon.pt(),muon.eta(),muon.phi(),muon.v4().E());



  TVector3 wlep_boost = wlep_t.BoostVector(); 

  neutrino_t.Boost(-wlep_boost);
  bjet_t.Boost(-wlep_boost);
  muon_t.Boost(-wlep_boost);

  double deltaphi = fabs(bjet_t.Phi() - neutrino_t.Phi());
  if(deltaphi > PI) deltaphi = 2* PI - deltaphi;

  double deltaeta = bjet_t.Eta() - neutrino_t.Eta();
  double deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);

  Hist("boost_phi_neutrino_blep" )->Fill(deltaphi,weight);
  Hist("boost_delR_neutrino_blep")->Fill(deltaR,weight);


  deltaphi = fabs(muon_t.Phi() - neutrino_t.Phi());
  if(deltaphi > PI) deltaphi = 2* PI - deltaphi;

  deltaeta = muon_t.Eta() - neutrino_t.Eta();
  deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);

  Hist("boost_delR_muon_nu")->Fill(deltaR,weight);
  Hist("boost_phi_muon_nu")->Fill(deltaphi,weight);



  TLorentzVector topboost;
  topboost.SetPtEtaPhiE(toplep.pt(),toplep.eta(),toplep.phi(),toplep.E());
  TLorentzVector neutrino_top;
  neutrino_top.SetPtEtaPhiE(hyp->neutrino_v4().pt(),hyp->neutrino_v4().eta(),hyp->neutrino_v4().phi(),hyp->neutrino_v4().E());
  TLorentzVector wlep_top;
  wlep_top.SetPtEtaPhiE(wlep.pt(),wlep.eta(),wlep.phi(),wlep.v4().E());
  TLorentzVector bjet_top;
  bjet_top.SetPtEtaPhiE(bjet.pt(),bjet.eta(),bjet.phi(),bjet.v4().E());
  TLorentzVector muon_top;
  muon_top.SetPtEtaPhiE(muon.pt(),muon.eta(),muon.phi(),muon.v4().E());

  TVector3 toplep_boost = topboost.BoostVector(); 
  

  if(toplep_boost.Mag2()>1) return;

  neutrino_top.Boost(-toplep_boost);
  bjet_top.Boost(-toplep_boost);
  muon_top.Boost(-toplep_boost);
  wlep_top.Boost(-toplep_boost);

  

  deltaphi = fabs(muon_top.Phi() - neutrino_top.Phi());
  if(deltaphi > PI) deltaphi = 2* PI - deltaphi;

  deltaeta = muon_top.Eta() - neutrino_top.Eta();
  deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);

  Hist("topboost_delR_muon_nu")->Fill(deltaR,weight);
  Hist("topboost_phi_muon_nu")->Fill(deltaphi,weight);



  deltaphi = fabs(wlep_top.Phi() - bjet_top.Phi());
  if(deltaphi > PI) deltaphi = 2* PI - deltaphi;

  deltaeta = wlep_top.Eta() - bjet_top.Eta();
  deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);

  Hist("topboost_phi_wlep_blep" )->Fill(deltaR,weight);
  Hist("topboost_delR_wlep_blep")->Fill(deltaphi,weight);

 

  deltaphi = fabs(neutrino_top.Phi() - bjet_top.Phi());
  if(deltaphi > PI) deltaphi = 2* PI - deltaphi;

  deltaeta = neutrino_top.Eta() - bjet_top.Eta();
  deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);

  Hist("topboost_phi_neutrino_blep" )->Fill(deltaphi,weight);
  Hist("topboost_delR_neutrino_blep")->Fill(deltaR,weight);


  
}

void RazorHists::Finish()
{
  // final calculations, like division and addition of certain histograms

}

