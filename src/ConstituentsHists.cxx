#include "include/ConstituentsHists.h"
#include "include/EventCalc.h"
#include "include/TopFitCalc.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include "TH3F.h"



using namespace std;


float IsoMuon( Muon mu, float deltaR ){
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

float relIsoParticle( Particle mu, float deltaR ){
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








ConstituentsHists::ConstituentsHists(const char* name,  HypothesisDiscriminator *discr) : BaseHists(name)
{
  // named default constructor
  m_discr = discr;
}

ConstituentsHists::~ConstituentsHists()
{
  // default destructor, does nothing
}

void ConstituentsHists::Init()
{
  


  // book all histograms here

  // leptons
  Book( TH1F("mu_size"  ,"Number of #mu" ,5,-.5,4.5));
  Book( TH1F("mu_pt_ly"  ," #mu p_{T}" ,50,0,900));
  Book( TH1F("mu_E"  ," #mu Energy" ,50,0,1000));
  Book( TH1F("mu_eta"  ," #eta #mu" ,50,-4,4));
  Book( TH1F("mu_phi"  ," #phi #mu" ,50,-M_PI,M_PI));
  Book( TH1F("mu_reliso"  ," relIso #mu" ,100,0,5));

  // jets
  Book( TH1F("mu_nextJet_pt_ly"  ," p_{T} nextJet to #mu" ,50,0,600));
  Book( TH1F("mu_nextJet_phi"  ,"#phi nextJet to #mu" ,50,-M_PI,M_PI));
  Book( TH1F("mu_nextJet_eta"  ,"#eta nextJet to #mu" ,50,-4,4));
  Book( TH1F("mu_nextJet_delR"  ,"#Delta R nextJet to #mu " ,50,0,6));

  Book( TH1F("mu_nextJet_area_num"  ,"Area/Number of Constituent fom nextJet to #mu " ,100,0,.2));
  Book( TH1F("neutrino_nextJet_area_num"  ,"Area/Number of Constituent fom nextJet to #nu " ,100,0,.2));

  // constituents
  Book( TH1F("con_size_ly"  ,"Number of Constituents" ,100,-1,100));
  Book( TH1F("con_pt_ly"  ,"Constituents  p_{T}" ,50,0,100));
  Book( TH1F("con_eta"  ," #eta Constituents" ,50,-4,4));
  Book( TH1F("con_phi"  ," #phi Constituents" ,50,-M_PI,M_PI));

  Book( TH1F("con_muondelphi"  ,"#Delta #phi muon-Constituents" ,100,0,M_PI/2));
  Book( TH1F("con_muondelR"  ," #Delta R muon-Constituents" ,100,0,2.5));

  Book( TH1F("con_neutrinodelphi"  ,"#Delta #phi #nu -Constituents" ,100,0,M_PI));
  Book( TH1F("con_neutrinodelR_ly"  ," #Delta R #nu -Constituents" ,100,0,5));
  Book( TH1F("filter"  ," ",100,40,900));

  
  Book( TH2F("pf_con_pt_mu_delR","",100,0,150,100,0,5));
  Book( TH2F("pf_con_pt_nu_delR","",100,0,150,100,0,5));


  Book( TH1F("nearcon","#Delta R(#mu,PF_{nearest}) ",100,0,0.4));
  Book( TH1F("nearcon_nomuon","#Delta R(#mu,PF_{nearest}) ",100,0,0.4));

  Book( TH1F("delR_sqrtE","#Delta R(#mu,PF_{con}})/#sqrt{E} ", 100,0,10));
  Book( TH1F("delR_E","#Delta R(#mu,PF_{con}})/E ", 100,0,10));
  Book( TH1F("delR_E_pf","#Delta R(#mu,PF_{con}})/E_{PFC} ", 100,0,10));
  Book( TH1F("delR_sqrtE_pf","#Delta R(#mu,PF_{con}})/E_{PFC} ", 100,0,10));

  Book( TH1F("delR_sqrtE_nomuon","#Delta R(#mu,PF_{con}})/#sqrt{E} Second ", 100,0,1));
  Book( TH1F("delR_E_nomuon","#Delta R(#mu,PF_{con}})/E Second", 100,0,1));
  Book( TH1F("delR_E_pf_nomuon","#Delta R(#mu,PF_{con}})/E_{PFC} Second", 100,0,1));
  Book( TH1F("delR_sqrtE_pf_nomuon","#Delta R(#mu,PF_{con}})/E_{PFC} Second", 100,0,1));

  Book( TH1F("relIsoMuon","relIsoMuon",100,0,5));
  Book( TH1F("relIsoMuon_sqrtE","relIsoMuon/#sqrt{E}",100,0,2));

  Book( TH2F("IsoMuon_2D","",100,0,500,100,0,1));
  Book( TH2F("relIsoMuon_2D","",100,0,10,100,0,1));
  Book( TH2F("delR_E_mu","",100,0,.01,100,0,600));
  Book( TH2F("delR_E_pfcon","",100,0,.01,100,0,100));
  Book( TH2F("delR_sqrtE_mu","",100,0,.01,100,0,600));
  Book( TH2F("delR_sqrtE_pfcon","",100,0,.01,100,0,10));

  Book( TH2F("IsoMuon_sqrtE","",100,0,1,100,0,1));
  Book( TH2F("IsoMuon_E","",100,0,1,100,0,1));   
  Book( TH2F("IsoMuon_sqrt","",100,0,1,100,0,1));   

  Book( TH2F("IsoMuon02_pt","",100,0,500,100,0,500)); 
  Book( TH2F("IsoMuon04_pt","",100,0,500,100,0,500)); 


  Book( TH1F("relIsoMuon_004","relIsoMuon #Delta R 0.04",100,0,5));
  Book( TH1F("relIsoMuon_05" ,"relIsoMuon #Delta R 0.5",100,0,5));
  Book( TH1F("relIsoMuon_01" ,"relIsoMuon #Delta R 0.1",100,0,5));
  Book( TH1F("relIsoMuon_008","relIsoMuon #Delta R 0.08",100,0,5));
  Book( TH1F("relIsoMuon_002","relIsoMuon #Delta R 0.02",100,0,5));

  Book( TH1F("relIsoSub_004","relIsoSub #Delta R 0.04",100,0,5));
  Book( TH1F("relIsoSub_05" ,"relIsoSub #Delta R 0.5",100,0,5));
  Book( TH1F("relIsoSub_01" ,"relIsoSub #Delta R 0.1",100,0,5));
  Book( TH1F("relIsoSub_008","relIsoSub #Delta R 0.08",100,0,5));
  Book( TH1F("relIsoSub_002","relIsoSub #Delta R 0.02",100,0,5));

  Book( TH1F("HT","HT",100,0,3000));


  Book( TH2F("relIsoMuon_04_pt","relIsoMuon #Delta R 0.04",100,0,5,100,0,500));
  Book( TH2F("relIsoMuon_05_pt" ,"relIsoMuon #Delta R 0.5",100,0,5,100,0,500));
  Book( TH2F("relIsoMuon_01_pt" ,"relIsoMuon #Delta R 0.1",100,0,5,100,0,500));
  Book( TH2F("relIsoMuon_008_pt","relIsoMuon #Delta R 0.08",100,0,5,100,0,500));
  Book( TH2F("relIsoMuon_002_pt","relIsoMuon #Delta R 0.02",100,0,5,100,0,500));

  Book( TH2F("IsoMuonHT_04_HT" ,"IsoMuonHT #Delta R 0.04",100,0,.5,100,350,5000));
  Book( TH2F("IsoMuonHT_05_HT" ,"IsoMuonHT #Delta R 0.5",100,0,.5,100,350,5000));
  Book( TH2F("IsoMuonHT_01_HT" ,"IsoMuonHT #Delta R 0.1",100,0,.5,100,350,5000));
  Book( TH2F("IsoMuonHT_008_HT","IsoMuonHT #Delta R 0.08",100,0,.5,100,350,5000));
  Book( TH2F("IsoMuonHT_002_HT","IsoMuonHT #Delta R 0.02",100,0,.5,100,350,5000));

  Book( TH2F("relIsoMuon_04_HT" ,"IsoMuonHT #Delta R 0.04",100,0,5,100,350,5000));
  Book( TH2F("relIsoMuon_05_HT" ,"IsoMuonHT #Delta R 0.5",100,0,5,100,350,5000));
  Book( TH2F("relIsoMuon_01_HT" ,"IsoMuonHT #Delta R 0.1",100,0,5,100,350,5000));
  Book( TH2F("relIsoMuon_008_HT","IsoMuonHT #Delta R 0.08",100,0,5,100,350,5000));
  Book( TH2F("relIsoMuon_002_HT","IsoMuonHT #Delta R 0.02",100,0,5,100,350,5000));

  //Book( TH3F("relIsoMuon_3D","",100,0,5,100,0,1,100,0,500));


}

void ConstituentsHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  TopFitCalc* fitcalc = TopFitCalc::Instance();


  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  LuminosityHandler* lumih = calc->GetLumiHandler();

  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  // important: get the event weight
  double weight = calc->GetWeight();

  std::vector<Muon>* muons = calc->GetMuons();
  
  std::vector<TopJet>* cajets = calc->GetCAJets();
  std::vector<Jet>* antikjets = calc->GetJets();

  Particle neutrino; 
  neutrino.set_v4(hyp->neutrino_v4());

  Particle wlep;
  wlep.set_v4(hyp->wlep_v4());

  double HT = calc->GetHT();

  Hist("HT")->Fill(HT,weight);

  std::vector< PFParticle >* IsoPFParticle = calc->GetIsoPFParticles();

  int Nmuons = muons->size();
  Hist("mu_size")->Fill(Nmuons, weight);
 
  if(cajets->size()==0) return;
  
  if(Nmuons>1 ||Nmuons<1) cout<< "not exactly one muon"<<endl;

  for (int i=0; i<Nmuons; ++i){
    Muon muon = muons->at(i);

    Hist("mu_E")->Fill(muon.v4().E(), weight);

    Jet* Mu_netxtJet = nextJet(&muon,antikjets);
    Jet* Nu_netxtJet = nextJet(&neutrino,antikjets);
    
 

    Hist("con_size_ly")->Fill(IsoPFParticle->size(), weight);
    Hist("mu_pt_ly")->Fill(muon.pt(), weight);
    Hist("mu_eta")->Fill(muon.eta(), weight);
    Hist("mu_phi")->Fill(muon.phi(), weight);
    Hist("mu_reliso")->Fill(muon.relIso(), weight);
    Hist("mu_phi")->Fill(muon.phi(), weight);

    Hist("mu_nextJet_pt_ly")->Fill(nextJet(&muon,antikjets)->pt(),weight);
    Hist("mu_nextJet_phi")->Fill(nextJet(&muon,antikjets)->phi(),weight);
    Hist("mu_nextJet_eta")->Fill(nextJet(&muon,antikjets)->eta(),weight);
    Hist("mu_nextJet_delR")->Fill(deltaRmin(&muon,antikjets),weight);


    Hist("mu_nextJet_area_num")->Fill(Mu_netxtJet->jetArea()/IsoPFParticle->size(),weight);
    Hist("neutrino_nextJet_area_num")->Fill(Nu_netxtJet->jetArea()/IsoPFParticle->size(),weight);

    double count =0;	

    for(unsigned int u = 0; u<IsoPFParticle->size();++u){
      Particle pfcon = IsoPFParticle->at(u);
      
      if(pfcon.deltaR(muon)<.5) count+=1;


      Hist("con_muondelR")->Fill(pfcon.deltaR(muon));
      Hist("con_neutrinodelR_ly")->Fill(pfcon.deltaR(neutrino));
      Hist("con_muondelphi")->Fill(pfcon.deltaPhi(muon));
      Hist("con_neutrinodelphi")->Fill(pfcon.deltaPhi(neutrino));
       
      Hist("pf_con_pt_mu_delR")->Fill(pfcon.pt(),pfcon.deltaR(muon));
      Hist("pf_con_pt_nu_delR")->Fill(pfcon.pt(),pfcon.deltaR(neutrino));


      Hist("con_pt_ly")->Fill(pfcon.pt(), weight);
      Hist("con_eta")->Fill(pfcon.eta(), weight);
      Hist("con_phi")->Fill(pfcon.phi(), weight);
    }

    if(count/((double)IsoPFParticle->size())>.8)
      Hist("filter")->Fill(muon.pt(),weight);
    
  
    double minDelR = 999;
    double pf_energy = -1;
    double minDelR_nomuon = 999;
    double pf_energy_nomuon = -1;

    double relIsoParticle004 =0;
    double relIsoParticle05  =0;
    double relIsoParticle01  =0;
    double relIsoParticle008 =0;
    double relIsoParticle002 =0;

    for(unsigned int u = 0; u<IsoPFParticle->size();++u){
      Particle pfcon = IsoPFParticle->at(u);

      relIsoParticle004+=relIsoParticle(pfcon,0.04);
      relIsoParticle05+=relIsoParticle(pfcon,0.5);
      relIsoParticle01+=relIsoParticle(pfcon,0.1);
      relIsoParticle008+=relIsoParticle(pfcon,0.08);
      relIsoParticle002+=relIsoParticle(pfcon,0.02);

      if(pfcon.deltaR(muon)<minDelR && pfcon.deltaR(muon)!=0){
	minDelR_nomuon = pfcon.deltaR(muon);
	pf_energy_nomuon = pfcon.v4().E();
      }
	
      
      if(pfcon.deltaR(muon)<minDelR){
	minDelR = pfcon.deltaR(muon);
	pf_energy = pfcon.v4().E();
      }
    }


    Hist("relIsoSub_004")->Fill(relIsoParticle004,weight);
    Hist("relIsoSub_05" )->Fill(relIsoParticle05 ,weight);
    Hist("relIsoSub_01" )->Fill(relIsoParticle01 ,weight);
    Hist("relIsoSub_008")->Fill(relIsoParticle008,weight);
    Hist("relIsoSub_002")->Fill(relIsoParticle002,weight);



    Hist("nearcon")->Fill(minDelR,weight);
    Hist("nearcon_nomuon")->Fill(minDelR_nomuon,weight);

    Hist("delR_sqrtE")->Fill(minDelR/sqrt(muon.v4().E())*1000,weight);
    Hist("delR_E")->Fill(minDelR/muon.v4().E()*1000,weight);
    Hist("delR_E_pf")->Fill(minDelR/pf_energy,weight);
    Hist("delR_sqrtE_pf")->Fill(minDelR/sqrt(pf_energy),weight);

    Hist("delR_sqrtE_nomuon")->Fill(minDelR_nomuon/sqrt(muon.v4().E()),weight);
    Hist("delR_E_nomuon")->Fill(minDelR_nomuon/muon.v4().E(),weight);
    Hist("delR_E_pf_nomuon")->Fill(minDelR_nomuon/pf_energy_nomuon,weight);
    Hist("delR_sqrtE_pf_nomuon")->Fill(minDelR_nomuon/sqrt(pf_energy_nomuon),weight);

    ((TH2F*)Hist("delR_E_mu"))->Fill(minDelR,muon.v4().E(),weight);
    ((TH2F*)Hist("delR_E_pfcon"))->Fill(minDelR,pf_energy,weight);
    ((TH2F*)Hist("delR_sqrtE_mu"))->Fill(minDelR,sqrt(muon.v4().E()),weight);
    ((TH2F*)Hist("delR_sqrtE_pfcon"))->Fill(minDelR,sqrt(pf_energy),weight);
   

  
    for(unsigned int i = 0; i<1000;++i){
      float p = (float)i;
      p = p*0.001;
      ((TH2F*) Hist("IsoMuon_2D"))->Fill(IsoMuon(muon,p),p,weight);
      ((TH2F*) Hist("relIsoMuon_2D"))->Fill(relIsoMuon(muon,p),p,weight);
      ((TH2F*) Hist("IsoMuon_sqrtE"))->Fill(IsoMuon(muon,p)/sqrt(muon.v4().E()),p,weight);
      ((TH2F*) Hist("IsoMuon_E"    ))->Fill(IsoMuon(muon,p)/muon.v4().E(),p,weight);
      ((TH2F*) Hist("IsoMuon_sqrt" ))->Fill(IsoMuon(muon,p)/sqrt(muon.v4().pt()),p,weight);
      //((TH3F*) Hist("relIsoMuon_3D"))->Fill(relIsoMuon(muon,p),p,muon.pt(),weight);

    }

    ((TH2F*) Hist("IsoMuon02_pt"))->Fill(IsoMuon(muon,0.2),muon.pt(),weight);
    ((TH2F*) Hist("IsoMuon04_pt"))->Fill(IsoMuon(muon,0.4),muon.pt(),weight);


    ((TH2F*) Hist("relIsoMuon_04_pt"))->Fill(relIsoMuon(muon,0.4),muon.pt(),weight);
    ((TH2F*) Hist("relIsoMuon_05_pt" ))->Fill(relIsoMuon(muon,0.5),muon.pt(),weight);
    ((TH2F*) Hist("relIsoMuon_01_pt"))->Fill(relIsoMuon(muon,0.1),muon.pt(),weight);
    ((TH2F*) Hist("relIsoMuon_008_pt"))->Fill(relIsoMuon(muon,0.08),muon.pt(),weight);
    ((TH2F*) Hist("relIsoMuon_002_pt"))->Fill(relIsoMuon(muon,0.02),muon.pt(),weight);

    ((TH2F*) Hist("IsoMuonHT_04_HT" ))->Fill(IsoMuon(muon,0.4)/HT,HT,weight);
    ((TH2F*) Hist("IsoMuonHT_05_HT" ))->Fill(IsoMuon(muon,0.5)/HT,HT,weight);
    ((TH2F*) Hist("IsoMuonHT_01_HT" ))->Fill(IsoMuon(muon,0.1)/HT,HT,weight);
    ((TH2F*) Hist("IsoMuonHT_008_HT"))->Fill(IsoMuon(muon,0.08)/HT,HT,weight);
    ((TH2F*) Hist("IsoMuonHT_002_HT"))->Fill(IsoMuon(muon,0.02)/HT,HT,weight);

    ((TH2F*) Hist("relIsoMuon_04_HT" ))->Fill(relIsoMuon(muon,0.4),HT,weight);
    ((TH2F*) Hist("relIsoMuon_05_HT" ))->Fill(relIsoMuon(muon,0.5),HT,weight);
    ((TH2F*) Hist("relIsoMuon_01_HT" ))->Fill(relIsoMuon(muon,0.1),HT,weight);
    ((TH2F*) Hist("relIsoMuon_008_HT"))->Fill(relIsoMuon(muon,0.08),HT,weight);
    ((TH2F*) Hist("relIsoMuon_002_HT"))->Fill(relIsoMuon(muon,0.02),HT,weight);

    Hist("relIsoMuon")->Fill(relIsoMuon(muon,0.4),weight);
    Hist("relIsoMuon_sqrtE")->Fill(relIsoMuon(muon,0.4)/sqrt(muon.v4().E()),weight);
    Hist("relIsoMuon_004")->Fill(relIsoMuon(muon,0.04),weight);
    Hist("relIsoMuon_05" )->Fill(relIsoMuon(muon,0.5),weight);
    Hist("relIsoMuon_01" )->Fill(relIsoMuon(muon,0.1),weight);
    Hist("relIsoMuon_008")->Fill(relIsoMuon(muon,0.08),weight);
    Hist("relIsoMuon_002")->Fill(relIsoMuon(muon,0.08),weight);

   

  }
}

void ConstituentsHists::Finish()
{
  // final calculations, like division and addition of certain histograms

}

