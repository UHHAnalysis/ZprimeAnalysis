#include "include/JetEffiHists.h"
#include "include/SelectionModules.h"
#include "include/EventCalc.h"
#include <iostream>

using namespace std;

JetEffiHists::JetEffiHists(const char* name, double  R_lep, double R_had, double pt_lep, double pt_had) : BaseHists(name)
{
  // named default constructor
  m_R_lep = R_lep;
  m_R_had = R_had;
  m_pt_lep= pt_lep;
  m_pt_had= pt_had;



}

JetEffiHists::~JetEffiHists()
{
  // default destructor, does nothing
}

void JetEffiHists::Init()
{





  // book all histograms here
  Book( TH1F( "gentophad_pt", "GenTop_{had} pt [GeV]", 50, 0, 2000 ) );
  Book( TH1F( "gen_mttbar", "M_{t\bar{t}} [GeV] ", 100, 0, 4500 ) );
  Book( TH1F( "gentop_matched_pt","pt", 100, 0, 2000 ) );
  Book( TH1F( "leading_jet_pt","pt of leading C/A jet",50, 0, 2000 ) );
  Book( TH1F( "leading_jet_eta","pt of leading C/A jet",50, -5, 5 ) );
  Book( TH1F( "leading_jet_phi","pt of leading C/A jet",50, -4, 4 ) );

  /* Book( TH1F( "genblep_pt", "GenB_{lep} pt ", 50, 0, 2000 ) );
  Book( TH1F( "cajets_tophad_matched_pt", "C/A-Jets matched Top_{had} pt ", 50, 0, 2000 ) );
  Book( TH1F( "cajets_blep_matched_pt", "C/A-Jets matched B_{lep} pt ", 50, 0, 2000 ) );
  Book( TH1F( "antijets_tophad_matched_pt", "anti-kt Jets matched Top_{had} pt ", 50, 0, 2000 ) );
  Book( TH1F( "antijets_blep_matched_pt", "anti-kt Jets matched b_{lep} pt ", 50, 0, 2000 ) );
  // Book( TH1F( "gentoplep_matched_pt", "GenTop_{lep} pt ", 50, 0, 2000 ) );
  Book( TH1F( "effi_tophad_ca_pt", "Effizienz C/A-Jet Top_{had} pt ", 50, 0, 2000 ) );
  Book( TH1F( "effi_tophad_anti_pt", "Effizienz anti-kt Jet Top_{had} pt ", 50, 0, 2000 ) );
  Book( TH1F( "effi_blep_ca_pt", "Effizienz C/A-Jet pt B_{lep}", 50, 0, 2000 ) );
  Book( TH1F( "effi_blep_anti_pt", "Effizienz anti-kt Jet B_{lep} pt", 50, 0, 2000 ) );
  Book( TH1F( "toplep_matcheble_pt", "Matchable T_{lep} pt", 50, 0, 2000 ) );
  Book( TH1F( "tophad_matcheble_pt", "Matchable T_{had} pt", 50, 0, 2000 ) );
  Book( TH1F( "toplep_matcheble_pt", "Matchable T_{lep} pt", 50, 0, 2000 ) );
  Book( TH1F( "Mttbar", "Matchable M_{tt}", 200, 0,5000) );
  Book( TH1F( "Mttbar_ly", "Matchable M_{tt}", 200, 0,5000) );
  Book( TH1F( "toptag_pt", "top tag pt",50,0,2000));
  Book( TH1F( "chi_pt", "top tag pt",50,0,2000));
  Book( TH1F( "effi_toptag_pt", "Effizienz top tag pt",50,0,2000));
  Book( TH1F( "effi_chi_pt", "Effizienz #chi^2 pt",50,0,2000));
  */

}

void JetEffiHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  bool IsRealData = calc->IsRealData();
  double weight = calc->GetWeight();	  

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  std::vector<Jet>* antikjets = calc->GetJets();	
  std::vector<TopJet>* cajets = calc->GetCAJets();	

  
  if(cajets->size()>0){
    sort(cajets->begin(),cajets->end(),HigherPt());
    Hist("leading_jet_pt")->Fill(cajets->at(0).pt(),weight);
    Hist("leading_jet_phi")->Fill(cajets->at(0).phi(),weight);
    Hist("leading_jet_eta")->Fill(cajets->at(0).eta(),weight);
  }


  if(IsRealData) return;
  if(calc->GetTTbarGen()->DecayChannel() != TTbarGen::e_muhad) return;

  GenParticle GenTopLep = calc->GetTTbarGen()->TopLep();
  GenParticle GenTopHad = calc->GetTTbarGen()->TopHad();
  GenParticle GenBLep   = calc->GetTTbarGen()->BLep();
  GenParticle GenWHad   = calc->GetTTbarGen()->WHad();


  Hist("gentophad_pt")->Fill(GenTopHad.pt(),weight);


  Hist("gen_mttbar")->Fill((GenTopHad.v4()+GenTopLep.v4()).M(),weight);



  double mjet = 0;
  int nsubjets= 0;
  double mmin = 0;

  for(unsigned int i = 0 ; i<cajets->size(); ++i){
    TopJet cajet = cajets->at(i);
    if(TopTag(cajet,mjet,nsubjets,mmin) && GenTopHad.deltaR(cajet)<0.7)
      Hist("gentop_matched_pt")->Fill(GenTopHad.pt(),weight);

  }

 


  //Hist("genblep_pt")->Fill(GenBLep.pt(),weight);

  /*
  double ca_tophad_pt=0;
  int ca_had_num =0;
  int* ca_had_pos = new int[cajets->size()];
	
  for(unsigned int i = 0 ; i<cajets->size(); ++i){
    TopJet cajet = cajets->at(i);
    if(GenTopHad.deltaR(cajet)<m_R_had){ 
      ca_had_num++;
      ca_had_pos[ca_had_num]=i;
      ca_tophad_pt += cajet.pt();
    }
  }
   
  if(GenTopHad.pt()> m_pt_had*ca_tophad_pt && m_pt_had*GenTopHad.pt()<ca_tophad_pt) Hist("cajets_tophad_matched_pt")->Fill(GenTopHad.pt(),weight);
  
  double anti_tophad_pt=0;
  int anti_had_num =0;
  int* anti_had_pos = new int[antikjets->size()];

  for(unsigned int i = 0 ; i<antikjets->size(); ++i){    
    Jet antijet = antikjets->at(i);
    if(GenTopHad.deltaR(antijet)<m_R_had){
      anti_had_num++;
      anti_had_pos[anti_had_num]=i;
      anti_tophad_pt += antijet.pt();
    }
  }

  if(GenTopHad.pt()> m_pt_had*anti_tophad_pt && m_pt_had*GenTopHad.pt()<anti_tophad_pt) Hist("antijets_tophad_matched_pt")->Fill(GenTopHad.pt(),weight);  
  
	
  double anti_blep_pt=0;
  int anti_blep_num =0;
  int* anti_blep_pos = new int[antikjets->size()];

  //cout<<antikjets->size()<<endl;

  for(unsigned int i = 0 ; i<antikjets->size(); ++i){
    Jet antijet = antikjets->at(i);
    if(GenBLep.deltaR(antijet)<m_R_lep){
      anti_blep_num++;
      anti_blep_pos[anti_had_num]=i;
      anti_blep_pt += antijet.pt();
    }
  }

  if(GenBLep.pt()> m_pt_lep*anti_blep_pt && m_pt_lep*GenBLep.pt()<anti_blep_pt) Hist("antijets_blep_matched_pt")->Fill(GenBLep.pt(),weight);

  double ca_blep_pt=0;
  int ca_lep_num =0;
  int* ca_blep_pos = new int[cajets->size()];

  for(unsigned int i = 0 ; i<cajets->size(); ++i){
    TopJet cajet = cajets->at(i);
    if(GenBLep.deltaR(cajet)<m_R_lep){
      ca_lep_num++;
      ca_blep_pos[ca_had_num]=i;
      ca_blep_pt += cajet.pt();
    }
  }

  if(GenBLep.pt()> m_pt_lep*ca_blep_pt && m_pt_lep*GenBLep.pt()<ca_blep_pt) Hist("cajets_blep_matched_pt")->Fill(GenBLep.pt(),weight);
  
  
  if((GenBLep.pt()> m_pt_lep*ca_blep_pt && m_pt_lep*GenBLep.pt()<ca_blep_pt) || (GenBLep.pt()> m_pt_lep*anti_blep_pt && m_pt_lep*GenBLep.pt()<anti_blep_pt)){
    if((GenTopHad.pt()> m_pt_had*anti_tophad_pt && m_pt_had*GenTopHad.pt()<anti_tophad_pt) || (GenTopHad.pt()> m_pt_had*ca_tophad_pt && m_pt_had*GenTopHad.pt()<ca_tophad_pt)){
      Hist("Mttbar")->Fill((GenTopLep.v4()+GenTopHad.v4()).M(),weight);
      Hist("Mttbar_ly")->Fill((GenTopLep.v4()+GenTopHad.v4()).M(),weight);
    }
  }

  if((GenTopHad.pt()> m_pt_had*anti_tophad_pt && m_pt_had*GenTopHad.pt()<anti_tophad_pt) || (GenTopHad.pt()> m_pt_had*ca_tophad_pt && m_pt_had*GenTopHad.pt()<ca_tophad_pt))Hist("tophad_matcheble_pt")->Fill(GenTopHad.pt(),weight);
  if((GenBLep.pt()> m_pt_lep*ca_blep_pt && m_pt_lep*GenBLep.pt()<ca_blep_pt) || (GenBLep.pt()> m_pt_lep*anti_blep_pt && m_pt_lep*GenBLep.pt()<anti_blep_pt)) Hist("toplep_matcheble_pt")->Fill(GenTopLep.pt(),weight);



  double mjet = 0;
  int nsubjets= 0;
  double mmin = 0;

  for(unsigned int i = 0 ; i<cajets->size(); ++i){
    TopJet cajet = cajets->at(i);
    if(TopTag(cajet,mjet,nsubjets,mmin) && GenTopHad.deltaR(cajet)<0.7)
      Hist("toptag_pt")->Fill(GenTopHad.pt(),weight);

  }

  bool notfound = true;

  for(unsigned int i = 0; i<antikjets->size(); ++i){
    Jet antijet = antikjets->at(i);
    if(pow(174-GenTopHad.v4().M(),2)+pow(181-GenTopLep.v4().M(),2)<10 && GenTopHad.deltaR(antijet)<0.7 && notfound){
      Hist("chi_pt")->Fill(GenTopHad.pt(),weight);
      notfound = false;
    }



  }


  */

}

void JetEffiHists::Finish()
{/*
  // final calculations, like division and addition of certain histograms`
  Hist("effi_tophad_ca_pt")->Divide(Hist("cajets_tophad_matched_pt"),Hist("gentophad_pt"));
  Hist("effi_tophad_anti_pt")->Divide(Hist("antijets_tophad_matched_pt"),Hist("gentophad_pt"));

  Hist("effi_blep_ca_pt")->Divide(Hist("cajets_blep_matched_pt"),Hist("genblep_pt"));
  Hist("effi_blep_anti_pt")->Divide(Hist("antijets_blep_matched_pt"),Hist("genblep_pt"));

  Hist("effi_toptag_pt")->Divide(Hist("toptag_pt"),Hist("gentophad_pt"));
  Hist("effi_chi_pt")->Divide(Hist("chi_pt"),Hist("gentophad_pt"));
 */


}

