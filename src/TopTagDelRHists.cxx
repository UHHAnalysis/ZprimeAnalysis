#include "include/TopTagDelRHists.h"
#include "include/EventCalc.h"
#include "include/SelectionModules.h"
#include "include/TopFitCalc.h"
#include "SFrameTools/include/SubJetTagger.h"
#include <iostream>
#include "TH3F.h"

using namespace std;

TopTagDelRHists::TopTagDelRHists(const char* name, HypothesisDiscriminator *discr) : BaseHists(name)
{
  // named default constructor
  m_discr = discr;
}

TopTagDelRHists::~TopTagDelRHists()
{
  // default destructor, does nothing
}

void TopTagDelRHists::Init()
{
  Book( TH1F("DelR_Gen_BLep_TopTag"  ,"#Delta R(B_{Lep,Gen},C/A_{TopTag})",10,0,5));
  Book( TH1F("DelR_Gen_TopLep_TopTag","#Delta R(Top_{Lep,Gen},C/A_{TopTag})",10,0,5));
  Book( TH1F("DelR_muon_TopTag"      ,"#Delta R(#mu,C/A_{TopTag})",10,0,5));

  Book( TH1F("DelPhi_Gen_BLep_TopTag"  ,"#Delta #Phi (B_{Lep,Gen},C/A_{TopTag})",10,-PI,PI));
  Book( TH1F("DelPhi_Gen_TopLep_TopTag","#Delta #Phi (Top_{Lep,Gen},C/A_{TopTag})",10,-PI,PI));
  Book( TH1F("DelPhi_muon_TopTag"      ,"#Delta #Phi (#mu,C/A_{TopTag})",10,-PI,PI));

  Book( TH1F("E_muon_TopTag"   ,"#Delta E(C/A_{TopTag},#mu)",100,0,300));
  Book( TH1F("Pt_muon_TopTag"  ,"#Delta p_{T}(C/A_{TopTag},#mu)",100,0,300));
  Book( TH1F("Mass_muon_TopTag","#Delta m(C/A_{TopTag},#mu)",100,0,300));
}

void TopTagDelRHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  TopFitCalc* fitcalc = TopFitCalc::Instance();
  double weight = calc -> GetWeight();

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  std::vector<TopJet>* cajets = calc->GetCAJets();
  std::vector<Jet>* antikjets = calc->GetJets();


  Muon muon  = bcc->muons->at(0);

  Particle neutrino;
  neutrino.set_v4(hyp->neutrino_v4());
  Particle wlep;
  wlep.set_v4(hyp->wlep_v4());


  LorentzVector toplep = hyp->toplep_v4();
  LorentzVector tophad = hyp->tophad_v4();

  GenParticle BLep_gen = calc->GetTTbarGen()->BLep();
  GenParticle TopLep_gen = calc->GetTTbarGen()->TopLep();
  GenParticle TopHad_gen = calc->GetTTbarGen()->TopHad();

  int cajet_tag_num=-1;

  CMSTopTagger toptag;
  toptag.SetTau32Cut();

  for(unsigned int m =0;m<cajets->size();++m)
    if(toptag.Tag(cajets->at(m))) cajet_tag_num=m;
  
  TopJet TopTagJet = cajets->at(cajet_tag_num); 
    
  Hist("DelR_Gen_BLep_TopTag"  )->Fill(BLep_gen.deltaR(TopTagJet),weight);
  Hist("DelR_Gen_TopLep_TopTag")->Fill(TopLep_gen.deltaR(TopTagJet),weight);
  Hist("DelR_muon_TopTag"      )->Fill(muon.deltaR(TopTagJet),weight);

  Hist("DelPhi_Gen_BLep_TopTag"  )->Fill(BLep_gen.deltaPhi(TopTagJet),weight);
  Hist("DelPhi_Gen_TopLep_TopTag")->Fill(TopLep_gen.deltaPhi(TopTagJet),weight);
  Hist("DelPhi_muon_TopTag"      )->Fill(muon.deltaPhi(TopTagJet),weight);


  Hist("E_muon_TopTag"   )->Fill((TopTagJet.v4()-muon.v4()).E(),weight);
  Hist("Pt_muon_TopTag"  )->Fill((TopTagJet.v4()-muon.v4()).pt(),weight);
  Hist("Mass_muon_TopTag")->Fill(sqrt((TopTagJet.v4()-muon.v4()).M2()),weight);


}

void TopTagDelRHists::Finish()
{
  // final calculations, like division and addition of certain histograms

}

