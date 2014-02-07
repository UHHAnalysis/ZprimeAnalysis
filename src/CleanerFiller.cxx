#include "include/CleanerFiller.h"
#include "include/EventCalc.h"
#include "include/TopFitCalc.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
#include "TH3F.h"



using namespace std;


CleanerFiller::CleanerFiller(const char* name,  HypothesisDiscriminator *discr, TTree* MetaOutputTree) : BaseHists(name), m_tree(MetaOutputTree)
{
  // named default constructor
  m_discr = discr;
}

CleanerFiller::~CleanerFiller()
{
  // default destructor, does nothing
}

void CleanerFiller::Init()
{
  
  //book output tree branches
  
  //m_tree->Branch("data");
 
  m_tree->Branch("toptag"    , &m_toptag     , "toptag/D");
  m_tree->Branch("chi2"      , &m_chi2       , "chi2/D");
  m_tree->Branch("mtt_clean" , &m_mtt_clean  , "mtt_clean/D");
  m_tree->Branch("heptoptag" , &m_heptoptag  , "heptoptag/D");
}

void CleanerFiller::Fill()
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

  LorentzVector toplep = hyp->toplep_v4();
  LorentzVector tophad = hyp->tophad_v4();

  std::vector<Muon>* muons = calc->GetMuons();
  
  std::vector<TopJet>* cajets = calc->GetCAJets();
  std::vector<Jet>* antikjets = calc->GetJets();
  std::vector<TopJet>* hepjets = bcc->toptagjets;


  m_mtt_clean = (toplep+tophad).isTimelike() ? (toplep+tophad).M() : -sqrt((toplep+tophad).M2());
  m_chi2 = hyp->discriminator("Chi2");
  
  m_toptag = 0; 
  m_heptoptag =0;


  double mjet = 0;
  double mmin = 0;
  int nsubjets = 0;

  for(unsigned int i =0; i<cajets->size(); ++i) 
     if(TopTag(cajets->at(i),mjet,nsubjets,mmin)) m_toptag += 1;	


  for(unsigned int i =0; i<hepjets->size();++i)
    if(variableHepTopTag(hepjets->at(i))) m_heptoptag += 1;


  m_tree->Fill();
  

 
}

void CleanerFiller::Finish()
{
  // final calculations, like division and addition of certain histograms

}

