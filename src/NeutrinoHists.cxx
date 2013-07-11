#include "include/NeutrinoHists.h"
#include "include/SelectionModules.h"
#include "include/TopFitCalc.h"

#include <iostream>

using namespace std;

NeutrinoHists::NeutrinoHists(const char* name, HypothesisDiscriminator *discr) : BaseHists(name)
{
  // named default constructor
   m_discr = discr;
}

NeutrinoHists::~NeutrinoHists()
{
  // default destructor, does nothing
}

void NeutrinoHists::Init()
{
  Book( TH1F("met_pT_ly" ,"p_{T} met [GeV]",100,0,1000));
  Book( TH1F("met_phi","#phi met",100,-PI,PI));
  Book( TH1F("met_eta","#eta met",100,-1,1));
  Book( TH1F("diff_met_neutrino_phi_ly","#Delta #phi met-Neutrino",100,0,1));

  Book( TH1F("neutrino_pT_ly" ,"p_{T} neutrino [GeV]",100,0,1000));
  Book( TH1F("neutrino_phi","#phi neutrino",100,-PI,PI));
  Book( TH1F("neutrino_eta","#eta neutrino",100,-3,3));
  Book( TH1F("neutrino_pz_ly" ,"p_{z} neutrino [GeV]",100,-2500,2500));
  Book( TH1F("neutrino_chi2" ,"neutrino #chi^{2}",100,0,4000));

  Book( TH1F("neutrino_std_pT_ly" ,"p_{T,img} neutrino  [GeV]",100,0,1000));
  Book( TH1F("neutrino_std_phi"   ,"#phi_{img} neutrino",100,-PI,PI));
  Book( TH1F("neutrino_std_eta"   ,"#eta_{img} neutrino",100,-PI,PI));
  Book( TH1F("neutrino_std_px"    ,"p_{x,img} neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino_std_py"    ,"p_{y,img} neutrino [GeV]",100,-1000,1000));	

  Book( TH1F("neutrino0_pT_ly" ,"p_{T,fit} neutrino [GeV]",100,0,1000));
  Book( TH1F("neutrino0_phi"   ,"#phi_{fit} neutrino",100,-PI,PI));
  Book( TH1F("neutrino0_eta"   ,"#eta_{fit} neutrino",100,-3,3));
  Book( TH1F("neutrino0_px"    ,"p_{x,fit} neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino0_py"    ,"p_{y,fit} neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("diff_met_neutrino0_phi","#Delta #phi met-1 Neutrino",100,0,1));


  Book( TH1F("neutrino1_MC_px"    ,"p_{x} MC neutrino [GeV]",100,-1000,1000));

  Book( TH1F("neutrino1_pT_ly" ,"p_{T} 1 out of 2 neutrino [GeV]",100,0,1000));
  Book( TH1F("neutrino1_phi","#phi 1 out of 2 neutrino",100,-PI,PI));
  Book( TH1F("neutrino1_eta","#eta 1 out of 2 neutrino",100,-3,3));
  Book( TH1F("neutrino1_px"    ,"p_{x} 1 1 out of 2  neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino1_py"    ,"p_{y} 1 1 out of 2  neutrino [GeV]",100,-1000,1000));	

  Book( TH1F("neutrino2_pT_ly" ,"p_{T} 2 out of 2  neutrino [GeV]",100,0,1000));
  Book( TH1F("neutrino2_phi","#phi 2 out of 2 neutrino",100,-PI,PI));
  Book( TH1F("neutrino2_eta","#eta 2 out of 2 neutrino",100,-3,3));
  Book( TH1F("neutrino2_px"    ,"p_{x} 2 out of 2 neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino2_py"    ,"p_{y} 2 out of 2 neutrino [GeV]",100,-1000,1000));	

  Book( TH1F("neutrino_MC_pt_ly" ,"p_{T} MC neutrino [GeV]",100,0,1000));
  Book( TH1F("neutrino_MC_eta","#eta MC neutrino",100,-PI,PI));
  Book( TH1F("neutrino_MC_phi","#phi MC neutrino",100,-PI,PI));
  Book( TH1F("neutrino_MC_E_ly"  ,"E MC neutrino",100,0,2800));
  Book( TH1F("neutrino_MC_px"    ,"p_{x} MC neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino_MC_py"    ,"p_{y} MC neutrino [GeV]",100,-1000,1000));	
  Book( TH1F("neutrino_MC_pz"    ,"p_{z} MC neutrino [GeV]",100,-1000,1000));	

  Book( TH1F("neutrino_diff_pt","Diff between the 2 Neutrinos",100,-500,500));

  Book( TH2F("neutrino_px_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino_py_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino_phi_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino_eta_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino_pt_reco_gen","", 60, -10, 150, 60, -5, 5));

  Book( TH2F("neutrino0_px_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino0_py_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino0_phi_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino0_eta_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino0_pt_reco_gen","", 60,-10, 150, 60, -5,5));

  Book( TH2F("neutrino_std_px_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino_std_py_reco_gen","", 60, -400, 400, 60, -400, 400));
  Book( TH2F("neutrino_std_phi_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino_std_eta_reco_gen","", 60, -4, 4, 60, -4, 4));
  Book( TH2F("neutrino_std_pt_reco_gen","", 60,-10, 150, 60, -5,5));

  /*
  Book( TH1D( "neutrino0_slice","Pt reco gen Slice",60,0,400));	
  Book( TH1D( "neutrino_std_slice","Pt reco gen Slice",60,0,400));	
  
  /*
  Book( TH1D( "pt_reco_gen_slice","Pt reco gen Slice",10,0,600));	
  Book( TH1D( "pt_fitreco_gen_slice","Pt fitreco gen Slice",10,0,600));

  Book( TH1D( "slice1_slice2","pt fitreco_gen / reco_gen",10,0,600));
  */

}


void NeutrinoHists::Fill()
{

  EventCalc* calc = EventCalc::Instance();
  double weight = calc -> GetWeight();

  TopFitCalc* fitcalc = TopFitCalc::Instance();

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  LorentzVector neutrino = hyp->neutrino_v4();
  std::vector<LorentzVector> neutrinos = fitcalc->NeutrinoFitPolar(calc->GetPrimaryLepton()->v4(),bcc->met->v4());
  std::vector<LorentzVector> neutrino_std = calc->NeutrinoReconstruction(calc->GetPrimaryLepton()->v4(),bcc->met->v4());
  

  Hist("met_pT_ly" )->Fill(bcc->met->v4().pt(),weight);
  Hist("met_eta"   )->Fill(bcc->met->v4().eta(),weight);
  Hist("met_phi"   )->Fill(bcc->met->v4().phi(),weight); 

  Hist("neutrino_pT_ly" )->Fill(neutrino.pt(),weight);
  Hist("neutrino_eta")->Fill(neutrino.eta(),weight);
  Hist("neutrino_phi")->Fill(neutrino.phi(),weight);  
  Hist("neutrino_pz_ly" )->Fill(neutrino.pz(),weight);
  Hist("diff_met_neutrino_phi_ly")->Fill(fitcalc->delPhi(neutrino.phi(),bcc->met->v4().phi()),weight); 
  
	
  if(calc->GetGenParticles() ){
    if(calc->GetTTbarGen()->DecayChannel()== TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad || calc->GetTTbarGen()->DecayChannel()== TTbarGen::e_tauhad){	     
      GenParticle GenNeutrino = calc->GetTTbarGen()->Neutrino();
	
      Hist("neutrino_MC_pt_ly"  )->Fill(GenNeutrino.v4().pt() ,weight);
      Hist("neutrino_MC_eta" )->Fill(GenNeutrino.v4().eta() ,weight);
      Hist("neutrino_MC_phi" )->Fill(GenNeutrino.v4().phi() ,weight);
      Hist("neutrino_MC_E_ly"   )->Fill(GenNeutrino.v4().E() ,weight);
      Hist("neutrino_MC_px"   )->Fill(GenNeutrino.v4().px() ,weight);
      Hist("neutrino_MC_py"   )->Fill(GenNeutrino.v4().py() ,weight);
      Hist("neutrino_MC_pz"   )->Fill(GenNeutrino.v4().pz() ,weight);
      
      Hist("neutrino_pt_reco_gen"  )->Fill(neutrino.pt() ,(neutrino.pt()-GenNeutrino.pt())/neutrino.pt());
      Hist("neutrino_px_reco_gen"  )->Fill(neutrino.px() , GenNeutrino.v4().px());
      Hist("neutrino_py_reco_gen"  )->Fill(neutrino.py() , GenNeutrino.v4().py());
      Hist("neutrino_phi_reco_gen" )->Fill(neutrino.phi(), neutrino.phi()-GenNeutrino.phi());
      Hist("neutrino_eta_reco_gen" )->Fill(neutrino.eta(), neutrino.eta()-GenNeutrino.eta());

      if(neutrinos.size()==1){
	Hist("neutrino0_px_reco_gen"  )->Fill(neutrinos.at(0).px() , GenNeutrino.v4().px());
	Hist("neutrino0_py_reco_gen"  )->Fill(neutrinos.at(0).py() , GenNeutrino.v4().py());
	Hist("neutrino0_phi_reco_gen" )->Fill(GenNeutrino.phi(), neutrinos.at(0).phi()-GenNeutrino.phi());
	Hist("neutrino0_eta_reco_gen" )->Fill(GenNeutrino.eta(), neutrinos.at(0).eta()-GenNeutrino.eta());
	Hist("neutrino0_pt_reco_gen"  )->Fill(GenNeutrino.pt() , (neutrinos.at(0).pt()-GenNeutrino.pt())/neutrinos.at(0).pt());

	Hist("neutrino_std_px_reco_gen"  )->Fill(neutrino_std.at(0).px() , GenNeutrino.v4().px());
	Hist("neutrino_std_py_reco_gen"  )->Fill(neutrino_std.at(0).py() , GenNeutrino.v4().py());
	Hist("neutrino_std_phi_reco_gen" )->Fill(GenNeutrino.phi(), neutrino_std.at(0).phi()-GenNeutrino.phi());
	Hist("neutrino_std_eta_reco_gen" )->Fill(GenNeutrino.eta(), neutrino_std.at(0).eta()-GenNeutrino.eta());
	Hist("neutrino_std_pt_reco_gen"  )->Fill(GenNeutrino.pt() , (neutrino_std.at(0).pt()-GenNeutrino.pt())/neutrinos.at(0).pt());

	Hist("neutrino1_MC_px"   )->Fill(GenNeutrino.v4().px() ,weight);
      }
    }
  }	
  if(neutrinos.size()== 1){
    Hist("neutrino_chi2")->Fill(pow(neutrinos.at(0).px()-bcc->met->v4().px(),2)+pow(neutrinos.at(0).py()-bcc->met->v4().py(),2),weight);
    
    Hist("neutrino0_pT_ly" )->Fill(neutrinos.at(0).pt(),weight);
    Hist("neutrino0_eta"   )->Fill(neutrinos.at(0).eta(),weight);
    Hist("neutrino0_phi"   )->Fill(neutrinos.at(0).phi(),weight); 
    Hist("neutrino0_px"    )->Fill(neutrinos.at(0).px(),weight); 
    Hist("neutrino0_py"    )->Fill(neutrinos.at(0).py(),weight); 
    
    Hist("diff_met_neutrino0_phi"    )->Fill(fitcalc->delPhi(neutrinos.at(0).phi(),bcc->met->v4().phi()),weight); 

    Hist("neutrino_diff_pt")->Fill(neutrino_std.at(0).pt()-neutrinos.at(0).pt(),weight);

    Hist("neutrino_std_pT_ly" )->Fill(neutrino_std.at(0).pt(),weight);
    Hist("neutrino_std_phi"   )->Fill(neutrino_std.at(0).phi(),weight);
    Hist("neutrino_std_eta"   )->Fill(neutrino_std.at(0).eta(),weight);
    Hist("neutrino_std_px"    )->Fill(neutrino_std.at(0).px(),weight); 
    Hist("neutrino_std_py"    )->Fill(neutrino_std.at(0).py(),weight);
      
      
  }
  else if(neutrinos.size() == 2){
    Hist("neutrino1_pT_ly" )->Fill(neutrinos.at(0).pt(),weight);
    Hist("neutrino1_eta"   )->Fill(neutrinos.at(0).eta(),weight);
    Hist("neutrino1_phi"   )->Fill(neutrinos.at(0).phi(),weight); 
    Hist("neutrino1_px"    )->Fill(neutrinos.at(0).px(),weight); 
    Hist("neutrino1_py"    )->Fill(neutrinos.at(0).py(),weight); 
    
    Hist("neutrino2_pT_ly" )->Fill(neutrinos.at(1).pt(),weight);
    Hist("neutrino2_eta"   )->Fill(neutrinos.at(1).eta(),weight);
    Hist("neutrino2_phi"   )->Fill(neutrinos.at(1).phi(),weight); 
    Hist("neutrino2_px"    )->Fill(neutrinos.at(1).px(),weight); 
    Hist("neutrino2_py"    )->Fill(neutrinos.at(1).py(),weight); 
  }



}    






void NeutrinoHists::Finish()
{
  
}
