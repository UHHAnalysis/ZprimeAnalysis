#include "include/DelRHists.h"
#include "include/EventCalc.h"
#include "include/SelectionModules.h"
#include "include/TopFitCalc.h"
#include <iostream>
#include "TH3F.h"

using namespace std;


double distance(LorentzVector v1, LorentzVector v2, double metrik , double R=1);
double jade(LorentzVector v1,LorentzVector v2, double mass);


double pTrel_Part( Particle p1, Particle p2)
{

 double ptrel=0;

 TVector3 p3(p1.v4().Px(),p1.v4().Py(),p1.v4().Pz());
 TVector3 p4(p2.v4().Px(),p2.v4().Py(),p2.v4().Pz());

 if(p3.Mag()!=0 && p4.Mag()!=0) {
 double sin_alpha = (p3.Cross(p4)).Mag()/p3.Mag()/p4.Mag();
 ptrel = p3.Mag()*sin_alpha;
 } else {
 std::cout << "something strange happend in the ptrel calculation: either lepton or jet momentum is 0" <<std::endl;
 }

 return ptrel;
 }

double pTrel_Part( LorentzVector p1, LorentzVector p2)
{

 double ptrel=0;

 TVector3 p3(p1.Px(),p1.Py(),p1.Pz());
 TVector3 p4(p2.Px(),p2.Py(),p2.Pz());

 if(p3.Mag()!=0 && p4.Mag()!=0) {
 double sin_alpha = (p3.Cross(p4)).Mag()/p3.Mag()/p4.Mag();
 ptrel = p3.Mag()*sin_alpha;
 } else {
 std::cout << "something strange happend in the ptrel calculation: either lepton or jet momentum is 0" <<std::endl;
 }

 return ptrel;
 }




double distance(LorentzVector v1, LorentzVector v2, double metrik , double R){
  if(metrik == 0)
    return deltaR(v1,v2)/R;
  else if(pow(v1.pt(),metrik) < pow(v2.pt(),2*metrik))
    return v1.pt()*deltaR(v1,v2)/R;
  else if(pow(v1.pt(),metrik) > pow(v2.pt(),2*metrik))
    return v2.pt()*deltaR(v1,v2)/R;
  else
    cerr<<"Could not measure the distance"<<endl;
}

double jade(LorentzVector v1,LorentzVector v2, double mass){

  TVector3 p1(v1.px(),v1.py(),v1.pz());
  TVector3 p2(v2.px(),v2.py(),v2.pz());
  
  return v1.E()*v2.E()* (1-cos(p1.Angle(p2)))/mass;

}

DelRHists::DelRHists(const char* name, HypothesisDiscriminator *discr) : BaseHists(name)
{
  // named default constructor
  m_discr = discr;
}

DelRHists::~DelRHists()
{
  // default destructor, does nothing
}

void DelRHists::Init()
{
  Book( TH1F("muon_jet_ly","#Delta R (Jet, #mu)",100,0,1.5));
  Book( TH1F("muon_jet","#Delta R (Jet, #mu)",100,0,1.5));

  Book( TH1F("muon_jet_long_ly","#Delta R (Jet, #mu)",100,0,5));
  Book( TH1F("muon_jet_long","#Delta R (Jet, #mu)",100,0,5));

  Book( TH1F("muon_bjet_delR_ly","#Delta R (b-Jet, #mu) ",100,0,1.5));
  Book( TH1F("muon_bjet_delR","#Delta R (b-Jet, #mu) ",100,0,1.5));

  Book( TH1F("muon_bjet_delR_long_ly","#Delta R (b-Jet, #mu) ",100,0,5));
  Book( TH1F("muon_bjet_delR_long","#Delta R (b-Jet, #mu) ",100,0,5));

  Book( TH1F("neutrino_bjet_delR_ly","#Delta R (b-Jet, #nu) ",100,0,1.5));
  Book( TH1F("neutrino_bjet_delR","#Delta R (b-Jet, #nu) ",100,0,1.5));

  Book( TH1F("neutrino_bjet_delR_long_ly","#Delta R (b-Jet, #nu) ",100,0,5));
  Book( TH1F("neutrino_bjet_delR_long","#Delta R (b-Jet, #nu) ",100,0,5));

  Book( TH1F("neutrino_jet_delR_ly","#Delta R (Jet, #nu) ",100,0,1.5));
  Book( TH1F("neutrino_jet_delR","#Delta R (Jet, #nu) ",100,0,1.5));

  Book( TH1F("neutrino_jet_delR_long_ly","#Delta R (Jet, #nu) ",100,0,5));
  Book( TH1F("neutrino_jet_delR_long","#Delta R (Jet, #nu) ",100,0,5));

  Book( TH1F("bjet_jet_delR_ly","#Delta R (b-Jet, Jet) ",100,0,1.5));
  Book( TH1F("bjet_jet_delR","#Delta R (b-Jet, Jet) ",100,0,1.5));

  Book( TH1F("bjet_jet_delR_long_ly","#Delta R (b-Jet, Jet) ",100,0,5));
  Book( TH1F("bjet_jet_delR_long","#Delta R (b-Jet, Jet) ",100,0,5));
 
  Book( TH1F("delRmax","#Delta Rmax ",100,0,7));
  Book( TH1F("delRmin","#Delta Rmin ",100,0,5));


  Book( TH1F("muon_bjet_pTrel"," p_{T,rel}(b-Jet,#mu) ",100,0,1000));
  Book( TH1F("neutrino_bjet_pTrel"," p_{T,rel}(b-Jet,#nu) ",100,0,1000));
  Book( TH1F("muon_neutrino_pTrel"," p_{T,rel}(#mu,#nu ) ",100,0,1000));

  Book( TH1F("muon_bjet_pTrel_ly"," p_{T,rel}(b-Jet,#mu) ",100,0,1000));
  Book( TH1F("neutrino_bjet_pTrel_ly"," p_{T,rel}(b-Jet,#nu) ",100,0,1000));
  Book( TH1F("muon_neutrino_pTrel_ly"," p_{T,rel}(#mu,#nu ) ",100,0,1000));

  Book( TH1F("deleta_top"," #Delta #eta (Top_{had},Top_{lep}) ",100,0,4));
  Book( TH1F("eta_tops"," #eta(Top_{had}+Top_{lep}) ",100,-7,7));
  Book( TH1F("deleta_top_ly"," #Delta #eta (Top_{had},Top_{lep}) ",100,0,4));
  Book( TH1F("eta_tops_ly"," #eta(Top_{had}+Top_{lep}) ",100,-7,7));


  Book( TH2F("deltaRmax_top_pt","",100,0,7,100,0,3000));
  Book( TH2F("deltaRmin_top_pt","",100,0,5,100,0,3000));

  Book( TH2F("deltaRmax_top_mass","",100,0,5,100,0,3000));
  Book( TH2F("deltaRmax_top_mtt","",100,0,5,100,0,5000));
  Book( TH2F("deltaRmax_top_pt_mtt","",100,0,5,100,0,5));

  Book( TH2F("deltaRmin_top_mass","",100,0,5,100,0,2000));


  Book( TH1F("DistanceMax","Distance max ",100,0,3000));
  Book( TH1F("DistanceMin","Distance min ",100,0,3000));
  Book( TH1F("DistanceMin_ly","Distance min ",100,0,3000));

  Book( TH2F("DistanceMax_top_pt","",100,0,2000,100,0,2000));
  Book( TH2F("DistanceMin_top_pt","",100,0,2000,100,0,2000));

  Book( TH2F("DistanceMax_top_mass","",100,0,3000,100,0,3000));
  Book( TH2F("DistanceMax_top_mtt","",100,0,3000,100,0,5000));
  Book( TH2F("DistanceMax_top_pt_mtt","",100,0,3000,100,0,5));

  Book( TH2F("DistanceMin_top_mass","",100,0,2000,100,0,2000));


  Book( TH1F("JadeMax","Jade max ",100,0,2000));
  Book( TH1F("JadeMin","Jade min ",100,0,2000));

  Book( TH2F("JadeMax_top_pt","",100,0,2000,100,0,3000));
  Book( TH2F("JadeMin_top_pt","",100,0,2000,100,0,3000));


  Book( TH2F("JadeMax_top_mass","",100,0,2000,100,0,3000));
  Book( TH2F("JadeMax_top_mtt","",100,0,2000,100,0,5000));
  Book( TH2F("JadeMax_top_pt_mtt","",100,0,2000,100,0,5));
 
  Book( TH2F("JadeMin_top_mass","",100,0,2000,100,0,3000));

  Book( TH1F("SumVec","SumVec p_{T} ",100,-1000,1000));
  
  Book( TH2F("muon_bjet_pTrel_delR","",100,0,500,100,0,4));
  Book( TH2F("muon_bjet_pTrel_muon_jet_delR","",100,0,500,100,0,4));
  Book( TH2F("muon_bjet_pTrel_bjet_jet_delR","",100,0,500,100,0,5));
  Book( TH2F("muon_bjet_pTrel_muon_neutrino_delR","",100,0,100,200,0,5));
  Book( TH2F("neutrino_bjet_pTrel_muon_bjet_delR","",100,0,500,200,0,5));
  Book( TH2F("neutrino_bjet_pTrel_neutrino_bjet_delR","",100,0,500,200,0,5));

  Book( TH2F("wlep_bjet_pTrel_delR"              ,"",100,0,500,200,0,5));
  Book( TH2F("wlep_neutrino_pTrel_delR"          ,"",100,0,500,200,0,5));
  Book( TH2F("wlep_muon_pTrel_delR"              ,"",100,0,500,200,0,5));
  Book( TH2F("wlep_neutrino_pTrel_wlep_muon_delR","",100,0,500,200,0,5));
  Book( TH2F("wlep_bjet_pTrel_wlep_muon_delR"    ,"",100,0,500,200,0,5));
  Book( TH2F("wlep_bjet_pTrel_wlep_neutrino_delR","",100,0,500,200,0,5));
  Book( TH2F("wlep_neutrino_pTrel_wlep_bjet_delR","",100,0,500,200,0,5));
  Book( TH2F("wlep_muon_pTrel_wlep_bjet_delR"    ,"",100,0,500,200,0,5));

  Book( TH2F("toplep_tophad_pTrel_delR"   ,"",100,0,100,200,0,5));
  Book( TH2F("tophad_muon_pTrel_delR"     ,"",100,0,100,200,0,5));
  Book( TH2F("tophad_neutrino_pTrel_delR" ,"",100,0,100,200,0,5));
  Book( TH2F("tophad_bjet_pTrel_delR"     ,"",100,0,100,200,0,5));


  Book( TH2F("razor_muon_included","",100,0,5000,100,0,1));
  Book( TH2F("razor_muon_excluded","",100,0,5000,100,0,1));

  Book( TH2F("razor2_muon_included","",100,0,5000,100,0,1));
  Book( TH2F("razor2_muon_excluded","",100,0,5000,100,0,1));


  Book( TH2F("razor_muon_included_blep","",100,0,5000,100,0,1));
  Book( TH2F("razor_muon_excluded_blep","",100,0,5000,100,0,1));

  Book( TH2F("razor2_muon_included_blep","",100,0,5000,100,0,1));
  Book( TH2F("razor2_muon_excluded_blep","",100,0,5000,100,0,1));

  Book( TH2F("chi2_pTrel_muon_bjet","",100,0,100,100,0,500));
  Book( TH2F("chi2_delR_muon_bjet" ,"",100,0,100,100,0,5)); 
  Book( TH2F("chi2_distance_max"    ,"",100,0,100,100,0,1000));

  Book( TH2F("chi2_delR_muon_neutrino" ,"",100,0,100,100,0,5));
  Book( TH2F("chi2_delR_bjet_neutrino" ,"",100,0,100,100,0,5)); 
  Book( TH2F("chi2_delR_muon_jet"      ,"",100,0,100,100,0,2));
  Book( TH2F("chi2_pTrel_muon_neutrino","",100,0,100,100,0,500));

  Book( TH2F("chi2_deltaRmax"  ,"",100,0,100,100,0,5));
  Book( TH2F("chi2_jade_max"   ,"",100,0,100,100,0,500));

  Book( TH1F("all_bjet_delR_ly","all #Delta R(b_{reco},b_{tag})",100,0,5));
  Book( TH1F("near_bjet_delR_ly","nearest #Delta R(b_{reco},b_{tag})",100,0,5));

  Book( TH1F("test","",100,0,1000));
  Book( TH1F("test2","",100,0,1000000));
  Book( TH1F("test_ly","",100,0,1000));
  Book( TH1F("test2_ly","",100,0,1000000));

  Book( TH2F("test2_pTrel_muon_bjet" ,"",100,0,500,100,0,1000000));
  Book( TH2F("test_pTrel_muon_bjet"  ,"",1000,0,500,1000,0,1000));

  // Book( TH3F("test_pTrel_muon_bjet_delR_muon_jet"  ,"",100,0,500,100,0,100,100,0,5));
  //Book( TH3F("test_pTrel_muon_bjet_delR_neutrino_jet"  ,"",100,0,500,100,0,100,100,0,5));


}

void DelRHists::Fill()
{
  // fill the histograms


  EventCalc* calc = EventCalc::Instance();
  TopFitCalc* fitcalc = TopFitCalc::Instance();
  double weight = calc -> GetWeight();

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
 

  double masstt = (toplep+tophad).isTimelike() ? (toplep+tophad).M() : -sqrt((toplep+tophad).M2()); 
  //double masstt = sqrt((toplep+tophad).M2()); 
  //if(!toplep.isTimelike()) return;	
  double toplepmass = (toplep).isTimelike() ? (toplep).M() : -sqrt((toplep).M2()); 

  
  
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

  double mrazorT = sqrt(met->pt()*(toplep.pt()+tophad.pt())- met3Vec*(plepT+phadT));
  double mrazor = sqrt(pow(toplep.E()+tophad.E(),2)-pow(toplep.Pz()+tophad.Pz(),2));

  //cout<<(toplep.pt()+tophad.pt())<<" "<<met->v4().E()<<endl;

  //cout<<mrazor<<" "<<mrazorT/mrazor<<" weight: "<<weight<<endl;


  ((TH2F*)Hist("razor_muon_included"))->Fill(mrazor,mrazorT/mrazor,weight);
  ((TH2F*)Hist("razor2_muon_included"))->Fill(mrazor,pow(mrazorT/mrazor,2),weight);
		       

  mrazorT = sqrt(met->pt()*(wlep.pt()+tophad.pt())- met3Vec*(pwlepT+phadT));
  mrazor = sqrt(pow(wlep.v4().E()+tophad.E(),2)-pow(wlep.v4().Pz()+tophad.Pz(),2));

  ((TH2F*)Hist("razor_muon_excluded"))->Fill(mrazor,mrazorT/mrazor,weight);
  ((TH2F*)Hist("razor2_muon_excluded"))->Fill(mrazor,pow(mrazorT/mrazor,2),weight);


  mrazorT = sqrt(met->pt()*((muon.v4()+bjet.v4()).pt()+tophad.pt())- met3Vec*(muonT+blepT+phadT));
  mrazor = sqrt(pow((muon.v4()+bjet.v4()).E()+tophad.E(),2)-pow((muon.v4()+bjet.v4()).Pz()+tophad.Pz(),2));

  ((TH2F*)Hist("razor_muon_included_blep"))->Fill(mrazor,mrazorT/mrazor,weight);
  ((TH2F*)Hist("razor2_muon_included_blep"))->Fill(mrazor,pow(mrazorT/mrazor,2),weight);
		       
  mrazorT = sqrt(met->pt()*(bjet.pt()+tophad.pt())- met3Vec*(blepT+phadT));
  mrazor = sqrt(pow(bjet.v4().E()+tophad.E(),2)-pow(bjet.v4().Pz()+tophad.Pz(),2));

  ((TH2F*)Hist("razor_muon_excluded_blep"))->Fill(mrazor,mrazorT/mrazor,weight);
  ((TH2F*)Hist("razor2_muon_excluded_blep"))->Fill(mrazor,pow(mrazorT/mrazor,2),weight);

  Hist("muon_jet_ly")->Fill(deltaRmin(&muon,antikjets),weight);
  Hist("muon_jet")->Fill(deltaRmin(&muon,antikjets),weight);

  Hist("muon_bjet_delR_ly")->Fill(deltaR(muon.v4(),bjet.v4()),weight);
  Hist("muon_bjet_delR")->Fill(deltaR(muon.v4(),bjet.v4()),weight);

  Hist("neutrino_bjet_delR_ly")->Fill(deltaR(neutrino.v4(),bjet.v4()),weight);
  Hist("neutrino_bjet_delR")->Fill(deltaR(neutrino.v4(),bjet.v4()),weight);

  Hist("neutrino_jet_delR_ly")->Fill(deltaRmin(&neutrino,antikjets),weight);
  Hist("neutrino_jet_delR")->Fill(deltaRmin(&neutrino,antikjets),weight);

  Hist("bjet_jet_delR_ly")->Fill(deltaRmin(&bjet,antikjets),weight);
  Hist("bjet_jet_delR")->Fill(deltaRmin(&bjet,antikjets),weight);

  Hist("muon_jet_long_ly")->Fill(deltaRmin(&muon,antikjets),weight);
  Hist("muon_jet_long")->Fill(deltaRmin(&muon,antikjets),weight);

  Hist("muon_bjet_delR_long_ly")->Fill(deltaR(muon.v4(),bjet.v4()),weight);
  Hist("muon_bjet_delR_long")->Fill(deltaR(muon.v4(),bjet.v4()),weight);

  Hist("neutrino_bjet_delR_long_ly")->Fill(deltaR(neutrino.v4(),bjet.v4()),weight);
  Hist("neutrino_bjet_delR_long")->Fill(deltaR(neutrino.v4(),bjet.v4()),weight);

  Hist("neutrino_jet_delR_long_ly")->Fill(deltaRmin(&neutrino,antikjets),weight);
  Hist("neutrino_jet_delR_long")->Fill(deltaRmin(&neutrino,antikjets),weight);

  Hist("bjet_jet_delR_long_ly")->Fill(deltaRmin(&bjet,antikjets),weight);
  Hist("bjet_jet_delR_long")->Fill(deltaRmin(&bjet,antikjets),weight);


  double deltaRmax = deltaR(muon.v4(),toplep);
  double deltaRminimum = deltaRmax;

  if(deltaRmax < deltaR(toplep,bjet.v4()))     deltaRmax = deltaR(toplep,bjet.v4());
  if(deltaRmax < deltaR(neutrino.v4(),toplep)) deltaRmax = deltaR(neutrino.v4(),toplep);
  //if(deltaRmax < deltaRmin(&neutrino,antikjets))   deltaRmax = deltaRmin(&neutrino,antikjets);
  //if(deltaRmax < deltaRmin(&bjet,antikjets))       deltaRmax = deltaRmin(&bjet,antikjets);

  if(deltaRminimum > deltaR(toplep,bjet.v4()))     deltaRminimum = deltaR(toplep,bjet.v4());
  if(deltaRminimum > deltaR(neutrino.v4(),toplep)) deltaRminimum = deltaR(neutrino.v4(),toplep);
  //if(deltaRminimum > deltaRmin(&neutrino,antikjets))   deltaRminimum = deltaRmin(&neutrino,antikjets);
  //if(deltaRminimum > deltaRmin(&bjet,antikjets))       deltaRminimum = deltaRmin(&bjet,antikjets);

  Hist("delRmax")->Fill(deltaRmax,weight);
  Hist("delRmin")->Fill(deltaRminimum,weight);
    
  Hist("deltaRmax_top_pt")->Fill(deltaRmax,toplep.pt());
  Hist("deltaRmin_top_pt")->Fill(deltaRminimum,toplep.pt());

  Hist("deltaRmax_top_mass")->Fill(deltaRmax,toplepmass);
  Hist("deltaRmax_top_mtt")->Fill(deltaRmax,(masstt));
  Hist("deltaRmax_top_pt_mtt")->Fill(deltaRmax,toplep.pt()/(masstt));
  Hist("deltaRmin_top_mass")->Fill(deltaRminimum,toplepmass);
 

  double distance_max = distance(muon.v4(),toplep,-1.,1.);
  double distance_minimum = distance_max;

  if(distance_max < distance(toplep,bjet.v4(),-1.,1.))     distance_max = distance(toplep,bjet.v4(),-1.,1.);
  if(distance_max < distance(neutrino.v4(),toplep,-1.,1.)) distance_max = distance(neutrino.v4(),toplep,-1.,1.);
  //if(distancemax < distancemin(&neutrino,antikjets))   distancemax = distancemin(&neutrino,antikjets);
  //if(distancemax < distancemin(&bjet,antikjets))       distancemax = distancemin(&bjet,antikjets);

  if(distance_minimum > distance(toplep,bjet.v4(),-1.,1.))     distance_minimum = distance(toplep,bjet.v4(),-1.,1.);
  if(distance_minimum > distance(neutrino.v4(),toplep,-1.,1.)) distance_minimum = distance(neutrino.v4(),toplep,-1.,1.);
  //if(distanceminimum > distancemin(&neutrino,antikjets))   distanceminimum = distancemin(&neutrino,antikjets);
  //if(distanceminimum > distancemin(&bjet,antikjets))       distanceminimum = distancemin(&bjet,antikjets);
    


  //cout<<distance_max<<" "<<distance_minimum<<endl;

  Hist("DistanceMax")->Fill(distance_max,weight);
  Hist("DistanceMin")->Fill(distance_minimum,weight);
  Hist("DistanceMin_ly")->Fill(distance_minimum,weight);

  Hist("DistanceMax_top_pt")->Fill(distance_max,toplep.pt());
  Hist("DistanceMin_top_pt")->Fill(distance_minimum,toplep.pt());

  Hist("DistanceMax_top_mass")->Fill(distance_max,toplepmass);
  Hist("DistanceMax_top_mtt")->Fill(distance_max,(masstt));
  Hist("DistanceMax_top_pt_mtt")->Fill(distance_max,toplep.pt()/(masstt));
  Hist("DistanceMin_top_mass")->Fill(distance_minimum,toplepmass);

  
  double jade_max = jade(muon.v4(),toplep,toplepmass);
  double jade_minimum = jade_max;

  if(jade_max < jade(toplep,bjet.v4(),toplepmass))     jade_max = jade(toplep,bjet.v4(),toplepmass);
  if(jade_max < jade(neutrino.v4(),toplep,toplepmass)) jade_max = jade(neutrino.v4(),toplep,toplepmass);
  //if(jademax < jademin(&neutrino,antikjets))   jademax = jademin(&neutrino,antikjets);
  //if(jademax < jademin(&bjet,antikjets))       jademax = jademin(&bjet,antikjets);

  if(jade_minimum > jade(toplep,bjet.v4(),toplepmass))     jade_minimum = jade(toplep,bjet.v4(),toplepmass);
  if(jade_minimum > jade(neutrino.v4(),toplep,toplepmass)) jade_minimum = jade(neutrino.v4(),toplep,toplepmass);
  //if(jademinimum > jademin(&neutrino,antikjets))   jademinimum = jademin(&neutrino,antikjets);
  //if(jademinimum > jademin(&bjet,antikjets))       jademinimum = jademin(&bjet,antikjets);
    


  //cout<<jade_max<<" "<<jade_minimum<<endl;

  Hist("JadeMax")->Fill(jade_max,weight);
  Hist("JadeMin")->Fill(jade_minimum,weight);

  Hist("JadeMax_top_pt")->Fill(jade_max,toplep.pt());
  Hist("JadeMin_top_pt")->Fill(jade_minimum,toplep.pt());

  Hist("JadeMax_top_mass")->Fill(jade_max,toplepmass);
  Hist("JadeMax_top_mtt")->Fill(jade_max,(masstt));
  Hist("JadeMax_top_pt_mtt")->Fill(jade_max,toplep.pt()/(masstt));
  Hist("JadeMin_top_mass")->Fill(jade_minimum,toplepmass);


  LorentzVector SumVec = toplep;


  for(int i =0; i < antikjets->size(); i++)
    if(deltaR(toplep,antikjets->at(i).v4())<1.6) SumVec -= antikjets->at(i).v4();

  Hist("SumVec")->Fill(SumVec.pt(),weight);

  Hist("eta_tops")->Fill((toplep+tophad).eta(),weight);
  Hist("deleta_top")->Fill(abs(toplep.eta()+tophad.eta()),weight);
  Hist("eta_tops_ly")->Fill((toplep+tophad).eta(),weight);
  Hist("deleta_top_ly")->Fill(abs(toplep.eta()+tophad.eta()),weight);

  Hist("muon_bjet_pTrel")->Fill(pTrel_Part(muon,bjet),weight);
  Hist("neutrino_bjet_pTrel")->Fill(pTrel_Part(neutrino,bjet),weight);
  Hist("muon_neutrino_pTrel")->Fill(pTrel_Part(muon,neutrino),weight);

  Hist("muon_bjet_pTrel_ly")->Fill(pTrel_Part(muon,bjet),weight);
  Hist("neutrino_bjet_pTrel_ly")->Fill(pTrel_Part(neutrino,bjet),weight);
  Hist("muon_neutrino_pTrel_ly")->Fill(pTrel_Part(muon,neutrino),weight);
  

  ((TH2F*)Hist("muon_bjet_pTrel_delR"))->Fill(pTrel_Part(muon,bjet),deltaR(muon.v4(),bjet.v4()),weight);
  ((TH2F*)Hist("muon_bjet_pTrel_muon_jet_delR"))->Fill(pTrel_Part(muon,bjet),deltaRmin(&muon,antikjets),weight);
  Hist("muon_bjet_pTrel_bjet_jet_delR")->Fill(pTrel_Part(muon,bjet),deltaRmin(&bjet,antikjets));
  Hist("muon_bjet_pTrel_muon_neutrino_delR")->Fill(pTrel_Part(muon,bjet),deltaR(muon.v4(),neutrino.v4()));
  Hist("neutrino_bjet_pTrel_muon_bjet_delR")->Fill(pTrel_Part(bjet,neutrino),deltaR(muon.v4(),bjet.v4()));
  Hist("neutrino_bjet_pTrel_neutrino_bjet_delR")->Fill(pTrel_Part(bjet,neutrino),deltaR(neutrino.v4(),bjet.v4()));


  Hist("wlep_bjet_pTrel_delR"                  )->Fill(pTrel_Part(wlep,bjet),deltaR(wlep.v4(),bjet.v4()));
  Hist("wlep_neutrino_pTrel_delR"              )->Fill(pTrel_Part(wlep,neutrino),deltaR(wlep.v4(),neutrino.v4()));
  Hist("wlep_muon_pTrel_delR"                  )->Fill(pTrel_Part(wlep,muon),deltaR(wlep.v4(),muon.v4()));
  Hist("wlep_neutrino_pTrel_wlep_muon_delR"    )->Fill(pTrel_Part(wlep,neutrino),deltaR(wlep.v4(),muon.v4()));
  Hist("wlep_bjet_pTrel_wlep_muon_delR"        )->Fill(pTrel_Part(wlep,bjet),deltaR(wlep.v4(),muon.v4()));
  Hist("wlep_bjet_pTrel_wlep_neutrino_delR"    )->Fill(pTrel_Part(wlep,bjet),deltaR(wlep.v4(),neutrino.v4()));
  Hist("wlep_neutrino_pTrel_wlep_bjet_delR"    )->Fill(pTrel_Part(wlep,neutrino),deltaR(wlep.v4(),bjet.v4()));
  Hist("wlep_muon_pTrel_wlep_bjet_delR"        )->Fill(pTrel_Part(wlep,muon),deltaR(wlep.v4(),bjet.v4()));

  Hist("toplep_tophad_pTrel_delR"              )->Fill(pTrel_Part(toplep,tophad),deltaR(toplep,tophad));
  Hist("tophad_muon_pTrel_delR"                )->Fill(pTrel_Part(muon.v4(),tophad),deltaR(muon.v4(),tophad));
  Hist("tophad_neutrino_pTrel_delR"            )->Fill(pTrel_Part(neutrino.v4(),tophad),deltaR(neutrino.v4(),tophad));
  Hist("tophad_bjet_pTrel_delR"                )->Fill(pTrel_Part(bjet.v4(),tophad),deltaR(bjet.v4(),tophad));

  bool nbtag = false;	  
  double delR_btag_bjet=100;	

  for(unsigned int i=0; i<bcc->jets->size(); ++i) {
    if(bcc->jets->at(i).btag_combinedSecondaryVertex()>0.244) nbtag = true;
    if(bcc->jets->at(i).btag_combinedSecondaryVertex()>0.679) nbtag= true;
    if(bcc->jets->at(i).btag_combinedSecondaryVertex()>0.898) nbtag= true;
    if(bcc->jets->at(i).btag_jetProbability()>0.275) nbtag= true;
    if(bcc->jets->at(i).btag_jetProbability()>0.545) nbtag= true;
    if(bcc->jets->at(i).btag_jetProbability()>0.790) nbtag= true;
  
    if(nbtag){
      Hist("all_bjet_delR_ly")->Fill(deltaR(bcc->jets->at(i).v4(),bjet.v4()),weight);
      if(delR_btag_bjet>deltaR(bcc->jets->at(i).v4(),bjet.v4())){ 
	delR_btag_bjet=deltaR(bcc->jets->at(i).v4(),bjet.v4());
  	Hist("near_bjet_delR_ly")->Fill(delR_btag_bjet,weight);
      }
    }

  }

  double mttbar_chi2 = 0;
  double exp_tophad_mass = 181;
  double exp_tophad_sig = 15;
  double exp_toplep_mass = 174;
  double exp_toplep_sig = 18;

  if(toplep.isTimelike())  
    mttbar_chi2 +=pow(toplepmass- exp_toplep_mass,2)/(exp_toplep_sig*exp_toplep_sig);
  else
    mttbar_chi2 +=pow(-sqrt(toplep.mass2())- exp_toplep_mass,2)/(exp_toplep_sig*exp_toplep_sig);
  if(tophad.isTimelike())
    mttbar_chi2 = pow(tophad.M()- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);
  else
    mttbar_chi2 = pow(-sqrt(tophad.mass2())- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);


  Hist("chi2_pTrel_muon_bjet")->Fill(mttbar_chi2,pTrel_Part(muon,bjet));
  Hist("chi2_delR_muon_bjet" )->Fill(mttbar_chi2,deltaR(muon.v4(),bjet.v4()));
  Hist("chi2_distance_max"    )->Fill(mttbar_chi2,distance_max);

  Hist("chi2_delR_muon_neutrino"  )->Fill(mttbar_chi2,deltaR(muon.v4(),neutrino.v4()));
  Hist("chi2_delR_bjet_neutrino"  )->Fill(mttbar_chi2,deltaR(bjet.v4(),neutrino.v4()));
  Hist("chi2_delR_muon_jet"       )->Fill(mttbar_chi2,deltaRmin(&muon,antikjets));
  Hist("chi2_pTrel_muon_neutrino" )->Fill(mttbar_chi2,pTrel(muon,neutrino));

  Hist("chi2_deltaRmax"    )->Fill(mttbar_chi2,deltaRmax);
  Hist("chi2_jade_max"     )->Fill(mttbar_chi2,jade_max);


  //TVector3 v1(toplep.px(),toplep.py(),toplep.pz());
  //TVector3 v2(tophad.px(),tophad.py(),tophad.pz());

  //double test = sqrt(pow(0.5*((toplep+tophad).M2() - exp_tophad_mass*exp_tophad_mass - exp_toplep_mass*exp_toplep_mass)+ v1*v2 - toplep.E()* tophad.E(),2));

  double test = sqrt(pow(0.5*(toplep.M2()+tophad.M2() - exp_tophad_mass*exp_tophad_mass - exp_toplep_mass*exp_toplep_mass),2))/1000;


 

  TVector3 v3(0,0,0);
  double lepE=0;
  TVector3 v4(0,0,0);
  double hadE =0;

  for(int i =0; i<antikjets->size(); ++i){
    TVector3 vjet(antikjets->at(i).v4().px(),antikjets->at(i).v4().py(),antikjets->at(i).v4().pz());

    if(deltaR(toplep,antikjets->at(i).v4())<1.5){
      v3 +=  vjet;
      lepE += antikjets->at(i).v4().E();
    }
    if(deltaR(tophad, antikjets->at(i).v4())<1.5){
      v4 +=   vjet;
      hadE += antikjets->at(i).v4().E();
    }

  }



  double test2 =  sqrt(pow(0.5*((toplep+tophad).M2() - exp_tophad_mass*exp_tophad_mass - exp_toplep_mass*exp_toplep_mass)+ v3*v4 - lepE* hadE,2));

  
  
  Hist("test")->Fill(test,weight);
  Hist("test2")->Fill(test2,weight);
  Hist("test_ly")->Fill(test,weight);
  Hist("test2_ly")->Fill(test2,weight);

  Hist("test2_pTrel_muon_bjet")->Fill(pTrel_Part(muon,bjet),test2);
  ((TH2F*)Hist("test_pTrel_muon_bjet" ))->Fill(pTrel_Part(muon,bjet),test,weight);

  //((TH3F*)Hist("test_pTrel_muon_bjet_delR_muon_jet"))->Fill(pTrel_Part(muon,bjet),test,deltaRmin(&muon,antikjets));
  //((TH3F*)Hist("test_pTrel_muon_bjet_delR_neutrino_jet"))->Fill(pTrel_Part(muon,bjet),test,deltaRmin(&neutrino,antikjets));


}

void DelRHists::Finish()
{
  // final calculations, like division and addition of certain histograms

}

