// Dear emacs, this is -*- c++ -*-

#include "include/TopFitCalc.h"
#include "TSystem.h"
#include <stdio.h>
#include <iostream>
#include <stdexcept>

TopFitCalc* TopFitCalc::m_instance = NULL;

TopFitCalc* TopFitCalc::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new TopFitCalc();
  }
  return m_instance;    
}



TopFitCalc::TopFitCalc() : m_logger("TopFitCalc")
{
  // constructor: initialise all variables

  
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;

  positiv = new TMinuit(5);

  //gSystem->Load("libMinuit");

  //Reset();



}

void TopFitCalc::Reset()
{



  b_Reconstruction = false;
}


TopFitCalc::~TopFitCalc()
{
  // default destructor
  
}


double DeltaPolarNeutrino(double PhiN, double metPx, double metPy, double PhiLep, double PtLep)
{
  using namespace std;

  double PyN;
  double PxN;
  const double mass_w = 80.399;
  
  double num = 10.e-7;	  
  	
  if(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN )) < num){
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*num);
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*num);
  }
  else{
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
  }

  return pow(PxN-metPx,2)+pow(PyN-metPy,2);

}

static void polarminuitfunc(int& nDim, double* gout, double& result, double par[], int flg){
  result = DeltaPolarNeutrino(par[0],par[1],par[2],par[3],par[4]);

}

std::vector<LorentzVector> TopFitCalc::NeutrinoFitPolar(const LorentzVector lepton, const LorentzVector met){


  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);
  
  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);
  
  const double mass_w = 80.399;
  
  
  
  double min = -2*PI;
  double max = 2*PI;
  double start = met.phi();
  double step = 10.e-5;
  
  double mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  double A = - (lepton_pT * lepton_pT);
  double B = mu * lepton.pz();
  double C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
  
  double discriminant = B * B - A * C;
  
  std::vector<LorentzVector> solutions;


  if (0 >= discriminant)
    {
    
    double resultPhi = 0;
    double error = 0;
    int ierflg;
    double* arg = new double[1];
      
    positiv->SetPrintLevel(-1); // -1 quiet, 0 normal, 1 verbose; Preset 0
    
    positiv->SetFCN(polarminuitfunc);
    
    positiv->DefineParameter(0,"PhiN",start, step,  min, max);
    positiv->DefineParameter(1,"metPx",met.px(),0,0,0);
    positiv->DefineParameter(2,"metPy",met.py(),0,0,0);
    positiv->DefineParameter(3,"PhiLep",lepton.phi(),0,0,0);
    positiv->DefineParameter(4,"PtLep",lepton.pt(),0,0,0);
    
    positiv->FixParameter(1);
    positiv->FixParameter(2);
    positiv->FixParameter(3);
    positiv->FixParameter(4);
    
    positiv->SetMaxIterations(500);
    
    
    arg[0]= 2;
    positiv->mnexcm("SET STR",arg,1,ierflg);
    
    positiv->Migrad();

    positiv->GetParameter(0,resultPhi,error);
    
    delete[] arg;
   
    if(resultPhi != resultPhi)
      {
	std::cerr << "neutrino phi is NAN " << std::endl;
      }      
   
    if(resultPhi > PI) resultPhi = resultPhi-2*PI;
    if(resultPhi < PI) resultPhi = resultPhi+2*PI;

 
    double PyN;
    double PxN;

    double num = 10.e-7;	  
    
    if(1-cos(delPhi(lepton.phi(), resultPhi)) < num){
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*num);
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*num);
    }
    else{
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*(1-cos(delPhi(lepton.phi(), resultPhi))));
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*(1-cos(delPhi(lepton.phi(), resultPhi))));
    }
  
    LorentzVectorXYZE neutrino_result(0,0,0,0);
    neutrino_result.SetPx(PxN);
    neutrino_result.SetPy(PyN);

    double pzfit =  lepton.pz()*neutrino_result.pt()/lepton.pt();
    
    LorentzVectorXYZE solution (0,0,0,0);
    solution.SetPx(PxN);
    solution.SetPy(PyN);
    solution.SetPz(pzfit);
    solution.SetE(toVector(solution).Mag());
      

    solutions.push_back(toPtEtaPhi(solution));

    }
  else
    {
      discriminant = sqrt(discriminant);
      
      LorentzVectorXYZE solution (0,0,0,0);
      solution.SetPx(met.Px());
      solution.SetPy(met.Py());
      solution.SetPz((-B - discriminant) / A);
      solution.SetE(toVector(solution).Mag());
      
      solutions.push_back(toPtEtaPhi(solution));
      
      LorentzVectorXYZE solution2 (0,0,0,0);
      solution2.SetPx(met.Px());
      solution2.SetPy(met.Py());
      solution2.SetPz((-B + discriminant) / A);
      solution2.SetE(toVector(solution2).Mag());
      
      solutions.push_back(toPtEtaPhi(solution2));
      
   }
  
  return solutions;
}


void TopFitCalc::CalculateSelection()
{

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* m_bcc = calc->GetBaseCycleContainer();

  //clear hypothesis list
  //m_bcc->recoHyps->clear();

  //find primary charged lepton
  Particle* lepton = &m_bcc->muons->at(0);

  std::vector<Jet>* antikjets = m_bcc->jets; 
  std::vector<TopJet>* cajets = m_bcc->topjets;

  //std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction(lepton->v4(),m_bcc->met->v4());
  std::vector<LorentzVector> neutrinos = NeutrinoFitPolar(lepton->v4(),m_bcc->met->v4());
 
  ReconstructionHypothesis hyp;

  hyp.set_lepton(*lepton);

  double cajets_pt   =  cajets->at(0).pt(); 
  int caposi_pt = 0;
  
  for(unsigned int m = 0; m<cajets->size(); ++m)
    {
      TopJet cajet = cajets->at(m);	 
      
      if(cajets_pt < cajet.pt())
	{
	  cajets_pt = cajet.pt();
	  caposi_pt = m;
	}
    }
  
 
  LorentzVector top_had = cajets->at(caposi_pt).v4();
  
  hyp.add_tophad_jet_index(caposi_pt);
  
  hyp.set_tophad_v4(top_had);

  for(unsigned int i = 0; i < neutrinos.size();i++){
  
    LorentzVector wboson_lep = lepton->v4()+neutrinos.at(i);
  
    for(unsigned int m=0; m<antikjets->size(); ++m){
      if( m_bcc->jets->at(m).pt()>50){
	LorentzVector top_lep = wboson_lep + m_bcc->jets->at(m).v4();
	hyp.set_blep_index(m);
	hyp.add_toplep_jet_index(m);
	hyp.set_blep_v4(m_bcc->jets->at(m).v4());

        hyp.set_neutrino_v4(neutrinos[i]);
	hyp.set_tophad_v4(top_had);
	hyp.set_toplep_v4(top_lep);

	m_bcc->recoHyps->push_back(hyp);
	}
    }
  }


}


void TopFitCalc::CalculateTopTag()
{

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* m_bcc = calc->GetBaseCycleContainer();
  m_bcc->recoHyps->clear();

  Particle* lepton = &m_bcc->muons->at(0);

  std::vector<Jet>* antikjets = m_bcc->jets; 
  std::vector<TopJet>* cajets = m_bcc->topjets; 
  

  //std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction(lepton->v4(),m_bcc->met->v4());
  std::vector<LorentzVector> neutrinos = NeutrinoFitPolar(lepton->v4(),m_bcc->met->v4());
 
  ReconstructionHypothesis hyp;

  hyp.set_lepton(*lepton);
   
  int best_neutrino = 0;

 
  int NCAJets = cajets->size();
  
  int caposi=-1;
  double mjet;
  int nsubjets;
  double mmin;	
  double deltaR_Lep_Tophad = .8;
  double deltaR_Jet_Tophad = 1.3;

  for(unsigned int m = 0; m<cajets->size(); ++m)
    {
      TopJet cajet = cajets->at(m);

      if(TopTag(cajet, mjet, nsubjets, mmin) && deltaR_Lep_Tophad < delR(cajet.v4(),lepton->v4())){
	caposi = m;
	//deltaR_Lep_Tophad = delR(cajet.v4(),lepton->v4());
	
	//if(caposi==-1) return;  
	LorentzVector top_had = cajet.v4();
 
	int n_jets = antikjets->size();
	if(n_jets>6) n_jets=6;
	int max_j = myPow(3, n_jets);
	

	for(unsigned int i = 0; i < neutrinos.size();++i){
	  
	  Particle wboson_lep;
	  wboson_lep.set_v4(lepton->v4()+neutrinos.at(i));
	  

	  for(int j=0; j<max_j; ++j){
	    LorentzVector top_lep(0,0,0,0);  
	    LorentzVector b_lep(0,0,0,0);  
	    int num = j;
    
	    for(unsigned int p=0; p<antikjets->size(); ++p){
	      if(delR(top_had,antikjets->at(p).v4())> deltaR_Jet_Tophad && num%3==0){
		b_lep = b_lep + antikjets->at(p).v4();
		top_lep = wboson_lep.v4() + b_lep;
		hyp.set_blep_index(p);
		hyp.set_blep_v4(b_lep);
		hyp.add_toplep_jet_index(p);
		hyp.add_tophad_jet_index(caposi);
		  
		hyp.set_neutrino_v4(neutrinos[i]);
		hyp.set_tophad_v4(top_had);
		hyp.set_toplep_v4(top_lep);
		  
		m_bcc->recoHyps->push_back(hyp);
	      }
	      num/=3;
	    }
	  }  
	}
      }
    }
}

