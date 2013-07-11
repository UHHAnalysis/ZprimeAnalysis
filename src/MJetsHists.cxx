#include "include/MJetsHists.h"
#include "include/SelectionModules.h"
#include "include/TopFitCalc.h"
#include "TH3F.h"
#include <iostream>
#include <cmath>
#include <stdio.h>

using namespace std;

MJetsHists::MJetsHists(const char* name, HypothesisDiscriminator *discr, int scalefactor) : BaseHists(name)
{
  // named default constructor
   m_discr = discr;
   scale = scalefactor;
}

MJetsHists::~MJetsHists()
{
  // default destructor, does nothing
}

void MJetsHists::Init()
{
  // book all histograms here
  Book( TH1F( "number","number of muons",5,-0.5,4.5));
  Book( TH1F( "pT","p_{T} muon [GeV]",1000/scale,0,900));
  Book( TH1F( "pT_ly","p_{T} muon[GeV]",1000/scale,0,900));
  Book( TH1F( "eta","#eta muon",1000/scale,-3,3));
  Book( TH1F( "eta_ly","#eta muon",1000/scale,-3,3));
  Book( TH1F( "phi","#phi muon",1000/scale,-PI,PI));
  Book( TH1F( "phi_ly","#phi muon",1000/scale,-PI,PI));
  Book( TH1F( "isolation","relIso muon",1000/scale,0,5));
  Book( TH1F( "isolation_ly","relIso muon",1000/scale,0,5));
  Book( TH1F("muon_pz_ly","p_{z} muon [GeV]",1000/scale,-1500,1500));

  Book( TH1F("muon_nearJet_pTrel","p_{T,rel} Jet-muon",1000/scale,0,1000));
  Book( TH1F("muon_nearJet_pTrel_ly","p_{T,rel} Jet-muon",1000/scale,0,1000));

  Book( TH1F("ak_pt","p_{t} AK-Jet [GeV]",1500/scale,0,1400));
  Book( TH1F("ak_pt1","p_{t} 1 leading AK-Jet [GeV]",1000/scale,0,1400));
  Book( TH1F("ak_pt2","p_{t} 2 leading AK-Jet [GeV]",1000/scale,0,1400));
  Book( TH1F("ak_pt3","p_{t} 3 leading AK-Jet [GeV]",1000/scale,0,1400));

  Book( TH1F("met_pT_ly" ,"p_{T} met [GeV]",1000/scale,0,1000));
  Book( TH1F("met_phi","#phi met",1000/scale,-PI,PI));

  Book( TH1F("HTmuon_pT_ly" ,"p_{T} H_{T,muon} [GeV]",1000/scale,0,1000));

  Book( TH1F("diff_met_neutrino_phi_ly","#Delta #phi met-Neutrino",1000/scale,0,1));

  Book( TH1F("neutrino_pT_ly" ,"p_{T} neutrino [GeV]",1000/scale,0,1000));
  Book( TH1F("neutrino_phi","#phi neutrino",1000/scale,-PI,PI));
  Book( TH1F("neutrino_eta","#eta neutrino",1000/scale,-3,3));
  Book( TH1F("neutrino_pz_ly" ,"p_{z} neutrino [GeV]",1000/scale,-2500,2500));

  Book( TH1F("wboson_pT_ly" ,"p_{T} W [GeV]",1000/scale,0,2000));
  Book( TH1F("wboson_phi","#phi W" ,1000/scale,-PI,PI));
  Book( TH1F("wboson_eta","#eta W" ,1000/scale,-PI,PI));
  Book( TH1F("wboson_pz_ly" ,"p_{z} W [GeV]",1000/scale,-2800,2800));

  Book( TH1F("wboson_MC_phi","#phi MC  W",1000/scale,-PI,PI));
  Book( TH1F("wboson_MC_pt_ly" ,"p_{T} MC W [GeV]" ,1000/scale,0,1000));
  Book( TH1F("wboson_MC_eta","#eta MC W" ,1000/scale,-PI,PI));
  Book( TH1F("wboson_MC_E_ly"  ,"E MC W [GeV]",1000/scale,0,2800));
  Book( TH1F("wboson_MC_mass_ly"  ,"MC W_{m} [GeV]",1000/scale,-20,200));

  Book( TH1F("wboson_E_ly"   ,"W Energy [GeV]",1000/scale,100,3000));
  Book( TH1F("wboson_mass_ly","W_{m} [GeV]"   ,1000/scale,60,200));
  Book( TH1F("wboson_mass","W_{m} [GeV]"   ,1000/scale,60,200));

  Book( TH1F("toplep_pt_ly"  ,"p_{T} Top_{lep} [GeV]",1000/scale,0,2000));
  Book( TH1F("toplep_phi"    ,"#phi Top_{lep}",1000/scale,-PI,PI));
  Book( TH1F("toplep_eta"    ,"#eta Top_{lep}",1000/scale,-PI,PI));
  Book( TH1F("toplep_E_ly"   ,"E Top_{lep}",1000/scale,0,2800));
  Book( TH1F("toplep_mass"   ,"M^{Top}_{lep} [GeV]",2000/scale,0,1000));
  Book( TH1F("toplep_mass_ly","M^{Top}_{lep} [GeV]",2000/scale,0,1000));

  Book( TH1F("toplep_MC_phi","#phi MC  Top_{lep}",1000/scale,-PI,PI));
  Book( TH1F("toplep_MC_pt_ly" ,"p_{T} MC Top_{lep} [GeV]" ,1000/scale,0,2000));
  Book( TH1F("toplep_MC_eta","#eta MC Top_{lep}" ,1000/scale,-PI,PI));
  Book( TH1F("toplep_MC_E_ly"  ,"E MC Top_{lep} [GeV]",1000/scale,0,2800));
  Book( TH1F("toplep_MC_mass_ly"  ,"MC Top^{lep}_{m} [GeV]",1000/scale,0,2000));

  Book( TH1F("deltaR_top_MC_reco_ly"     ,"#Delta R MC-Reco  Top_{lep}",1000/scale,0,6));
  Book( TH1F("deltaR_wboson_MC_reco_ly"  ,"#Delta R MC-Reco  W_{lep}"  ,1000/scale,0,6));
  Book( TH1F("deltaR_neutrino_Mc_reco_ly","#Delta R MC-Reco  Neutrino" ,1000/scale,0,6));

  Book( TH1F("near_phi_delphi","nearest Jet to W in #phi with #Delta #phi ",1000/scale,-PI,PI));
  Book( TH1F("near_phi_delphi_ly","nearest Jet to W in #phi with #Delta #phi ",1000/scale,-PI,PI));

  Book( TH1F("near_phi_delR","nearest Jet to W in #phi with #Delta R  ",1000/scale,-PI,PI));
  Book( TH1F("near_phi_delR_ly","nearest Jet to W in #phi with #Delta R  ",1000/scale,-PI,PI));

  Book( TH1F("near_delR_delR","nearest Jet to W in #Delta R with #Delta R",1000/scale,0,7));
  Book( TH1F("near_delR_delR_ly","nearest Jet to W in #Delta R with #Delta R",1000/scale,0,7));

  Book( TH1F("delphi_delphi","#Delta #phi W-Jets with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_delphi_ly","#Delta #phi W-Jets with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_delR","#Delta #phi W-Jets with #Delta R",1000/scale,0,PI));
  Book( TH1F("delphi_delR_ly","#Delta #phi W-Jets with #Delta R",1000/scale,0,PI));
  
  Book( TH1F("number_Jets","number of anti k_{T}-Jets",14,1.5,15.5));
  Book( TH1F("number_Jets_ly","number of anti k_{T}-Jets",14,1.5,15.5));
  
  Book( TH1F("deltaR_wboson_jet_ly","#Delta R W-anti-kT Jet ",1000/scale,0,2*PI));
  Book( TH1F("deltaR_wboson_jet",   "#Delta R W-anti-kT Jet",1000/scale,0,2*PI));

  Book( TH1F("number_CAJets","number of CAJets",9,-0.5,8.5));
  Book( TH1F("number_CAJets_ly","number of CAJets",9,-0.5,8.5));

  Book( TH1F("ca_pt","p_{t} CA-Jet [GeV]",1500/scale,0,1400));
  Book( TH1F("ca_pt1","p_{t} 1 leading CA-Jet [GeV]",1000/scale,0,1400));
  Book( TH1F("ca_pt2","p_{t} 2 leading CA-Jet [GeV]",1000/scale,0,1400));
  Book( TH1F("ca_pt3","p_{t} 3 leading CA-Jet [GeV]",1000/scale,0,1400));
  Book( TH1F("ca_pt4","p_{t} 4 leading CA-Jet [GeV]",1000/scale,0,1400));

  Book( TH1F("ca_phi","#phi CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_phi1","#phi 1 leading CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_phi2","#phi 2 leading CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_phi3","#phi 3 leading CA-Jet",200/scale,-PI,PI));
  Book( TH1F("ca_phi4","#phi 4 leading CA-Jet",200/scale,-PI,PI));

  Book( TH1F("ca_eta" ,"#eta CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_eta1","#eta 1 leading CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_eta2","#eta 2 leading CA-Jet",1000/scale,-PI,PI));
  Book( TH1F("ca_eta3","#eta 3 leading CA-Jet",200/scale,-PI,PI));
  Book( TH1F("ca_eta4","#eta 4 leading CA-Jet",200/scale,-PI,PI));

  Book( TH1F("ca_px","p_{x} CA-Jet [GeV]",1000/scale,-1000,1000));
  Book( TH1F("ca_py","p_{y} CA-Jet [GeV]",1000/scale,-1000,1000));
  Book( TH1F("ca_pz","p_{z} CA-Jet [GeV]",1000/scale,-20,20));

  Book( TH1F("number_ca_subjets","number of CA-Subjets",9,-0.5,8.5));

  Book( TH1F("ca_subjet_pt1","p_{t} 1 leading CA-Subjet [GeV]",1000/scale,0,1000));
  Book( TH1F("ca_subjet_pt2","p_{t} 2 leading CA-Subjet [GeV]",1000/scale,0,500));
  Book( TH1F("ca_subjet_pt3","p_{t} 3 leading CA-Subjet [GeV]",1000/scale,0,200));
  Book( TH1F("ca_subjet_pt4","p_{t} 4 leading CA-Subjet [GeV]",1000/scale,0,200));


  Book( TH1F("ca_subjet_phi1","#phi 1 leading CA-Subjet",1000/scale,-PI,PI));
  Book( TH1F("ca_subjet_phi2","#phi 2 leading CA-Subjet",1000/scale,-PI,PI));
  Book( TH1F("ca_subjet_phi3","#phi 3 leading CA-Subjet",200/scale,-PI,PI));
  Book( TH1F("ca_subjet_phi4","#phi 4 leading CA-Subjet",200/scale,-PI,PI));


  Book( TH1F("ca_subjet_eta1","#eta 1 leading CA-Subjet",1000/scale,-PI,PI));
  Book( TH1F("ca_subjet_eta2","#eta 2 leading CA-Subjet",1000/scale,-PI,PI));
  Book( TH1F("ca_subjet_eta3","#eta 3 leading CA-Subjet",200/scale,-PI,PI));
  Book( TH1F("ca_subjet_eta4","#eta 4 leading CA-Subjet",200/scale,-PI,PI));


  Book( TH1F("delphi_ca_delphi","#Delta #phi CAJets to -T_{lep} with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_ca_delphi_ly","#Delta #phi CAJets to -T_{lep} with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_ca_delR","#Delta #phi CAJets to -T_{lep} with #Delta R",1000/scale,0,PI));
  Book( TH1F("delphi_ca_delR_ly","#Delta #phi #mu-CAJets with #Delta R",1000/scale,0,PI));
  Book( TH1F("near_ca_phi_delphi_ly","nearest CAJet to -T_{lep} in #phi with #Delta #phi ",1000/scale,0,PI));
  Book( TH1F("near_ca_phi_delphi","nearest CAJet to -T_{lep} in #phi with #Delta #phi ",1000/scale,0,PI));
  Book( TH1F("near_ca_phi_delR","nearest CAJet to -T_{lep} in #phi with #Delta R  ",1000/scale,0,PI));
  Book( TH1F("near_ca_phi_delR_ly","nearest CAJet to -T_{lep} in #phi with #Delta R  ",1000/scale,0,PI));
 

  Book( TH1F("number_antik_withoutca","number of anti k_{T}-Jets without used ones",10,-0.5,9.5));
  Book( TH1F("delphi_antik_delphi","#Delta #phi anti k_{T}-Jets to -T_{lep} with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_antik_delphi_ly","#Delta #phi anti k_{T}-Jets to -T_{lep} with #Delta #phi",1000/scale,0,PI));
  Book( TH1F("delphi_antik_delR","#Delta #phi anti k_{T}-Jets to -T_{lep} with #Delta R",1000/scale,0,PI));
  Book( TH1F("delphi_antik_delR_ly","#Delta #phi anti k_{T}-Jets to -T_{lep} with #Delta R",1000/scale,0,PI));
  Book( TH1F("near_antik_phi_delphi","nearest anti k_{T}-Jet to -T_{lep} in #phi with #Delta #phi ",1000/scale,-PI,PI));
  Book( TH1F("near_antik_phi_delphi_ly","nearest anti k_{T}-Jet to -T_{lep} in #phi with #Delta #phi ",1000/scale,-PI,PI));
  Book( TH1F("near_antik_phi_delR","nearest anti k_{T}-Jet to -T_{lep} in #phi with #Delta R  ",1000/scale,-PI,PI));
  Book( TH1F("near_antik_phi_delR_ly","nearest anti k_{T}-Jet to -T_{lep} in #phi with #Delta R  ",1000/scale,-PI,PI));

  Book( TH1F("topmass_had_ly","Top^{had}_{m} [GeV]",2000/scale,0,500));
  Book( TH1F("topmass_had","Top^{had}_{m} [GeV]",2000/scale,0,500));
  Book( TH1F("tophad_phi"    ,"#phi Top^{had}",1000/scale,-PI,PI));
  Book( TH1F("tophad_pt_ly"  ,"p_{T} Top^{had} [GeV]" ,1000/scale,0,2000));
  Book( TH1F("tophad_eta"    ,"#eta Top^{had}" ,1000/scale,-PI,PI));
  Book( TH1F("tophad_E_ly"   ,"E Top^{had} [GeV]",1000/scale,0,2800));

  Book( TH1F("tophad_MC_phi","#phi MC  Top^{had}",1000/scale,-PI,PI));
  Book( TH1F("tophad_MC_pt_ly" ,"p_{T} MC Top^{had} [GeV]" ,1000/scale,0,2000));
  Book( TH1F("tophad_MC_eta","#eta MC Top^{had}" ,1000/scale,-PI,PI));
  Book( TH1F("tophad_MC_E_ly"  ,"E MC Top^{had} [GeV]",1000/scale,0,2800));
  Book( TH1F("tophad_MC_mass_ly","Top^{had}_{m} [GeV]",1000/scale,0,600));

  Book( TH1F("tophad_subjet1_phi"    ,"#phi Top^{had} 1 leading Subjet",1000/scale,-PI,PI));
  Book( TH1F("tophad_subjet1_pt_ly"  ,"p_{T} Top^{had} 1 leading Subjet [GeV]" ,1000/scale,0,1400));
  Book( TH1F("tophad_subjet1_eta"    ,"#eta Top^{had} 1 leading Subjet" ,1000/scale,-PI,PI));
  Book( TH1F("tophad_subjet1_mass"   ,"Mass Top^{had} 1 leading Subjet [GeV]" ,1000/scale,0,400));

  Book( TH1F("tophad_subjet2_phi"    ,"#phi Top^{had} 2 leading Subjet",1000/scale,-PI,PI));
  Book( TH1F("tophad_subjet2_pt_ly"  ,"p_{T} Top^{had} 2 leading Subjet [GeV]" ,1000/scale,0,800));
  Book( TH1F("tophad_subjet2_eta"    ,"#eta Top^{had} 2 leading Subjet" ,1000/scale,-PI,PI));
  Book( TH1F("tophad_subjet2_mass"   ,"Mass Top^{had} 2 leading Subjet [GeV]" ,1000/scale,0,300));

  Book( TH1F("tophad_subjet3_phi"    ,"#phi Top^{had} 3 leading Subjet",700/scale,-PI,PI));
  Book( TH1F("tophad_subjet3_pt_ly"  ,"p_{T} Top^{had} 3 leading Subjet" ,700/scale,0,500));
  Book( TH1F("tophad_subjet3_eta"    ,"#eta Top^{had} 3 leading Subjet" ,700/scale,-PI,PI));
  Book( TH1F("tophad_subjet3_mass"   ,"Mass Top^{had} 3 leading Subjet" ,700/scale,0,100));
 
  Book( TH1F("tophad_subjet4_phi"    ,"#phi Top^{had} 4 leading Subjet",200/scale,-PI,PI));
  Book( TH1F("tophad_subjet4_pt_ly"  ,"p_{T} Top^{had} 4 leading Subjet [GeV]" ,200/scale,0,200));
  Book( TH1F("tophad_subjet4_eta"    ,"#eta Top^{had} 4 leading Subjet" ,200/scale,-PI,PI));
  Book( TH1F("tophad_subjet4_mass"   ,"Mass Top^{had} 4 leading Subjet [GeV]" ,200/scale,0,100));

  Book( TH1F("tophad_subjet_all_mass","Mass Top^{had} Sum of all Subjets [GeV]" ,2000/scale,0,500));
  Book( TH1F("tophad_subjet_all_mass_ly","Mass Top^{had} Sum of all Subjets [GeV]" ,2000/scale,0,500));
  Book( TH1F("tophad_subjet_all_mass_short","Mass Top^{had} Sum of all Subjets [GeV]" ,40,140,250));

  Book( TH1F("toptag_min_pt","Min() of 2 Subjets [GeV]" ,40,0,200));

  Book( TH1F("tophad_delR_antikt"          ," Top^{had} nearest anti-kt #Delta R" ,500/scale,0,10));
  Book( TH1F("tophad_delR_antikt_ly"       ," Top^{had} nearest anti-kt #Delta R" ,500/scale,0,10));
  Book( TH1F("tophad_delR_antikt_clean"    ," Top^{had} nearest anti-kt with distance 0.8 #Delta R " ,500/scale,0,10));
  Book( TH1F("tophad_delR_antikt_clean_ly" ," Top^{had} nearest anti-kt with distance 0.8 #Delta R" ,500/scale,0,10));

  Book( TH1F("tophad_MC_phi","#phi MC  Top^{had}",1000/scale,-PI,PI));
  Book( TH1F("tophad_MC_pt_ly" ,"p_{T} MC Top^{had} [GeV]" ,1000/scale,0,2000));
  Book( TH1F("tophad_MC_eta","#eta MC Top^{had}" ,1000/scale,-PI,PI));
  Book( TH1F("tophad_MC_E_ly"  ,"E MC Top^{had} [GeV]",1000/scale,0,2800));
  Book( TH1F("tophad_MC_mass_ly","Top^{had}_{m} [GeV]",1000/scale,0,600));
 
  Book( TH1F("Mttbar_ly","M_{tt} [GeV]",2000/scale,0,5000));
  Book( TH1F("Mttbar","M_{tt} [GeV]",2000/scale,0,5000));//geaendert fuer die untersuchung des binnings
  Book( TH1F("Mttbar_gen","M_{tt, gen} [GeV]",2000/scale,0,5000));

  Book( TH1F("Mttbar_short","M_{tt} [GeV]",2000/scale,2000,5000));

  Book( TH1F("Ptttbar_ly","P_{T tt} [GeV]",2000/scale,0,2000));
  Book( TH1F("Mttbar_chi2","M_{tt}#chi^2 [GeV]",1000/scale,0,150));
  Book( TH1F("M_ttbar_gen_rec","#frac{M_{tt,gen}-M_{tt,rec}}{M_{tt,gen}}",1000/scale,-2,2));
  Book( TH1F("toplep_pt_gen_rec","Top_{lep} #frac{p_{T,gen}-p_{T,reco}}{p_{T,gen}}", 1000/scale,-5,5));  

  Book( TH1F( "pt_toplep_norm", "P_{T}^{top,lep} [GeV/c]", 100, -5,5) );
  Book( TH1F( "pt_tophad_norm", "P_{T}^{top,had} [GeV/c]", 100, -5,5) );

  Book( TH1F( "eta_toplep_norm", "P_{T}^{top,lep} [GeV/c]", 100, -5,5 ) );
  Book( TH1F( "eta_tophad_norm", "P_{T}^{top,had} [GeV/c]", 100, -5,5 ) );

  Book( TH1F( "phi_toplep_norm", "P_{T}^{top,lep} [GeV/c]", 100, -5,5 ) );
  Book( TH1F( "phi_tophad_norm", "P_{T}^{top,had} [GeV/c]", 100, -5,5 ) );



  Book( TH1F( "Cut_Variable", "M_{Top} [GeV/c]", 100, -5000,5000 ) );



  //Book( TH1F("simple_mttbar","simple m_{t#bar{t}}",100,0,2500));

  //Book( TH1F("percentage_delphi_delR_ly","percentage #Delta #phi/#Delta R ",100,-PI,PI));


  Book( TH2F("M_ttbar_rec_vs_M_ttbar_gen","M_{t#bar{t}}^{rec} [GeV/c^{2}] vs M_{t#bar{t}}^{gen} [GeV/c^{2}]",600/scale,0,4000,600/scale,-10,10));
  Book( TH2F( "Pt_ttbar_rec_vs_Pt_ttbar_gen", "P_{T,t#bar{t}}^{rec} [GeV/c] vs P_{T,t#bar{t}}^{gen} [GeV/c]", 600/scale, 0, 600 ,60, -200, 200));

  Book( TH2F("wbosonlep_pt_reco_gen","",600/scale,-10,1000,60,-10,2));
  Book( TH2F("wbosonlep_phi_reco_gen","",600/scale,-PI,PI,60,-7,7));
  Book( TH2F("wbosonlep_eta_reco_gen","",600/scale,-4,4,60,-4,4));
 
  Book( TH2F("wbosonhad_pt_reco_gen","",1000/scale,0,1000,100,0,1000));
	     				    
  Book( TH2F("bjetlep_pt_reco_gen"   ,"",600/scale,0,1000,60,-6,2));
  Book( TH2F("bjetlep_phi_reco_gen"  ,"",600/scale,-PI,PI,60,-7,7));
  Book( TH2F("bjetlep_eta_reco_gen"  ,"",600/scale,-4,4  ,60,-6,6));


  Book( TH2F("bjethad_pt_reco_gen"  ,"",1000/scale,0,1000,100,0,1000));
	     				    

  Book( TH2F("toplep_pt_reco_gen"    ,"",600/scale,0,1500,60,-2,2));
  Book( TH2F("tophad_pt_reco_gen"    ,"",600/scale,0,1500,60,-2,2));
  Book( TH2F("toplep_phi_reco_gen"   ,"",600/scale,-PI,PI,60,-7,7));
  Book( TH2F("tophad_phi_reco_gen"   ,"",600/scale,-PI,PI,60,-7,7));
  Book( TH2F("toplep_eta_reco_gen"   ,"",600/scale,-4,4  ,60,-4,4));
  Book( TH2F("tophad_eta_reco_gen"   ,"",600/scale,-4,4  ,60,-4,4));
  Book( TH2F("toplep_m_reco_gen"     ,"",600/scale,0,2000,60,-20,30));
  Book( TH2F("tophad_m_reco_gen"     ,"",600/scale,0,2000,60,-20,30));

  Book( TH2F("cajet_muon_pt","",1000/scale,0,3000,100,0,3000));

  Book( TH2F("delR_Iso","",100,0,5,100,0,2));
  Book( TH2F("delR_pTrel","",100,0,5,1000,0,1000));
  Book( TH2F("Iso_pTrel","",100,0,5,100,0,100));
  Book( TH2F("delR_pTnearJet","",100,0,1,100,0,1000));
  Book( TH2F("Iso_pTnearJet","",100,0,20,100,0,1000));
  Book( TH2F("muon_bjetpt","",100,0,1000,100,0,1000));
  Book( TH2F("delR_pTmuon","",100,0,5,100,0,400));

  Book( TH2F("neutrino_delR_pTrel","",100,0,2,100,0,300));
  Book( TH2F("neutrino_delR_pTnearJet","",100,0,2,100,0,1000));
  Book( TH2F("neutrino_bjet_pT","",100,0,1000,100,0,1000));
  Book( TH2F("neutrino_delR_pT","",100,0,5,100,0,1000));

  // Book( TH3F("neutrino_delR_pTrel_eta","",100,0,5,100,0,300,100,-4,4));
  //Book( TH3F("neutrino_delR_pTrel_phi","",100,0,5,100,0,300,100,-4,4));

  Book( TH2F("muon_neutrino_delR_pTrel","",100,0,2,100,0,300));
  Book( TH2F("neutrino_muon_delR_pTrel","",100,0,2,100,0,300));
}


void MJetsHists::Fill()
{

 


   // important: get the event weight
   
 

  

  EventCalc* calc = EventCalc::Instance();
  TopFitCalc* fitcalc = TopFitCalc::Instance();
  double weight = calc -> GetWeight();

  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();




  ReconstructionHypothesis* hyp = m_discr->GetBestHypothesis();

  std::vector<TopJet>* cajets = calc->GetCAJets();
  int NCAJets = cajets->size();

  double mttbar_gen = 0;
  double ptttbar_gen = 0;
  double topmasshad = 0;



  int NMuons = bcc->muons->size();
  Hist("number")-> Fill(NMuons,weight);
 
  std::vector<Jet>* antikjets = calc->GetJets();

  for(unsigned int i=0; i< bcc->muons->size(); ++i)
    {
      Muon muon  = bcc->muons->at(i);
     
      


      double tag_mjet = 0;
      int tag_nsubjets = 0;
      double tag_mmin = 0;

      for(unsigned int i = 0; i<cajets->size(); ++i){
	if(TopTag(cajets->at(i),tag_mjet,tag_nsubjets,tag_mmin)) 
	  Hist("toptag_min_pt")->Fill(tag_mmin,weight);
      }



     
      LorentzVector neutrino = hyp->neutrino_v4();

      Particle neutrino_particle;
      neutrino_particle.set_v4(neutrino);

      std::vector<LorentzVector> neutrinos = fitcalc->NeutrinoFitPolar(muon.v4(),bcc->met->v4());
      std::vector<LorentzVector> neutrino_std = calc->NeutrinoReconstruction(muon.v4(),bcc->met->v4());

      //if(neutrinos.at(0).phi() < 3.14 && neutrinos.at(0).phi() > -3.14 ) return; 
      //if(neutrinos.at(0).phi() > 3.0 || neutrinos.at(0).phi() < -3. ) return; 
      //if(pow(neutrinos.at(0).px()-bcc->met->v4().px(),2)+pow(neutrinos.at(0).py()-bcc->met->v4().py(),2)>2500) return;
      /*
      if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_tauhad){
	if(calc->GetTTbarGen()->Neutrino().v4().px()-bcc->met->v4().px()>20 || calc->GetTTbarGen()->Neutrino().v4().py()-bcc->met->v4().py()> 20 ) return;
      }
      else
	return;
      */

	//cout<<"neutrino phi "<<neutrinos.at(0).phi()<<" lepton phi "<< muon.phi()<<" first "<< cos(neutrinos.at(0).phi())/(1-cos(muon.phi()-neutrinos.at(0).phi()))<<" second "<<sin(neutrinos.at(0).phi())/(1-cos(muon.phi()-neutrinos.at(0).phi()))<<" chi2 1 "<< pow(neutrinos.at(0).px()-bcc->met->v4().px(),2)<<" chi2 2 "<<pow(neutrinos.at(0).py()-bcc->met->v4().py(),2) <<endl;


	// cout<<neutrino_std.size()<<endl;
	
	
	
      Hist("number_CAJets_ly")->Fill(NCAJets,weight);
      Hist("number_CAJets"   )->Fill(NCAJets,weight);
      
      Hist("met_pT_ly" )->Fill(bcc->met->v4().pt(),weight);
      Hist("met_phi"   )->Fill(bcc->met->v4().phi(),weight); 

      Hist("HTmuon_pT_ly" )->Fill(bcc->met->v4().pt()+muon.pt(),weight);
    
      // Muon Control Hists similar to the ones in Muons

      Hist("pT"             )->Fill(muon.pt(),weight);
      Hist("pT_ly"          )->Fill(muon.pt(),weight);
      Hist("eta"            )->Fill(muon.eta(),weight);
      Hist("eta_ly"         )->Fill(muon.eta(),weight);
      Hist("phi"            )->Fill(muon.phi(),weight);
      Hist("phi_ly"         )->Fill(muon.phi(),weight);
      Hist("isolation"      )->Fill(muon.relIso(),weight);
      Hist("isolation_ly"   )->Fill(muon.relIso(),weight);

      Hist("muon_nearJet_pTrel")->Fill(pTrel(&muon,antikjets),weight);
      Hist("muon_nearJet_pTrel_ly")->Fill(pTrel(&muon,antikjets),weight);
  
      //Hist("TwoDCut")->Fill(deltaRmin(&muon,antikjets),pTrel(&muon,antikjets ));

      Hist("delR_Iso")->Fill(deltaRmin(&muon,antikjets),muon.relIso());
      ((TH2F*)Hist("delR_pTrel"))->Fill(deltaRmin(&muon,antikjets),pTrel(&muon,antikjets),weight);
      Hist("Iso_pTrel")->Fill(muon.relIso(),pTrel(&muon,antikjets));
      Hist("Iso_pTnearJet")->Fill(muon.relIso(),nextJet(&muon,antikjets)->pt());
  
      Hist("neutrino_delR_pTrel")->Fill(deltaRmin(&neutrino_particle,antikjets),pTrel(&neutrino_particle,antikjets));

      Hist("muon_neutrino_delR_pTrel")->Fill(deltaRmin(&muon,antikjets),pTrel(&neutrino_particle,antikjets));
      Hist("neutrino_muon_delR_pTrel")->Fill(deltaRmin(&neutrino_particle,antikjets),pTrel(&muon,antikjets));

      //((TH3F*) Hist("neutrino_delR_pTrel_eta"))->Fill(deltaRmin(&neutrino_particle,antikjets),pTrel(&neutrino_particle,antikjets),neutrino_particle.eta());
      //((TH3F*) Hist("neutrino_delR_pTrel_phi"))->Fill(deltaRmin(&neutrino_particle,antikjets),pTrel(&neutrino_particle,antikjets),neutrino_particle.phi());

      Hist("neutrino_delR_pTnearJet")->Fill(deltaRmin(&neutrino_particle,antikjets),nextJet(&neutrino_particle,antikjets)->pt());
      Hist("neutrino_delR_pT")->Fill(deltaRmin(&neutrino_particle,antikjets),neutrino_particle.pt());



      double m_phi = muon.phi();    


      Hist("muon_pz_ly")->Fill(muon.v4().pz(),weight);
      
      
      
      int NJets = antikjets->size();

      Hist("number_Jets_ly")->Fill(NJets,weight);
      Hist("number_Jets"   )->Fill(NJets,weight);
      
      if(antikjets)
	{
	  Particle wboson;
	  wboson.set_v4(hyp->wlep_v4());
	  

	  //Neutrino Control Hists
	  
	  Hist("neutrino_pT_ly" )->Fill(neutrino.pt(),weight);
	  Hist("neutrino_eta")->Fill(neutrino.eta(),weight);
	  Hist("neutrino_phi")->Fill(neutrino.phi(),weight);  
	  Hist("neutrino_pz_ly" )->Fill(neutrino.pz(),weight);
	  Hist("diff_met_neutrino_phi_ly")->Fill(fitcalc->delPhi(neutrino.phi(),bcc->met->v4().phi()),weight); 
	  
	  //WBoson Hists 

	  Hist("wboson_pT_ly" )->Fill(wboson.v4().pt(),weight);
	  Hist("wboson_phi")->Fill(wboson.v4().phi(),weight);
	  Hist("wboson_eta")->Fill(wboson.v4().eta(),weight);  
	  Hist("wboson_pz_ly" )->Fill(wboson.v4().pz(),weight);

	  Hist("wboson_E_ly"     )->Fill(wboson.v4().E(),weight);
	  Hist("wboson_mass_ly"  )->Fill(sqrt(wboson.v4().M2()),weight);
	  Hist("wboson_mass"     )->Fill(sqrt(wboson.v4().M2()),weight);

	  double mjets_delphi = wboson.deltaPhi(antikjets->at(0));
	  double mjets_delR = wboson.deltaR(antikjets->at(0)); 
	  int bjetposi_delphi = 0;
	  unsigned int bjetposi_delR = 0;

	

	  for(unsigned int j=0; j<antikjets->size(); ++j)
	    {
	      Jet jet = antikjets->at(j);
	      
 
	      if(mjets_delphi > wboson.deltaPhi(jet)) 
		{
		  mjets_delphi = wboson.deltaPhi(jet);
		  bjetposi_delphi = j;
		}

	      if(mjets_delR > wboson.deltaR(jet))
		{
		  mjets_delR = wboson.deltaR(jet);
		  bjetposi_delR = j;
		}	    


	    }



	  Hist("near_phi_delphi"     )->Fill(antikjets->at(bjetposi_delphi).phi(),weight);
	  Hist("near_phi_delR"       )->Fill(antikjets->at(bjetposi_delR).phi(),weight);
	 
	  Hist("near_delR_delR"      )->Fill( mjets_delR,weight);
	  Hist("near_delR_delR_ly"   )->Fill( mjets_delR,weight);


	  Hist("delphi_delphi"       )->Fill(mjets_delphi,weight);
	  Hist("delphi_delR"         )->Fill(wboson.deltaPhi(antikjets->at(bjetposi_delR)),weight);
	  
	  Hist("near_phi_delphi_ly"  )->Fill(antikjets->at(bjetposi_delphi).phi(),weight);
	  Hist("near_phi_delR_ly"    )->Fill(antikjets->at(bjetposi_delR).phi(),weight);
	  
	  Hist("delphi_delphi_ly"    )->Fill(mjets_delphi,weight);
	  Hist("delphi_delR_ly"      )->Fill(wboson.deltaPhi(antikjets->at(bjetposi_delR)),weight);
	 
	  Hist("deltaR_wboson_jet_ly")->Fill(deltaR(wboson.v4(),antikjets->at(bjetposi_delR).v4()),weight);
	  Hist("deltaR_wboson_jet"   )->Fill(deltaR(wboson.v4(),antikjets->at(bjetposi_delR).v4()),weight);



	  Particle bjet;
	  bjet.set_v4(hyp->blep_v4()); 
		      
		     
	  Hist("muon_bjetpt")->Fill(muon.pt(),bjet.pt());

	  Hist("neutrino_bjet_pT")->Fill(neutrino_particle.pt(),bjet.pt());

 
	  Particle top_lep; 
	  top_lep.set_v4(hyp->toplep_v4());


	  Hist("toplep_pt_ly"   )->Fill(top_lep.v4().pt(),weight);
	  Hist("toplep_phi"  )->Fill(top_lep.v4().phi(),weight);
	  Hist("toplep_eta"  )->Fill(top_lep.v4().eta(),weight);
	  Hist("toplep_E_ly"    )->Fill(top_lep.v4().E(),weight);

	  Hist("toplep_mass"    )->Fill(sqrt(top_lep.v4().M2()),weight);
	  Hist("toplep_mass_ly" )->Fill(sqrt(top_lep.v4().M2()),weight);



	  if(!calc->IsRealData()){
	    if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_tauhad){
	      GenParticle GenNeutrino = calc->GetTTbarGen()->Neutrino();
	      GenParticle GenB_lep = calc->GetTTbarGen()->BLep();
	      GenParticle GenW_lep = calc->GetTTbarGen()->WLep();
	      GenParticle GenTop_lep = calc->GetTTbarGen()->TopLep();
	      GenParticle GenTop_had = calc->GetTTbarGen()->TopHad();
      		     
	      Hist("wboson_MC_phi"   )->Fill(GenW_lep.v4().phi(),weight);
	      Hist("wboson_MC_pt_ly"    )->Fill(GenW_lep.v4().pt() ,weight);
	      Hist("wboson_MC_eta"   )->Fill(GenW_lep.v4().eta(),weight);
	      Hist("wboson_MC_E_ly"     )->Fill(GenW_lep.v4().E() ,weight);
	      Hist("wboson_MC_mass_ly"  )->Fill(GenW_lep.v4().isTimelike()? GenW_lep.v4().M(): -sqrt(GenW_lep.v4().M2()) ,weight);
	      
	      Hist("toplep_MC_phi"   )->Fill(GenTop_lep.v4().phi(),weight);
	      Hist("toplep_MC_pt_ly"    )->Fill(GenTop_lep.v4().pt() ,weight);
	      Hist("toplep_MC_eta"   )->Fill(GenTop_lep.v4().eta(),weight);
	      Hist("toplep_MC_E_ly"     )->Fill(GenTop_lep.v4().E() ,weight);
	      Hist("toplep_MC_mass_ly"  )->Fill(GenTop_lep.v4().isTimelike()?GenTop_lep.v4().M():-sqrt(GenTop_lep.v4().M2()) ,weight);

	      Hist("tophad_MC_phi"   )->Fill(GenTop_had.v4().phi(),weight);
	      Hist("tophad_MC_pt_ly"    )->Fill(GenTop_had.v4().pt() ,weight);
	      Hist("tophad_MC_eta"   )->Fill(GenTop_had.v4().eta(),weight);
	      Hist("tophad_MC_E_ly"     )->Fill(GenTop_had.v4().E() ,weight);
	      Hist("tophad_MC_mass_ly"  )->Fill(GenTop_had.v4().isTimelike()?GenTop_had.v4().M():-sqrt(GenTop_had.v4().M2()) ,weight);

	      Hist("deltaR_top_MC_reco_ly"      )->Fill(deltaR(top_lep.v4(),GenTop_lep.v4()));
	      Hist("deltaR_wboson_MC_reco_ly"   )->Fill(deltaR(wboson.v4(),GenW_lep.v4()));
	      Hist("deltaR_neutrino_Mc_reco_ly" )->Fill(deltaR(neutrino,GenNeutrino.v4()));

	      Hist("bjetlep_pt_reco_gen" )->Fill(bjet.pt() ,(bjet.pt()-GenB_lep.pt())/bjet.pt());
	      Hist("bjetlep_phi_reco_gen")->Fill(bjet.phi(),bjet.phi()-GenB_lep.phi());
	      Hist("bjetlep_eta_reco_gen")->Fill(bjet.eta(),bjet.eta()-GenB_lep.eta());

	      Hist("wbosonlep_pt_reco_gen" )->Fill(wboson.pt() ,(wboson.pt()- GenW_lep.pt())/wboson.pt());
	      Hist("wbosonlep_phi_reco_gen")->Fill(wboson.phi(),wboson.phi()- GenW_lep.phi());
	      Hist("wbosonlep_eta_reco_gen")->Fill(wboson.eta(),wboson.eta()- GenW_lep.eta());

	      Hist("toplep_m_reco_gen"  )->Fill(sqrt(GenTop_lep.v4().M2()),(sqrt(top_lep.v4().M2())-sqrt(GenTop_lep.v4().M2()))/sqrt(GenTop_lep.v4().M2())); 
	      Hist("toplep_pt_reco_gen" )->Fill(GenTop_lep.pt(),(top_lep.pt()-GenTop_lep.pt())/GenTop_lep.pt()); 
	      Hist("toplep_phi_reco_gen")->Fill(GenTop_lep.phi(),top_lep.phi()-GenTop_lep.phi()); 
	      Hist("toplep_eta_reco_gen")->Fill(GenTop_lep.eta(),top_lep.eta()-GenTop_lep.eta()); 

	      Hist("pt_toplep_norm")->Fill((calc->GetTTbarGen()->TopLep().v4().pt()-hyp->toplep_v4().pt())/calc->GetTTbarGen()->TopLep().v4().pt(),weight);
	      Hist("pt_tophad_norm")->Fill((calc->GetTTbarGen()->TopHad().v4().pt()-hyp->tophad_v4().pt())/calc->GetTTbarGen()->TopHad().v4().pt(),weight);

	      Hist("phi_toplep_norm")->Fill((calc->GetTTbarGen()->TopLep().v4()-hyp->toplep_v4()).phi(),weight);
	      Hist("phi_tophad_norm")->Fill((calc->GetTTbarGen()->TopHad().v4()-hyp->tophad_v4()).phi(),weight);
	      
	      Hist("eta_toplep_norm")->Fill((calc->GetTTbarGen()->TopLep().v4()-hyp->toplep_v4()).eta(),weight);
	      Hist("eta_tophad_norm")->Fill((calc->GetTTbarGen()->TopHad().v4()-hyp->tophad_v4()).eta(),weight);


	      Hist("Mttbar_gen")->Fill((GenTop_lep.v4()+GenTop_had.v4()).isTimelike()?(GenTop_lep.v4()+GenTop_had.v4()).M():-sqrt((GenTop_lep.v4()+GenTop_had.v4()).M2()),weight);

	    }
	    else{
      //this is an misreconstructed event, there should be a plot like number ?
      //Hist("missreco_number")->Fill(1,weight);
      
	    } 
	  }

          //LorentzVector MJets_delphi = top_lep.v4();//antikjets->at(bjetposi_delphi).v4()+muon.v4();
	  //LorentzVector MJets_delR   = top_lep.v4();

	  double phisum_delphi = -top_lep.phi();
	  //double phisum_delR   = -MJets_delR.phi();

	  
	  Particle antiMJets_delR;
	  antiMJets_delR.set_pt(-1*top_lep.pt());
	  antiMJets_delR.set_eta(-1*top_lep.eta());
	  antiMJets_delR.set_phi(-1*top_lep.phi());
	 
	  if(NCAJets>0)
	    {

	      double cajets_delphi =  fitcalc->delPhi(phisum_delphi, cajets->at(0).phi());
	      double cajets_delR   =  antiMJets_delR.deltaR(cajets->at(0));

	      int caposi_delphi = 0;
	      int caposi_delR = 0;

	      for(unsigned int m = 0; m<cajets->size(); ++m)
		{
		  Jet cajet = cajets->at(m);
		  //if( muon.phi()>=cajet.phi()-.5 && muon.phi()<=cajet.phi()+.5 ) Hist("cajet_muon_pt")->Fill(muon.pt(),cajet.pt());//cajet.set_v4(cajet.v4()-muon.v4());		 

		  if(cajets_delphi > fitcalc->delPhi(phisum_delphi, cajets->at(m).phi()))
		    {
		      cajets_delphi = fitcalc->delPhi(phisum_delphi, cajets->at(m).phi());
		      caposi_delphi = m;
		    }
		 
		  if(cajets_delR > antiMJets_delR.deltaR(cajet))
		    {
		      cajets_delR = antiMJets_delR.deltaR(cajet);
		      caposi_delR = m;
		    }
		}


	      double mass_w = 80.399;

	      Particle top_had;
	      top_had.set_v4(hyp->tophad_v4()); 	      
	      
	      

	      

	      Hist("Cut_Variable")->Fill(muon.v4().E()*bjet.v4().E()+bjet.v4().E()*neutrino.E()-muon.v4().px()*bjet.v4().px()+bjet.v4().px()*neutrino.px()-muon.v4().py()*bjet.v4().py()+bjet.v4().py()*neutrino.py()-muon.v4().pz()*bjet.v4().pz()+bjet.v4().pz()*neutrino.pz()-12071.76,weight);
	      
	      int caposi_subjets = -1;
	      double mjet_subjets;
	      int nsubjets_subjets;
	      double mmin_subjets;	
	      double deltaR_Lep_Tophad = .8;
	      
	      for(unsigned int m = 0; m<cajets->size(); ++m)
		{
		  TopJet cajet_subjets = cajets->at(m);

		  if(TopTag(cajet_subjets, mjet_subjets, nsubjets_subjets, mmin_subjets) && deltaR_Lep_Tophad < fitcalc->delR(cajet_subjets.v4(),muon.v4())){
		    caposi_subjets = m;
		    deltaR_Lep_Tophad = fitcalc->delR(cajet_subjets.v4(),muon.v4());
		  }
		}
	      

	      //cout<<hyp->tophad_jets_indices().at(0)<<" "<<caposi_subjets<<endl;
	      

	      if(caposi_subjets != -1){
		TopJet topjet = cajets->at(caposi_subjets);
		
		
		//if(TopTag(topjet, mjet, nsubjets, mmin)) cerr<<"TopJet without Toptag choosed"<<endl;
		
		
		std::vector<Particle> topsub = topjet.subjets();
		sort(topsub.begin(),topsub.end(),HigherPt());
		LorentzVector sumsubv4 (0,0,0,0);
		
		for(unsigned int y = 0; y < topsub.size(); y++){
		  
		  TString hname = TString::Format("tophad_subjet%d_pt_ly", y+1);
		  Hist(hname)->Fill(topsub.at(y).pt(),weight);
		  TString hname_eta = TString::Format("tophad_subjet%d_eta", y+1);
		  Hist(hname_eta)->Fill(topsub.at(y).eta(),weight);
		  TString hname_phi = TString::Format("tophad_subjet%d_phi", y+1);
		  Hist(hname_phi)->Fill(topsub.at(y).phi(),weight);
		  TString hname_mass = TString::Format("tophad_subjet%d_mass", y+1);
		  double submass = topsub.at(y).v4().isTimelike() ?  topsub.at(y).v4().M() : -sqrt(topsub.at(y).v4().M2());  
		  Hist(hname_mass)->Fill(submass,weight);
		  sumsubv4 += topsub.at(y).v4();
		}
		
		double sumsubmass = sumsubv4.isTimelike() ? sumsubv4.M() : -sqrt(sumsubv4.M2());
		Hist("tophad_subjet_all_mass")->Fill(sumsubmass,weight);
		Hist("tophad_subjet_all_mass_ly")->Fill(sumsubmass,weight);
		Hist("tophad_subjet_all_mass_short")->Fill(sumsubmass,weight);
		
	      }
	      
	      double mttbar_rec = 0;
	      double mttbar_chi2 = 0;
	      double exp_tophad_mass = 181;
	      double exp_tophad_sig = 15;
	      double exp_toplep_mass = 174;
	      double exp_toplep_sig = 18;
	      /*
	      if(sumsubv4.isTimelike())
		mttbar_chi2 = pow(sumsubv4.M()- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);
	      else
		mttbar_chi2 = pow(-sqrt(sumsubv4.mass2())- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);
	      */
	      if(top_lep.v4().isTimelike())  
		mttbar_chi2 +=pow(top_lep.v4().M()- exp_toplep_mass,2)/(exp_toplep_sig*exp_toplep_sig);
	      else
		mttbar_chi2 +=pow(-sqrt(top_lep.v4().mass2())- exp_toplep_mass,2)/(exp_toplep_sig*exp_toplep_sig);
	      /*
	      if((sumsubv4+top_lep.v4()).isTimelike())
		mttbar_rec  = (sumsubv4+top_lep.v4()).M();
	      else
		mttbar_rec = -sqrt( (sumsubv4+top_lep.v4()).mass2());
	      */
	      
	      if(top_had.v4().isTimelike())
		mttbar_chi2 = pow(top_had.v4().M()- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);
	      else
		mttbar_chi2 = pow(-sqrt(top_had.v4().mass2())- exp_tophad_mass,2)/(exp_tophad_sig*exp_tophad_sig);
 
	      if((top_had.v4()+top_lep.v4()).isTimelike())
		mttbar_rec  = (top_had.v4()+top_lep.v4()).M();
	      else
		mttbar_rec = -sqrt( (top_had.v4()+top_lep.v4()).mass2());
	      


	      //if(mttbar_chi2>10) return;

	      Hist("Mttbar_ly"      )->Fill(mttbar_rec,weight);
	      Hist("Mttbar"         )->Fill(mttbar_rec,weight);
	      Hist("Mttbar_short"   )->Fill(mttbar_rec,weight);
	      Hist("Ptttbar_ly"     )->Fill((top_had.v4()+top_lep.v4()).pt(),weight);
	      Hist("Mttbar_chi2"    )->Fill(mttbar_chi2,weight);


	      
	      double tophad_antikt_delR_clean = 9999;
	      double tophad_antikt_delR       = 9999;


	      for(unsigned int p=0; p<antikjets->size(); ++p)
		{
		  Jet antiktjet = antikjets->at(p);
		  
		  if(top_had.deltaR(antiktjet)> 0.8 && bjetposi_delR !=  p && top_had.deltaR(antiktjet)<tophad_antikt_delR_clean)
		    {
		      tophad_antikt_delR_clean = top_had.deltaR(antiktjet);
		    }
		  if(bjetposi_delR !=  p && top_had.deltaR(antiktjet)<tophad_antikt_delR)
		    {
		      tophad_antikt_delR = top_had.deltaR(antiktjet);
		    }
		}


	      if(tophad_antikt_delR_clean < 9999) Hist("tophad_delR_antikt_clean"  )->Fill(tophad_antikt_delR_clean,weight);
	      if(tophad_antikt_delR       < 9999) Hist("tophad_delR_antikt"        )->Fill(tophad_antikt_delR,weight);
	      if(tophad_antikt_delR_clean < 9999) Hist("tophad_delR_antikt_clean_ly"  )->Fill(tophad_antikt_delR_clean,weight);
	      if(tophad_antikt_delR       < 9999) Hist("tophad_delR_antikt_ly"        )->Fill(tophad_antikt_delR,weight);

	      if(calc->GetGenParticles()){
		mttbar_gen =  (calc->GetTTbarGen()->Top().v4() + calc->GetTTbarGen()->Antitop().v4()).isTimelike()? ( calc->GetTTbarGen()->Top().v4() + calc->GetTTbarGen()->Antitop().v4()).M():-sqrt(( calc->GetTTbarGen()->Top().v4() + calc->GetTTbarGen()->Antitop().v4()).M2());
		ptttbar_gen = ( calc->GetTTbarGen()->Top().v4() + calc->GetTTbarGen()->Antitop().v4()).Pt();



		Hist("M_ttbar_gen_rec"  )->Fill((mttbar_gen-mttbar_rec)/mttbar_gen,weight);
                //Hist("toplep_pt_gen_rec")->Fill((calc->GetTTbarGen()->TopLep().pt()-top_lep.pt())/calc->GetTTbarGen()->TopLep().pt(),weight);
		Hist("M_ttbar_rec_vs_M_ttbar_gen"  )->Fill(mttbar_gen,(mttbar_rec-mttbar_gen)/mttbar_gen);
		Hist("Pt_ttbar_rec_vs_Pt_ttbar_gen")->Fill(ptttbar_gen,(((top_had.v4()+top_lep.v4()).pt())-ptttbar_gen)/ptttbar_gen);
	      }
                 
	      if(!bcc->isRealData){
	       if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_tauhad){
	         GenParticle GenTopHad = calc->GetTTbarGen()->TopHad();
 	         Hist("tophad_m_reco_gen")->Fill(sqrt(GenTopHad.v4().M2()),(sqrt(top_had.v4().M2())-sqrt(GenTopHad.v4().M2()))/sqrt(GenTopHad.v4().M2()));
 	         Hist("tophad_pt_reco_gen")->Fill(GenTopHad.pt(),(top_had.pt()-GenTopHad.pt())/GenTopHad.pt());
 	         Hist("tophad_phi_reco_gen")->Fill(GenTopHad.phi(),top_had.phi()-GenTopHad.phi());
 	         Hist("tophad_eta_reco_gen")->Fill(GenTopHad.eta(),top_had.eta()-GenTopHad.eta());
	       }
	      }

	      //if(!bcc->isRealData) Hist("ca_pz")->Fill(cajets->at(caposi_delR).genjet_v4().pz(),weight);
	      Hist("ca_pz")->Fill(cajets->at(caposi_delR).genjet_v4().pz(),weight);

	      Hist("near_ca_phi_delphi")->Fill(cajets->at(caposi_delphi).phi(),weight);
	      Hist("near_ca_phi_delR"  )->Fill(cajets->at(caposi_delR).phi(),weight);
	      Hist("delphi_ca_delphi"  )->Fill(cajets_delphi,weight);
	      Hist("delphi_ca_delR"    )->Fill(cajets->at(caposi_delR).deltaPhi(antiMJets_delR),weight);

	      Hist("near_ca_phi_delphi_ly")->Fill(cajets->at(caposi_delphi).phi(),weight);
	      Hist("near_ca_phi_delR_ly"  )->Fill(cajets->at(caposi_delR).phi(),weight);
	      Hist("delphi_ca_delphi_ly"  )->Fill(cajets_delphi,weight);
	      Hist("delphi_ca_delR_ly"    )->Fill(cajets->at(caposi_delR).deltaPhi(antiMJets_delR),weight);

	           
	      double minusmjets_delphi = 999.;//antiMJets_delR.deltaPhi(antikjets->at(0));
	      double minusmjets_delR = 999.;//antiMJets_delR.deltaR(antikjets->at(0));
	      
	      int jetposi_delphi = -1;
	      int jetposi_delR   = -1;
	      int antiktnum      =  0;



	      for(unsigned int p=0; p<antikjets->size(); ++p)
		{
		  Jet minusjet = antikjets->at(p);
		  
		  if(cajets->at(caposi_delR).deltaR(minusjet)> 0.8 && bjetposi_delR !=  p)
		    {
		      antiktnum++; 
		      
		      if(minusmjets_delphi> antiMJets_delR.deltaPhi( minusjet))
			{
			  minusmjets_delphi = antiMJets_delR.deltaPhi(minusjet);
			  jetposi_delphi = p;
			}
		  
		      if(minusmjets_delR > antiMJets_delR.deltaR(minusjet))
			{
			  minusmjets_delR = antiMJets_delR.deltaR(minusjet);
			  jetposi_delR = p;
			}	     
		    }
		}

	    
	      Hist("number_antik_withoutca")->Fill(antiktnum,weight);

	      if(jetposi_delphi>-1) Hist("delphi_antik_delphi"     )->Fill(minusmjets_delphi,weight) ;
	      if(jetposi_delphi>-1) Hist("delphi_antik_delphi_ly"  )->Fill(minusmjets_delphi,weight) ;
	      if(jetposi_delR  >-1) Hist("delphi_antik_delR"       )->Fill(antikjets->at(jetposi_delR).deltaPhi(antiMJets_delR),weight);
	      if(jetposi_delR  >-1) Hist("delphi_antik_delR_ly"    )->Fill(antikjets->at(jetposi_delR).deltaPhi(antiMJets_delR),weight);
	      if(jetposi_delphi>-1) Hist("near_antik_phi_delphi"   )->Fill(antikjets->at(jetposi_delphi).phi(),weight);
	      if(jetposi_delphi>-1) Hist("near_antik_phi_delphi_ly")->Fill(antikjets->at(jetposi_delphi).phi(),weight);
	      if(jetposi_delR  >-1) Hist("near_antik_phi_delR"     )->Fill(antikjets->at(jetposi_delR).phi(),weight);
	      if(jetposi_delR  >-1) Hist("near_antik_phi_delR_ly"  )->Fill(antikjets->at(jetposi_delR).phi(),weight);

	      Hist("tophad_phi"      )->Fill(top_had.v4().phi(),weight);
	      Hist("tophad_pt_ly"    )->Fill(top_had.v4().pt() ,weight);
	      Hist("tophad_eta"      )->Fill(top_had.v4().eta(),weight);
	      Hist("tophad_E_ly"     )->Fill(top_had.v4().E() ,weight);
	      Hist("topmass_had_ly"  )->Fill(top_had.v4().isTimelike()? top_had.v4().M() : -sqrt(top_had.v4().M2()),weight);
	      Hist("topmass_had"  )->Fill(top_had.v4().isTimelike()? top_had.v4().M() : -sqrt(top_had.v4().M2()),weight);
	    }
	 }
	    

      }
  
  sort(antikjets->begin(),antikjets->end(),HigherPt());
  
  for (unsigned int i = 0; i<3; i++)
    {
      if(unsigned(antikjets->size())>i){
	TString hname = TString::Format("ak_pt%d", i+1);
	Hist(hname)->Fill(antikjets->at(i).pt(),weight);
      }
    }
  
  for(unsigned int i = 0; i < antikjets->size(); i++)
    Hist("ak_pt")->Fill(antikjets->at(i).pt(),weight);
  

  sort(cajets->begin(),cajets->end(),HigherPt());
      
  for (unsigned int i = 0; i<4; i++)
    {
      if( unsigned(NCAJets) > i)
	{
	  
	  TopJet CAJet =  cajets->at(i); 
	  if(bcc->muons->at(0).phi() == CAJet.phi()) CAJet.set_v4(cajets->at(i).v4()-bcc->muons->at(0).v4()); 
	  
	  TString hname = TString::Format("ca_pt%d", i+1);
	  Hist(hname)->Fill(CAJet.pt(),weight);
	  TString hname_eta = TString::Format("ca_eta%d", i+1);
	  Hist(hname_eta)->Fill(CAJet.eta(),weight);
	  TString hname_phi = TString::Format("ca_phi%d", i+1);
	  Hist(hname_phi)->Fill(CAJet.phi(),weight);
	    
	}
    }
  

  for(int i = 0; i < NCAJets; i++){
    
    Hist("ca_pt")->Fill(cajets->at(i).pt(),weight);
    Hist("ca_phi")->Fill(cajets->at(i).phi(),weight);
    Hist("ca_eta")->Fill(cajets->at(i).eta(),weight);
    Hist("ca_px")->Fill(cajets->at(i).v4().px(),weight);
    Hist("ca_py")->Fill(cajets->at(i).v4().py(),weight);
    //Hist("ca_pz")->Fill(cajets->at(i).pz(),weight);
    

  }
  
  for(unsigned int p = 0; p < cajets->size();p++){
    std::vector<Particle> casubjet = cajets->at(p).subjets();
    Hist("number_ca_subjets")->Fill(casubjet.size(),weight);
    sort(casubjet.begin(),casubjet.end(),HigherPt());
    
    for(unsigned int m = 0; m < casubjet.size(); m++){
      
      TString hname = TString::Format("ca_subjet_pt%d", m+1);
      Hist(hname)->Fill(casubjet.at(m).pt(),weight);
      TString hname_eta = TString::Format("ca_subjet_eta%d", m+1);
      Hist(hname_eta)->Fill(casubjet.at(m).eta(),weight);
      TString hname_phi = TString::Format("ca_subjet_phi%d", m+1);
      Hist(hname_phi)->Fill(casubjet.at(m).phi(),weight);
            
    }	     

  }



}

void MJetsHists::Finish()
{
 
}

