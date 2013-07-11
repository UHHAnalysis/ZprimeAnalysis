// $Id: ZPrimeTopTaggerCycle.cxx,v 1.0 2012/11/09 10:06:24 gonvaq Exp $

#include <iostream>

using namespace std;

// Local include(s):
//#include "JetEffiHists.h"
#include "include/ZPrimeTopTaggerCycle.h"
#include "SelectionModules.h"
#include "include/ExampleHists.h"
#include "HypothesisHists.h"
#include "JetHists.h"
#include "MuonHists.h"
#include "TopJetHists.h"
#include "EventHists.h"
#include "include/TopEffiHists.h"
#include "include/MJetsHists.h"
#include "include/NeutrinoHists.h"
#include "TopFitCalc.h"
#include "include/DelRHists.h"



ClassImp( ZPrimeTopTaggerCycle );

ZPrimeTopTaggerCycle::ZPrimeTopTaggerCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(250.);
  m_mttgencut = false;
  DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
}

ZPrimeTopTaggerCycle::~ZPrimeTopTaggerCycle()
{
  // destructor
}

void ZPrimeTopTaggerCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void ZPrimeTopTaggerCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void ZPrimeTopTaggerCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{

   // declaration of histograms and selections

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

  // -------------------- set up the selections ---------------------------

  Selection* BSel = new Selection( "BSelection");
  BSel->addSelectionModule(new NBTagSelection(1)); //at least one b tag

  Selection* NoBSel = new Selection( "NoBSelection");
  NoBSel->addSelectionModule(new NBTagSelection(0,0)); //no b tags

  Selection* MuonSel = new Selection("MuonSelection");
  MuonSel->addSelectionModule(new NMuonSelection (1, 1));


  Selection * TwoDMuonSel = new Selection("TwoDSelection");
  TwoDMuonSel->addSelectionModule(new TwoDCutMuon(.5,25));

  static Selection* mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
  if (((id.GetVersion() == "TTbar_0to700") || (id.GetVersion() == "TTbar") )  ) {
    m_logger << INFO << "Applying mttbar generator cut from 0 to 700 GeV." << SLogger::endmsg;
    mttbar_gen_selection->addSelectionModule(new MttbarGenCut(0,700));
    mttbar_gen_selection->EnableSelection();
  }
  else {
    m_logger << INFO << "Disabling mttbar generator cut." << SLogger::endmsg;
    mttbar_gen_selection->DisableSelection();
  }
  /*  
  Selection* HCALlaser = new Selection("HCAL_laser_events");
  HCALlaser->addSelectionModule(new EventFilterSelection(m_filter_file.c_str()) );
  */




  Selection* TopSel = new Selection("TopSelection");
  //DO NOT use trigger selection in PROOF mode at the moment
  //TopSel->addSelectionModule(new TriggerSelection("HLT_PFJet320_v"));
  TopSel->addSelectionModule(new NJetSelection (1, int_infinity(), 150., 2.5));
  TopSel->addSelectionModule(new NJetSelection (2, int_infinity(), 50.));
  TopSel->addSelectionModule(new METCut(20));
  TopSel->addSelectionModule(new HTlepCut(150));
  TopSel->addSelectionModule(new NTopJetSelection(1,int_infinity(),0,2.5));
  //TopSel->addSelectionModule(new NTopTagSelection(1,int_infinity()));
  TopSel->addSelectionModule(new CAAntiktJetSelection(1,2,1,int_infinity(),int_infinity()));
  //TopSel->addSelectionModule(new TopTagAntiktJetSelection(1,2,0.8,int_infinity(),int_infinity()));
  //TopSel->addSelectionModule(new METCut(50));
  
  Selection* TopTagSel = new Selection("TopTagSelection");
  TopTagSel->addSelectionModule(new NTopTagSelection(1));
  TopTagSel->addSelectionModule(new TopTagAntiktJetSelection(1,2,0.8,int_infinity(),int_infinity()));
  //Selection* TTreco = new Selection("TTreco");

  m_chi2discr = new Chi2Discriminator();
  m_tagchi2discr = new Chi2Discriminator();
 
  //RegisterSelection(TTreco);

  //RegisterSelection(HCALlaser);
  RegisterSelection(mttbar_gen_selection);
  RegisterSelection(BSel);
  RegisterSelection(NoBSel);
  RegisterSelection(TopSel);
  RegisterSelection(TopTagSel);
  RegisterSelection(MuonSel);	     
  RegisterSelection(TwoDMuonSel);

  // ---------------- set up the histogram collections --------------------
  
  // histograms without any cuts
  RegisterHistCollection( new ExampleHists("NoCuts") );

  //static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();

  //RegisterHistCollection( new HypothesisHists("TopTag_NoCuts", m_chi2discr ) );

  //histograms with and without b tagging
  //RegisterHistCollection( new ExampleHists("BTag") );
  //RegisterHistCollection( new ExampleHists("NoBTag") );


  /*
  RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
  RegisterHistCollection( new HypothesisHists("Chi2_NoBTag", m_chi2discr ) );

  // histograms after the top selection
  RegisterHistCollection( new ExampleHists("TopSel") );
  RegisterHistCollection( new HypothesisHists("Chi2_TopSel", m_chi2discr ) );

  
  //histograms to see if everything works
  RegisterHistCollection( new JetHists("Jets"));
  */
  RegisterHistCollection( new MuonHists("Muons"));
  RegisterHistCollection( new MuonHists("TagMuons"));
  /*
  RegisterHistCollection( new TopJetHists("TopJets"));
  RegisterHistCollection( new EventHists("Events"));

  //histograms for efficiency
  RegisterHistCollection( new TopEffiHists("TopEffi"));
  */

  //histograms for the effi calculation of jets
  //RegisterHistCollection( new JetEffiHists("JetEffiHists",0.7,0.7,0,0));

  //histograms with muons and jets
  RegisterHistCollection( new MJetsHists("MJetsHists", m_chi2discr,10 ));
  RegisterHistCollection( new MJetsHists("MJetsHistsCut", m_chi2discr,10 ));
  //RegisterHistCollection( new MJetsHists("MJetsHistsBTag", m_chi2discr,10 ));
  //RegisterHistCollection( new MJetsHists("MJetsHistsNoBTag", m_chi2discr,10 ));

  //histograms with muons and jets and a TopTag
  RegisterHistCollection( new MJetsHists("TagMJetsHists", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_neutrino", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHistsBTag_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHistsNoBTag_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR05", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR03", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR02", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR015", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR01", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR005", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagMJetsHists_delR001", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("TagHists_MuonSel_old", m_chi2discr,10));

  RegisterHistCollection( new MJetsHists("NoTagMJetsHists", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_neutrino", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHistsBTag_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHistsNoBTag_neutrino_muon", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR05", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR03", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR02", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR015", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR01", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR005", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagMJetsHists_delR001", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("NoTagHists_MuonSel_old", m_chi2discr,10));



  RegisterHistCollection( new MJetsHists("MJetsHists_delR05", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR03", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR02", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR015", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR01", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR005", m_chi2discr,10));
  RegisterHistCollection( new MJetsHists("MJetsHists_delR001", m_chi2discr,10));

  RegisterHistCollection( new MJetsHists("TagMJetsHistsCut", m_chi2discr,10));
  
  RegisterHistCollection( new MJetsHists("TagHists_Scaled", m_chi2discr,40));
  /*
  RegisterHistCollection( new MJetsHists("TagHists_ScaledBTag", m_chi2discr,40));
  RegisterHistCollection( new MJetsHists("TagHists_ScaledNoBTag", m_chi2discr,40));
  */
  RegisterHistCollection( new MJetsHists("TagHists_ScaledMuonSel", m_chi2discr,40));
  RegisterHistCollection( new MJetsHists("TagHists_ScaledMuonSel_old", m_chi2discr,40));
  RegisterHistCollection( new MJetsHists("TagHists_ScaledNoMuonSel", m_chi2discr,40));


  RegisterHistCollection( new DelRHists("DelRHists", m_chi2discr));

  /*
  //static HypothesisDiscriminator* discr = new HypothesisDiscriminator();
  RegisterHistCollection( new HypothesisHists("TTrecon", m_chi2discr));


  // histograms about the hypothesis
  RegisterHistCollection( new HypothesisHists("Chi2", m_chi2discr) );


  //with toptag
  RegisterHistCollection( new HypothesisHists("Chi2TopTag", m_chi2discr) );
  */
  //control for neutrinos
  RegisterHistCollection( new NeutrinoHists("NeutrinoControlHists", m_chi2discr) );

  // important: initialise histogram collections after their definition
  InitHistos();

  /*
  m_bp_chi2 = new HypothesisStatistics("b.p. vs. Chi2");
  m_bp_sumdr = new HypothesisStatistics("b.p. vs. SumDR");
  m_cm_chi2 = new HypothesisStatistics("matched vs. Chi2");
  m_cm_sumdr = new HypothesisStatistics("matched vs. SumDR");
  m_cm_bp = new HypothesisStatistics("matched vs. b.p.");
  */  



  return;

}

void ZPrimeTopTaggerCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );

  /*
  m_bp_chi2->PrintStatistics();
  m_bp_sumdr->PrintStatistics();
  m_cm_chi2->PrintStatistics();
  m_cm_sumdr->PrintStatistics();
  m_cm_bp->PrintStatistics();
  */


  return;

}

void ZPrimeTopTaggerCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  return;

}

void ZPrimeTopTaggerCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );

  // get the selections

  //static Selection* TwoDMuonSel = GetSelection("TwoDSelection");
  static Selection* BSel = GetSelection("BSelection");
  static Selection* MuonSel = GetSelection("MuonSelection");
  static Selection* NoBSel = GetSelection("NoBSelection");
  static Selection* TopSel = GetSelection("TopSelection");
  //static Selection* TTreco = GetSelection("TTreco");
  static Selection* TopTagSel = GetSelection("TopTagSelection");
  static Selection* TwoDMuon = GetSelection("TwoDSelection");
  //static Selection* HCALlaser = GetSelection("HCAL_laser_events");

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
 
  TopFitCalc* topfit = TopFitCalc::Instance();

  //if(!TwoDMuonSel->passSelection())  throw SError( SError::SkipEvent );
 
  /*
  if (calc->IsRealData()){
    if (!HCALlaser->passSelection()) throw SError( SError::SkipEvent );
  }
  */

 
  if(!MuonSel->passSelection())  throw SError( SError::SkipEvent );
  static Selection* mttbar_gen_selection = GetSelection("Mttbar_Gen_Selection");

  //if(calc->GetMuons()->at(0).relIso() > 0.4) throw SError( SError::SkipEvent );	


  if(!mttbar_gen_selection->passSelection())  throw SError( SError::SkipEvent );

  //BaseHists* JetEffiHists = GetHistCollection("JetEffiHists");	
  //JetEffiHists->Fill();



  //do reconstruction here
  bcc->recoHyps->clear();    

  if(!TopSel->passSelection())  throw SError( SError::SkipEvent );



  //if(!TopTagSel->passSelection()) throw SError( SError::SkipEvent );
  
  //topfit->FillTTbarTagHypotheses();


  //MonteCarlo Eff

  //if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_tauhad)
  //if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad){}
    //  topfit->MCSelection();
  //else  
  //throw SError( SError::SkipEvent );

  //----------------------
  //these are the more recent lines befor the monte carlo test
  //

  bool notoptag = 0;




  if(TopTagSel->passSelection()){
    topfit->CalculateTopTag();
  }
  else{
    topfit->CalculateSelection();
    //calc->FillHighMassTTbarHypotheses();
    notoptag = 1;
   

  }






  //-------------------------
 

  // topfit->CalculateEvent();
  //calc->FillHighMassTTbarHypotheses();
  
  if(!bcc->recoHyps || bcc->recoHyps->size()==0) throw SError( SError::SkipEvent );
  //m_tagchi2discr->FillDiscriminatorValues();
  //if(!topfit->Match(m_tagchi2discr->GetBestHypothesis(),0.7,0.7,0,0))throw SError( SError::SkipEvent );
  //bcc->recoHyps->clear();    
  
  


  //topfit->MCSelection();
  m_chi2discr->FillDiscriminatorValues();

  ReconstructionHypothesis *hyp = m_chi2discr->GetBestHypothesis();
   
  // get the histogram collections
  /*
  BaseHists* HistsNoCuts = GetHistCollection("NoCuts");
  
  BaseHists* HistsBTag = GetHistCollection("BTag");
  BaseHists* HistsNoBTag = GetHistCollection("NoBTag");
  
  BaseHists* HistsTopSel = GetHistCollection("TopSel");
  
  BaseHists* Chi2_HistsNoCuts = GetHistCollection("Chi2_NoCuts");
  BaseHists* Chi2_HistsBTag = GetHistCollection("Chi2_BTag");
  BaseHists* Chi2_HistsNoBTag = GetHistCollection("Chi2_NoBTag");
  BaseHists* Chi2_HistsTopSel = GetHistCollection("Chi2_TopSel");

 
  BaseHists* Chi2Hists = GetHistCollection("Chi2");

  
  BaseHists* Chi2TopTagHists = GetHistCollection("Chi2TopTag");


  BaseHists* JetsHists = GetHistCollection("Jets");
  */
  BaseHists* MuonsHists = GetHistCollection("Muons");
  BaseHists* TagMuonsHists = GetHistCollection("TagMuons");
  /*
  BaseHists* TopJetsHists = GetHistCollection("TopJets");
  BaseHists* EventHists = GetHistCollection("Events");

  BaseHists* TopEffiHists = GetHistCollection("TopEffi");
  */
  BaseHists* MJets = GetHistCollection("MJetsHists");
  BaseHists* MJetsCut = GetHistCollection("MJetsHistsCut");
  /*
  BaseHists* MJetsBTag = GetHistCollection("MJetsHistsBTag");
  BaseHists* MJetsNoBTag = GetHistCollection("MJetsHistsNoBTag");
  */
  BaseHists* TagMJets = GetHistCollection("TagMJetsHists");
  BaseHists* TagMJets_neutrino = GetHistCollection("TagMJetsHists_neutrino");
  BaseHists* TagMJets_neutrino_muon = GetHistCollection("TagMJetsHists_neutrino_muon");
  BaseHists* TagMJetsBTag_neutrino_muon = GetHistCollection("TagMJetsHistsBTag_neutrino_muon");
  BaseHists* TagMJetsNoBTag_neutrino_muon = GetHistCollection("TagMJetsHistsNoBTag_neutrino_muon");
  BaseHists* TagMJets_delR05 =  GetHistCollection("TagMJetsHists_delR05");
  BaseHists* TagMJets_delR03 =  GetHistCollection("TagMJetsHists_delR03");
  BaseHists* TagMJets_delR02 =  GetHistCollection("TagMJetsHists_delR02");
  BaseHists* TagMJets_delR015 = GetHistCollection("TagMJetsHists_delR015");
  BaseHists* TagMJets_delR01 =  GetHistCollection("TagMJetsHists_delR01");
  BaseHists* TagMJets_delR005 = GetHistCollection("TagMJetsHists_delR005");
  BaseHists* TagMJets_delR001 = GetHistCollection("TagMJetsHists_delR001");
  BaseHists* TagMuonSel_old = GetHistCollection("TagHists_MuonSel_old");

  BaseHists* NoTagMJets = GetHistCollection("NoTagMJetsHists");
  BaseHists* NoTagMJets_neutrino = GetHistCollection("NoTagMJetsHists_neutrino");
  BaseHists* NoTagMJets_neutrino_muon = GetHistCollection("NoTagMJetsHists_neutrino_muon");
  BaseHists* NoTagMJetsBTag_neutrino_muon = GetHistCollection("NoTagMJetsHistsBTag_neutrino_muon");
  BaseHists* NoTagMJetsNoBTag_neutrino_muon = GetHistCollection("NoTagMJetsHistsNoBTag_neutrino_muon");
  BaseHists* NoTagMJets_delR05 =  GetHistCollection("NoTagMJetsHists_delR05");
  BaseHists* NoTagMJets_delR03 =  GetHistCollection("NoTagMJetsHists_delR03");
  BaseHists* NoTagMJets_delR02 =  GetHistCollection("NoTagMJetsHists_delR02");
  BaseHists* NoTagMJets_delR015 = GetHistCollection("NoTagMJetsHists_delR015");
  BaseHists* NoTagMJets_delR01 =  GetHistCollection("NoTagMJetsHists_delR01");
  BaseHists* NoTagMJets_delR005 = GetHistCollection("NoTagMJetsHists_delR005");
  BaseHists* NoTagMJets_delR001 = GetHistCollection("NoTagMJetsHists_delR001");
  BaseHists* NoTagMuonSel_old = GetHistCollection("NoTagHists_MuonSel_old");


  BaseHists* MJets_delR05 =  GetHistCollection("MJetsHists_delR05");
  BaseHists* MJets_delR03 =  GetHistCollection("MJetsHists_delR03");
  BaseHists* MJets_delR02 =  GetHistCollection("MJetsHists_delR02");
  BaseHists* MJets_delR015 = GetHistCollection("MJetsHists_delR015");
  BaseHists* MJets_delR01 =  GetHistCollection("MJetsHists_delR01");
  BaseHists* MJets_delR005 = GetHistCollection("MJetsHists_delR005");
  BaseHists* MJets_delR001 = GetHistCollection("MJetsHists_delR001");

  BaseHists* TagMJetsCut = GetHistCollection("TagMJetsHistsCut");
  //BaseHists* TagMJetsBTag = GetHistCollection("TagMJetsHistsBTag");
  // BaseHists* TagMJetsNoBTag = GetHistCollection("TagMJetsHistsNoBTag");

  BaseHists* TagScaled = GetHistCollection("TagHists_Scaled");

  /*
  BaseHists* TagScaledBTag = GetHistCollection("TagHists_ScaledBTag");
  BaseHists* TagScaledNoBTag = GetHistCollection("TagHists_ScaledNoBTag");
  */
  BaseHists* TagScaledMuonSel = GetHistCollection("TagHists_ScaledMuonSel");
  BaseHists* TagScaledMuonSel_old = GetHistCollection("TagHists_ScaledMuonSel_old");
  BaseHists* TagScaledNoMuonSel = GetHistCollection("TagHists_ScaledNoMuonSel");

  /*
  BaseHists* HypothesisHists = GetHistCollection("TTrecon");
  */
  BaseHists* NeutrinoHists = GetHistCollection("NeutrinoControlHists");

  BaseHists* DeltaRHists = GetHistCollection("DelRHists");


  Particle neutrino;
  neutrino.set_v4(hyp->neutrino_v4());
  Muon muon = bcc->muons->at(0);

 
  // start the analysis
  
  /*HistsNoCuts->Fill();
  Chi2_HistsNoCuts->Fill();
  JetsHists->Fill();
  */
  MuonsHists->Fill();
  /*
  TopJetsHists->Fill();
  EventHists->Fill();

  HypothesisHists->Fill();

  TopEffiHists->Fill();
  */

  DeltaRHists->Fill();
  MJets->Fill();
 
  
  if(deltaRmin(&(bcc->muons->at(0)),  bcc->jets)>0.5)MJets_delR05->Fill();
  if(deltaRmin(&(bcc->muons->at(0)),  bcc->jets)>0.3)MJets_delR03->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)),  bcc->jets)>0.2)MJets_delR02->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.15)MJets_delR015->Fill();
  if(deltaRmin(&(bcc->muons->at(0)),  bcc->jets)>0.1)MJets_delR01->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.05)MJets_delR005->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.01)MJets_delR001->Fill();



  NeutrinoHists->Fill();
  //HistsNoBTag->Fill();
 

  if(bcc->muons->at(0).relIso()<0.4 && deltaRmin(&(bcc->muons->at(0)), bcc->jets)>.1 ){
   
  }
  else{
   MJetsCut->Fill();
  }
 
  if(notoptag ==1)NoTagMJets->Fill();
  if(TwoDMuon->passSelection() && notoptag ==1)NoTagMuonSel_old->Fill();

  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.5 && notoptag ==1)NoTagMJets_delR05->Fill();

  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.3 && notoptag ==1)NoTagMJets_delR03->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.2 && notoptag ==1)NoTagMJets_delR02->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.15 && notoptag ==1)NoTagMJets_delR015->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.1 && notoptag ==1)NoTagMJets_delR01->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.05 && notoptag ==1)NoTagMJets_delR005->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.01 && notoptag ==1)NoTagMJets_delR001->Fill();
  if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*6.25 + pow(pTrel(&neutrino, bcc->jets)-100,2)*0.0000694) > 1 && notoptag ==1 && (hyp->toplep_v4()+hyp->tophad_v4()).M()>500 )NoTagMJets_neutrino->Fill();
  if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*1.5625 + pow(pTrel(&muon, bcc->jets),2)*.00308641975) > 1 && notoptag ==1 && (hyp->toplep_v4()+hyp->tophad_v4()).M()>500){
    NoTagMJets_neutrino_muon->Fill();
    if(BSel->passSelection())
      NoTagMJetsBTag_neutrino_muon->Fill();
    else
      NoTagMJetsNoBTag_neutrino_muon->Fill();
  }
 
  /*
  if(BSel->passSelection()){
    HistsBTag->Fill();
    MJetsBTag->Fill();
    //Chi2_HistsBTag->Fill();
  }
  if(NoBSel->passSelection()){
    HistsNoBTag->Fill(); 
    MJetsNoBTag->Fill(); 
    // Chi2_HistsNoBTag->Fill();
  }
  */
  //if(!TopTagSel->passSelection())  throw SError( SError::SkipEvent );
  
  //  if(TopTagSel->passSelection()) 


 
  if(!TopTagSel->passSelection()) throw SError( SError::SkipEvent );
  //bcc->recoHyps->clear();
  //topfit->CalculateTopTag();
  //m_TagChi2discr->FillDiscriminatorValues();

  TagMJets->Fill();


  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.5)TagMJets_delR05->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.3)TagMJets_delR03->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.2)TagMJets_delR02->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.15)TagMJets_delR015->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.1)TagMJets_delR01->Fill(); 
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.05)TagMJets_delR005->Fill();
  if(deltaRmin(&(bcc->muons->at(0)), bcc->jets)>0.01)TagMJets_delR001->Fill();


  if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*1.5625 + pow(pTrel(&muon, bcc->jets),2)*.00308641975) > 1 && (hyp->toplep_v4()+hyp->tophad_v4()).M()>500){
    TagMJets_neutrino_muon->Fill();
    if(BSel->passSelection())
      TagMJetsBTag_neutrino_muon->Fill();
    else
      TagMJetsNoBTag_neutrino_muon->Fill();
  }

  TagScaled->Fill();
  
  /*
  if(BSel->passSelection()){
    TagScaledBTag->Fill();
  }
  if(NoBSel->passSelection()){
    TagScaledNoBTag->Fill();
  }
  
  if(TwoDMuonSel->passSelection()){
    TagScaledMuonSel->Fill();
  }
  else{
    TagScaledNoMuonSel->Fill();
  }
  */
  
  TagMuonsHists->Fill();

  //if(deltaRmin(&(bcc->muons->at(0)), bcc->jets) > 150/bcc->muons->at(0).pt() || pTrel(&(bcc->muons->at(0)), bcc->jets)> 25){///bcc->muons->at(0).pt()

  //  if(bcc->muons->at(0).relIso()<0.4 && deltaRmin(&(bcc->muons->at(0)), bcc->jets)>.1 ){

  // if(deltaRmin(&neutrino, bcc->jets)< 0.42+pow(pTrel(&neutrino, bcc->jets)-65,2)){
  //if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*6.92520775623 + pow(pTrel(&neutrino, bcc->jets)-80,2)*0.0000510204) > 1){
  if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*6.25 + pow(pTrel(&neutrino, bcc->jets)-100,2)*0.0000694) > 1 && (hyp->toplep_v4()+hyp->tophad_v4()).M()>500){
    TagScaledMuonSel->Fill();
    TagMJets_neutrino->Fill();
  }
  else{
    TagScaledNoMuonSel->Fill();
    TagMJetsCut->Fill();
  }

  if( (pow(deltaRmin(&neutrino, bcc->jets)-.8,2)*1.5625 + pow(pTrel(&muon, bcc->jets),2)*.00308641975) > 1 && notoptag ==1)NoTagMJets_neutrino_muon->Fill();
  if(TwoDMuon->passSelection()){
    TagScaledMuonSel_old->Fill();
    TagMuonSel_old->Fill();
    }

  
  /*
  Chi2TopTagHists->Fill();
  TagMJets->Fill();
  HistsTopSel->Fill();
  Chi2_HistsTopSel->Fill();
  */
 
  
  
  return;
  
}

