
// $Id: TreeCycle.cxx,v 1.10 2012/12/07 14:21:51 peiffer Exp $


#include <iostream>

using namespace std;

// Local include(s):
#include "include/TreeCycle.h"

ClassImp( TreeCycle );

TreeCycle::TreeCycle()
   : AnalysisCycle() {

  // constructor, declare additional variables that should be 
  // obtained from the steering-xml file
  
  // set the integrated luminosity per bin for the lumi-yield control plots
  SetIntLumiPerBin(500.);

  m_mttgencut = false;
  DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
}

TreeCycle::~TreeCycle() 
{
  // destructor
}

void TreeCycle::BeginCycle() throw( SError ) 
{
  // Start of the job, general set-up and definition of 
  // objects are done here

  // Important: first call BeginCycle of base class
  AnalysisCycle::BeginCycle();

  return;

}

void TreeCycle::EndCycle() throw( SError ) 
{
  // clean-up, info messages and final calculations after the analysis

  
  // call the base cycle class for all standard methods
  AnalysisCycle::EndCycle();

  return;

}

void TreeCycle::BeginInputData( const SInputData& id ) throw( SError ) 
{
  // declaration of histograms and selections

  // Important: first call BeginInputData of base class
  AnalysisCycle::BeginInputData( id );

  Selection* TopTagSel = new Selection("TopTagSelection");
  TopTagSel->addSelectionModule(new NTopTagSelection(1));
  TopTagSel->addSelectionModule(new TopTagOverlapSelection());


  static Selection* mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
 
  if ((id.GetVersion() == "TTbar") ){
    m_logger << INFO << "Applying mttbar generator cut from 0 to 700 GeV." << SLogger::endmsg;
    mttbar_gen_selection->addSelectionModule(new MttbarGenCut(0,700));
    mttbar_gen_selection->EnableSelection();
  }
  else {
    m_logger << INFO << "Disabling mttbar generator cut." << SLogger::endmsg;
    mttbar_gen_selection->DisableSelection();
  }
  
  m_chi2discr = new Chi2Discriminator();


  RegisterSelection(mttbar_gen_selection);
  RegisterSelection(TopTagSel);

  // -------------------- set up the selections ---------------------------

  RegisterHistCollection( new MYTMVATree("Tree", m_chi2discr ) );

  // important: initialise histogram collections after their definition
  InitHistos();

  return;

}

void TreeCycle::EndInputData( const SInputData& id ) throw( SError ) 
{
  AnalysisCycle::EndInputData( id );
  return;

}

void TreeCycle::BeginInputFile( const SInputData& id ) throw( SError ) 
{
  // Connect all variables from the Ntuple file with the ones needed for the analysis
  // The variables are commonly stored in the BaseCycleContaincer

  // important: call to base function to connect all variables to Ntuples from the input tree
  AnalysisCycle::BeginInputFile( id );

  return;

}

void TreeCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError ) 
{
  // this is the most important part: here the full analysis happens
  // user should implement selections, filling of histograms and results

  // first step: call Execute event of base class to perform basic consistency checks
  // also, the good-run selection is performed there and the calculator is reset
  AnalysisCycle::ExecuteEvent( id, weight );


  static Selection* mttbar_gen_selection = GetSelection("Mttbar_Gen_Selection");
  if(!mttbar_gen_selection->passSelection()) ClearEvent();

  EventCalc* calc = EventCalc::Instance();
  TopFitCalc* topfit = TopFitCalc::Instance();

  calc->FillHighMassTTbarHypotheses();

  /*
  if(TopTagSel->passSelection()){
    topfit->CalculateTopTag();
  }
  else{
    calc->FillHighMassTTbarHypotheses();
  }
  */
  m_chi2discr->FillDiscriminatorValues();
  

  BaseHists* Tree = GetHistCollection("Tree");

  Tree->Fill();

  ClearEvent();
  return;
  
}


