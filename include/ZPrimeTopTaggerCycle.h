// Dear emacs, this is -*- c++ -*-
// $Id: ZPrimeTopTaggerCycle.h,v 1.0 2012/06/05 14:45:24 gonvaq Exp $
#ifndef ZPrimeTopTaggerCycle_H
#define ZPrimeTopTaggerCycle_H

// SFrame include(s):
#include "include/AnalysisCycle.h"
#include "HypothesisDiscriminator.h"
#include "HypothesisHists.h"
#include "include/SelectionModules.h"
#include "HypothesisStatistics.h"

#include "EventHists.h"
#include "JetHists.h"
#include "ElectronHists.h"
#include "MuonHists.h"
#include "TauHists.h"
#include "TopJetHists.h"




/**
 *   @short Example of an analysis cycle
 *
 *          This is an example of an analysis cycle. It can be used
 *          as a template for writing your own analysis. Also should
 *          be used for quick cross checks of the system setup.
 *
 *  @author Roman Kogler
 *  @version $Revision: 1.1 $
 */

class ZPrimeTopTaggerCycle : public AnalysisCycle {

public:
  /// Default constructor
  ZPrimeTopTaggerCycle();
  /// Default destructor
  ~ZPrimeTopTaggerCycle();

  /// Function called at the beginning of the cycle
  void BeginCycle() throw( SError );
  /// Function called at the end of the cycle
  void EndCycle() throw( SError );

  /// Function called at the beginning of a new input data
  void BeginInputData( const SInputData& ) throw( SError );
  /// Function called after finishing to process an input data
  void EndInputData  ( const SInputData& ) throw( SError );

  /// Function called after opening each new input file
  void BeginInputFile( const SInputData& ) throw( SError );

  /// Function called for every event
  void ExecuteEvent( const SInputData&, Double_t ) throw( SError );

private:
  //
  // Put all your private variables here
  //
  
  Chi2Discriminator* m_chi2discr;
  Chi2Discriminator* m_tagchi2discr;
  BestPossibleDiscriminator* m_bpdiscr;
  SumDeltaRDiscriminator* m_sumdrdiscr;
  CorrectMatchDiscriminator* m_cmdiscr;


  HypothesisStatistics* m_bp_chi2;
  HypothesisStatistics* m_bp_sumdr;
  HypothesisStatistics* m_cm_chi2;
  HypothesisStatistics* m_cm_sumdr; 
  HypothesisStatistics* m_cm_bp; 


  bool m_mttgencut;



  // Macro adding the functions for dictionary generation
  ClassDef( ZPrimeTopTaggerCycle, 0 );

}; // class ZPrimeTopTaggerCycle

#endif // ZPrimeTopTaggerCycle_H

