// Dear emacs, this is -*- c++ -*-
// $Id: QCDCycle.h,v 1.0 2012/06/05 14:45:24 gonvaq Exp $
#ifndef QCDCycle_H
#define QCDCycle_H

// SFrame include(s):
#include "include/AnalysisCycle.h"
#include "HypothesisDiscriminator.h"
#include "HypothesisHists.h"
#include "include/SelectionModules.h"
#include "HypothesisStatistics.h"
#include "Cleaner.h"
#include "JetCorrectionUncertainty.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"


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
 *  @author Daniel Gonzalez
 *  @version $Revision: 1.1 $
 */

class QCDCycle : public AnalysisCycle {

public:
  /// Default constructor
  QCDCycle();
  /// Default destructor
  ~QCDCycle();

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

  /// Fill control histograms
  void FillControlHistos(TString postfix="");


private:
  //
  // Put all your private variables here
  //
  
  int m_Nbtags_max;
  int m_Nbtags_min;  


  TTree* m_tree;

  Chi2Discriminator* m_chi2discr;
  Chi2Discriminator* m_tagchi2discr;
  BestPossibleDiscriminator* m_bpdiscr;
  SumDeltaRDiscriminator* m_sumdrdiscr;
  CorrectMatchDiscriminator* m_cmdiscr;

  E_BtagType m_btagtype;

  std::string m_dobsf;
  BTaggingScaleFactors* m_bsf;
  std::string m_filter_file;

  std::string m_Electron_Or_Muon_Selection;

  bool m_mttgencut;

  JetCorrectionUncertainty* m_jes_unc;
  Cleaner* m_cleaner;

  // Macro adding the functions for dictionary generation
  ClassDef( QCDCycle, 0 );

}; // class QCDCycle

#endif // QCDCycle_H

