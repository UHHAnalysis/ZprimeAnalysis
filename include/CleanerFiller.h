#ifndef CleanerFiller_H
#define CleanerFiller_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"
#include "HypothesisDiscriminator.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   class for booking and filling histograms to discriminate 
 *   QCD from TTbar with the constituents of the jets
 *
 *  
 *   
 *   @version $Revision: 0.1 $
 */

class CleanerFiller : public BaseHists {

public:
   /// Named constructor
  CleanerFiller(const char* name, HypothesisDiscriminator *discr, TTree* MetaOutputTree);

   /// Default destructor
   ~CleanerFiller();

   void Init();

   void Fill();

   void Finish();

private:
   TTree* m_tree;
   
   Double_t m_toptag;
   Double_t m_chi2;	
   Double_t m_mtt_clean;	
   Double_t m_heptoptag;

   HypothesisDiscriminator* m_discr;

}; // class CleanerFiller


#endif // CleanerFiller_H
