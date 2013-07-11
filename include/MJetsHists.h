#ifndef MJetsHists_H
#define MJetsHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"
#include "HypothesisDiscriminator.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   Class for booking and filling muon+jets histograms
 *   for interessting Z' analysis
 *   
 *   @version $Revision: 1.2 $
 */

class MJetsHists : public BaseHists {

public:
   /// Named constructor
  MJetsHists(const char* name, HypothesisDiscriminator *discr, int scalefactor);

   /// Default destructor
   ~MJetsHists();

   void Init();

   void Fill();

   void Finish();

private:

   HypothesisDiscriminator* m_discr;
   int scale;
}; // class MJetsHists


#endif // MJetsHists_H


