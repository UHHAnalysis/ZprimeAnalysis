#ifndef RazorHists_H
#define RazorHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"
#include "include/HypothesisDiscriminator.h"

/**
 *   Razor class for booking and filling histograms
 *

 *   @version $Revision: 1.2 $
 */

class RazorHists : public BaseHists {

public:
   /// Named constructor
   RazorHists(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~RazorHists();

   void Init();

   void Fill();

   void Finish();

private:
   HypothesisDiscriminator *m_discr;


}; // class RazorHists


#endif // RazorHists_H
