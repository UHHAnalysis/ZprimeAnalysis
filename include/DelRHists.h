#ifndef DelRHists_H
#define DelRHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"
#include "HypothesisDiscriminator.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   DelR class for booking and filling histograms
 *
 *   @version $Revision: 1.2 $
 */

class DelRHists : public BaseHists {

public:
   /// Named constructor
  DelRHists(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~DelRHists();

   void Init();

   void Fill();

   void Finish();

private:
   HypothesisDiscriminator* m_discr;

}; // class DelRHists


#endif // DelRHists_H
