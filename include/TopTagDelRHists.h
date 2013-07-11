#ifndef TopTagDelRHists_H
#define TopTagDelRHists_H

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

class TopTagDelRHists : public BaseHists {

public:
   /// Named constructor
  TopTagDelRHists(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~TopTagDelRHists();

   void Init();

   void Fill();

   void Finish();

private:
   HypothesisDiscriminator* m_discr;

}; // class TopTagDelRHists


#endif // TopTagDelRHists_H
