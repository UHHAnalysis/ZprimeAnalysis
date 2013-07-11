#ifndef TopEffiHists_H
#define TopEffiHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   Class for booking and filling Top Tagging efficiency histograms
 *
 *   
 *   @version $Revision: 1.1 $
 */

class TopEffiHists : public BaseHists {

public:
   /// Named constructor
   TopEffiHists(const char* name);

   /// Default destructor
   ~TopEffiHists();

   void Init();

   void Fill();

   void Finish();

private:

}; // class TopEffiHists


#endif // TopEffiHists_H


