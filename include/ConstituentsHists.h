#ifndef ConstituentsHists_H
#define ConstituentsHists_H

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

class ConstituentsHists : public BaseHists {

public:
   /// Named constructor
  ConstituentsHists(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~ConstituentsHists();

   void Init();

   void Fill();

   void Finish();

private:
 
  

   HypothesisDiscriminator* m_discr;

}; // class ConstituentsHists


#endif // ConstituentsHists_H
