#ifndef NeutrinoHists_H
#define NeutrinoHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"
#include "HypothesisDiscriminator.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"



class NeutrinoHists : public BaseHists {

public:
   /// Named constructor
   NeutrinoHists(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~NeutrinoHists();

   void Init();

   void Fill();

   void Finish();

private:

 HypothesisDiscriminator* m_discr;

}; // class NeutrinoHists


#endif // NeutrinoHists_H
