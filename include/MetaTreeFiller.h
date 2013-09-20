#ifndef MetaTreeFiller_H
#define MetaTreeFiller_H

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

class MetaTreeFiller : public BaseHists {

public:
   /// Named constructor
  MetaTreeFiller(const char* name, HypothesisDiscriminator *discr, TTree* MetaOutputTree);

   /// Default destructor
   ~MetaTreeFiller();

   void Init();

   void Fill();

   void Finish();

private:
   TTree* m_tree;
   
   Double_t m_weight;
   Double_t m_pT_mu;
   Double_t m_HT; 
   Double_t m_2D; 
   Double_t m_mtt;
   Double_t m_met_pt; 
   Double_t m_Jet_pt_max;

   Double_t m_Iso05; 
   Double_t m_Iso04; 

   Double_t m_Iso02;
   Double_t m_Iso018;
   Double_t m_Iso016;
   Double_t m_Iso014;
   Double_t m_Iso012;
   Double_t m_Iso01; 
   Double_t m_Iso008;
   Double_t m_Iso006;
   Double_t m_Iso004;
   Double_t m_Iso002;


   HypothesisDiscriminator* m_discr;

}; // class MetaTreeFiller


#endif // MetaTreeFiller_H
