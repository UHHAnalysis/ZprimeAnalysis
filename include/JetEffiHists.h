#ifndef JetEffiHists_H
#define JetEffiHists_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include "TH2.h"

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"

/**
 *   JetEffi class for booking and filling histograms
 *
 *   This class books and fills a collection of histograms.
 *   It should have a unique name, such that the histograms
 *   of multiple instances of this class are ordered in the
 *   output file. 
 *   Always sort your histograms and used methods topically.
 *   This example collection can be used for data and reconstructed
 *   MC events.
 *   
 *   @version $Revision: 1.2 $
 */

class JetEffiHists : public BaseHists {

public:
   /// Named constructor
  JetEffiHists(const char* name, double R_lep = 0, double R_had = 0, double pt_lep = 0, double pt_had = 0);

   /// Default destructor
   ~JetEffiHists();

   void Init();

   void Fill();

   void Finish();

private:
   double m_R_lep;
   double m_R_had;
   double m_pt_lep; 
   double m_pt_had;


}; // class JetEffiHists


#endif // JetEffiHists_H
