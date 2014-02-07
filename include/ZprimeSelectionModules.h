#ifndef ZprimeSelectionModules_H
#define ZprimeSelectionModules_H

#include "SFrameAnalysis/include/SelectionModules.h"


class IsoConeSelection : public SelectionModule {
 public:
  IsoConeSelection(TString type="mu", double a=29.1356, double b=164.383, double c=0.023111, double iso_num=0.2);
  ~IsoConeSelection(){};

  virtual bool pass(BaseCycleContainer*);
  virtual std::string description();
  
 private:
 TString m_type;
 double m_a;
 double m_b;
 double m_c;
 double m_iso;

};






#endif
