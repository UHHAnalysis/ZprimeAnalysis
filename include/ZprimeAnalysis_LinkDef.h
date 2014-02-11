// Dear emacs, this is -*- c++ -*-
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;

// Add the declarations of your cycles, and any other classes for which you
// want to generate a dictionary, here. The usual format is:
//
//#pragma link C++ class TMVATreeFiller;

#pragma link C++ class JetLeptonCleanerCycle+;
#pragma link C++ class ZprimePreSelectionCycle+;
#pragma link C++ class ZprimeSelectionCycle+;
#pragma link C++ class ZprimePostSelectionCycle+;
#pragma link C++ class ZprimeEleTrigCycle+;
#pragma link C++ class ZprimeSelectionTrigTestCycle+;
#pragma link C++ class ZprimeSelectionDataTrigTestCycle+;
#pragma link C++ class ZprimeEleTrigCycle+;
#pragma link C++ class DileptonPreSelectionCycle+;


#pragma link C++ class TopFitCalc+;

#pragma link C++ class TopTagTMVACycle+;
#pragma link C++ class GenTTbarCycle+;

#pragma link C++ class DelRHists+;
#pragma link C++ class ConstituentsHists+;
#pragma link C++ class MetaTreeFiller+;
#pragma link C++ class CleanerFiller+;
#pragma link C++ class QCDCycle+;
#pragma link C++ class RazorHists+;
#pragma link C++ class TopTagDelRHists+;
#pragma link C++ class ZPrimeTopTaggerCycle+;
#pragma link C++ class MJetsHists+;
#pragma link C++ class NeutrinoHists+;
#pragma link C++ class TopEffiHists+;
#pragma link C++ class ExampleHists+;

#endif // __CINT__
