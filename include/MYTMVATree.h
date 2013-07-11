#ifndef MYTMVATree_H
#define MYTMVATree_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>

// Local include(s):
#include "include/BaseHists.h"
#include "include/BaseCycleContainer.h"
#include "include/HypothesisDiscriminator.h"
#include "include/ReconstructionHypothesis.h"   

/**
 *   Example class for booking and filling histograms
 *
 *   This class books and fills a collection of histograms.
 *   It should have a unique name, such that the histograms
 *   of multiple instances of this class are ordered in the
 *   output file. 
 *   Always sort your histograms and used methods topically.
 *   This example collection can be used for data and reconstructed
 *   MC events.
 *   
 *   @version $Revision: 1.1 $
 */

class MYTMVATree : public BaseHists {

public:
   /// Named constructor
   MYTMVATree(const char* name, HypothesisDiscriminator *discr);

   /// Default destructor
   ~MYTMVATree();

   void Init();

   void Fill();

   void Finish();

   void FillTopJetProperties(TopJet topjet, GenParticle topquark);
   void FillJetProperties(Jet jet);
   
   void ClearVariables();

   void ClearCutVariables();
   void ClearJetVariables();


   Double_t RecGenMatch(TopJet tj, GenParticle gen);

private:
   
   HypothesisDiscriminator* m_discr;


   TTree* m_tree;
   TTree* m_treeJet;
   TTree* m_treeCut;
   
   // variables for the TMVA training   
   Double_t m_pt;       // px of the CA-fat jet	 
   Double_t m_phi;       // py of the CA-fat jet	 
   Double_t m_eta;       // pz of the CA-fat jet	 
   Double_t m_mass;     // mass of the CA-fat jet

   Double_t m_genpx;    // px of the top quark (gen)
   Double_t m_genpy;    // py of the top quark (gen)
   Double_t m_genpz;    // pz of the top quark (gen)

   Int_t m_npv;         // number of reconstructed primary vertices
   Double_t m_weight;   // event weight

   Int_t m_nsubs;       // number of sub-jets
   Bool_t m_HEPTopTag;  // tag from HEP Top Tagger
   Bool_t m_CMSTopTag;  // tag from CMS Top Tagger

   Double_t m_sub1_px;       // px of subjet 1
   Double_t m_sub1_py;       // py of subjet 1
   Double_t m_sub1_pz;       // pz of subjet 1
   Double_t m_sub1_mass;     // mass of subjet 1

   Double_t m_sub2_px;       // px of subjet 2
   Double_t m_sub2_py;       // py of subjet 2
   Double_t m_sub2_pz;       // pz of subjet 2
   Double_t m_sub2_mass;     // mass of subjet 2

   Double_t m_sub3_px;       // px of subjet 3
   Double_t m_sub3_py;       // py of subjet 3
   Double_t m_sub3_pz;       // pz of subjet 3
   Double_t m_sub3_mass;     // mass of subjet 3

   Double_t m_sub4_px;       // px of subjet 4
   Double_t m_sub4_py;       // py of subjet 4
   Double_t m_sub4_pz;       // pz of subjet 4
   Double_t m_sub4_mass;     // mass of subjet 4

   Double_t m_m12;      // mass of subjets 1+2
   Double_t m_m13;      // mass of subjets 1+3
   Double_t m_m23;      // mass of subjets 2+3

   Double_t m_pruned_mass;     // mass of the pruned CA-fat jet
   Double_t m_pruned_m12;      // mass of pruned subjets 1+2
   Double_t m_pruned_m13;      // mass of pruned subjets 1+3
   Double_t m_pruned_m23;      // mass of pruned subjets 2+3

   
   Double_t m_Jet_pt ;
   Double_t m_Jet_eta ;
   Double_t m_Jet_phi ;
   Double_t m_Jet_E;

   Double_t m_pTrel_muon_jet;
   Double_t m_pTrel_muon_neutrino;
   Double_t m_pTrel_muon_bjet;
   Double_t m_pTrel_neutrino_bjet;
   Double_t m_pTrel_neutrino_jet;

   Double_t m_deltaR_muon_jet;
   Double_t m_deltaR_neutrino_jet;
   Double_t m_deltaR_muon_bjet;
   Double_t m_deltaR_neutrino_bjet;

   Double_t m_razor_mrt_w_tophad;
   Double_t m_razor_mr_w_tophad;
   Double_t m_razor_mrt_toplep_tophad;
   Double_t m_razor_mr_toplep_tophad;

   Double_t m_razor_mrt_w_tophad_addmiss;    
   Double_t m_razor_mr_w_tophad_addmiss;    

   Double_t m_mass2;
   Double_t m_muon_pt;
   Double_t m_muon_phi;
   Double_t m_muon_eta;
   Double_t m_beta;

   Double_t m_mtt;


}; // class MYTMVATree


#endif // MYTMVATree_H
