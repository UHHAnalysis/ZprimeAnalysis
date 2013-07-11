#include "include/TopEffiHists.h"
#include "include/SelectionModules.h"
#include "include/EventCalc.h"
#include "TH2.h"

#include <iostream>

using namespace std;

TopEffiHists::TopEffiHists(const char* name) : BaseHists(name)
{
  // named default constructor
  
}

TopEffiHists::~TopEffiHists()
{
  // default destructor, does nothing
}

void TopEffiHists::Init()
{
  // book all histograms here
 
  Book( TH1F( "Mjet", "m_{jet}", 100,0,500 ) );
  Book( TH1F( "Mmin", "m_{min}", 100,0,200 ) );
  Book( TH2F( "Mjet_Mmin", "m_{jet} : m_{min}", 100,0,500,100,0,200 ) );
  Book( TH1F( "Nsubjet", "N^{subjet}", 10,0,10 ) );
  Book( TH1F( "pT","pT",100,0,2000));
  Book( TH1F( "top_pT","Top p_{T}",50,0,2000));
  Book( TH1F( "pT_tag","p_{T} Tag",50,0,2000));
  Book( TH1F( "Eta","#eta",50,-5,5));
  Book( TH1F( "top_Eta","#eta Top",50,-4,4));
  Book( TH1F( "Eta_tag","#eta Tag",50,-4,4));
  Book( TH1F( "Phi","#phi",50,-PI,PI));
  Book( TH1F( "top_Phi","#phi Top",50,-PI,PI));
  Book( TH1F( "Phi_tag","#phi Tag",50,-PI,PI));
  Book( TH2F( "ToppT_JetpT","Top p_{T}, Jet p_{T}",50,0,2000,100,0,2000) );
  Book( TH1F( "Mjet_W", " W m_{Jet}", 50,0,200 ) );
  Book( TH1F( "mTopjet", "m_{Top}", 50,0,500 ) );
  Book( TH1F( "r", "#Delta r", 50,0,5 ) );
  Book( TH1F( "r_topjet", " #Delta r < 0.8", 50,0,5 ) );
  Book( TH1F( "Efficiency_pT", "Top p_{T}", 50,0,2000 ) );
  Book( TH1F( "Efficiency_Eta", "#eta tag/top", 50,-4,4 ) );
  Book( TH1F( "Efficiency_Phi", "#phi tag/top", 50,-PI,PI ) );
}


void TopEffiHists::Fill()
{
   // important: get the event weight
  EventCalc* calc = EventCalc::Instance();
  double weight = calc -> GetWeight();
  
  if (calc->IsRealData()) return;
   
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();



  
  for (unsigned int i =0; i<bcc->topjets->size(); ++i)
    {
      TopJet topjet =  bcc->topjets->at(i);

      double mmin=0;
      double mjet=0;
      int nsubjets=0;

      bool NTopTag = TopTag(topjet,mjet,nsubjets,mmin);

      Hist( "Mjet" )->Fill( mjet, weight );
      if(nsubjets>=3)
       {
	 Hist( "Mmin" )->Fill( mmin, weight );
	 TH2F* h2 = (TH2F*)Hist("Mjet_Mmin");
	 h2->Fill(mjet,mmin,weight);
       }

      Hist( "Nsubjet" )->Fill( nsubjets, weight );
     
      Hist("pT")->Fill(topjet.pt(), weight);    
      Hist("Eta")->Fill(topjet.eta(), weight);
      Hist("Phi")->Fill(topjet.phi(), weight);
      

   
      for(unsigned int j=0; j<bcc->genparticles->size(); ++j)
	{
	  GenParticle genp = bcc->genparticles->at(j);

	  if (genp.pdgId() == 24 || genp.pdgId() == -24 )
	    Hist("Mjet_W")->Fill(mjet,weight);

	  Hist("r")->Fill(topjet.deltaR(genp),weight);

	  if (genp.pdgId() == 6 || genp.pdgId() == -6 )
	    {
	      double r = topjet.deltaR(genp);
	      if (r < 0.8)
		{
		 
		  Hist("top_pT")->Fill(genp.pt(),weight);
		  Hist("top_Eta")->Fill(genp.eta(), weight);
		  Hist("top_Phi")->Fill(genp.phi(), weight);
		  Hist("r_topjet")->Fill(topjet.deltaR(genp),weight);
		  Hist("mTopjet")->Fill(mjet,weight);


		  if (NTopTag)
		    {
		      Hist("pT_tag")->Fill(genp.pt(), weight);
		      Hist("Eta_tag")->Fill(genp.eta(), weight);
		      Hist("Phi_tag")->Fill(genp.phi(), weight);
		      TH2F* h3 = (TH2F*)Hist("ToppT_JetpT");
		      h3->Fill(genp.pt(), topjet.pt(),weight);
		    }
		}
	    }
	} 
    }
}

void TopEffiHists::Finish()
{
  // final calculations, like division and addition of certain histograms
  Hist("Efficiency_pT") ->Divide( Hist("pT_tag"), Hist("top_pT"));
  Hist("Efficiency_Eta")->Divide( Hist("Eta_tag"), Hist("top_Eta"));
  Hist("Efficiency_Phi")->Divide( Hist("Phi_tag"), Hist("top_Phi"));
}

