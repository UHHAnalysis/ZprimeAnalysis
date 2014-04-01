// $Id: ZprimePostSelectionCycle.cxx,v 1.12 2013/04/30 07:51:59 peiffer Exp $

using namespace std;

// Local include(s):
#include "include/ZprimePostSelectionCycle.h"
#include "include/EventFilterFromListStandAlone.h"

ClassImp( ZprimePostSelectionCycle );

ZprimePostSelectionCycle::ZprimePostSelectionCycle()
    : AnalysisCycle()
{
    // constructor, declare additional variables that should be
    // obtained from the steering-xml file

    m_dobsf = "None";
    m_dotsf = "None";
    m_mttgencut = false;
    m_flavor_selection = "None";
    m_writeeventlist = false;

    // steerable properties for making qcd (pre) selection
    DeclareProperty( "Electron_Or_Muon_Selection", m_Electron_Or_Muon_Selection );
    DeclareProperty( "BTaggingScaleFactors", m_dobsf );
    DeclareProperty( "TopTaggingScaleFactors", m_dotsf );
    DeclareProperty( "ApplyMttbarGenCut", m_mttgencut );
    DeclareProperty( "ApplyFlavorSelection", m_flavor_selection );
    DeclareProperty( "EventFilterFile", m_filter_file );
    DeclareProperty( "WriteEventList", m_writeeventlist);

    // set the integrated luminosity per bin for the lumi-yield control plots
    SetIntLumiPerBin(500.);

    // set the btagging operating point
    m_btagtype = e_CSVT;
    x_btagtype = e_CSVL;

    // apply SF for the "Ele30 OR PFJet320" trigger (electron channel)
    m_applyEleORJetTriggerSF = false;
    DeclareProperty( "applyEleORJetTriggerSF", m_applyEleORJetTriggerSF);
}

ZprimePostSelectionCycle::~ZprimePostSelectionCycle()
{
    // destructor
}

void ZprimePostSelectionCycle::BeginCycle() throw( SError )
{
    // Start of the job, general set-up and definition of
    // objects are done here

    // Important: first call BeginCycle of base class
    AnalysisCycle::BeginCycle();

    return;

}

void ZprimePostSelectionCycle::EndCycle() throw( SError )
{
    // clean-up, info messages and final calculations after the analysis


    // call the base cycle class for all standard methods
    AnalysisCycle::EndCycle();

    return;

}

void ZprimePostSelectionCycle::BeginInputData( const SInputData& id ) throw( SError )
{
    // declaration of histograms and selections

    // Important: first call BeginInputData of base class
    AnalysisCycle::BeginInputData( id );

    // -------------------- set up the selections ---------------------------

    bool doEle=false;
    bool doMu=false;

    if(m_Electron_Or_Muon_Selection=="Electrons" || m_Electron_Or_Muon_Selection=="Electron" || m_Electron_Or_Muon_Selection=="Ele" || m_Electron_Or_Muon_Selection=="ELE") {
        doEle=true;
    } else if(m_Electron_Or_Muon_Selection=="Muon" || m_Electron_Or_Muon_Selection=="Muons" || m_Electron_Or_Muon_Selection=="Mu" || m_Electron_Or_Muon_Selection=="MU") {
        doMu=true;
    } else {
        m_logger << ERROR << "Electron_Or_Muon_Selection is not defined in your xml config file --- should be either `ELE` or `MU`" << SLogger::endmsg;
    }

    // Chi2 reconstruction and discriminant
    static Chi2Discriminator* m_chi2discr = new Chi2Discriminator();

    // event filter for HCAL laser events
    // NOTE: the file is missing, therefore the HCAL laser selection is commented out
    //Selection* HCALlaser = new Selection("HCAL_laser_events");
    //HCALlaser->addSelectionModule(new EventFilterSelection(m_filter_file.c_str()) );

    // Leading jet selection
    Selection* LeadingJetSelection = new Selection("LeadingJetSelection");
    LeadingJetSelection->addSelectionModule(new NJetSelection(1,int_infinity(),150,2.5));

    // Kinematic selection
    Selection* KinematicSelection = new Selection("KinematicSelection");
    KinematicSelection->addSelectionModule(new METCut(50));
    if (doEle)
        KinematicSelection->addSelectionModule(new HypothesisLeptopPtCut( m_chi2discr, 140.0, double_infinity()));

    static Selection* mttbar_gen_selection = new Selection("Mttbar_Gen_Selection");
    if ( m_mttgencut && ((id.GetVersion() == "TTbar_0to700") || (id.GetVersion() == "TTbar") )  ) {
      m_logger << INFO << "Applying mttbar generator cut from 0 to 700 GeV." << SLogger::endmsg;
      mttbar_gen_selection->addSelectionModule(new MttbarGenCut(0,700));
      mttbar_gen_selection->EnableSelection();
    } else {
      m_logger << INFO << "Disabling mttbar generator cut." << SLogger::endmsg;
      mttbar_gen_selection->DisableSelection();
    }

    std::transform(
        m_flavor_selection.begin(), m_flavor_selection.end(), m_flavor_selection.begin(), ::tolower
    );
    if (m_flavor_selection == "bflavor") {
        m_logger << INFO << "Applying b flavor selection" << SLogger::endmsg;
        KinematicSelection->addSelectionModule(new EventFlavorSelection(e_BFlavor));
    } else if (m_flavor_selection == "cflavor") {
        m_logger << INFO << "Applying c flavor selection" << SLogger::endmsg;
        KinematicSelection->addSelectionModule(new EventFlavorSelection(e_CFlavor));
    } else if (m_flavor_selection == "lflavor") {
        m_logger << INFO << "Applying l flavor selection" << SLogger::endmsg;
        KinematicSelection->addSelectionModule(new EventFlavorSelection(e_LFlavor));
    } else if (m_flavor_selection != "none") {
        m_logger << ERROR << "Unknown ApplyFlavorSelection option --- should be either `BFlavor`, `CFlavor` or `LFlavor`" << SLogger::endmsg;
    }

    Selection* Chi2Seletion40 = new Selection("Chi2Selection40");
    Chi2Seletion40->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 40));

    Selection* Chi2Seletion10 = new Selection("Chi2Selection10");
    Chi2Seletion10->addSelectionModule(new HypothesisDiscriminatorCut( m_chi2discr, -1*double_infinity(), 10));

    Selection* BTagSelection = new Selection("BTagSelection");
    BTagSelection->addSelectionModule(new NBTagSelection(1,int_infinity(),m_btagtype)); //at least one b tag

    Selection* NoBTagSelection = new Selection("NoBTagSelection");
    NoBTagSelection->addSelectionModule(new NBTagSelection(0,0,m_btagtype)); //no b tags

    Selection* BTagAntiktJetSelection = new Selection("BTagAntiktJetSelection");
    BTagAntiktJetSelection->addSelectionModule(new NBTagAntiktJetSelection(1,int_infinity(),m_btagtype,1.3)); //at least one AK5 b jet

    Selection* TopTagSelection = new Selection("TopTagSelection");
    TopTagSelection->addSelectionModule(new NTopJetSelection(1,int_infinity(),350,2.5));// top jet
    TopTagSelection->addSelectionModule(new NTopTagSelection(1,int_infinity())); //top tag
    TopTagSelection->addSelectionModule(new TopTagOverlapSelection()); //top tag

    Selection* CMSSubBTagNsubjSelection = new Selection ("CMSSubBTagNsubjSelection");
    CMSSubBTagNsubjSelection->addSelectionModule(new NCMSSubBTagSelection(1,int_infinity(),1,int_infinity(),x_btagtype,0.7)); // with Nsubjettiness cut

    Selection* CMSSubBTagSelection = new Selection ("CMSSubBTagSelection");
    CMSSubBTagSelection->addSelectionModule(new NCMSSubBTagSelection(1,int_infinity(),1,int_infinity(),x_btagtype,int_infinity())); // without Nsubjettiness

    Selection* SumBTags0Selection = new Selection ("SumBTags0Selection");
    SumBTags0Selection->addSelectionModule(new NSumBTagsSelection(0,0,x_btagtype)); //Sum B-Tags

    Selection* SumBTags1Selection = new Selection ("SumBTags1Selection");
    SumBTags1Selection->addSelectionModule(new NSumBTagsSelection(1,1,x_btagtype)); //Sum B-Tags

    Selection* SumBTags2Selection = new Selection ("SumBTags2Selection");
    SumBTags2Selection->addSelectionModule(new NSumBTagsSelection(2,int_infinity(),x_btagtype)); //Sum B-Tags

    RegisterSelection(mttbar_gen_selection);
    //RegisterSelection(HCALlaser);
    RegisterSelection(LeadingJetSelection);
    RegisterSelection(KinematicSelection);
    RegisterSelection(TopTagSelection);
    RegisterSelection(Chi2Seletion40);
    RegisterSelection(Chi2Seletion10);
    RegisterSelection(BTagSelection);
    RegisterSelection(NoBTagSelection);
    RegisterSelection(CMSSubBTagNsubjSelection);
    RegisterSelection(CMSSubBTagSelection);
    RegisterSelection(BTagAntiktJetSelection);
    RegisterSelection(SumBTags0Selection);
    RegisterSelection(SumBTags1Selection);
    RegisterSelection(SumBTags2Selection);

    // ---------------- set up the histogram collections --------------------

    // histograms without any cuts
    RegisterHistCollection( new HypothesisHists("Chi2_Presel", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_Presel") );
    RegisterHistCollection( new JetHists("Jets_Presel") );
    RegisterHistCollection( new ElectronHists("Electron_Presel") );
    RegisterHistCollection( new MuonHists("Muon_Presel") );
    RegisterHistCollection( new TauHists("Tau_Presel") );
    RegisterHistCollection( new TopJetHists("TopJets_Presel") );

    // histograms with leading jet, ....
    RegisterHistCollection( new HypothesisHists("Chi2_LJetsel", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_LJetsel") );
    RegisterHistCollection( new JetHists("Jets_LJetsel") );
    RegisterHistCollection( new ElectronHists("Electron_LJetsel") );
    RegisterHistCollection( new MuonHists("Muon_LJetsel") );
    RegisterHistCollection( new TauHists("Tau_LJetsel") );
    RegisterHistCollection( new TopJetHists("TopJets_LJetsel") );

    // histogram with kinematic cuts
    RegisterHistCollection( new HypothesisHists("Chi2_Kinesel", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_Kinesel") );
    RegisterHistCollection( new JetHists("Jets_Kinesel") );
    RegisterHistCollection( new ElectronHists("Electron_Kinesel") );
    RegisterHistCollection( new MuonHists("Muon_Kinesel") );
    RegisterHistCollection( new TauHists("Tau_Kinesel") );
    RegisterHistCollection( new TopJetHists("TopJets_Kinesel") );

    // histograms with chi2 cut
    RegisterHistCollection( new HypothesisHists("Chi2_Chi2sel", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_Chi2sel") );
    RegisterHistCollection( new JetHists("Jets_Chi2sel") );
    RegisterHistCollection( new ElectronHists("Electron_Chi2sel") );
    RegisterHistCollection( new MuonHists("Muon_Chi2sel") );
    RegisterHistCollection( new TauHists("Tau_Chi2sel") );
    RegisterHistCollection( new TopJetHists("TopJets_Chi2sel") );
    RegisterHistCollection( new BTagEffHists("BTagEff_Chi2sel", m_btagtype) );

    // histograms with Btag and NoBtag and Chi2
    RegisterHistCollection( new HypothesisHists("Chi2_BTag", m_chi2discr ) );
    RegisterHistCollection( new HypothesisHists("Chi2_NoBTag", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_BTag") );
    RegisterHistCollection( new JetHists("Jets_BTag") );
    RegisterHistCollection( new ElectronHists("Electron_BTag") );
    RegisterHistCollection( new MuonHists("Muon_BTag") );
    RegisterHistCollection( new TauHists("Tau_BTag") );
    RegisterHistCollection( new TopJetHists("TopJets_BTag") );

    RegisterHistCollection( new EventHists("Event_NoBTag") );
    RegisterHistCollection( new JetHists("Jets_NoBTag") );
    RegisterHistCollection( new ElectronHists("Electron_NoBTag") );
    RegisterHistCollection( new MuonHists("Muon_NoBTag") );
    RegisterHistCollection( new TauHists("Tau_NoBTag") );
    RegisterHistCollection( new TopJetHists("TopJets_NoBTag") );

    // histograms after the top selection
    RegisterHistCollection( new HypothesisHists("Chi2_TopTag", m_chi2discr ) );

    RegisterHistCollection( new EventHists("Event_TopTag") );
    RegisterHistCollection( new JetHists("Jets_TopTag") );
    RegisterHistCollection( new ElectronHists("Electron_TopTag") );
    RegisterHistCollection( new MuonHists("Muon_TopTag") );
    RegisterHistCollection( new TauHists("Tau_TopTag") );
    RegisterHistCollection( new TopJetHists("TopJets_TopTag") );

    //histograms of Top-Tags with SumBTags

    // No-TopTag (+0,1,>=2 SumBTag)
    RegisterHistCollection( new HypothesisHists("Chi2_NoTopTagSumBTag0",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_NoTopTagSumBTag0") );
    RegisterHistCollection( new JetHists("Jets_NoTopTagSumBTag0") );
    RegisterHistCollection( new ElectronHists("Electron_NoTopTagSumBTag0") );
    RegisterHistCollection( new MuonHists("Muon_NoTopTagSumBTag0") );
    RegisterHistCollection( new TauHists("Tau_NoTopTagSumBTag0") );
    RegisterHistCollection( new TopJetHists("TopJets_NoTopTagSumBTag0") );

    RegisterHistCollection( new HypothesisHists("Chi2_NoTopTagSumBTag1",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_NoTopTagSumBTag1") );
    RegisterHistCollection( new JetHists("Jets_NoTopTagSumBTag1") );
    RegisterHistCollection( new ElectronHists("Electron_NoTopTagSumBTag1") );
    RegisterHistCollection( new MuonHists("Muon_NoTopTagSumBTag1") );
    RegisterHistCollection( new TauHists("Tau_NoTopTagSumBTag1") );
    RegisterHistCollection( new TopJetHists("TopJets_NoTopTagSumBTag1") );

    RegisterHistCollection( new HypothesisHists("Chi2_NoTopTagSumBTag2",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_NoTopTagSumBTag2") );
    RegisterHistCollection( new JetHists("Jets_NoTopTagSumBTag2") );
    RegisterHistCollection( new ElectronHists("Electron_NoTopTagSumBTag2") );
    RegisterHistCollection( new MuonHists("Muon_NoTopTagSumBTag2") );
    RegisterHistCollection( new TauHists("Tau_NoTopTagSumBTag2") );
    RegisterHistCollection( new TopJetHists("TopJets_NoTopTagSumBTag2") );

    // TopTag (+0,1,>=2 SumBTag)
    RegisterHistCollection( new HypothesisHists("Chi2_TopTagSumBTag0",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_TopTagSumBTag0") );
    RegisterHistCollection( new JetHists("Jets_TopTagSumBTag0") );
    RegisterHistCollection( new ElectronHists("Electron_TopTagSumBTag0") );
    RegisterHistCollection( new MuonHists("Muon_TopTagSumBTag0") );
    RegisterHistCollection( new TauHists("Tau_TopTagSumBTag0") );
    RegisterHistCollection( new TopJetHists("TopJets_TopTagSumBTag0") );

    RegisterHistCollection( new HypothesisHists("Chi2_TopTagSumBTag1",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_TopTagSumBTag1") );
    RegisterHistCollection( new JetHists("Jets_TopTagSumBTag1") );
    RegisterHistCollection( new ElectronHists("Electron_TopTagSumBTag1") );
    RegisterHistCollection( new MuonHists("Muon_TopTagSumBTag1") );
    RegisterHistCollection( new TauHists("Tau_TopTagSumBTag1") );
    RegisterHistCollection( new TopJetHists("TopJets_TopTagSumBTag1") );

    RegisterHistCollection( new HypothesisHists("Chi2_TopTagSumBTag2",m_chi2discr));
    RegisterHistCollection( new EventHists("Event_TopTagSumBTag2") );
    RegisterHistCollection( new JetHists("Jets_TopTagSumBTag2") );
    RegisterHistCollection( new ElectronHists("Electron_TopTagSumBTag2") );
    RegisterHistCollection( new MuonHists("Muon_TopTagSumBTag2") );
    RegisterHistCollection( new TauHists("Tau_TopTagSumBTag2") );
    RegisterHistCollection( new TopJetHists("TopJets_TopTagSumBTag2") );

    // important: initialise histogram collections after their definition
    InitHistos();

    // Data-MC b-tagging reweighting
    m_bsf = NULL;
    std::transform(m_dobsf.begin(), m_dobsf.end(), m_dobsf.begin(), ::tolower);
    if(m_dobsf != "none") {
        E_SystShift sys_bjets = e_Default;
        E_SystShift sys_ljets = e_Default;
        if (m_dobsf == "default") {
            m_logger << INFO << "Applying btagging scale factor" << SLogger::endmsg;
        } else if (m_dobsf == "up-bjets") {
            m_logger << INFO << "Applying btagging up scale factor for b-jets" << SLogger::endmsg;
            sys_bjets = e_Up;
        } else if (m_dobsf == "down-bjets") {
            m_logger << INFO << "Applying btagging down scale factor for b-jets" << SLogger::endmsg;
            sys_bjets = e_Down;
        } else if (m_dobsf == "up-ljets") {
            m_logger << INFO << "Applying btagging up scale factor for l-jets" << SLogger::endmsg;
            sys_ljets = e_Up;
        } else if (m_dobsf == "down-ljets") {
            m_logger << INFO << "Applying btagging down scale factor for l-jets" << SLogger::endmsg;
            sys_ljets = e_Down;
        }
        else
            m_logger << ERROR << "Unknown BTaggingScaleFactors option, default option is applied --- should be either `Default`, `Up-bjets`, `Down-bjets`, `Up-ljets`, or `Down-ljets`" << SLogger::endmsg;
        if(doEle)
            m_bsf = new BTaggingScaleFactors(x_btagtype, e_Electron, sys_bjets, sys_ljets);
        else if(doMu)
            m_bsf = new BTaggingScaleFactors(x_btagtype, e_Muon, sys_bjets, sys_ljets);
    }

    m_tsf = NULL;
    std::transform(m_dotsf.begin(), m_dotsf.end(), m_dotsf.begin(), ::tolower);
    if(m_dotsf != "none") {
        E_SystShift sys_toptag = e_Default;
        E_SystShift sys_mistag = e_Default;
        if (m_dotsf == "default") {
            m_logger << INFO << "Applying toptagging scale factor" << SLogger::endmsg;
        } else if (m_dotsf == "up-mistag") {
            m_logger << INFO << "Applying topmistagging up scale factor" << SLogger::endmsg;
            sys_mistag = e_Up;
        } else if (m_dotsf == "down-mistag") {
            m_logger << INFO << "Applying topmistagging down scale factor" << SLogger::endmsg;
            sys_mistag = e_Down;
        } else if (m_dotsf == "up-toptag") {
            m_logger << INFO << "Applying toptagging up scale factor" << SLogger::endmsg;
            sys_toptag = e_Up;
        } else if (m_dotsf == "down-toptag") {
            m_logger << INFO << "Applying toptagging down scale factor" << SLogger::endmsg;
            sys_toptag = e_Down;
        }
        else
            m_logger << ERROR << "Unknown TopTaggingScaleFactors option, default option is applied --- should be either `Default`, `Up`, or `Down`" << SLogger::endmsg;

        m_tsf = new TopTaggingScaleFactors(sys_toptag, sys_mistag);
    }


    if(m_writeeventlist)
      m_eventlist.open( id.GetVersion()+"_eventlist.txt" );

    return;
}

void ZprimePostSelectionCycle::EndInputData( const SInputData& id ) throw( SError )
{
    AnalysisCycle::EndInputData( id );
    if(m_writeeventlist)
      m_eventlist.close();
    return;
}

void ZprimePostSelectionCycle::BeginInputFile( const SInputData& id ) throw( SError )
{
    // Connect all variables from the Ntuple file with the ones needed for the analysis
    // The variables are commonly stored in the BaseCycleContaincer

    // important: call to base function to connect all variables to Ntuples from the input tree
    AnalysisCycle::BeginInputFile( id );

    return;
}

void ZprimePostSelectionCycle::ExecuteEvent( const SInputData& id, Double_t weight) throw( SError )
{
    // this is the most important part: here the full analysis happens
    // user should implement selections, filling of histograms and results

    // first step: call Execute event of base class to perform basic consistency checks
    // also, the good-run selection is performed there and the calculator is reset
    AnalysisCycle::ExecuteEvent( id, weight );

    EventCalc* calc = EventCalc::Instance();

    // get the selections
    //static Selection* HCALlaser = GetSelection("HCAL_laser_events");
    static Selection* LeadingJetSelection = GetSelection("LeadingJetSelection");
    static Selection* KinematicSelection = GetSelection("KinematicSelection");
    static Selection* TopTagSelection = GetSelection("TopTagSelection");
    static Selection* Chi2Selection40 = GetSelection("Chi2Selection40");
    static Selection* Chi2Selection10 = GetSelection("Chi2Selection10");
    static Selection* BTagSelection = GetSelection("BTagSelection");
    static Selection* NoBTagSelection = GetSelection("NoBTagSelection");
    static Selection* CMSSubBTagNsubjSelection = GetSelection("CMSSubBTagNsubjSelection");
    static Selection* CMSSubBTagSelection = GetSelection("CMSSubBTagSelection");
    static Selection* SumBTags0Selection = GetSelection("SumBTags0Selection");
    static Selection* SumBTags1Selection = GetSelection("SumBTags1Selection");
    static Selection* SumBTags2Selection = GetSelection("SumBTags2Selection");

    // get the histogram collections
    BaseHists* Chi2_HistsPresel = GetHistCollection("Chi2_Presel");
    BaseHists* Chi2_HistsLJetsel = GetHistCollection("Chi2_LJetsel");
    BaseHists* Chi2_HistsKinesel = GetHistCollection("Chi2_Kinesel");
    BaseHists* Chi2_HistsChi2sel = GetHistCollection("Chi2_Chi2sel");
    BaseHists* Chi2_HistsBTag = GetHistCollection("Chi2_BTag");
    BaseHists* Chi2_HistsNoBTag = GetHistCollection("Chi2_NoBTag");
    BaseHists* Chi2_HistsTopTag = GetHistCollection("Chi2_TopTag");
    BaseHists* Chi2_HistsNoTopTagSumBTag0 = GetHistCollection("Chi2_NoTopTagSumBTag0");
    BaseHists* Chi2_HistsNoTopTagSumBTag1 = GetHistCollection("Chi2_NoTopTagSumBTag1");
    BaseHists* Chi2_HistsNoTopTagSumBTag2 = GetHistCollection("Chi2_NoTopTagSumBTag2");
    BaseHists* Chi2_HistsTopTagSumBTag0 = GetHistCollection("Chi2_TopTagSumBTag0");
    BaseHists* Chi2_HistsTopTagSumBTag1 = GetHistCollection("Chi2_TopTagSumBTag1");
    BaseHists* Chi2_HistsTopTagSumBTag2 = GetHistCollection("Chi2_TopTagSumBTag2");
    BaseHists* BTagEff_HistsChi2sel = GetHistCollection("BTagEff_Chi2sel");


    // reject laser events only for data
    //if (calc->IsRealData()){
    //  if (!HCALlaser->passSelection()) throw SError( SError::SkipEvent );
    //}

    // cut out events from the inclusive ttbar sample to avoid double counting
    static Selection* mttbar_gen_selection = GetSelection("Mttbar_Gen_Selection");
    if(!mttbar_gen_selection->passSelection())  throw SError( SError::SkipEvent );

    // b tagging scale factor
    if(m_bsf && m_addGenInfo) {
        calc->ProduceWeight(m_bsf->GetWeight());
    }

    // top tagging scale factor
    if(m_tsf && m_addGenInfo) {
        calc->ProduceWeight(m_tsf->GetWeight());
    }

    // Ele30_OR_PFJet320 trigger Scale Factor
    if(m_applyEleORJetTriggerSF && !calc->IsRealData()) calc->ProduceWeight( m_lsf->GetElectronORJetTrigWeight() );

    if(calc->GetJets()->size()>=12) {
        std::cout << "run: " << calc->GetRunNum() << "   lb: " << calc->GetLumiBlock() << "  event: " << calc->GetEventNum() << "   N(jets): " << calc->GetJets()->size() << std::endl;
    }

    Chi2_HistsPresel->Fill();
    FillControlHistos("_Presel");

    if(!LeadingJetSelection->passSelection())  throw SError( SError::SkipEvent );

    Chi2_HistsLJetsel->Fill();
    FillControlHistos("_LJetsel");

    if(!KinematicSelection->passSelection())  throw SError( SError::SkipEvent );

    Chi2_HistsKinesel->Fill();
    FillControlHistos("_Kinesel");

    if(!Chi2Selection40->passSelection()) throw SError( SError::SkipEvent );

    Chi2_HistsChi2sel->Fill();
    if(m_addGenInfo) BTagEff_HistsChi2sel->Fill();
    FillControlHistos("_Chi2sel");

    // BTag-NoBTag categories: do a chi2 selection of 10 for comparison with published analysis
    if(Chi2Selection10->passSelection()){
      if(BTagSelection->passSelection()) {
        Chi2_HistsBTag->Fill();
        FillControlHistos("_BTag");
      }
      if(NoBTagSelection->passSelection()) {
        Chi2_HistsNoBTag->Fill();
        FillControlHistos("_NoBTag");
      }
    }

    // NoTopTag categories
    if(!TopTagSelection->passSelection()){

      if(SumBTags0Selection->passSelection()){
        Chi2_HistsNoTopTagSumBTag0->Fill();
        FillControlHistos("_NoTopTagSumBTag0");
      }
      if(SumBTags1Selection->passSelection()){
        Chi2_HistsNoTopTagSumBTag1->Fill();
        FillControlHistos("_NoTopTagSumBTag1");
      }
      if(SumBTags2Selection->passSelection()){
        Chi2_HistsNoTopTagSumBTag2->Fill();
        FillControlHistos("_NoTopTagSumBTag2");
      }
    }

    // TopTag categories
    if(TopTagSelection->passSelection()){

      Chi2_HistsTopTag->Fill();
      FillControlHistos("_TopTag");

      if(SumBTags0Selection->passSelection()){
        Chi2_HistsTopTagSumBTag0->Fill();
        FillControlHistos("_TopTagSumBTag0");
      }
      if(SumBTags1Selection->passSelection()){
        Chi2_HistsTopTagSumBTag1->Fill();
        FillControlHistos("_TopTagSumBTag1");
      }
      if(SumBTags2Selection->passSelection()){
        Chi2_HistsTopTagSumBTag2->Fill();
        FillControlHistos("_TopTagSumBTag2");
      }
    }

    if(m_writeeventlist){
      if(id.GetType()=="DATA" || id.GetType()=="Data" || id.GetType()=="data" )
	m_eventlist << calc->GetRunNum() << ":" << calc->GetLumiBlock() << ":" << calc->GetEventNum() << std::endl;
      else
	//don't fill the random run number produced by LumiHandler for MC samples
	m_eventlist << "1:" << calc->GetLumiBlock() << ":" << calc->GetEventNum() << std::endl;
    }

    return;
}


void ZprimePostSelectionCycle::FillControlHistos(TString postfix)
{
    // fill some control histograms, need to be defined in BeginInputData

    BaseHists* eventhists = GetHistCollection((std::string)("Event"+postfix));
    BaseHists* jethists = GetHistCollection((std::string)("Jets"+postfix));
    BaseHists* elehists = GetHistCollection((std::string)("Electron"+postfix));
    BaseHists* muonhists = GetHistCollection((std::string)("Muon"+postfix));
    BaseHists* tauhists = GetHistCollection((std::string)("Tau"+postfix));
    BaseHists* topjethists = GetHistCollection((std::string)("TopJets"+postfix));

    eventhists->Fill();
    jethists->Fill();
    elehists->Fill();
    muonhists->Fill();
    tauhists->Fill();
    topjethists->Fill();
}
