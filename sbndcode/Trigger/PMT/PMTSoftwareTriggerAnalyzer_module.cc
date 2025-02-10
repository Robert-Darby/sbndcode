////////////////////////////////////////////////////////////////////////
// Class:       PMTSoftwareTriggerAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        PMTSoftwareTriggerAnalyzer_module.cc
//
// Generated at Wed Feb  5 05:07:08 2025 by Robert Darby using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h" // Find associations as pointers
#include "canvas/Persistency/Common/FindOneP.h"

#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"

#include "TTree.h"

namespace sbnd
{
  class PMTSoftwareTriggerAnalyzer;
}

class sbnd::PMTSoftwareTriggerAnalyzer : public art::EDAnalyzer
{
public:
  explicit PMTSoftwareTriggerAnalyzer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTSoftwareTriggerAnalyzer(PMTSoftwareTriggerAnalyzer const &) = delete;
  PMTSoftwareTriggerAnalyzer(PMTSoftwareTriggerAnalyzer &&) = delete;
  PMTSoftwareTriggerAnalyzer &operator=(PMTSoftwareTriggerAnalyzer const &) = delete;
  PMTSoftwareTriggerAnalyzer &operator=(PMTSoftwareTriggerAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  // Declare member data here.
  TTree *fTree;
  int fRun;
  int fSubrun;
  int fEvent;
  int fNAboveThreshold;
  int fTrigTs;
  double fPromptPE;
  double fPrelimPE;
  double fPeakPE;
  double fPeakTime;
  double fMichelPeakTime;
  double fMichelPeakPE;

  const std::string fTriggerLabel;
};

sbnd::PMTSoftwareTriggerAnalyzer::PMTSoftwareTriggerAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      fTriggerLabel(p.get<std::string>("TriggerLabel", "pmtmetricproducermichel")) // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::PMTSoftwareTriggerAnalyzer::analyze(art::Event const &e)
{
  // Implementation of required member function here.
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  art::Handle<std::vector<sbnd::trigger::pmtSoftwareTrigger>> triggerHandle;
  std::vector<art::Ptr<sbnd::trigger::pmtSoftwareTrigger>> triggerVect;
  if (e.getByLabel(fTriggerLabel, triggerHandle)) // Make sure artHandle is from module
    art::fill_ptr_vector(triggerVect, triggerHandle);

  if (!triggerHandle)
  {
    mf::LogWarning("PMTTriggerAnalyzer") << "No PMTSoftwareTrigger data found in event.";
    return;
  }

  for (auto const &trigger : triggerVect)
  {
    fNAboveThreshold = trigger->nAboveThreshold;
    fTrigTs = trigger->trig_ts;
    fPromptPE = trigger->promptPE;
    fPrelimPE = trigger->prelimPE;
    fPeakPE = trigger->peakPE;
    fPeakTime = trigger->peaktime;
    fMichelPeakTime = trigger->michelpeaktime;
    fMichelPeakPE = trigger->michelpeakPE;
    fTree->Fill();
  }
}

void sbnd::PMTSoftwareTriggerAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TriggerTree", "PMT Trigger Data");
  fTree->Branch("Run", &fRun);
  fTree->Branch("Subrun", &fSubrun);
  fTree->Branch("Event", &fEvent);
  fTree->Branch("NAboveThreshold", &fNAboveThreshold, "NAboveThreshold/I");
  fTree->Branch("TrigTs", &fTrigTs, "TrigTs/I");
  fTree->Branch("PromptPE", &fPromptPE, "PromptPE/D");
  fTree->Branch("PrelimPE", &fPrelimPE, "PrelimPE/D");
  fTree->Branch("PeakPE", &fPeakPE, "PeakPE/D");
  fTree->Branch("PeakTime", &fPeakTime, "PeakTime/D");
  fTree->Branch("MichelPeakTime", &fMichelPeakTime);
  fTree->Branch("MichelPeakPE", &fMichelPeakPE);
}

DEFINE_ART_MODULE(sbnd::PMTSoftwareTriggerAnalyzer)
