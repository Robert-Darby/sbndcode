// Include necessary headers
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Include necessary ROOT headers
#include "TTree.h"

// Include sbnd StoppingMuonTrigger
#include "sbnobj/SBND/LateLight/LateLight.h"

namespace sbnd
{

  class LateLightAnalyser;
}

class sbnd::LateLightAnalyser : public art::EDAnalyzer
{
public:
  explicit LateLightAnalyser(fhicl::ParameterSet const &p);
  void analyze(art::Event const &e) override;

  LateLightAnalyser(LateLightAnalyser const &) = delete;
  LateLightAnalyser(LateLightAnalyser &&) = delete;
  LateLightAnalyser &operator=(LateLightAnalyser const &) = delete;
  LateLightAnalyser &operator=(LateLightAnalyser &&) = delete;

private:
  // Variables to store in the TTree
  float fMuonTime;
  float fMichelTime;
  int fCRTPlane;
  int fG4ID;
  int fG4PDG;
  bool fHasOpFlash;

  // Variables for LateLightTail information
  std::vector<int> fChannel;
  std::vector<float> fStartTime;
  std::vector<std::vector<float>> fWaveforms;

  // Event identification numbers
  unsigned int fEvent;
  unsigned int fSubRun;
  unsigned int fRun;

  // Tree
  TTree *fTree;

  // Input tag for StoppingMuonTrigger
  art::InputTag fStoppingMuonTriggerTag;
};

// Constructor
sbnd::LateLightAnalyser::LateLightAnalyser(fhicl::ParameterSet const &p)
    : EDAnalyzer(p),
      fStoppingMuonTriggerTag(p.get<std::string>("StoppingMuonTriggerTag"))
{
  // Get access to the TFileService to create the TTree
  art::ServiceHandle<art::TFileService> tfs;

  // Create the TTree
  fTree = tfs->make<TTree>("LateLightTree", "Tree storing StoppingMuonTrigger information");

  // Setup branches for the tree
  fTree->Branch("MuonTime", &fMuonTime, "MuonTime/F");
  fTree->Branch("MichelTime", &fMichelTime, "MichelTime/F");
  fTree->Branch("CRTPlane", &fCRTPlane, "CRTPlane/I");
  fTree->Branch("G4ID", &fG4ID, "G4ID/I");
  fTree->Branch("G4PDG", &fG4PDG, "G4PDG/I");
  fTree->Branch("HasOpFlash", &fHasOpFlash, "HasOpFlash/O");

  // Branches for the LateLightTails
  fTree->Branch("Channel", &fChannel);
  fTree->Branch("StartTime", &fStartTime);
  fTree->Branch("Waveforms", &fWaveforms);

  // Branches for the event identification numbers
  fTree->Branch("Event", &fEvent, "Event/i");
  fTree->Branch("SubRun", &fSubRun, "SubRun/i");
  fTree->Branch("Run", &fRun, "Run/i");
}

// Analyze method, called once per event
void sbnd::LateLightAnalyser::analyze(art::Event const &e)
{
  // Get the event identification numbers
  fEvent = e.event();
  fSubRun = e.subRun();
  fRun = e.run();

  // Get the StoppingMuonTrigger from the event
  art::Handle<std::vector<sbnd::StoppingMuonTrigger>> handle;
  e.getByLabel(fStoppingMuonTriggerTag, handle);
  std::cout << (*handle).size() << "\n";
  if (!handle.isValid())
  {
    // Handle the case where the product is not found
    mf::LogError("LateLightAnalyser") << "StoppingMuonTrigger not found in event " << e.event();
    return;
  }

  // Access the data from StoppingMuonTrigger
  for (const sbnd::StoppingMuonTrigger &trigger : *handle)
  {

    // Fill the main trigger information
    fMuonTime = trigger.MuonTime;
    fMichelTime = trigger.MichelTime;
    fCRTPlane = trigger.CRTPlane;
    fG4ID = trigger.G4ID;
    fG4PDG = trigger.G4PDG;
    fHasOpFlash = trigger.HasOpFlash;

    // Clear vectors for LateLightTail information
    fChannel.clear();
    fStartTime.clear();
    fWaveforms.clear();

    // Fill the LateLightTail information
    std::cout << trigger.LateLightTails.size() << "\n";
    for (const auto &tail : trigger.LateLightTails)
    {
      fChannel.push_back(tail.Channel);
      fStartTime.push_back(tail.StartTime);
      fWaveforms.push_back(tail.Waveform); // Store the waveform as a vector of floats
    }

    // Fill the tree
    fTree->Fill();
  }
}

// Define the module as a plugin
DEFINE_ART_MODULE(sbnd::LateLightAnalyser)
