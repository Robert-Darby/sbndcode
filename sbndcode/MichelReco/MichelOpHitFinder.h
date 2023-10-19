#ifndef sbndcode_MichelReco_MichelOpHitFinder_h
#define sbndcode_MichelReco_MichelOpHitFinder_h
////////////////////////////////////////////////////////////////////////
// Class:       MichelOpHitFinder
// Plugin Type: producer (Unknown Unknown)
// File:        MichelOpHitFinder.h
//
// Generated at Wed Oct 18 11:51:42 2023 by Robert Darby using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Provenance/ProductID.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

namespace sbnd {
  class MichelOpHitFinder;
}


class sbnd::MichelOpHitFinder : public art::EDProducer {
public:
  explicit MichelOpHitFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelOpHitFinder(MichelOpHitFinder const&) = delete;
  MichelOpHitFinder(MichelOpHitFinder&&) = delete;
  MichelOpHitFinder& operator=(MichelOpHitFinder const&) = delete;
  MichelOpHitFinder& operator=(MichelOpHitFinder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Member functions
  void initTree(void);
  void resetTree(void);

  // Declare member data here.

  // FHICL Parameters
  const std::string fOpDetWaveformLabel;
  const std::vector<std::string> fOpDetTypes;
  const std::vector<float> fMinADC;
  const float fMuonDelay;
  const bool fPositivePolarity;
  const bool fMakeTree;
  float fSamplingFreq;

  // Tree Variables
  TTree* _mophf_tree;
  int _run, _sub, _evt;

  std::vector<double> _muon_times;
};

#endif /* sbndcode_MichelReco_MichelOpHitFinder_h */
