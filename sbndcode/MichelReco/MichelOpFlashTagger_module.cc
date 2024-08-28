////////////////////////////////////////////////////////////////////////
// Class:       MichelOpFlashTagger
// Plugin Type: producer (Unknown Unknown)
// File:        MichelOpFlashTagger_module.cc
//
// Generated at Mon Aug 12 05:58:45 2024 by Robert Darby using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace sbnd {
  class MichelOpFlashTagger;
}


class sbnd::MichelOpFlashTagger : public art::EDProducer {
public:
  explicit MichelOpFlashTagger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelOpFlashTagger(MichelOpFlashTagger const&) = delete;
  MichelOpFlashTagger(MichelOpFlashTagger&&) = delete;
  MichelOpFlashTagger& operator=(MichelOpFlashTagger const&) = delete;
  MichelOpFlashTagger& operator=(MichelOpFlashTagger&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // FHICL Params
  const std::string fRawWaveformLabel;
  const float fWindowLength;
  const float fMichelOpHitThreshold;
  const float fMuonOpHitThreshold;
};


sbnd::MichelOpFlashTagger::MichelOpFlashTagger(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
  , fRawWaveformLabel(p.get<std::string>("RawWaveformLabel"))
  , fWindowLength(p.get<float>("WindowLength"))
  , fMichelOpHitThreshold(p.get<float>("MichelOpHitThreshold"))
  , fMuonOpHitThreshold(p.get<float>("MuonOpHitThreshold"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::MichelOpFlashTagger::produce(art::Event& e)
{
  // Uses deconvolved waveforms
}

void sbnd::MichelOpFlashTagger::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::MichelOpFlashTagger::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::MichelOpFlashTagger)
