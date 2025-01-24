////////////////////////////////////////////////////////////////////////
// Class:       StoppingMuonFilter
// Plugin Type: filter (Unknown Unknown)
// File:        StoppingMuonFilter_module.cc
//
// Generated at Fri Oct  4 09:34:49 2024 by Robert Darby using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>

namespace sbnd {
  class StoppingMuonFilter;
}

class sbnd::StoppingMuonFilter : public art::EDFilter {
public:
  explicit StoppingMuonFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoppingMuonFilter(StoppingMuonFilter const&) = delete;
  StoppingMuonFilter(StoppingMuonFilter&&) = delete;
  StoppingMuonFilter& operator=(StoppingMuonFilter const&) = delete;
  StoppingMuonFilter& operator=(StoppingMuonFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  art::InputTag fMCParticleTag; // Input tag for MCParticle collection

  double fMinX, fMaxX;
  double fMinY, fMaxY;
  double fMinZ, fMaxZ;
};

sbnd::StoppingMuonFilter::StoppingMuonFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fMCParticleTag(p.get<std::string>("MCParticleTag", "largeant")),
    fMinX(p.get<double>("MinX", -200.0)),
    fMaxX(p.get<double>("MaxX", 200.0)),
    fMinY(p.get<double>("MinY", -200.0)),
    fMaxY(p.get<double>("MaxY", 200.0)),
    fMinZ(p.get<double>("MinZ", 0.0)),
    fMaxZ(p.get<double>("MaxZ", 500.0))
{
  // Call appropriate produces<>() functions here if needed.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool sbnd::StoppingMuonFilter::filter(art::Event& e)
{
  // Retrieve the MCParticle collection from the event
  auto const& mcParticleHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
  auto const& mcParticles = *mcParticleHandle;

  // Loop through each MCParticle
  for (auto const& particle : mcParticles) {
    // Check if the particle is a muon (PDG code ±13)
    if (std::abs(particle.PdgCode()) == 13) {
      // Get the stopping point coordinates
      auto const& endX = particle.EndX();
      auto const& endY = particle.EndY();
      auto const& endZ = particle.EndZ();

      // Check if the stopping point is within the specified boundaries
      if (fMinX < endX && endX < fMaxX &&
          fMinY < endY && endY < fMaxY &&
          fMinZ < endZ && endZ < fMaxZ &&
          particle.EndProcess()=="Decay") {
        // If a muon satisfies the condition, accept the event
        return true;
      }
    }
  }

  // If no muons satisfy the condition, reject the event
  return false;
}

void sbnd::StoppingMuonFilter::beginJob()
{
  // Optional: Code to execute before starting event processing.
  mf::LogInfo("StoppingMuonFilter") << "Starting StoppingMuonFilter job...";
}

void sbnd::StoppingMuonFilter::endJob()
{
  // Optional: Code to execute after finishing event processing.
  mf::LogInfo("StoppingMuonFilter") << "Stopping StoppingMuonFilter job.";
}

DEFINE_ART_MODULE(sbnd::StoppingMuonFilter)
