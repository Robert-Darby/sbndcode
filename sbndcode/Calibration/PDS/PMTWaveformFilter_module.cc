////////////////////////////////////////////////////////////////////////
// Class:       PMTWaveformFilter
// Plugin Type: filter (Unknown Unknown)
// File:        PMTWaveformFilter_module.cc
//
// Generated at Tue Aug 27 05:31:24 2024 by Robert Darby using cetskelgen
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

#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>
#include <vector>
#include <string>

namespace sbnd {
  class PMTWaveformFilter;
}


class sbnd::PMTWaveformFilter : public art::EDFilter {
public:
  explicit PMTWaveformFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTWaveformFilter(PMTWaveformFilter const&) = delete;
  PMTWaveformFilter(PMTWaveformFilter&&) = delete;
  PMTWaveformFilter& operator=(PMTWaveformFilter const&) = delete;
  PMTWaveformFilter& operator=(PMTWaveformFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  const std::string fWaveformInstanceName;
};


sbnd::PMTWaveformFilter::PMTWaveformFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
  ,fWaveformInstanceName(p.get<std::string>("WaveformInstanceName", "PMTChannels"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool sbnd::PMTWaveformFilter::filter(art::Event& e)
{
   std::vector<art::Handle<std::vector<raw::OpDetWaveform>>> waveformHandles = e.getMany<std::vector<raw::OpDetWaveform>>();
   std::cout << "\n\n\n Module Initialised \n";

   std::string wv_name = fWaveformInstanceName;
   waveformHandles.erase(
     std::remove_if(waveformHandles.begin(), waveformHandles.end(),
            [&wv_name](const art::Handle<std::vector<raw::OpDetWaveform>>& handle) {
                if (handle.isValid()) {
                    auto const& prov = handle.provenance();
                    return prov->productInstanceName() != wv_name;}
                return false;}),
     waveformHandles.end());

  if(waveformHandles[0]->size()==0 || waveformHandles.empty()) return false;
  return true;
}

void sbnd::PMTWaveformFilter::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::PMTWaveformFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::PMTWaveformFilter)
