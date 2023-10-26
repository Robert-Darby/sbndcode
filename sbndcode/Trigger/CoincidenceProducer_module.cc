////////////////////////////////////////////////////////////////////////
// Class:       CoincidenceProducer
// Plugin Type: producer (Unknown Unknown)
// File:        CoincidenceProducer_module.cc
//
// Generated at Wed Oct 25 14:55:24 2023 by Robert Darby using cetskelgen
// from  version .
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

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"

#include "sbnobj/SBND/Trigger/Coincidence.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include<map>
#include<vector>
#include <memory>

namespace sbnd {
  class CoincidenceProducer;
}


class sbnd::CoincidenceProducer : public art::EDProducer {
public:
  explicit CoincidenceProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CoincidenceProducer(CoincidenceProducer const&) = delete;
  CoincidenceProducer(CoincidenceProducer&&) = delete;
  CoincidenceProducer& operator=(CoincidenceProducer const&) = delete;
  CoincidenceProducer& operator=(CoincidenceProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void getCRTTimeStamps(
    artdaq::Fragment& frag,
    std::multimap<double, unsigned>& fragTimeStamps);
  void getCAEN1730FragmentTimeStamp(
    const artdaq::Fragment &frag,
    std::multimap<double, unsigned>& fragTimeStamps);
  void getCoincidence(
    const std::multimap<double, unsigned>& fragTimeStamps,
    std::unique_ptr<std::vector<sbnd::Coincidence>>& coincidence_v);

  // FHICL Params
  const double fCoincidenceWindow;
  const float fCRTTriggerOffset;
  const float fCRTClockFreq;
  const double fCAENTriggerOffset;
  const bool fVerbose;
};

sbnd::CoincidenceProducer::CoincidenceProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  , fCoincidenceWindow(p.get<double>("CoincidenceWindow", 0.24))
  , fCRTTriggerOffset(p.get<float>("CRTTriggerOffset", 1700000))
  , fCRTClockFreq(p.get<float>("CRTClockFreq", 1000.))
  , fCAENTriggerOffset(p.get<double>("CAENTriggerOffset", 0.5))
  , fVerbose(p.get<bool>("Verbose", false))
{
  produces<std::vector<sbnd::Coincidence>>();
}

void sbnd::CoincidenceProducer::produce(art::Event& e)
{
  std::unique_ptr<std::vector<sbnd::Coincidence>>
    coincidence_v(new std::vector<sbnd::Coincidence>);

  std::multimap<double, unsigned> fragTimeStamps;
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  // loop over fragment handles
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      // container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() == sbndaq::detail::FragmentType::BERNCRTV2){
          for (size_t ii = 0; ii < contf.block_count(); ++ii) getCRTTimeStamps(*contf[ii].get(), fragTimeStamps);
        } // Check fragment type
      } // Fragment loop
    } // if CRT Fragment
    else {
      // normal fragment
      for (auto frag : *handle){
        if (frag.type()==sbndaq::detail::FragmentType::BERNCRTV2) {
          getCRTTimeStamps(frag, fragTimeStamps);
        }
        else if (frag.type()==sbndaq::detail::FragmentType::CAENV1730) {
          getCAEN1730FragmentTimeStamp(frag, fragTimeStamps);
        } //if is pmt frag
      }
    } // Fragment loop
  } // Fragment handle loop
  getCoincidence(fragTimeStamps, coincidence_v);
  for(auto coinc = coincidence_v->begin(); coinc != coincidence_v->end(); coinc++) std::cout << "Coincidence found at: " << coinc->TimeStamp << std::endl;
  e.put(std::move(coincidence_v));  
}

void sbnd::CoincidenceProducer::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::CoincidenceProducer::endJob()
{
  // Implementation of optional member function here.
}

void sbnd::CoincidenceProducer::getCRTTimeStamps(
  artdaq::Fragment& frag,
  std::multimap<double, unsigned>& fragTimeStamps)
{
  // use  fragment ID to get plane information
  sbndaq::BernCRTFragmentV2 bern_fragment(frag);
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  auto thisone = frag.fragmentID();  uint plane = (thisone & 0x0700) >> 8;
  if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}
  // Check if this fragment is considered in the coincidence check
  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
    // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
    auto thisflag = bevt->flags;
    if (thisflag & 0x2 && !(thisflag & 0xC) ) {
      auto thistime = (bevt->ts1 - fCRTTriggerOffset) / fCRTClockFreq;
      fragTimeStamps.insert(std::pair<float, size_t>(thistime, plane));
    }
  }
} // CoincidenceProducer::getCRTTimeStamps

void sbnd::CoincidenceProducer::getCAEN1730FragmentTimeStamp(
  const artdaq::Fragment &frag,
  std::multimap<double, unsigned>& fragTimeStamps)
{
  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp
  double timestamp = (double)md->timeStampNSec;
  timestamp -= fCAENTriggerOffset*1e9; timestamp /= 1000;
  fragTimeStamps.insert(std::pair<float, size_t>(timestamp, 8));
}

void sbnd::CoincidenceProducer::getCoincidence(
  const std::multimap<double, unsigned>& fragTimeStamps,
  std::unique_ptr<std::vector<sbnd::Coincidence>>& coincidence_v)
{
  // Loop over chronological CRT/PDS tiggers
  for(auto it = fragTimeStamps.begin(); it != fragTimeStamps.end(); it++) {
    if(fVerbose) std::cout << "Plane: " << it->second << " Time: " << it->first << std::endl;
    std::vector<unsigned> coinc_frags; coinc_frags.push_back(it->second);
    std::vector<bool> frag_present(8, false);
    frag_present[it->second] = true;
    // Loop over all triggers in specified coincidence window
    float window_end = it->first + fCoincidenceWindow;
    int frags_in_window = 0;
    for(auto jt = it; jt != fragTimeStamps.end(); jt++) {
      frags_in_window++;
      if(jt->first > window_end) break;
      if(jt->second == it->second) continue;
      coinc_frags.push_back(jt->second);
    }
    if(coinc_frags.size() > 1) {
      std::sort(coinc_frags.begin(), coinc_frags.end());
      coincidence_v->emplace_back(sbnd::Coincidence(it->first, coinc_frags));
      std::advance(it, frags_in_window);
    }
  } // Looop over timestamps
}

DEFINE_ART_MODULE(sbnd::CoincidenceProducer)
