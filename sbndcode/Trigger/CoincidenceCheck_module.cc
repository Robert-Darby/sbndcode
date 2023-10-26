////////////////////////////////////////////////////////////////////////
// Class:       CoincidenceCheck
// Plugin Type: filter (Unknown Unknown)
// File:        CoincidenceCheck_module.cc
//
// Generated at Thu Oct 19 15:06:07 2023 by Robert Darby using cetskelgen
// from  version .
//
// Thise module checks for coincidence between CRT/PMT trigger times
// within a given window. Can be used for CRT and PDS or CRT only.
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


#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"


#include <memory>

class CoincidenceCheck;


class CoincidenceCheck : public art::EDFilter {
public:
  explicit CoincidenceCheck(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CoincidenceCheck(CoincidenceCheck const&) = delete;
  CoincidenceCheck(CoincidenceCheck&&) = delete;
  CoincidenceCheck& operator=(CoincidenceCheck const&) = delete;
  CoincidenceCheck& operator=(CoincidenceCheck&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void getCRTTimeStamps(
    artdaq::Fragment& frag,
    std::multimap<double, size_t>& fragTimeStamps);
  void getCAEN1730FragmentTimeStamp(
    const artdaq::Fragment &frag,
    std::multimap<double, size_t>& fragTimeStamps);
  bool getCoincidence(
    const std::multimap<double, size_t>& fragTimeStamps);

  // FHICL Params
  const double fCoincidenceWindow;
  const std::vector<size_t> fCRTPlanes;
  const float fCRTTriggerOffset;
  const float fCRTClockFreq;
  const std::string fCRTLogic;
  const bool fCRTOnly;
  const double fCAENTriggerOffset;
  const bool fVerbose;

  unsigned fNCoincidentFrags;
  float fSamplingFreq;
};


CoincidenceCheck::CoincidenceCheck(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  , fCoincidenceWindow(p.get<double>("CoincidenceWindow", 0.24))
  , fCRTPlanes(p.get<std::vector<size_t>>("CRTPlanes"))
  , fCRTTriggerOffset(p.get<float>("CRTTriggerOffset", 1700000))
  , fCRTClockFreq(p.get<float>("CRTClockFreq", 1000.))
  , fCRTLogic(p.get<std::string>("CRTLogic", "AND"))
  , fCRTOnly(p.get<bool>("CRTOnly", false))
  , fCAENTriggerOffset(p.get<double>("CAENTriggerOffset", 0.5))
  , fVerbose(p.get<bool>("Verbose", false))
{
  fNCoincidentFrags = std::count(fCRTPlanes.begin(), fCRTPlanes.end(), 0) +
                      std::count(fCRTPlanes.begin(), fCRTPlanes.end(), 1);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSamplingFreq = clockData.OpticalClock().Frequency(); // MHz
}

bool CoincidenceCheck::filter(art::Event& e)
{
  std::multimap<double, size_t> fragTimeStamps;
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
        else if (!fCRTOnly && frag.type()==sbndaq::detail::FragmentType::CAENV1730) {
          getCAEN1730FragmentTimeStamp(frag, fragTimeStamps);
        } //if is pmt frag
      }
    } // Fragment loop
  } // Fragment handle loop 
  return getCoincidence(fragTimeStamps);
}

void CoincidenceCheck::beginJob()
{
  // Implementation of optional member function here.
}

void CoincidenceCheck::endJob()
{
  // Implementation of optional member function here.
}

void CoincidenceCheck::getCRTTimeStamps(
  artdaq::Fragment& frag,
  std::multimap<double, size_t>& fragTimeStamps)
{
  // use  fragment ID to get plane information
  sbndaq::BernCRTFragmentV2 bern_fragment(frag);
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  auto thisone = frag.fragmentID();  uint plane = (thisone & 0x0700) >> 8;
  if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}
  // Check if this fragment is considered in the coincidence check
  if(fCRTPlanes[plane] == 2) return;
  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
    // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
    auto thisflag = bevt->flags;
    if (thisflag & 0x2 && !(thisflag & 0xC) ) {
      auto thistime = (bevt->ts1 - fCRTTriggerOffset) / fCRTClockFreq;
      fragTimeStamps.insert(std::pair<float, size_t>(thistime, plane));
    }
  }
} // CoincidenceCheck::getCRTTimeStamps

void CoincidenceCheck::getCAEN1730FragmentTimeStamp(
  const artdaq::Fragment &frag,
  std::multimap<double, size_t>& fragTimeStamps)
{
  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp
  double timestamp = (double)md->timeStampNSec;
  timestamp -= fCAENTriggerOffset*1e9; timestamp /= 1000;
  fragTimeStamps.insert(std::pair<float, size_t>(timestamp, 8));
}

bool CoincidenceCheck::getCoincidence(
  const std::multimap<double, size_t>& fragTimeStamps)
{
  // Loop over chronological CRT/PDS tiggers
  for(auto it = fragTimeStamps.begin(); it != fragTimeStamps.end(); it++) {
    if(fVerbose) std::cout << "Plane: " << it->second << " Time: " << it->first << std::endl;
    unsigned count_coincident_frags = 0;
    std::vector<bool> frag_present(8, false);
    frag_present[it->second] = true;
    // Loop over all triggers in specified coincidence window
    float window_end = it->first + fCoincidenceWindow;
    for(auto jt = it; jt != fragTimeStamps.end(); jt++) {
      frag_present[jt->second] = true;
      if(jt->first > window_end) break;
    }
    if(!fCRTOnly && !frag_present[7]) continue;
    // Check for coincidence between specified CRT planes
    for(size_t fid = 0; fid < 8; fid++) {
      if((fCRTPlanes[fid] == 2) ||
         (fCRTPlanes[fid] != frag_present[fid])) continue;
      else {count_coincident_frags++;}
    }
    if(count_coincident_frags == fNCoincidentFrags) {
      std::cout << "Found coincidence at: " << it->first << std::endl;
      return true;
    }
  } // Looop over timestamps
  return false;
}

DEFINE_ART_MODULE(CoincidenceCheck)
