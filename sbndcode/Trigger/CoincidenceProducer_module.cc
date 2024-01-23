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
#include "art_root_io/TFileService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"

#include "sbnobj/SBND/Trigger/Coincidence.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include<map>
#include<vector>
#include <memory>

#include "TH1.h"

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
    const art::Event& e,
    const artdaq::Fragment &frag,
    std::multimap<double, unsigned>& fragTimeStamps,
    std::map<double, std::vector<std::vector<uint16_t>>>& fWvfmsMap);
  std::map<int, unsigned> getMultiplicities(
    std::map<double, std::vector<std::vector<uint16_t>>>& fWvfmsMap);
  void getCoincidence(
    const std::multimap<double, unsigned>& fragTimeStamps,
    std::map<int, unsigned>& multMap,
    std::unique_ptr<std::vector<sbnd::Coincidence>>& coincidence_v);

  // Tree params
  std::stringstream histname;

  // FHICL Params
  const double fCoincidenceWindow;
  const float fCRTTriggerOffset;
  const float fCRTClockFreq;
  const double fCAENTriggerOffset;
  const std::vector<float> fPMTThresholds;
  const unsigned fWvfmLength;
  const bool fMakeHists;
  const bool fVerbose;

  art::ServiceHandle<art::TFileService> tfs;
  opdet::sbndPDMapAlg pdMap; // photon detector map
};

sbnd::CoincidenceProducer::CoincidenceProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  , fCoincidenceWindow(p.get<double>("CoincidenceWindow", 0.24))
  , fCRTTriggerOffset(p.get<float>("CRTTriggerOffset", 1700000))
  , fCRTClockFreq(p.get<float>("CRTClockFreq", 1000.))
  , fCAENTriggerOffset(p.get<double>("CAENTriggerOffset", 0.5))
  , fPMTThresholds(p.get<std::vector<float>>("PMTThresholds"))
  , fWvfmLength(p.get<unsigned>("WvfmLength", 5120))
  , fMakeHists(p.get<bool>("MakeHists", false))
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
  std::map<double, std::vector<std::vector<uint16_t>>> fWvfmsMap;

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
          getCAEN1730FragmentTimeStamp(e, frag, fragTimeStamps, fWvfmsMap);
        } //if is pmt frag
      }
    } // Fragment loop
  } // Fragment handle loop
  auto multMap = getMultiplicities(fWvfmsMap);
  getCoincidence(fragTimeStamps, multMap, coincidence_v);
  if(fVerbose) {
    for(auto coinc = coincidence_v->begin(); coinc != coincidence_v->end(); coinc++) {
      std::cout << "Coincidence found at: " << coinc->TimeStamp << " Multiplicity: " << coinc->Multiplicity<< std::endl;
      for(auto& plane : coinc->Planes) std::cout << plane << " ";
      std::cout << std::endl;
    }
  }
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
  const art::Event& e,
  const artdaq::Fragment &frag,
  std::multimap<double, unsigned>& fragTimeStamps,
  std::map<double, std::vector<std::vector<uint16_t>>>& fWvfmsMap)
{
  // get fragment metadata
  sbndaq::CAENV1730Fragment bb(frag);
  auto const* md = bb.Metadata();

  // access timestamp
  double timestamp = (double)md->timeStampNSec;
  timestamp -= fCAENTriggerOffset*1e9; timestamp /= 1000; timestamp -= 1510.;
  fragTimeStamps.insert(std::pair<float, size_t>(timestamp, 7));
  if(timestamp > 0) return;

  // Add CEAN data to map
  // access fragment ID; index of fragment out of set of 8 fragments
  int fragId = static_cast<int>(frag.fragmentID());

  // access waveforms in fragment and save
  const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes()
                 + sizeof(sbndaq::CAENV1730EventHeader));
  const uint16_t* value_ptr =  data_begin;
  uint16_t value = 0;

  // channel offset
  size_t nChannels = 15; // 15 pmts per fragment
  size_t ch_offset = 0;

  // loop over channels
  if(fWvfmsMap.find(timestamp)==fWvfmsMap.end()) fWvfmsMap[timestamp].resize(120);
  for (size_t i_ch = 0; i_ch < nChannels; ++i_ch){
    fWvfmsMap[timestamp][i_ch + nChannels*fragId].resize(fWvfmLength);
    ch_offset = (size_t)(i_ch * fWvfmLength);
    //--loop over waveform samples
    histname.str(std::string());
    std::string prefix = (timestamp < 0) ? "m" : "";
    histname << "frag_" << e.id().event() << "_" << prefix << abs(int(timestamp)) << "_" << i_ch + nChannels*fragId;
    TH1D *fraghist = tfs->make<TH1D>(histname.str().c_str(), "Fragment", fWvfmLength, 0, fWvfmLength);
    for(size_t i_t = 0; i_t < fWvfmLength; ++i_t) {
      value_ptr = data_begin + ch_offset + i_t; // pointer arithmetic
      value = *(value_ptr);
      fWvfmsMap[timestamp][i_ch + nChannels*fragId][i_t] = value;
      fraghist->SetBinContent(i_t + 1, value);
    } //--end loop samples
  } //--end loop channels
}

std::map<int, unsigned> sbnd::CoincidenceProducer::getMultiplicities(
  std::map<double, std::vector<std::vector<uint16_t>>>& fWvfmsMap)
{
  std::map<int, unsigned> multMap; multMap[-9999] = 9999;
  for(auto mult_it = fWvfmsMap.begin(); mult_it != fWvfmsMap.end(); mult_it++) {
    std::vector<unsigned> multProf(fWvfmLength, 0);
    auto& wvfmsVec = mult_it->second;
    for(unsigned i_ch=0; i_ch<wvfmsVec.size(); i_ch++) {
      auto wvfm = wvfmsVec[i_ch];
      std::string opdetType = pdMap.pdType(i_ch);
      double threshold = (opdetType=="pmt_uncoated") ? fPMTThresholds[0] : fPMTThresholds[1];
      for(unsigned i = 0; i < fWvfmLength; i++) {
        auto value = wvfm[i];
        if(value < threshold) multProf[i]++;
      } // Loop over individual waveform
    } // Loop over all waveforms
    multMap.insert(std::pair<int, unsigned>(int(mult_it->first), *max_element(multProf.begin(), multProf.end())));
  } // Loop over fWvfmsMap
//  for(auto mult_it = multMap.begin(); mult_it != multMap.end(); mult_it++) {
//    std::cout << "Time: " << mult_it->first << " Multiplicity: " << mult_it->second << std::endl;
//  }
  return multMap;
}

void sbnd::CoincidenceProducer::getCoincidence(
  const std::multimap<double, unsigned>& fragTimeStamps,
  std::map<int, unsigned>& multMap,
  std::unique_ptr<std::vector<sbnd::Coincidence>>& coincidence_v)
{
  for(auto it=multMap.begin(); it!=multMap.end();it++)
    std::cout << "Time: " << it->first << " Multiplicity " << it->second << std::endl;
  // Loop over chronological CRT/PDS tiggers
  for(auto it = fragTimeStamps.begin(); it != fragTimeStamps.end(); it++) {
    int pmt_time = (it->second==7) ? int(it->first) : -9999;
    std::vector<unsigned> coinc_frags; coinc_frags.push_back(it->second);
    std::vector<bool> frag_present(8, false);
    frag_present[it->second] = true;
    // Loop over all triggers in specified coincidence window
    float window_end = it->first + fCoincidenceWindow;
    int frags_in_window = 0;
    for(auto jt = it; jt != fragTimeStamps.end(); jt++) {
      if(jt->first > window_end) break;
      frags_in_window++;
      if(jt->second == it->second) continue;
      if(jt->second==7) pmt_time = int(jt->first);
      coinc_frags.push_back(jt->second);
    }
    if(coinc_frags.size() > 1) {
      std::sort(coinc_frags.begin(), coinc_frags.end());
      unsigned multiplicity = multMap[pmt_time];
      std::cout << "Time: " << it->first << " PMT Time: " << pmt_time << " Multiplicity: " << multiplicity << std::endl;
      coincidence_v->emplace_back(sbnd::Coincidence(it->first, coinc_frags, multiplicity));
      std::advance(it, frags_in_window);
    }
  } // Looop over timestamps
}

DEFINE_ART_MODULE(sbnd::CoincidenceProducer)
