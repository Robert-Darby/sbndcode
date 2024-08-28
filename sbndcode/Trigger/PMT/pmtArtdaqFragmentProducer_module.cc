////////////////////////////////////////////////////////////////////////
// Class:       pmtArtdaqFragmentProducer
// Plugin Type: producer (Unknown Unknown)
// File:        pmtArtdaqFragmentProducer_module.cc
//
// Generated at Mon Feb  7 14:18:34 2022 by Patrick Green using cetskelgen
// from  version .

// Module to convert simulated PMT waveforms into artdaq fragment format
// Input: PMT waveforms and output from hardware trigger producer (pmtTriggerProducer_module.cc)
// For each hardware trigger, creates and saves CAEN1730 art-daq fragments for each PMT waveform. 
// Time range: -1us to +9us from hardware trigger time. 8 fragments per trigger, each containing 15 PMT waveforms (channels 0-14)
// Sets trigger time to 0.5 seconds +- PMT trigger time; and adds step function for beam time (relative to trigger) (channel 15)
// Output: adds std::vector<artdaq::Fragment> containing the fragments from each event
// To do -- change timerange offset and trigger time offset to be fhicl parameters + any other hardcoded pieces  
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"

// artdaq includes
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"

#include"sbndaq-artdaq-core/Obj/SBND/Coincidence.hh"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// SBN/SBND includes
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

// ROOT includes
#include "TFile.h"
#include "TH1D.h"

// random numbers
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <ctime>

namespace sbnd {
  namespace trigger {
    class pmtArtdaqFragmentProducer;
  }
}


class sbnd::trigger::pmtArtdaqFragmentProducer : public art::EDProducer {
public:
  explicit pmtArtdaqFragmentProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pmtArtdaqFragmentProducer(pmtArtdaqFragmentProducer const&) = delete;
  pmtArtdaqFragmentProducer(pmtArtdaqFragmentProducer&&) = delete;
  pmtArtdaqFragmentProducer& operator=(pmtArtdaqFragmentProducer const&) = delete;
  pmtArtdaqFragmentProducer& operator=(pmtArtdaqFragmentProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void getCRTTimeStamps(
    art::Event& e,
    std::map<int, unsigned>& coincTimeStamps);
  double fMinStartTime;
  double fMaxEndTime;


  // Declare member data here.

  // fhicl parameters
  std::string fInputModuleNameWvfm;
  std::string fWaveformInstanceName;
  std::string fInputModuleNameTrigger;
  int fBaseline; // baseline in simulation, default 8000 ADC (for expanding waveforms only, when not fully simulated)
  int fMultiplicityThreshold; // number of PMT pairs in hardware trigger to pass
  double fBeamWindowStart;
  double fBeamWindowLength;
  uint32_t nChannelsFrag;
  uint32_t wfm_length; // ~10us, 2ns tick 
  int fCoincidenceWindow;
  float fCRTTriggerOffset;
  float fCRTClockFreq;
  std::string fCRTFragLabel;
  std::vector<std::vector<unsigned>> fCRTCoincidence; 
  bool fVerbose;

  // event information
  int fRun, fSubrun;
  art::EventNumber_t fEvent;

  // PD information
  opdet::sbndPDMapAlg pdMap; // photon detector map
  std::vector<unsigned int> channelList; 

  // waveforms
  std::vector<std::vector<short>> wvf_channel;
   
  // sampling rate
  double fSampling;

  // services
  art::ServiceHandle<art::TFileService> tfs;

  // random numbers
  CLHEP::HepRandomEngine&       fTriggerTimeEngine;

};


sbnd::trigger::pmtArtdaqFragmentProducer::pmtArtdaqFragmentProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fInputModuleNameWvfm(p.get<std::string>("InputModuleNameWvfm")),
  fWaveformInstanceName(p.get<std::string>("WaveformInstanceName", "PMTChannels")),
  fInputModuleNameTrigger(p.get<std::string>("InputModuleNameTrigger")),
  fBaseline(p.get<int>("Baseline",8000)),
  fMultiplicityThreshold(p.get<int>("MultiplicityThreshold")),
  fBeamWindowStart(p.get<double>("BeamWindowStart", 1510.)),
  fBeamWindowLength(p.get<double>("BeamWindowLength", 1.6)),
  nChannelsFrag(p.get<double>("nChannelsFrag", 15)),
  wfm_length(p.get<double>("WfmLength", 5120)),
  fCoincidenceWindow(p.get<int>("CoincidenceWindow", 1000)),
  fCRTTriggerOffset(p.get<float>("CRTTriggerOffset", 1700000)),
  fCRTClockFreq(p.get<float>("CRTClockFreq", 1000.)),
  fCRTFragLabel(p.get<std::string>("CRTFragLabel")),
  fCRTCoincidence(p.get<std::vector<std::vector<unsigned>>>("CRTCoincidence")),
  fVerbose(p.get<bool>("Verbose", false)),
    fTriggerTimeEngine(art::ServiceHandle<rndm::NuRandomService>{}->registerAndSeedEngine(
                         createEngine(0, "HepJamesRandom", "trigger"), "HepJamesRandom", "trigger", p, "SeedTriggerTime"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<sbnd::trigger::Coincidence> >(); 
  produces< std::vector<artdaq::Fragment> >();

  // get clock
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSampling = clockData.OpticalClock().Frequency(); // MHz

  // build PD map and channel list

  auto subsetCondition = [](auto const& i)->bool { return i["pd_type"] == "pmt_coated" || i["pd_type"] == "pmt_uncoated"; };
  auto pmtMap = pdMap.getCollectionFromCondition(subsetCondition);
  if (fVerbose) std::cout << "Number of PDs selected: \t" << pmtMap.size() << "\n";
  for(auto const& i:pmtMap){
    channelList.push_back(i["channel"]);
  }
   
}

void sbnd::trigger::pmtArtdaqFragmentProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  // event information
  fRun = e.run();
  fSubrun = e.subRun();
  fEvent = e.id().event();

  if (fVerbose) std::cout << "Processing Run: " << fRun << ", Subrun: " << fSubrun << ", Event: " << fEvent << std::endl;

  // access PMT waveforms and hardware trigger information
  art::Handle< std::vector< raw::OpDetWaveform > > wvfmHandle;
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

  wvfmHandle = waveformHandles[0];

  art::Handle< std::vector< sbnd::comm::pmtTrigger > > triggerHandle;
  e.getByLabel(fInputModuleNameTrigger, triggerHandle);

  if(!wvfmHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any waveforms, input must contain OpDetWaveforms." << "\n";
  }
  if(!triggerHandle.isValid()) {
    throw art::Exception(art::errors::Configuration)
          << "Could not find any PMT hardware trigger object, must run PMT trigger producer before running this module." << "\n";
  }

  // create empty vectors to hold waveforms for each channel
  fMinStartTime = std::numeric_limits<double>::max();//in us
  fMaxEndTime = -std::numeric_limits<double>::max();//in us

  for(auto const& wvf : (*wvfmHandle)) {
    double fChNumber = wvf.ChannelNumber();
    std::string opdetType = pdMap.pdType(fChNumber);
    // only look at pmts
    if (opdetType != "pmt_coated" && opdetType != "pmt_uncoated") continue;
      if (wvf.TimeStamp() < fMinStartTime){ fMinStartTime = wvf.TimeStamp(); }
      auto end_time = (double)wvf.size() / fSampling + wvf.TimeStamp();
      if ((double)end_time > fMaxEndTime){ fMaxEndTime = (double)end_time;}
  }

  if (fVerbose){std::cout<<"MinStartTime: "<<fMinStartTime<<" MaxEndTime: "<<fMaxEndTime<<std::endl;}

  std::vector<short> wvf_0; wvf_0.reserve((int)(3020*1e6/2)); 
  for (double i = fMinStartTime; i<fMaxEndTime+(1./fSampling); i+=(1./fSampling)){
    wvf_0.push_back(fBaseline);
  }
  for (size_t i = 0; i < channelList.size(); i++){
    wvf_channel.push_back(wvf_0);
  }

  // counters
  int num_pmt_wvf = 0;
  
  // loop through waveform handles
  size_t wvf_id = -1;
  size_t hist_id = -1;
  for(auto const& wvf : (*wvfmHandle)) {
    wvf_id++;
    double fChNumber = wvf.ChannelNumber();
    std::string opdetType = pdMap.pdType(fChNumber);
    
    // only look at pmts
    if (opdetType != "pmt_coated" && opdetType != "pmt_uncoated") continue;
    
    num_pmt_wvf++;

    double fStartTime = wvf.TimeStamp(); // in us
    double fEndTime = double(wvf.size()) / fSampling + fStartTime; // in us

    // create full waveform
    std::vector<short> wvf_full; wvf_full.reserve((short)(3020*1e6/2));

    if (fStartTime > fMinStartTime){
      for (double i = fStartTime-fMinStartTime; i>0.; i-=(1./fSampling)){
        wvf_full.push_back(fBaseline);
      }
    }

    for(unsigned int i = 0; i < wvf.size(); i++) {
      wvf_full.push_back(wvf[i]);
    }

    if (fEndTime < fMaxEndTime){
      for (double i = fMaxEndTime-fEndTime; i>0.; i-=(1./fSampling)){
        wvf_full.push_back(fBaseline);
      }
    }

    // combine waveform with any other waveforms from same channel
    int i_ch = -1.;
    auto ich = std::find(channelList.begin(), channelList.end(), fChNumber);
    if (ich != channelList.end()){
      i_ch = ich - channelList.begin();
    }
    if (wvf_channel.at(i_ch).size() < wvf_full.size()){
      if (fVerbose) std::cout<<"Full waveform -- Previous Channel" << fChNumber <<" Size: "<<wvf_channel.at(i_ch).size()<<"New Channel" << fChNumber <<" Size: "<<wvf_full.size()<<std::endl;
      for(unsigned int i = wvf_channel.at(i_ch).size(); i < wvf_full.size(); i++) {
        wvf_channel.at(i_ch).push_back(fBaseline);
      }
    }
    for(unsigned int i = 0; i < wvf_full.size(); i++) {
       wvf_channel.at(i_ch)[i] += (wvf_full[i] - fBaseline);
    } 

    hist_id++;

    wvfmHandle.clear();
  } // waveform handle loop

  std::map<int, unsigned> CRTTimes;
  getCRTTimeStamps(e, CRTTimes);

  // access hardware trigger information
  std::vector<size_t> triggerIndex; 
  for(auto const& trigger : (*triggerHandle)) {
    for (size_t idx = 0; idx < trigger.numPassed.size(); idx++) {
      if (trigger.numPassed[idx] >= fMultiplicityThreshold) {
        // save index
        triggerIndex.push_back(idx);
        // skip ahead by 5120 samples (downsampled by 4) before checking again
        idx += 1280;
      }
    }
  } // trigger handle loop

  // Only make fragments with CRT coincident triggers
  for(auto i_pmt : triggerIndex) {CRTTimes.insert(std::pair<int, unsigned>(i_pmt, 7));}
  for(auto coinc : CRTTimes) std::cout << coinc.second << " " << coinc.first << "\n";
  std::map<int, unsigned> coincTimes;
  unsigned id = 0;
  std::unique_ptr<std::vector<sbnd::trigger::Coincidence>> coincidence_v = std::make_unique<std::vector<sbnd::trigger::Coincidence>>();
  for(auto it = CRTTimes.begin(); it != CRTTimes.end(); it++) {
    for(auto jt = it; jt!=CRTTimes.end(); jt++) {
      if(jt==it) continue;
      if(jt->first > it->first + (fCoincidenceWindow / 8)) break;
      if((it->second == 7 && jt->second != 7)) {
        coincTimes.insert(std::pair<int, unsigned>(it->first, jt->second*10 + it->second));
        coincidence_v->push_back(sbnd::trigger::Coincidence(id,
                                                        jt->second*10 + it->second,
							(float)it->first / 125. -1510.,
                                                        (float)it->first / 125. -1510.));
        id++;
      } else if(it->second != 7 && jt->second == 7) {
        coincTimes.insert(std::pair<int, unsigned>(jt->first, it->second*10 + jt->second));
        coincidence_v->push_back(sbnd::trigger::Coincidence(id,
                                                        it->second*10 + jt->second,
                                                        (float)it->first / 125. -1510.,
                                                        (float)jt->first / 125. -1510.));
        id++;
      }
    }
  }

  e.put(std::move(coincidence_v));
  if (fVerbose) std::cout << "Number of PMT hardware triggers found: " << triggerIndex.size() << std::endl; 
  
  // fragments vector
  std::unique_ptr<std::vector<artdaq::Fragment>> vecFrag = std::make_unique<std::vector<artdaq::Fragment>>();

  // set properties of fragment that are common to event
  uint32_t nChannelsTotal = channelList.size();
  //uint32_t nChannelsFrag = 15;
  uint32_t nFrag = (uint32_t)nChannelsTotal/nChannelsFrag;
  //uint32_t wfm_length = 5120; // ~10us, 2ns tick

  // fragment properties
  uint32_t sequenceIDVal = fEvent;

  // create and populate common metadata
  sbndaq::CAENV1730FragmentMetadata metadata;
  metadata.nChannels = nChannelsFrag + 1; // 15 PMT channels + final channel to store beam window / trigger information
  metadata.nSamples = wfm_length;
  
  // fragment handle properties
  uint32_t eventCounterVal = fEvent;
  uint32_t boardIDVal = 0;
  uint32_t triggerTimeTagVal = (uint32_t)CLHEP::RandFlat::shoot(&fTriggerTimeEngine, 0, 1e9);
  uint32_t eventSizeVal = ((wfm_length * (nChannelsFrag+1)) * sizeof(uint16_t) + sizeof(sbndaq::CAENV1730EventHeader)) / sizeof(uint32_t);
  
  // loop over PMT hardware triggers
  for (auto coinc : coincTimes) {

    // index in full waveform, 2ns tick
    size_t trigIdx = coinc.first *4;
    // size_t startIdx = trigIdx-500); // -1us
    size_t startIdx = trigIdx-500;

    // determine and set timestamp for particular trigger
    // double triggerTime = fMinStartTime + wvfIdx*0.008; // in us
    double triggerTime = coinc.first *0.008; // in us
    double timestampVal = 0.5 + (triggerTime*1e-6); // in seconds // std::time(nullptr); // current time
    metadata.timeStampSec = (uint32_t)timestampVal;
    metadata.timeStampNSec = (uint32_t)(timestampVal*1e9);

    // create fragments to hold waveforms, set properties and populate
    // 15 PMTs stored per fragment, 120/15 = 8 fragments per trigger
    for (size_t counter = 0; counter < nFrag; counter++) {
      
      // create fragment
      const auto fragment_datasize_bytes = metadata.ExpectedDataSize();
      uint32_t fragmentIDVal = counter;
      auto fragment_uptr = artdaq::Fragment::FragmentBytes(fragment_datasize_bytes, sequenceIDVal, fragmentIDVal, sbndaq::detail::FragmentType::CAENV1730, metadata); // unique pointer
      fragment_uptr->setTimestamp(timestampVal);

      // populate fragment header
      auto header_ptr = reinterpret_cast<sbndaq::CAENV1730EventHeader*>(fragment_uptr->dataBeginBytes());
        
      header_ptr->eventCounter = eventCounterVal;
      header_ptr->boardID = boardIDVal;
      header_ptr->triggerTimeTag = triggerTimeTagVal;  // ns // set timetag as random value for event
      header_ptr->eventSize = eventSizeVal;
        
      // populate waveforms  
      // populate fragment with waveform
      uint16_t* data_begin = reinterpret_cast<uint16_t*>(fragment_uptr->dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
      uint16_t* value_ptr = data_begin;
      uint16_t value = 0;
      size_t ch_offset = 0;

      // loop over channels
      for (size_t i_ch = 0; i_ch < nChannelsFrag; i_ch++) {
        ch_offset = (size_t)(i_ch*wfm_length);
        // loop over waveform
        for (size_t i_t = 0; i_t < wfm_length; i_t++) {
          // set value
          value = wvf_channel[counter*nChannelsFrag + i_ch][startIdx+i_t];
          value_ptr = data_begin + ch_offset + i_t;
          *value_ptr = value;
        }
      }

      // create add beam window trigger waveform
      ch_offset = (size_t)(nChannelsFrag*wfm_length);
      size_t beamStartIdx = 0;
      size_t beamEndIdx = beamStartIdx + fBeamWindowLength*1000/2;
      // loop over waveform
      for (size_t i_t = 0; i_t < wfm_length; i_t++) {
        // set value
        if (startIdx + i_t >= beamStartIdx && startIdx + i_t <= beamEndIdx) value = 1;
        else value = 0;
        value_ptr = data_begin + ch_offset + i_t;
        *value_ptr = value;
      }

      // add fragment to vector
      vecFrag->push_back(*fragment_uptr);  
    }
  }

  if(fVerbose) std::cout << "Fragments written: " << vecFrag->size() << std::endl;

  // add fragments to event
  e.put(std::move(vecFrag));

  // clear variables
  wvf_channel.clear();
  wvf_channel.shrink_to_fit();
}

void sbnd::trigger::pmtArtdaqFragmentProducer::getCRTTimeStamps(
  art::Event& e,
  std::map<int, unsigned>& coincTimeStamps)
{
  std::map<int, unsigned> fragTimeStamps;

  // Access CRT Fragments
  art::Handle< std::vector<artdaq::Fragment> > crtHandle;
  std::vector< art::Ptr<artdaq::Fragment> > crtVect;
  if(e.getByLabel(fCRTFragLabel, crtHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(crtVect, crtHandle);

  // use  fragment ID to get plane information
  for(auto& frag : crtVect) {
    sbndaq::BernCRTFragmentV2 bern_fragment(*(frag));
    const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
    auto thisone = frag->fragmentID();  uint plane = (thisone & 0x0700) >> 8;
    if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}
    // Check if this fragment is considered in the coincidence check
    for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
      sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
      // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
      auto thisflag = bevt->flags;
      if (thisflag & 0x2 && !(thisflag & 0xC) ) {
        auto thistime = (((bevt->ts1 - fCRTTriggerOffset) / fCRTClockFreq) - fMinStartTime)*1000.;
        std::cout << (bevt->ts1 - fCRTTriggerOffset) / fCRTClockFreq + (fMinStartTime)*1000.<< "\n";
        if(thistime >= 0.) fragTimeStamps.insert(std::pair<int, unsigned>((int)(thistime / 8.), plane));
      }
    }
  }

  // Only select times with coincident CRTs from fhicl
  for(auto hit_it=fragTimeStamps.begin(); hit_it != fragTimeStamps.end(); hit_it++) {
    if(std::distance(fragTimeStamps.begin(), hit_it) >= (int)fragTimeStamps.size()) break;
    int n_hits=0;
    std::vector<unsigned> cur_coinc(7); cur_coinc[hit_it->second] = 1;
    for(auto hit_jt = hit_it; hit_jt != fragTimeStamps.end(); hit_jt++) {
      if(hit_jt->first > hit_it->first + fCoincidenceWindow) break;
      if(hit_it==hit_jt)continue;
      n_hits++;
      cur_coinc[hit_jt->second] = 1;
    }
    int coinc_int = 0;
    for(unsigned i=0; i<7; i++) {
      if(cur_coinc[i]==1) {coinc_int *= 10; coinc_int += i;}
    }
    for(auto& coinc_vec : fCRTCoincidence) {
      if(coinc_vec==cur_coinc) {
        coincTimeStamps.insert(std::pair<int, unsigned>(hit_it->first, coinc_int));
        break;
      }
    }
    std::advance(hit_it, n_hits);
  }

}

DEFINE_ART_MODULE(sbnd::trigger::pmtArtdaqFragmentProducer)
