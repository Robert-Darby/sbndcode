////////////////////////////////////////////////////////////////////////
// Class:       MichelWaveformProducer
// Plugin Type: producer (Unknown Unknown)
// File:        MichelWaveformProducer_module.cc
//
// Generated at Wed Mar 12 06:21:32 2025 by Robert Darby using cetskelgen
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

#include "sbndcode/MichelRecoUtil/MichelRecoUtil.hh"
#include "canvas/Persistency/Common/Assns.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "sbnobj/SBND/Trigger/MichelTag.hh"
#include "art_root_io/TFileService.h"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

#include <memory>

namespace sbnd
{
  class MichelWaveformProducer;
}

class sbnd::MichelWaveformProducer : public art::EDProducer
{
public:
  explicit MichelWaveformProducer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelWaveformProducer(MichelWaveformProducer const &) = delete;
  MichelWaveformProducer(MichelWaveformProducer &&) = delete;
  MichelWaveformProducer &operator=(MichelWaveformProducer const &) = delete;
  MichelWaveformProducer &operator=(MichelWaveformProducer &&) = delete;

  // Required functions.
  void produce(art::Event &e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  int fEventID;
  int fChannel;
  float muonThreshold, michelThreshold;
  float wvfStartTime, wvfEndTime;
  float fMuonPeakTime, fMichelPeakTime, fMuonPeakAmp, fMichelPeakAmp;
  int fTriggerG4ID;
  bool fHasMichelPeak;
  size_t nSamples;

  // Constants. These shouldn't need to be changed in the fhicl
  const std::vector<int> fPair1, fPair2, fUnpaired;
  const std::vector<std::string> fPDs;
  const float fDaphneFreq;
  const float fCAENFreq;
  const float fPMTReadoutDelay, fARAReadoutDelay;

  // For finding peaks in summed waveform
  const std::vector<std::string> fOpFlashLabels;
  const std::string fFinderInputLabel;
  const std::vector<std::string> fFinderOpTypes;
  const std::string fPMTTriggerLabel;
  const int fMinPMTMultiplicity;
  const float fFinderPolarity;
  const float fFindMuonThreshold, fFindMichelPeak;
  const float fMuonMichelMaxRatio;
  const std::string fCRTLabel;
  const sbnd::crt::CRTGeoAlg fCRTGeo;
  const float fCRTOffset, fCRTClockSpeed;
  const float fCRTCoincidence;
  const int fRunningAvgSampleWidth;

  // For finding peaks in individual waveforms
  const bool fProduceWaveforms;
  const art::InputTag fInputLabel;
  const art::InputTag fMichelTagLabel;
  const std::vector<std::string> fInputOpTypes;
  const std::vector<float> fInputPolarity;
  const std::vector<float> fMuonThresholds, fMichelThresholds;
  const float fPeakCoincWindow;
  const float fTimeBinWidth;
  const float fStartTimeAfterMuon;
  const float fEndTimeAfterMuon;

  // MC Info
  const bool fUseMC;
  const std::string fMCTruthLabel;
  const std::string fSEDLabel;
  const std::string fSEDOutLabel;
  const bool fOnlyStopping;
  std::string fAuxDetHitLabel;
  std::vector<std::string> fAuxDetHitTags;

  // Misc
  const bool fRequireMichel;
  const float fMinLifetime;
  const float fMaxLifetime;
  const bool fSaveHists;
  const bool fMakeTree;
  const std::vector<int> fSaveEvents;
  const bool fBinWaveform;
  const bool fVerbose;

  const opdet::sbndPDMapAlg fPDMap;

  void findPeakPairs(
      const std::vector<std::pair<size_t, float>> &peakIndices,
      const std::vector<MichelTag> mtags,
      std::vector<std::pair<size_t, size_t>> &peakPairs,
      int channel);

  // Declare member data here.
};

sbnd::MichelWaveformProducer::MichelWaveformProducer(fhicl::ParameterSet const &p)
    : EDProducer{p},
      // More initializers here.
      fPair1(p.get<std::vector<int>>("Pair1")),
      fPair2(p.get<std::vector<int>>("Pair2")),
      fUnpaired(p.get<std::vector<int>>("Unpaired")),
      fPDs(p.get<std::vector<std::string>>("PDs")),
      fDaphneFreq(p.get<float>("DaphneFreq")),
      fCAENFreq(p.get<float>("CAENFreq")),
      fPMTReadoutDelay(p.get<float>("PMTReadoutDelay", 0.135)),
      fARAReadoutDelay(p.get<float>("XARAPUCAReadoutDelay", 0.0)),

      // For finding peaks in summed waveform
      fOpFlashLabels(p.get<std::vector<std::string>>("OpFlashLabels")),
      fFinderInputLabel(p.get<std::string>("FinderInputLabel")),
      fFinderOpTypes(p.get<std::vector<std::string>>("FinderOpTypes")),
      fPMTTriggerLabel(p.get<std::string>("PMTTriggerLabel")),
      fMinPMTMultiplicity(p.get<int>("MinPMTMultiplicity", 10)),
      fFinderPolarity(p.get<float>("FinderPolarity", 1.)),
      fFindMuonThreshold(p.get<float>("FindMuonThreshold", 5000.)),
      fFindMichelPeak(p.get<float>("FindMichelPeak", 1000.)),
      fMuonMichelMaxRatio(p.get<float>("MuonMichelMaxRatio", 1.)),
      fCRTLabel(p.get<std::string>("CRTLabel")),
      fCRTGeo(p.get<fhicl::ParameterSet>("CRTGeoParams")),
      fCRTOffset(p.get<float>("CRTOffset", 1700000.)),
      fCRTClockSpeed(p.get<float>("CRTClockspeed", 1.)),
      fCRTCoincidence(p.get<float>("CRTCoincidence")),
      fRunningAvgSampleWidth(p.get<int>("RunningSumSampleWidth", 20)),

      // For finding peaks in individual waveforms
      fProduceWaveforms(p.get<bool>("ProduceWaveforms", false)),
      fInputLabel(p.get<art::InputTag>("InputLabel", "opdecopmt")),
      fMichelTagLabel(p.get<art::InputTag>("MichelTagLabel")),
      fInputOpTypes(p.get<std::vector<std::string>>("InputOpTypes")),
      fInputPolarity(p.get<std::vector<float>>("InputPolarity")),
      fMuonThresholds(p.get<std::vector<float>>("MuonThresholds")),
      fMichelThresholds(p.get<std::vector<float>>("MichelThresholds")),
      fPeakCoincWindow(p.get<float>("PeakCoincWindow", 80.)),
      fTimeBinWidth(p.get<float>("TimeBinWidth")),
      fStartTimeAfterMuon(p.get<float>("StartTimeAfterMuon")),
      fEndTimeAfterMuon(p.get<float>("EndTimeAfterMuon")),

      // MC Info
      fUseMC(p.get<bool>("UseMC", false)),
      fMCTruthLabel(p.get<std::string>("MCTruthLabel")),
      fSEDLabel(p.get<std::string>("SEDLabel", "ionandscint")),
      fSEDOutLabel(p.get<std::string>("SEDOutLabel", "ionandscintout")),
      fOnlyStopping(p.get<bool>("OnlyStopping", true)),
      fAuxDetHitLabel(p.get<std::string>("AuxDetHitLabel")),
      fAuxDetHitTags(p.get<std::vector<std::string>>("AuxDetHitTags")),

      // Misc
      fRequireMichel(p.get<bool>("RequireMichel", true)),
      fMinLifetime(p.get<float>("MinLifetime", 2000.0f)),
      fMaxLifetime(p.get<float>("MaxLifetime", 10000.)),
      fSaveHists(p.get<bool>("SaveHists", false)),
      fMakeTree(p.get<bool>("MakeTree", false)),
      fSaveEvents(p.get<std::vector<int>>("SaveEvents")),
      fBinWaveform(p.get<bool>("BinWaveforms", false)),
      fVerbose(p.get<bool>("Verbose", false))

{
  produces<std::vector<raw::OpDetWaveform>>();
}

void sbnd::MichelWaveformProducer::produce(art::Event &e)
{
  std::unique_ptr<std::vector<raw::OpDetWaveform>>
      michelwvfms_v(std::make_unique<std::vector<raw::OpDetWaveform>>());

  fEventID = e.id().event();

  auto wvfmHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(fInputLabel);
  auto mtagHandle = e.getHandle<std::vector<sbnd::MichelTag>>(fMichelTagLabel);

  for (const auto &waveform : *wvfmHandle)
  {
    fChannel = waveform.ChannelNumber();
    auto opType = fPDMap.pdType(fChannel);
    if (std::find(fInputOpTypes.begin(), fInputOpTypes.end(), opType) == fInputOpTypes.end())
      continue;
    auto conversionFactor = (fPDMap.electronicsType(fChannel) == "daphne") ? fDaphneFreq : fCAENFreq;
    unsigned opTypeIdx = std::distance(fPDs.begin(), std::find(fPDs.begin(), fPDs.end(), opType));
    auto polarity = fInputPolarity[opTypeIdx];
    muonThreshold = fMuonThresholds.at(opTypeIdx);
    michelThreshold = fMichelThresholds.at(opTypeIdx);

    const std::vector<short> rawWaveformData = waveform.Waveform();
    std::vector<float> waveformData(rawWaveformData.size(), 0.);
    std::transform(
        rawWaveformData.begin(), rawWaveformData.end(), waveformData.begin(),
        [polarity](short adc)
        { return polarity * (float)adc; });
    auto correctedWaveform = MichelRecoUtil::subtractBaseline(waveformData, conversionFactor);
    // Apply rolling sum to the waveform
    auto smooth_wvf = MichelRecoUtil::CalcRunningAvg(fRunningAvgSampleWidth, correctedWaveform);
    auto rollingSumWaveform = MichelRecoUtil::applyRollingSum(correctedWaveform);

    nSamples = waveformData.size();

    // Find if a muon and michel time is within the waveform, skip if not
    wvfEndTime = waveform.TimeStamp() + nSamples * conversionFactor / 1000.;
    wvfStartTime = waveform.TimeStamp();

    std::vector<sbnd::MichelTag> mtagsInWaveform;
    for (const auto &mtag : *mtagHandle)
    {
      if (mtag.MuonTime > wvfStartTime && mtag.MuonTime < wvfEndTime)
      {
        mtagsInWaveform.push_back(mtag);
      }
    }
    if (mtagsInWaveform.empty())
      continue;

    // Subtract baseline

    // Find peaks
    std::vector<std::pair<size_t, float>> peakIndices;
    MichelRecoUtil::findPeaks(rollingSumWaveform, correctedWaveform, michelThreshold, peakIndices);

    // Find pairs of peaks
    std::vector<std::pair<size_t, size_t>> peakPairs;
    findPeakPairs(peakIndices, mtagsInWaveform, peakPairs, fChannel);

    if (peakPairs.empty())
      continue;

    // Process each pair
    fMuonPeakTime = -9999.;
    fMichelPeakTime = -9999.;
    for (const auto &[muonPeakIndex, michelPeakIndex] : peakPairs)
    {
      fHasMichelPeak = (michelPeakIndex < nSamples);
      fMuonPeakTime = ((float)muonPeakIndex * conversionFactor / 1000.) + (float)wvfStartTime;
      fMichelPeakTime = (fHasMichelPeak) ? ((float)michelPeakIndex * conversionFactor / 1000.) + (float)wvfStartTime : -9999.;
      fMuonPeakAmp = correctedWaveform[muonPeakIndex];
      fMichelPeakAmp = correctedWaveform[michelPeakIndex];

      if (fHasMichelPeak)
      {
        unsigned save_start_bin, save_end_bin;
        float save_start_time;
        if (fMichelPeakTime - fMuonPeakTime > 0.1)
        {
          save_start_bin = michelPeakIndex - (unsigned)(100 / conversionFactor);
        }
        else
        {
          auto min_part = std::min_element(correctedWaveform.begin() + muonPeakIndex, correctedWaveform.begin() + michelPeakIndex);
          save_start_bin = std::distance(correctedWaveform.begin(), min_part);
        }
        save_end_bin = (fMichelPeakTime + 2. > wvfEndTime) ? nSamples - 1 : michelPeakIndex + (unsigned)(2000 / conversionFactor);
        save_start_time = ((float)save_start_bin * conversionFactor / 1000.) + (float)wvfStartTime;
        std::vector<short unsigned> save_wvf;
        for (unsigned i_adc = save_start_bin; i_adc < save_end_bin; i_adc++)
          save_wvf.push_back((short unsigned)correctedWaveform[i_adc]);
        michelwvfms_v->push_back(raw::OpDetWaveform(save_start_time, fChannel, save_wvf));
      } // If michel peak

      if (fTriggerG4ID > 0 && fSaveHists && !peakPairs.empty() &&
          std::find(fSaveEvents.begin(), fSaveEvents.end(), fEventID) != fSaveEvents.end() &&
          fMichelPeakTime - fMuonPeakTime > 0.1)
      {
        MichelRecoUtil::ConvertWaveformToHistogram(waveform, fMuonPeakTime, fMichelPeakTime, conversionFactor, fEventID);
      } // Save histograms
    } // Loop over peak pairs
  } // Loop over waveforms

  e.put(std::move(michelwvfms_v));
}

void sbnd::MichelWaveformProducer::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::MichelWaveformProducer::findPeakPairs(
    const std::vector<std::pair<size_t, float>> &peakIndices,
    const std::vector<MichelTag> mtags,
    std::vector<std::pair<size_t, size_t>> &peakPairs,
    int channel = 6)
{
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  std::vector<std::pair<size_t, size_t>> muon_windows, michel_windows;
  size_t mu_coincidence_window, michel_coincidence_window;
  std::vector<float> muon_times, michel_times;

  mu_coincidence_window = (size_t)(fPeakCoincWindow / conversionFactor);
  michel_coincidence_window = mu_coincidence_window;
  for (unsigned i = 0; i < mtags.size(); i++)
  {
    muon_times.push_back((mtags[i].MuonTime - wvfStartTime) * 1000.);
    michel_times.push_back((mtags[i].MichelTime - wvfStartTime) * 1000.);
  }
  for (unsigned i_mu = 0; i_mu < muon_times.size(); i_mu++)
  {
    size_t mu_bin = (size_t)(muon_times[i_mu] / conversionFactor);
    size_t mu_win_start = (mu_bin < mu_coincidence_window) ? 0 : mu_bin - mu_coincidence_window;
    size_t mu_win_end = (mu_bin + mu_coincidence_window > nSamples) ? nSamples - 1 : mu_bin + mu_coincidence_window;
    muon_windows.push_back(std::make_pair(mu_win_start, mu_win_end));
    size_t michel_bin = (size_t)(michel_times[i_mu] / conversionFactor);
    size_t michel_win_start = (michel_bin < michel_coincidence_window) ? 0 : michel_bin - michel_coincidence_window;
    size_t michel_win_end = (michel_bin + michel_coincidence_window > nSamples) ? nSamples - 1 : michel_bin + michel_coincidence_window;
    michel_windows.push_back(std::make_pair(michel_win_start, michel_win_end));
  }

  for (unsigned i_mu = 0; i_mu < muon_windows.size(); i_mu++)
  {
    float mu_peak = -std::numeric_limits<float>::max();
    float michel_peak = -std::numeric_limits<float>::max();
    float mu_peak2 = -std::numeric_limits<float>::max();
    float michel_peak2 = -std::numeric_limits<float>::max();
    size_t mu_peak_idx = -1;
    size_t michel_peak_idx = -1;
    size_t michel_peak2_idx = -1;
    for (unsigned i_peak = 0; i_peak < peakIndices.size(); i_peak++)
    {
      if (peakIndices[i_peak].first > michel_windows[i_mu].second)
        break;

      if (peakIndices[i_peak].first < muon_windows[i_mu].first)
        continue;

      // If new tallest peak in muon window
      if (peakIndices[i_peak].first > muon_windows[i_mu].first && peakIndices[i_peak].first < muon_windows[i_mu].second &&
          peakIndices[i_peak].second > mu_peak2 && peakIndices[i_peak].second > muonThreshold)
      {
        if (peakIndices[i_peak].second > mu_peak)
        {
          mu_peak2 = mu_peak;
          mu_peak = peakIndices[i_peak].second;
          mu_peak_idx = peakIndices[i_peak].first;
        }
        else
        {
          mu_peak2 = peakIndices[i_peak].second;
        } // If one of the two biggest peaks in muon search window
      }
      // If new tallest peak in michel window
      if (peakIndices[i_peak].first > michel_windows[i_mu].first && peakIndices[i_peak].first < michel_windows[i_mu].second &&
          peakIndices[i_peak].second > michel_peak2 && peakIndices[i_peak].second > michelThreshold &&
          (peakIndices[i_peak].second < fMuonMichelMaxRatio * mu_peak))
      {
        if (peakIndices[i_peak].second > michel_peak)
        {
          michel_peak2 = michel_peak;
          michel_peak2_idx = michel_peak_idx;
          michel_peak = peakIndices[i_peak].second;
          michel_peak_idx = peakIndices[i_peak].first;
        }
        else
        {
          michel_peak2 = peakIndices[i_peak].second;
          michel_peak2_idx = peakIndices[i_peak].first;
        }
      }
    } // Loop to find michelon peak
    if (mu_peak_idx == michel_peak_idx)
    {
      michel_peak_idx = michel_peak2_idx;
      michel_peak = michel_peak2;
    }
    if (mu_peak > muonThreshold &&
        ((fRequireMichel && michel_peak > michelThreshold) || !fRequireMichel))
    {
      if (mu_peak_idx > michel_peak_idx && michel_peak > 0)
      {
        auto temp = mu_peak_idx;
        mu_peak_idx = michel_peak_idx;
        michel_peak_idx = temp;
      }
      peakPairs.push_back(std::make_pair(mu_peak_idx, michel_peak_idx));
    }
  } // Loop over candidate times
} // findPeakPairs

DEFINE_ART_MODULE(sbnd::MichelWaveformProducer)
