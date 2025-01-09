////////////////////////////////////////////////////////////////////////
// Class:       MaMOpHT
// Plugin Type: producer (art::EDProducer)
// File:        MaMOpHT_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <limits>

namespace sbnd {
  class MaMOpHT;
}

class sbnd::MaMOpHT : public art::EDProducer {
public:
  explicit MaMOpHT(fhicl::ParameterSet const& p);
  
  MaMOpHT(MaMOpHT const&) = delete;
  MaMOpHT(MaMOpHT&&) = delete;
  MaMOpHT& operator=(MaMOpHT const&) = delete;
  MaMOpHT& operator=(MaMOpHT&&) = delete;

  void produce(art::Event& e) override;

private:
  float CalculateBaselineAndSubtract(std::vector<short>& waveform) const;
  void ApplyRollingWindow(const std::vector<short>& waveform, std::map<int, short>& peaks) const;
  void SelectPeaks(const std::map<int, short>& peaks, std::vector<std::pair<float, short>>& selectedPeaks) const;

  const art::InputTag fInputTag;
  const short fThreshold;
  const size_t fMinPeakWidth;
  const size_t fMaxPeakSeparation;
  const float fPeakAmplitudeRatio;

  float startTime, endTime;
};

sbnd::MaMOpHT::MaMOpHT(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fInputTag(p.get<std::string>("InputTag")),
    fThreshold(p.get<short>("Threshold", 1000)),
    fMinPeakWidth(p.get<size_t>("MinPeakWidth", 5)),
    fMaxPeakSeparation(p.get<size_t>("MaxPeakSeparation", 4500)),
    fPeakAmplitudeRatio(p.get<float>("PeakAmplitudeRatio", 0.01))
{
  produces<std::vector<std::pair<float, short>>>();
}

void sbnd::MaMOpHT::produce(art::Event& e) {
  auto waveforms = e.getValidHandle<std::vector<raw::OpDetWaveform>>(fInputTag);
  startTime = std::numeric_limits<float>::max();
  endTime = std::numeric_limits<float>::lowest();

  // First loop: Determine the bounds of the cumulative waveform
  for (auto const& waveform : *waveforms) {
    float waveformStart = waveform.TimeStamp();
    float waveformEnd = waveformStart + waveform.Waveform().size() / 500.0;  // Example frequency of 500 Hz

    startTime = std::min(startTime, waveformStart);
    endTime = std::max(endTime, waveformEnd);
  }

  // Initialize the cumulative waveform with appropriate size based on the calculated bounds
  size_t waveformSize = static_cast<size_t>((endTime - startTime) * 500);  // Example frequency of 500 Hz
  std::vector<short> cumulativeWaveform(waveformSize, 0);

  // Second loop: Process and accumulate the waveforms
  for (auto const& waveform : *waveforms) {
    std::vector<short> wf = waveform.Waveform();
    float baseline = CalculateBaselineAndSubtract(wf);
    unsigned int startBin = static_cast<unsigned int>((waveform.TimeStamp() - startTime) * 500);

    for (size_t i = 0; i < wf.size(); ++i) {
      if (startBin + i < cumulativeWaveform.size()) {
        cumulativeWaveform[startBin + i] += (baseline - wf[i]);
      }
    }
  }

  std::map<int, short> peaks;
  ApplyRollingWindow(cumulativeWaveform, peaks);
  std::vector<std::pair<float, short>> selectedPeaks;
  SelectPeaks(peaks, selectedPeaks);

  e.put(std::make_unique<std::vector<std::pair<float, short>>>(selectedPeaks));

  // Log the start and end time of the cumulative waveform
  mf::LogInfo("WaveformTiming") << "Cumulative Waveform Start Time: " << startTime
                                << ", End Time: " << endTime;
}

float sbnd::MaMOpHT::CalculateBaselineAndSubtract(std::vector<short>& waveform) const {
  float baseline = std::accumulate(waveform.begin(), waveform.begin() + 200, 0.0) / 200;
  for (auto& bin : waveform) {
    bin -= static_cast<short>(baseline);
  }
  return baseline;
}

void sbnd::MaMOpHT::ApplyRollingWindow(const std::vector<short>& waveform, std::map<int, short>& peaks) const {
  int rollingSum = 0;
  int peakStart = -1;

  for (size_t i = 1; i < waveform.size(); ++i) {
    if (waveform[i] > waveform[i - 1]) {
      rollingSum += (waveform[i] - waveform[i - 1]);
    } else {
      rollingSum = 0;
    }

    if (rollingSum > fThreshold) {
      if (peakStart == -1) {
        peakStart = i;
      }
      if (i - peakStart >= fMinPeakWidth) {
        peaks[peakStart] = waveform[i];
      }
    } else {
      peakStart = -1;
    }
  }
}

void sbnd::MaMOpHT::SelectPeaks(const std::map<int, short>& peaks, std::vector<std::pair<float, short>>& selectedPeaks) const {
  for (auto it = peaks.begin(); it != peaks.end(); ++it) {
    auto next_it = std::next(it);
    if (next_it != peaks.end() && (next_it->first - it->first <= (int)fMaxPeakSeparation)) {
      if (next_it->second >= fPeakAmplitudeRatio * it->second) {
        float timestamp = it->first / 500.0 + startTime;  // Convert bin number back to time
        selectedPeaks.push_back({timestamp, it->second});
      }
    }
  }
}

DEFINE_ART_MODULE(sbnd::MaMOpHT)
