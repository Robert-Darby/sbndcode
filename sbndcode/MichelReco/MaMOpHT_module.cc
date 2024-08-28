////////////////////////////////////////////////////////////////////////
// Class:       MaMOpHT
// Plugin Type: producer (Unknown Unknown)
// File:        MaMOpHT_module.cc
//
// Generated at Tue Aug 13 08:17:12 2024 by Robert Darby using cetskelgen
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
#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>

namespace sbnd {
  class MaMOpHT;
}

class sbnd::MaMOpHT : public art::EDProducer {
public:
  explicit MaMOpHT(fhicl::ParameterSet const& p);
  
  // Plugins should not be copied or assigned.
  MaMOpHT(MaMOpHT const&) = delete;
  MaMOpHT(MaMOpHT&&) = delete;
  MaMOpHT& operator=(MaMOpHT const&) = delete;
  MaMOpHT& operator=(MaMOpHT&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  float CalculateBaselineAndSubtract(std::vector<short>& waveform) const;
  void ApplyRollingWindow(const std::vector<short>& waveform, std::map<int, short>& peaks) const;
  void SelectPeaks(const std::map<int, short>& peaks, std::vector<std::pair<int, int>>& selectedPeaks) const;

  const art::InputTag fInputTag;
  const short fThreshold;
  const size_t fMinPeakWidth;
  const size_t fMaxPeakSeparation;
  const float fPeakAmplitudeRatio;
};

sbnd::MaMOpHT::MaMOpHT(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fInputTag(p.get<std::string>("InputTag")),
    fThreshold(p.get<short>("Threshold", 1000)),
    fMinPeakWidth(p.get<size_t>("MinPeakWidth", 5)),
    fMaxPeakSeparation(p.get<size_t>("MaxPeakSeparation", 4500)),
    fPeakAmplitudeRatio(p.get<float>("PeakAmplitudeRatio", 0.01))
{
  produces<std::vector<std::pair<int, int>>>();
}

void sbnd::MaMOpHT::produce(art::Event& e) {
  auto waveforms = e.getValidHandle<std::vector<raw::OpDetWaveform>>(fInputTag);

  std::vector<short> cumulativeWaveform(3020 * 500, 0);

  for (auto const& waveform : *waveforms) {
    std::vector<short> wf = waveform.Waveform();
    float baseline = CalculateBaselineAndSubtract(wf);
    unsigned int startBin = std::max(0, static_cast<int>((waveform.TimeStamp() + 1510) * 500));
    std::cout << waveform.TimeStamp() << "     " << startBin << "\n";

    for (size_t i = 0; i < wf.size(); ++i) {
      if (startBin + i < cumulativeWaveform.size()) {
        cumulativeWaveform[startBin + i] += (baseline - wf[i]);
      }
    }
  }

  std::map<int, short> peaks;
  ApplyRollingWindow(cumulativeWaveform, peaks);
  std::vector<std::pair<int, int>> selectedPeaks;
  SelectPeaks(peaks, selectedPeaks);

  e.put(std::make_unique<std::vector<std::pair<int, int>>>(selectedPeaks));
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

void sbnd::MaMOpHT::SelectPeaks(const std::map<int, short>& peaks, std::vector<std::pair<int, int>>& selectedPeaks) const {
  for (auto it = peaks.begin(); it != peaks.end(); ++it) {
    auto next_it = std::next(it);
    if (next_it != peaks.end() && ((unsigned)next_it->first - (unsigned)it->first <= fMaxPeakSeparation)) {
      if (next_it->second >= fPeakAmplitudeRatio * it->second) {
        selectedPeaks.push_back(*it);
        selectedPeaks.push_back(*next_it);
      }
    }
  }
}

DEFINE_ART_MODULE(sbnd::MaMOpHT)

