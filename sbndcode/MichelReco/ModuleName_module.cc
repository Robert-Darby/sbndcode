////////////////////////////////////////////////////////////////////////
// Class:       LateLightAnalyser
// Module Type: analyzer
// File:        LateLightAnalyser_module.cc
//
// Generated at <timestamp> by [Your Name] using artmod
// from cetpkgsupport v1_14_00.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RawData/OpDetWaveform.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/ToolMacros.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

namespace sbnd {
  class LateLightAnalyser;
}

class sbnd::LateLightAnalyser : public art::EDAnalyzer {
public:
  explicit LateLightAnalyser(fhicl::ParameterSet const& p);
  void analyze(art::Event const& e) override;

  // Operators
  LateLightAnalyser(LateLightAnalyser const&) = delete;
  LateLightAnalyser(LateLightAnalyser&&) = delete;
  LateLightAnalyser& operator=(LateLightAnalyser const&) = delete;
  LateLightAnalyser& operator=(LateLightAnalyser&&) = delete;

private:
  const std::vector<art::InputTag> fInputLabels;
  const double fThreshold;
  const double fTimeBinWidth;
  const double fDaphneFreq;
  const double fCAENFreq;
  const double fStartTimeAfterMuon;
  const double fEndTimeAfterMuon;

  const opdet::sbndPDMapAlg fPDMap;

  TTree* fTree;

  int fChannel;
  int fEventID;
  int fRun;
  int fSubRun;
  double fMuonPeakTime;
  double fMichelPeakTime;
  double fMuonFitParam0;
  double fMuonFitParam1;
  double fMichelFitParam0;
  double fMichelFitParam1;
  std::vector<double> fBinnedWaveformMuon;
  std::vector<double> fBinnedWaveformMichel;

  void findPeaks(const std::vector<double>& waveform, std::vector<std::pair<size_t, double>>& peakIndices);
  void findPeakPairs(const std::vector<std::pair<size_t, double>>& peakIndices, std::vector<std::pair<size_t, size_t>>& peakPairs, int channel);
  std::pair<std::vector<double>, std::vector<double>> binAndFit(const std::vector<double>& data, size_t peakIndex, int channel);
};

// Constructor
sbnd::LateLightAnalyser::LateLightAnalyser(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fInputLabels(p.get<std::vector<art::InputTag>>("InputLabels")),
    fThreshold(p.get<double>("Threshold", 100.)),
    fTimeBinWidth(p.get<double>("TimeBinWidth")),
    fDaphneFreq(p.get<double>("DaphneFreq")),
    fCAENFreq(p.get<double>("CAENFreq")),
    fStartTimeAfterMuon(p.get<double>("StartTimeAfterMuon")),
    fEndTimeAfterMuon(p.get<double>("EndTimeAfterMuon")),
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("MyTree", "Waveform Analysis Tree");
  fTree->Branch("Channel", &fChannel, "Channel/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");
  fTree->Branch("Run", &fRun, "Run/I");
  fTree->Branch("SubRun", &fSubRun, "SubRun/I");
  fTree->Branch("MuonPeakTime", &fMuonPeakTime, "MuonPeakTime/D");
  fTree->Branch("MichelPeakTime", &fMichelPeakTime, "MichelPeakTime/D");
  fTree->Branch("MuonFitParam0", &fMuonFitParam0, "MuonFitParam0/D");
  fTree->Branch("MuonFitParam1", &fMuonFitParam1, "MuonFitParam1/D");
  fTree->Branch("MichelFitParam0", &fMichelFitParam0, "MichelFitParam0/D");
  fTree->Branch("MichelFitParam1", &fMichelFitParam1, "MichelFitParam1/D");
  fTree->Branch("BinnedWaveformMuon", &fBinnedWaveformMuon);
  fTree->Branch("BinnedWaveformMichel", &fBinnedWaveformMichel);
}

void sbnd::LateLightAnalyser::analyze(art::Event const& e) {
  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  for (const auto& label : fInputLabels) {
    auto handle = e.getHandle<std::vector<raw::OpDetWaveform>>(label);

    if (!handle) continue;

    for (const auto& waveform : *handle) {
      fChannel = waveform.ChannelNumber();
      auto opType = fPDMap.pdType(fChannel);
      auto conversionFactor = (fPDMap.electronicsType(fChannel) == "daphne") ? fDaphneFreq : fCAENFreq;

      const auto& waveformData = waveform.Waveform();
      const size_t nSamples = waveformData.size();

      // Calculate baseline (mean of the first 100 ns)
      size_t baselineSampleCount = static_cast<size_t>(100 / conversionFactor);
      double baseline = std::accumulate(waveformData.begin(), waveformData.begin() + baselineSampleCount, 0.0) / baselineSampleCount;

      // Subtract baseline from waveform
      std::vector<double> correctedWaveform(nSamples);
      std::transform(waveformData.begin(), waveformData.end(), correctedWaveform.begin(),
                     [baseline](double val) { return val - baseline; });

      // Apply rolling sum to the waveform
      std::vector<double> rollingSumWaveform(nSamples, 0.0);
      for (size_t i = 1; i < nSamples; ++i) {
        if (correctedWaveform[i] > correctedWaveform[i - 1]) {
          rollingSumWaveform[i] = rollingSumWaveform[i - 1] + (correctedWaveform[i] - correctedWaveform[i - 1]);
        } else {
          rollingSumWaveform[i] = 0.0;
        }
      }

      // Find peaks
      std::vector<std::pair<size_t, double>> peakIndices;
      findPeaks(correctedWaveform, peakIndices);

      // Find pairs of peaks
      std::vector<std::pair<size_t, size_t>> peakPairs;
      findPeakPairs(peakIndices, peakPairs, fChannel);

      // Process each pair
      for (const auto& [muonPeakIndex, michelPeakIndex] : peakPairs) {
        fMuonPeakTime = muonPeakIndex * conversionFactor;
        fMichelPeakTime = michelPeakIndex * conversionFactor;

        // Muon peak binning and fitting
        std::tie(fBinnedWaveformMuon, std::tie(fMuonFitParam0, fMuonFitParam1)) = binAndFit(correctedWaveform, muonPeakIndex, fChannel);

        // Subtract extrapolated fit from Michel peak
        double extrapolatedValueAtMichel = fMuonFitParam0 * std::exp(-(michelPeakIndex - muonPeakIndex) * conversionFactor / fMuonFitParam1);
        for (size_t i = michelPeakIndex; i < michelPeakIndex + static_cast<size_t>(fEndTimeAfterMuon / conversionFactor); ++i) {
          correctedWaveform[i] -= extrapolatedValueAtMichel;
        }

        // Michel peak binning and fitting
        std::tie(fBinnedWaveformMichel, std::tie(fMichelFitParam0, fMichelFitParam1)) = binAndFit(correctedWaveform, michelPeakIndex, fChannel);

        // Fill the tree
        fTree->Fill();
      }
    }
  }
}

void sbnd::LateLightAnalyser::findPeaks(const std::vector<double>& waveform, std::vector<std::pair<size_t, double>>& peakIndices) {
  const size_t nSamples = waveform.size();

  for (size_t i = 1; i < nSamples - 1; ++i) {
    if (waveform[i] > fThreshold && waveform[i] > waveform[i - 1] && waveform[i] > waveform[i + 1]) {
      double peak = waveform[i];
      peakIndices.push_back(std::make_pair(i, peak));
    }
  }
}

void sbnd::LateLightAnalyser::findPeakPairs(
  const std::vector<std::pair<size_t, double>>& peakIndices,
  std::vector<std::pair<size_t, size_t>>& peakPairs,
  int channel)
{
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  size_t minIntervalSamples = static_cast<size_t>(2000 / conversionFactor); // 2 us
  size_t maxIntervalSamples = static_cast<size_t>(10000 / conversionFactor); // 10 us
  size_t muonPreInterval = static_cast<size_t>(10000 / conversionFactor);

  for (size_t i = 0; i < peakIndices.size(); ++i) {
    if (i == 0 || (peakIndices[i].first - peakIndices[i - 1].first) > muonPreInterval) {
      double max_peak = std::numeric_limits<double>::min();
      size_t max_peak_idx = std::numeric_limits<size_t>::max();
      for (size_t j = i + 1; j < peakIndices.size(); ++j) {
        size_t interval = peakIndices[j].first - peakIndices[i].first;
        if (interval >= maxIntervalSamples) { i += j - i; break; }
        if (interval >= minIntervalSamples) {
          if (peakIndices[j].second > max_peak) {
            max_peak = peakIndices[j].second;
            max_peak_idx = peakIndices[j].first;
          }
        }
      }
      if (max_peak != std::numeric_limits<double>::min()) {
        peakPairs.emplace_back(peakIndices[i].first, max_peak_idx);
      }
    }
  }
}

std::pair<std::vector<double>, std::pair<double, double>> sbnd::LateLightAnalyser::binAndFit(const std::vector<double>& data, size_t peakIndex, int channel) {
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  size_t startBin = peakIndex + static_cast<size_t>(fStartTimeAfterMuon / conversionFactor);
  size_t endBin = peakIndex + static_cast<size_t>(fEndTimeAfterMuon / conversionFactor);
  size_t binWidthSamples = static_cast<size_t>(fTimeBinWidth / conversionFactor);
  size_t nBins = (endBin - startBin) / binWidthSamples;

  std::vector<double> bins(nBins);
  std::vector<double> binValues(nBins);

  for (size_t i = 0; i < nBins; ++i) {
    double binValue = std::accumulate(data.begin() + startBin + i * binWidthSamples,
                                      data.begin() + startBin + (i + 1) * binWidthSamples, 0.0) / binWidthSamples;
    bins[i] = i * fTimeBinWidth;
    binValues[i] = binValue;
  }

  // Fit an exponential decay
  TH1D hist("hist", "hist", nBins, 0, nBins * fTimeBinWidth);
  for (size_t i = 0; i < nBins; ++i) {
    hist.SetBinContent(i + 1, binValues[i]);
  }

  TF1 fit("fit", "[0]*exp(-x/[1])", 0, nBins * fTimeBinWidth);
  hist.Fit(&fit, "Q");

  double fitParam0 = fit.GetParameter(0);
  double fitParam1 = fit.GetParameter(1);

  return { binValues, { fitParam0, fitParam1 } };
}

DEFINE_ART_MODULE(sbnd::LateLightAnalyser)

