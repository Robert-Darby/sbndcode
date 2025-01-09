////////////////////////////////////////////////////////////////////////
// Class:       LateLightProducer
// Module Type: producer
// File:        LateLightProducer_module.cc
//
// Generated at <timestamp> by [Your Name] using artmod
// from cetpkgsupport v1_14_00.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RawData/OpDetWaveform.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h" // Find associations as pointers
#include "canvas/Persistency/Common/FindOneP.h"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "sbnobj/SBND/LateLight/LateLight.h"
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

namespace sbnd
{
  class LateLightProducer;
}

class sbnd::LateLightProducer : public art::EDProducer
{
public:
  explicit LateLightProducer(fhicl::ParameterSet const &p);
  void produce(art::Event &e) override;

  // Operators
  LateLightProducer(LateLightProducer const &) = delete;
  LateLightProducer(LateLightProducer &&) = delete;
  LateLightProducer &operator=(LateLightProducer const &) = delete;
  LateLightProducer &operator=(LateLightProducer &&) = delete;

  struct MCMuon
  {
    int G4ID;
    float Time;
    int MichelG4ID;
    float MichelTime, MichelEnergy;
    float FlashTime, MichelFlashTime;
    bool HasFlash, Stopping, EntersTPC;
  };
  struct CRTHit
  {
    int Plane;
    float Time;
    unsigned G4ID;
    int PDG;
    float AuxTime;
    bool HasFlash;
  };

private:
  // Constants. These shouldn't need to be changed in the fhicl
  const std::vector<int> fPair1, fPair2, fUnpaired;
  const std::vector<std::string> fPDs;
  const float fDaphneFreq;
  const float fCAENFreq;
  const float fPMTReadoutDelay, fARAReadoutDelay;

  // For finding peaks in summed waveform
  const std::string fFinderInputLabel;
  const std::string fPMTTriggerLabel;
  const int fMinPMTMultiplicity;
  const float fFinderPolarity;
  const float fFindMuonThreshold, fFindMichelPeak;
  const std::string fCRTLabel;
  const float fCRTOffset, fCRTClockSpeed;
  const float fCRTCoincidence;

  // For finding peaks in individual waveforms
  const std::vector<art::InputTag> fInputLabels;
  const std::vector<float> fInputPolarity;
  const std::vector<float> fMuonThresholds, fMichelThresholds;
  const float fPeakCoincWindow;
  const float fTimeBinWidth;
  const float fStartTimeAfterMuon;
  const float fEndTimeAfterMuon;

  // MC Info
  const bool fUseMC;
  const std::string fMCTruthLabel;
  const bool fOnlyStopping;
  std::string fAuxDetHitLabel;
  std::vector<std::string> fAuxDetHitTags;

  // Misc
  const bool fRequireMichel;
  const float fMinLifetime;
  const float fMaxLifetime;
  const bool fSaveHists;
  const std::vector<int> fSaveEvents;
  const bool fBinWaveform;
  const bool fVerbose;

  const opdet::sbndPDMapAlg fPDMap;
  std::vector<MCMuon> muon_tuple_vect;
  std::vector<std::pair<float, float>> opMuons;
  std::vector<CRTHit> crtHits, wvfCRTHits;
  std::vector<float> muonTimes, michelTimes;
  float wvfStartTime, wvfEndTime;
  float muonThreshold, michelThreshold;
  unsigned nSamples;

  // Other variables shared between different methods.
  geo::GeometryCore const *fGeometryService;
  sbnd::crt::CRTGeoAlg fCrtGeo;
  std::map<int, std::vector<std::pair<int, float>>> muonMultCoatMap, michelMultCoatMap, muonMultUncoatMap, michelMultUncoatMap;

  TTree *fTree;

  int fChannel;
  int fEventID;
  int fRun;
  int fSubRun;
  int fOpChannel, fOpChannelType;
  float fMuonPeakTime;
  float fMuonPeakAmp;
  int fMichelID;
  float fMichelTime, fMichelEnergy;
  float fMichelPeakAmp;
  bool fHasMichelPeak;
  float fMichelPeakTime;
  float fMuonFitParam0;
  float fMuonFitParam1;
  float fMichelFitParam0;
  float fMichelFitParam1;
  std::vector<float> fBinnedWaveformMuon;
  std::vector<float> fBinnedWaveformMichel;
  std::vector<float> fRawWaveform;
  int fMCMuonG4ID;
  float fMCMuonTime;
  bool fMCMuonStopping;
  bool fMCMuonEntersTPC;
  bool fMCMuonHasFlash;
  float fMuonFlashTime, fMichelFlashTime;

  TTree *fTriggerTree;
  int fTriggerG4ID, fTriggerPDG, fTriggerPlane;
  bool fTriggerHasFlash;
  std::vector<float> fMuonMultCoat, fMichelMultCoat, fMuonMultUncoat, fMichelMultUncoat;

  int findChannelPair(int opChannel);
  void findMCMuons(const art::Event &e);
  std::vector<float> subtractBaseline(const std::vector<float> &waveform, const float conversionFactor);
  std::vector<float> applyRollingSum(const std::vector<float> &waveform);
  void findPeaks(const std::vector<float> &seazrchWaveform, const std::vector<float> &waveform, std::vector<std::pair<size_t, float>> &peakIndices);
  void findPeakPairs(const std::vector<std::pair<size_t, float>> &peakIndices, const bool use_opmuons, std::vector<std::pair<size_t, size_t>> &peakPairs, int channel);
  void findCRTTimes(const art::Event &e);
  std::vector<std::pair<float, float>> findOpMuons(const art::Event &e, std::unique_ptr<std::vector<sbnd::StoppingMuonTrigger>> &stoppingmu_v);
  void addPeakToMap(std::map<int, std::vector<std::pair<int, float>>> &multMap, int g4id, int opChannel, int pair_channel, float peakAmp, bool requirePositive);
  std::pair<std::vector<float>, std::vector<float>> binAndFit(const std::vector<float> &data, size_t peakIndex, int channel);
  void ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, int eventNumber);
};

// Constructor
sbnd::LateLightProducer::LateLightProducer(fhicl::ParameterSet const &p)
    : EDProducer{p},
      // Constants: These shouldn't need to be changed in the fhicl
      fPair1(p.get<std::vector<int>>("Pair1")),
      fPair2(p.get<std::vector<int>>("Pair2")),
      fUnpaired(p.get<std::vector<int>>("Unpaired")),
      fPDs(p.get<std::vector<std::string>>("PDs")),
      fDaphneFreq(p.get<float>("DaphneFreq")),
      fCAENFreq(p.get<float>("CAENFreq")),
      fPMTReadoutDelay(p.get<float>("PMTReadoutDelay", 0.135)),
      fARAReadoutDelay(p.get<float>("XARAPUCAReadoutDelay", 0.0)),

      // For finding peaks in summed waveform
      fFinderInputLabel(p.get<std::string>("FinderInputLabel")),
      fPMTTriggerLabel(p.get<std::string>("PMTTriggerLabel")),
      fMinPMTMultiplicity(p.get<int>("MinPMTMultiplicity", 10)),
      fFinderPolarity(p.get<float>("FinderPolarity", 1.)),
      fFindMuonThreshold(p.get<float>("FindMuonThreshold", 5000.)),
      fFindMichelPeak(p.get<float>("FindMichelPeak", 1000.)),
      fCRTLabel(p.get<std::string>("CRTLabel")),
      fCRTOffset(p.get<float>("CRTOffset", 1700000.)),
      fCRTClockSpeed(p.get<float>("CRTClockspeed", 1.)),
      fCRTCoincidence(p.get<float>("CRTCoincidence")),

      // For finding peaks in individual waveforms
      fInputLabels(p.get<std::vector<art::InputTag>>("InputLabels")),
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
      fOnlyStopping(p.get<bool>("OnlyStopping", true)),
      fAuxDetHitLabel(p.get<std::string>("AuxDetHitLabel")),
      fAuxDetHitTags(p.get<std::vector<std::string>>("AuxDetHitTags")),

      // Misc
      fRequireMichel(p.get<bool>("RequireMichel", true)),
      fMinLifetime(p.get<float>("MinLifetime", 2000.0f)),
      fMaxLifetime(p.get<float>("MaxLifetime", 10000.)),
      fSaveHists(p.get<bool>("SaveHists", false)),
      fSaveEvents(p.get<std::vector<int>>("SaveEvents")),
      fBinWaveform(p.get<bool>("BinWaveforms", false)),
      fVerbose(p.get<bool>("Verbose", false))
{
  produces<std::vector<sbnd::StoppingMuonTrigger>>();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TailTree", "Waveform Producis Tree");
  fTree->Branch("Channel", &fChannel, "Channel/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");
  fTree->Branch("Run", &fRun, "Run/I");
  fTree->Branch("SubRun", &fSubRun, "SubRun/I");
  fTree->Branch("OpChannel", &fOpChannel);
  fTree->Branch("OpChannelType", &fOpChannelType);
  fTree->Branch("MuonPeakTime", &fMuonPeakTime);
  fTree->Branch("MichelPeakTime", &fMichelPeakTime);
  fTree->Branch("MuonPeakAmp", &fMuonPeakAmp);
  fTree->Branch("MichelPeakAmp", &fMichelPeakAmp);
  fTree->Branch("HasMichelPeak", &fHasMichelPeak);
  fTree->Branch("MuonFitParam0", &fMuonFitParam0, "MuonFitParam0/F");
  fTree->Branch("MuonFitParam1", &fMuonFitParam1, "MuonFitParam1/F");
  fTree->Branch("MichelFitParam0", &fMichelFitParam0, "MichelFitParam0/F");
  fTree->Branch("MichelFitParam1", &fMichelFitParam1, "MichelFitParam1/F");
  fTree->Branch("BinnedWaveformMuon", &fBinnedWaveformMuon);
  fTree->Branch("raw_waveform", &fRawWaveform);
  fTree->Branch("BinnedWaveformMichel", &fBinnedWaveformMichel);
  fTree->Branch("MCMuonG4ID", &fMCMuonG4ID);
  fTree->Branch("MCMuonTime", &fMCMuonTime);
  fTree->Branch("MCMuonStopping", &fMCMuonStopping);
  fTree->Branch("MCMUonEntersTPC", &fMCMuonEntersTPC);

  // Set branches
  fTriggerTree = tfs->make<TTree>("TriggerTree", "Waveform Producers Tree");
  fTriggerTree->Branch("Run", &fRun);
  fTriggerTree->Branch("Subrun", &fSubRun);
  fTriggerTree->Branch("Event", &fEventID);
  fTriggerTree->Branch("G4ID", &fMCMuonG4ID);
  fTriggerTree->Branch("mich_g4id", &fMichelID);
  fTriggerTree->Branch("michel_energy", &fMichelEnergy);
  fTriggerTree->Branch("michel_time", &fMichelTime);
  fTriggerTree->Branch("is_stopping", &fMCMuonStopping);
  fTriggerTree->Branch("Enters_tpc", &fMCMuonEntersTPC);
  fTriggerTree->Branch("HasFlash", &fMCMuonHasFlash);
  fTriggerTree->Branch("MuonFlashTime", &fMuonFlashTime);
  fTriggerTree->Branch("MichelFlashTime", &fMichelFlashTime);
  fTriggerTree->Branch("Timestamp", &fMCMuonTime);
  fTriggerTree->Branch("CRTG4ID", &fTriggerG4ID);
  fTriggerTree->Branch("CRTPDG", &fTriggerPDG);
  fTriggerTree->Branch("CRTPlane", &fTriggerPlane);
  fTriggerTree->Branch("MuonMultCoat", &fMuonMultCoat);
  fTriggerTree->Branch("MichelMultCoat", &fMichelMultCoat);
  fTriggerTree->Branch("MuonMultUncoat", &fMuonMultUncoat);
  fTriggerTree->Branch("MichelMultUncoat", &fMichelMultUncoat);
}

void sbnd::LateLightProducer::produce(art::Event &e)
{
  std::unique_ptr<std::vector<sbnd::StoppingMuonTrigger>>
      stoppingmu_v(new std::vector<sbnd::StoppingMuonTrigger>);

  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();
  muon_tuple_vect.clear();

  muonMultCoatMap.clear();
  muonMultUncoatMap.clear();
  michelMultCoatMap.clear();
  michelMultUncoatMap.clear();
  fMuonMultCoat.clear();
  fMuonMultUncoat.clear();
  fMichelMultCoat.clear();
  fMichelMultUncoat.clear();

  if (fUseMC)
  {
    findMCMuons(e);
    for (auto &muon : muon_tuple_vect)
    {
      std::vector<std::pair<int, float>> blank_vec;
      muonMultCoatMap.insert(std::make_pair(muon.G4ID, blank_vec));
      michelMultCoatMap.insert(std::make_pair(muon.G4ID, blank_vec));
      muonMultUncoatMap.insert(std::make_pair(muon.G4ID, blank_vec));
      michelMultUncoatMap.insert(std::make_pair(muon.G4ID, blank_vec));
    }
  }

  opMuons.clear();
  opMuons = findOpMuons(e, stoppingmu_v);
  if (fVerbose)
  {
    std::cout << "OpMuons: \n";
    for (auto muonmichel : opMuons)
    {
      std::cout << "   " << muonmichel.first << "   " << muonmichel.second << "\n";
    }
  }
  if (fVerbose)
  {
    std::cout << "MC Muons\n";
    for (auto &muon : muon_tuple_vect)
    {
      std::cout << "    " << muon.G4ID << "   " << muon.Time << "   " << muon.MichelTime;
      if (muon.HasFlash)
        std::cout << "    TRIGGERED \n";
      else
        std::cout << "\n";
    }
  }

  for (size_t label_idx = 0; label_idx < fInputLabels.size(); label_idx++)
  {
    auto label = fInputLabels[label_idx];
    auto polarity = fInputPolarity[label_idx];
    auto handle = e.getHandle<std::vector<raw::OpDetWaveform>>(label);

    if (!handle)
      continue;

    for (const auto &waveform : *handle)
    {
      fChannel = waveform.ChannelNumber();
      int pair_channel = findChannelPair(fChannel);
      auto opType = fPDMap.pdType(fChannel);
      auto conversionFactor = (fPDMap.electronicsType(fChannel) == "daphne") ? fDaphneFreq : fCAENFreq;
      float readoutDelay = (opType == "pmt_coated" || "pmt_uncoated") ? fPMTReadoutDelay : fARAReadoutDelay;
      unsigned opTypeIdx = std::distance(fPDs.begin(), std::find(fPDs.begin(), fPDs.end(), opType));
      muonThreshold = fMuonThresholds.at(opTypeIdx);
      michelThreshold = fMichelThresholds.at(opTypeIdx);

      const std::vector<short> rawWaveformData = waveform.Waveform();
      std::vector<float> waveformData(rawWaveformData.size(), 0.);
      std::transform(
          rawWaveformData.begin(), rawWaveformData.end(), waveformData.begin(),
          [](short adc)
          { return (float)adc; });
      nSamples = waveformData.size();

      // Find if a muon and michel time is within the waveform, skip if not
      wvfEndTime = waveform.TimeStamp() - readoutDelay + nSamples * conversionFactor / 1000.;
      wvfStartTime = waveform.TimeStamp() - readoutDelay;
      muonTimes.clear();
      michelTimes.clear();

      for (auto muonmichel : opMuons)
      {
        if (muonmichel.first > wvfStartTime && muonmichel.first < wvfEndTime)
        {
          muonTimes.push_back(muonmichel.first);
          michelTimes.push_back(muonmichel.second);
        }
      }
      if (muonTimes.empty())
        continue;
      wvfCRTHits.clear();
      for (auto crthit : crtHits)
      {
        if (crthit.Time > wvfStartTime && crthit.Time < wvfEndTime)
          wvfCRTHits.push_back(crthit);
      }

      // Subtract baseline
      auto correctedWaveform = subtractBaseline(waveformData, conversionFactor);
      for (auto &adc : correctedWaveform)
        adc *= polarity;
      // Apply rolling sum to the waveform
      auto rollingSumWaveform = applyRollingSum(correctedWaveform);

      // Find peaks
      std::vector<std::pair<size_t, float>> peakIndices;
      findPeaks(rollingSumWaveform, correctedWaveform, peakIndices);

      // Find pairs of peaks
      std::vector<std::pair<size_t, size_t>> peakPairs;
      findPeakPairs(peakIndices, true, peakPairs, fChannel);

      if (peakPairs.empty())
        continue;

      // Process each pair
      fMuonPeakTime = -9999.;
      fMichelPeakTime = -9999.;
      for (const auto &[muonPeakIndex, michelPeakIndex] : peakPairs)
      {
        fHasMichelPeak = (michelPeakIndex < nSamples);
        fMuonPeakTime = ((float)muonPeakIndex * conversionFactor / 1000.) + (float)wvfStartTime ;
        fMichelPeakTime = (fHasMichelPeak) ? ((float)michelPeakIndex * conversionFactor / 1000.) + (float)wvfStartTime: -9999.;
        fMuonPeakAmp = correctedWaveform[muonPeakIndex];
        fMichelPeakAmp = correctedWaveform[michelPeakIndex];

        fRawWaveform.clear();
        if (fHasMichelPeak && fMichelPeakTime - fMuonPeakTime > 2.)
        {
          int startBin = (muonPeakIndex < 50) ? 0 : muonPeakIndex - 50;
          int endBin = (muonPeakIndex + 4500 > nSamples) ? nSamples : muonPeakIndex + 4500;
          for (unsigned i = (unsigned)startBin; i < (unsigned)endBin; i++)
            fRawWaveform.push_back(correctedWaveform[i]);
        }
        sbnd::LateLightTail latelighttail;
        latelighttail.Channel = fOpChannel;
        latelighttail.Waveform = fRawWaveform;
        for (auto &stoppingmu : *stoppingmu_v)
        {
          std::cout << stoppingmu.MuonTime << "   " << fMuonPeakTime << "\n";
          if (abs(stoppingmu.MuonTime - fMuonPeakTime) < 0.05)
          {
            std::cout << "Added latelighttail \n";
            stoppingmu.LateLightTails.push_back(latelighttail);
            break;
          }
        }

        fMCMuonG4ID = -9999;
        fMCMuonTime = std::numeric_limits<float>::max();
        fMCMuonStopping = false;
        fMCMuonEntersTPC = false;
        fMCMuonHasFlash = false;
        fMuonFlashTime = -9999.;
        fMichelFlashTime = -9999.;
        fTriggerG4ID = -1;
        if (fUseMC)
        {
          for (auto &crthit : wvfCRTHits)
          {
            if (abs(fMuonPeakTime - crthit.Time) * 1000. < fCRTCoincidence)
            {
              fTriggerG4ID = crthit.G4ID;
              fTriggerPDG = crthit.PDG;
              break;
            }
          }
          if (abs(fTriggerPDG) == 13)
          {
            for (auto &mcmuon : muon_tuple_vect)
            {
              if (abs(mcmuon.Time - fMuonPeakTime) < 2.4)
              {
                fMCMuonG4ID = mcmuon.G4ID;
                fMCMuonTime = mcmuon.Time;
                fMCMuonStopping = mcmuon.Stopping;
                fMCMuonEntersTPC = mcmuon.EntersTPC;
                fMCMuonHasFlash = mcmuon.HasFlash;
                fMuonFlashTime = mcmuon.FlashTime;
                fMichelFlashTime = mcmuon.MichelFlashTime;
              }
            }
          } // If flash came from a muon
        } // If fUseMC

        fOpChannel = fChannel;
        if (opType == "pmt_coated")
          fOpChannelType = 0;
        else if (opType == "pmt_uncoated")
          fOpChannelType = 1;
        else
          continue;

        // Muon peak binning and fitting
        std::vector<float> fMuonParams(2, 0.);
        if (fHasMichelPeak)
        {
          std::tie(fBinnedWaveformMuon, fMuonParams) = binAndFit(correctedWaveform, muonPeakIndex, fChannel);

          // Subtract extrapolated fit from Michel peak
          fMuonFitParam0 = fMuonParams[0];
          fMuonFitParam1 = fMuonParams[1];
          for (size_t i = michelPeakIndex; i < michelPeakIndex + static_cast<size_t>(fEndTimeAfterMuon / conversionFactor); ++i)
          {
            float extrapolatedValueAtMichel = fMuonFitParam0 * std::exp(-(i - muonPeakIndex) * conversionFactor / fMuonFitParam1);
            correctedWaveform[i] -= extrapolatedValueAtMichel;
          }
        }

        // Michel peak binning and fitting
        std::vector<float> fMichelParams;
        std::tie(fBinnedWaveformMichel, fMichelParams) = binAndFit(correctedWaveform, michelPeakIndex, fChannel);
        fMichelFitParam0 = fMichelParams[0];
        fMichelFitParam1 = fMichelParams[1];
        //        std::cout << fMichelFitParam0 << "    " << fMichelFitParam1 << "\n";

        // Fill the tree
        if (fTriggerG4ID < 0)
          continue;
        if (opType == "pmt_coated")
        {
          addPeakToMap(muonMultCoatMap, fTriggerG4ID, fOpChannel, pair_channel, fMuonPeakAmp, false);
          addPeakToMap(michelMultCoatMap, fTriggerG4ID, fOpChannel, pair_channel, fMichelPeakAmp, true);
        }
        else if (opType == "pmt_uncoated")
        {
          addPeakToMap(muonMultUncoatMap, fTriggerG4ID, fOpChannel, pair_channel, fMuonPeakAmp, false);
          addPeakToMap(michelMultUncoatMap, fTriggerG4ID, fOpChannel, pair_channel, fMichelPeakAmp, true);
        } // If PMT uncoated
        fRawWaveform.clear();
        if (fHasMichelPeak && fMichelPeakTime - fMuonPeakTime > 2.)
        {
          int startBin = (muonPeakIndex < 50) ? 0 : muonPeakIndex - 50;
          int endBin = (muonPeakIndex + 4500 > nSamples) ? nSamples : muonPeakIndex + 4500;
          for (unsigned i = (unsigned)startBin; i < (unsigned)endBin; i++)
            fRawWaveform.push_back(correctedWaveform[i]);
        }
        fTree->Fill();
        if (fTriggerG4ID > 0 && fSaveHists && !peakPairs.empty() &&
            std::find(fSaveEvents.begin(), fSaveEvents.end(), fEventID) != fSaveEvents.end() &&
            fMichelPeakTime - fMuonPeakTime > 2.)
          ConvertWaveformToHistogram(waveform, fEventID);
      } // Loop over waveform handles
    } // Loop over peak pairs

  } // Loop over labels

  // Fill Trigger Tree
  for (auto it = muonMultCoatMap.begin(); it != muonMultCoatMap.end(); ++it)
  {
    auto g4id = it->first;
    fMuonMultCoat.clear();
    fMuonMultUncoat.clear();
    fMichelMultCoat.clear();
    fMichelMultUncoat.clear();

    if (muonMultCoatMap.find(g4id) != muonMultCoatMap.end())
    {
      for (auto muon_opamp : muonMultCoatMap[g4id])
        fMuonMultCoat.push_back(muon_opamp.second);
    }
    if (michelMultCoatMap.find(g4id) != michelMultCoatMap.end())
    {
      for (auto michel_opamp : michelMultCoatMap[g4id])
        fMichelMultCoat.push_back(michel_opamp.second);
    }
    if (muonMultUncoatMap.find(g4id) != muonMultUncoatMap.end())
    {
      for (auto muon_opamp : muonMultUncoatMap[g4id])
        fMuonMultUncoat.push_back(muon_opamp.second);
    }
    if (michelMultUncoatMap.find(g4id) != michelMultUncoatMap.end())
    {
      for (auto michel_opamp : michelMultUncoatMap[g4id])
        fMichelMultUncoat.push_back(michel_opamp.second);
    }

    if (fVerbose)
    {
      std::cout << "\nFound " << fMuonMultCoat.size() << " coated muon peaks for " << g4id << "\n";
      std::cout << "Found " << fMichelMultCoat.size() << " coated michel peaks for " << g4id << "\n";
      std::cout << "Found " << fMuonMultUncoat.size() << " uncoated muon peaks for " << g4id << "\n";
      std::cout << "Found " << fMichelMultUncoat.size() << " uncoated michel peaks for " << g4id << "\n";
    }

    std::sort(fMuonMultCoat.begin(), fMuonMultCoat.end());
    std::sort(fMichelMultCoat.begin(), fMichelMultCoat.end());
    std::sort(fMuonMultUncoat.begin(), fMuonMultUncoat.end());
    std::sort(fMichelMultUncoat.begin(), fMichelMultUncoat.end());

    fMichelMultCoat.erase(std::remove_if(fMichelMultCoat.begin(), fMichelMultCoat.end(),
                                         [](float adc)
                                         { return (adc < 100); }),
                          fMichelMultCoat.end());
    fMichelMultUncoat.erase(std::remove_if(fMichelMultUncoat.begin(), fMichelMultUncoat.end(),
                                           [](float adc)
                                           { return (adc < 100); }),
                            fMichelMultUncoat.end());

    // Fill MC info if its a muon boi
    fMuonFlashTime = -9999.;
    fMichelFlashTime = -9999.;
    fMCMuonG4ID = -1;
    fMCMuonHasFlash = false;
    fMCMuonStopping = false;
    fMCMuonEntersTPC = false;
    fMichelEnergy = -9999.;
    fMichelTime = -9999.;
    fMichelID = false;
    fMCMuonTime = -9999.;
    for (auto &muon_tuple : muon_tuple_vect)
    {
      if (muon_tuple.G4ID == g4id)
      {
        fMCMuonG4ID = g4id;
        fMCMuonTime = muon_tuple.Time;
        fMichelID = muon_tuple.MichelG4ID;
        fMichelTime = muon_tuple.MichelTime;
        fMichelEnergy = muon_tuple.MichelEnergy;
        fMCMuonStopping = muon_tuple.Stopping;
        fMCMuonEntersTPC = muon_tuple.EntersTPC;
        fMCMuonHasFlash = muon_tuple.HasFlash;
        fMuonFlashTime = muon_tuple.FlashTime;
        fMichelFlashTime = muon_tuple.MichelFlashTime;
      }
    }

    // Fill CRT Info
    fTriggerG4ID = -1;
    fTriggerPDG = -1;
    fTriggerPlane = -1;
    for (auto &crthit : crtHits)
    {
      if (crthit.G4ID == (unsigned)g4id)
      {
        fTriggerG4ID = crthit.G4ID;
        fTriggerPDG = crthit.PDG;
        fTriggerPlane = crthit.Plane;
        fTriggerHasFlash = crthit.HasFlash;
      }
    }

    fTriggerTree->Fill();
  } // Loop over MC muon_tuple_vect
  e.put(std::move(stoppingmu_v));
}

int sbnd::LateLightProducer::findChannelPair(int opChannel)
{
  int pair_channel = -1;
  int pair_idx;
  auto pair_it = std::find(fPair1.begin(), fPair1.end(), opChannel); // Fixing typo fOpChannel to opChannel
  if (std::find(fUnpaired.begin(), fUnpaired.end(), opChannel) != fUnpaired.end())
    return opChannel;
  if (pair_it != fPair1.end())
  {
    pair_idx = std::distance(fPair1.begin(), pair_it);
    pair_channel = fPair2[pair_idx];
  }
  else
  {
    pair_it = std::find(fPair2.begin(), fPair2.end(), opChannel);
    pair_idx = std::distance(fPair2.begin(), pair_it);
    pair_channel = fPair1[pair_idx];
  }
  return pair_channel; // Moved outside the else block to fix error
}

void sbnd::LateLightProducer::findMCMuons(const art::Event &e)
{
  art::Handle<std::vector<simb::MCParticle>> mctruthHandle;
  std::vector<art::Ptr<simb::MCParticle>> mctruthVect;
  if (e.getByLabel(fMCTruthLabel, mctruthHandle)) // Make sure artHandle is from module
    art::fill_ptr_vector(mctruthVect, mctruthHandle);

  for (const auto &mcp : mctruthVect)
  {
    if (abs(mcp->PdgCode()) != 13 || abs(mcp->T()) / 1000. > 1510.) // ||
                                                                    // abs(mcp->EndX()) > 400. ||
                                                                    // abs(mcp->EndY()) > 400. ||
                                                                    // mcp->EndZ() > 700 || mcp->EndZ() < -200.)
      continue;
    auto mu_id = mcp->TrackId();
    bool mu_decayintpc = (abs(mcp->EndX()) < 200. &&
                          abs(mcp->EndY()) < 200. &&
                          mcp->EndZ() < 500. && mcp->EndZ() > 0.);
    if (fOnlyStopping && !mu_decayintpc)
      continue;
    float mu_time = std::numeric_limits<float>::max();
    bool mu_enterstpc = false;
    for (unsigned i = 0; i < mcp->NumberTrajectoryPoints(); i++)
    {
      if (abs(mcp->Position(i).X()) < 200. &&
          abs(mcp->Position(i).Y()) < 200. &&
          mcp->Position(i).Z() < 500. && mcp->Position(i).Z() > 0.)
      {
        mu_time = mcp->T(i) / 1000.;
        mu_enterstpc = true;
        break;
      }
    }
    int michel_id = -1;
    float michel_time = -9999.;
    float michel_energy = -9999.;
    for (auto &mcp2 : mctruthVect)
    {
      if (mcp2->Mother() != mcp->TrackId() ||
          (mcp2->Position() - mcp->EndPosition()).Mag() > 10. ||
          abs(mcp2->PdgCode()) != 11)
        continue;
      michel_id = mcp2->TrackId();
      michel_time = mcp2->T() / 1000.;
      michel_energy = mcp2->E() * 1000.;
    }
    // if (michel_energy < 0.)
    //   continue;
    if (!mu_enterstpc)
    {
      mu_time = mcp->T() / 1000.;
    }
    if (michel_time - mu_time < fMinLifetime / 1000.)
    {
      continue;
    }
    MCMuon mcmuon;
    mcmuon.G4ID = mu_id;
    mcmuon.Time = mu_time;
    mcmuon.MichelG4ID = michel_id;
    mcmuon.MichelTime = michel_time;
    mcmuon.MichelEnergy = michel_energy;
    mcmuon.Stopping = mu_decayintpc;
    mcmuon.EntersTPC = mu_enterstpc;
    mcmuon.HasFlash = false;
    mcmuon.FlashTime = -9999.;
    mcmuon.MichelFlashTime = -9999.;
    muon_tuple_vect.push_back(mcmuon);
  }
  std::sort(muon_tuple_vect.begin(), muon_tuple_vect.end(),
            [](const MCMuon &muon1, const MCMuon &muon2)
            { return (muon2.Time > muon1.Time); });
}

std::vector<float> sbnd::LateLightProducer::subtractBaseline(const std::vector<float> &waveform, const float conversionFactor)
{
  // Calculate baseline (mean of the first 100 ns)
  size_t baselineSampleCount = static_cast<size_t>(100 / conversionFactor);
  float baseline = std::accumulate(waveform.begin(), waveform.begin() + baselineSampleCount, 0.0f) / baselineSampleCount;

  // Subtract baseline from waveform
  std::vector<float> correctedWaveform(waveform.size());
  std::transform(waveform.begin(), waveform.end(), correctedWaveform.begin(),
                 [baseline](float val)
                 { return val - baseline; });

  return correctedWaveform;
}

std::vector<float> sbnd::LateLightProducer::applyRollingSum(const std::vector<float> &waveform)
{
  std::vector<float> rollingSum(waveform.size(), 0.);
  for (size_t i = 1; i < waveform.size(); ++i)
  {
    if (waveform[i] > waveform[i - 1])
    {
      rollingSum[i] = rollingSum[i - 1] + (waveform[i] - waveform[i - 1]);
    }
    else
    {
      rollingSum[i] = 0.0;
    }
  }
  return rollingSum;
}

void sbnd::LateLightProducer::findPeaks(
    const std::vector<float> &searchWaveform,
    const std::vector<float> &waveform,
    std::vector<std::pair<size_t, float>> &peakIndices)
{
  const size_t nSamples = waveform.size();

  for (size_t i = 1; i < nSamples - 1; ++i)
  {
    if (searchWaveform[i] > michelThreshold && searchWaveform[i] > searchWaveform[i - 1] && searchWaveform[i] > searchWaveform[i + 1])
    {
      float peak = searchWaveform[i];
      peakIndices.push_back(std::make_pair(i, peak));
    }
  }
  std::sort(peakIndices.begin(), peakIndices.end(),
            [](std::pair<size_t, float> peak1, std::pair<size_t, float> peak2)
            { return peak2.first > peak1.first; });
}

void sbnd::LateLightProducer::findPeakPairs(
    const std::vector<std::pair<size_t, float>> &peakIndices,
    const bool use_opmuons,
    std::vector<std::pair<size_t, size_t>> &peakPairs,
    int channel = 6)
{
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  std::vector<std::pair<size_t, size_t>> muon_windows, michel_windows;
  size_t mu_coincidence_window, michel_coincidence_window;
  std::vector<float> muon_times, michel_times;
  if (!use_opmuons)
  {
    mu_coincidence_window = (size_t)(fCRTCoincidence / conversionFactor);
    michel_coincidence_window = (size_t)((fMaxLifetime - fMinLifetime) / (2. * conversionFactor));
    for (auto crthit : wvfCRTHits)
    {
      muon_times.push_back((crthit.Time - wvfStartTime) * 1000.);
      michel_times.push_back((1000. * (crthit.Time - wvfStartTime) + fMinLifetime + (fMaxLifetime - fMinLifetime) / 2.));
    }
  }
  else
  {
    mu_coincidence_window = (size_t)(fPeakCoincWindow / conversionFactor);
    michel_coincidence_window = mu_coincidence_window;
    for (unsigned i = 0; i < muonTimes.size(); i++)
    {
      muon_times.push_back((muonTimes[i] - wvfStartTime) * 1000.);
      michel_times.push_back((michelTimes[i] - wvfStartTime) * 1000.);
    }
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
    size_t mu_peak2_idx = -1;
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
          mu_peak2_idx = mu_peak_idx;
          if (false)
            std::cout << mu_peak2_idx << "\n";
          mu_peak = peakIndices[i_peak].second;
          mu_peak_idx = peakIndices[i_peak].first;
        }
        else
        {
          mu_peak2 = peakIndices[i_peak].second;
          mu_peak2_idx = peakIndices[i_peak].first;
        } // If one of the two biggest peaks in muon search window
      }
      // If new tallest peak in michel window
      if (peakIndices[i_peak].first > michel_windows[i_mu].first && peakIndices[i_peak].first < michel_windows[i_mu].second &&
          peakIndices[i_peak].second > michel_peak2 && peakIndices[i_peak].second > michelThreshold)
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
      peakPairs.push_back(std::make_pair(mu_peak_idx, michel_peak_idx));
  } // Loop over candidate times
} // findPeakPairs

void sbnd::LateLightProducer::findCRTTimes(const art::Event &e)
{
  crtHits.clear();
  auto febData = e.getValidHandle<std::vector<sbnd::crt::FEBData>>(fCRTLabel);
  std::vector<std::pair<float, unsigned>> valid_hits;

  for (const auto &febdata : *febData)
  {
    auto mac5 = febdata.Mac5();
    auto plane = sbnd::crt::CRTCommonUtils::GetTaggerEnum(fCrtGeo.GetTaggerName(fCrtGeo.ChannelToStripName(mac5 * 32)));
    auto ts1 = ((float)febdata.Ts1() - fCRTOffset) / fCRTClockSpeed;
    auto flag = febdata.Flags();
    if (flag != 3 || plane > 6)
      continue;
    valid_hits.push_back(std::make_pair((float)ts1, plane));
  }
  std::sort(valid_hits.begin(), valid_hits.end(),
            [](std::pair<float, unsigned> hit1, std::pair<float, unsigned> hit2)
            { return (hit2.first > hit1.first); });

  for (unsigned i_hit = 0; i_hit < valid_hits.size(); i_hit++)
  {
    bool top_high = (valid_hits[i_hit].second == 6);
    bool top_low = (valid_hits[i_hit].second == 5);
    bool veto = false;
    if (valid_hits[i_hit].second == 0)
      veto = true;
    unsigned diff = 0;
    for (unsigned j_hit = i_hit + 1; j_hit < valid_hits.size(); j_hit++)
    {
      if (valid_hits[j_hit].first - valid_hits[i_hit].first > 240. ||
          j_hit == valid_hits.size() - 1)
      {
        diff = j_hit - i_hit;
        break;
      }
      if (valid_hits[j_hit].second == valid_hits[i_hit].second)
        continue;
      if (valid_hits[j_hit].second == 6)
        top_high = true;
      if (valid_hits[j_hit].second == 5)
        top_low = true;
      if (valid_hits[j_hit].second == 0 ||
          (valid_hits[i_hit].second != valid_hits[j_hit].second &&
           (valid_hits[i_hit].second < 5 || valid_hits[j_hit].second < 5)))
        veto = true;
    }
    if (!veto)
    {
      unsigned plane = valid_hits[i_hit].second;
      if (top_high && top_low)
        plane = 56;
      CRTHit crthit;
      crthit.Plane = plane;
      crthit.Time = valid_hits[i_hit].first / 1000.;
      crthit.G4ID = -1;
      crthit.PDG = 9999;
      crthit.HasFlash = false;
      crtHits.push_back(crthit);
    }
    i_hit += diff;
  } // Find valid CRT hits

  if (fUseMC)
  {
    for (auto &crthit : crtHits)
    {
      unsigned g4id = (unsigned)std::numeric_limits<int>::max();
      float mintime = std::numeric_limits<float>::max();
      float auxtime = -9999.;
      auto auxdethithandles = e.getMany<std::vector<sim::AuxDetHit>>();
      for (auto &auxdethithandle : auxdethithandles)
      {
        for (auto &auxdethit : *auxdethithandle)
        {
          if (abs((auxdethit.GetEntryT() / 1000.) - crthit.Time) < mintime)
          {
            mintime = abs((auxdethit.GetEntryT() / 1000.) - crthit.Time);
            auxtime = auxdethit.GetEntryT() / 1000.;
            g4id = auxdethit.GetTrackID();
          }
        }
      }
      crthit.G4ID = g4id;
      crthit.AuxTime = auxtime;
    } // Loop over crtHits
    art::Handle<std::vector<simb::MCParticle>> mctruthHandle;
    std::vector<art::Ptr<simb::MCParticle>> mctruthVect;
    if (e.getByLabel(fMCTruthLabel, mctruthHandle))
    { // Make sure artHandle is from module
      art::fill_ptr_vector(mctruthVect, mctruthHandle);
    }

    for (auto &crthit : crtHits)
    {
      for (auto &mcp : mctruthVect)
      {
        if (mcp->TrackId() == (int)crthit.G4ID)
          crthit.PDG = mcp->PdgCode();
      }
    }
  } // If fUseMC
}

std::vector<std::pair<float, float>> sbnd::LateLightProducer::findOpMuons(
    const art::Event &e,
    std::unique_ptr<std::vector<sbnd::StoppingMuonTrigger>> &stoppingmu_v)
{
  findCRTTimes(e);

  muonThreshold = fFindMuonThreshold;
  michelThreshold = fFindMichelPeak;
  auto waveforms = e.getValidHandle<std::vector<raw::OpDetWaveform>>(fFinderInputLabel);
  float startTime = std::numeric_limits<float>::max();
  float endTime = -std::numeric_limits<float>::max();
  float readoutDelay = 0.;
  // First loop: Determine the bounds of the cumulative waveform
  for (auto const &waveform : *waveforms)
  {
    float waveformStart = waveform.TimeStamp();
    float waveformEnd = waveformStart + waveform.Waveform().size() / 500.0; // Example frequency of 500 Hz
    auto opChannel = waveform.ChannelNumber();
    auto channelType = fPDMap.pdType(opChannel);
    if (channelType == "pmt_coated" || channelType == "pmt_uncoated")
      readoutDelay = fPMTReadoutDelay;
    else
      readoutDelay = fARAReadoutDelay;
   

    startTime = std::min(startTime, waveformStart);
    endTime = std::max(endTime, waveformEnd);
  }
  startTime -= readoutDelay;
  endTime -= readoutDelay;
  if (fVerbose)
    std::cout << "StartTime: " << startTime << "    End Time: " << endTime << "\n";
  wvfStartTime = startTime;
  wvfEndTime = endTime;
  wvfCRTHits.clear();
  for (auto crthit : crtHits)
  {
    if (crthit.Time < endTime && crthit.Time > startTime)
    {
      wvfCRTHits.push_back(crthit);
      if (fUseMC)
      {
        std::vector<std::pair<int, float>> blank_vec;
        if (muonMultCoatMap.find(crthit.PDG) == muonMultCoatMap.end())
        {
          muonMultCoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          muonMultUncoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          michelMultCoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          michelMultUncoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
        }
      }
    }
  }
  // Initialize the cumulative waveform with appropriate size based on the calculated bounds
  size_t waveformSize = static_cast<size_t>((endTime - startTime) * 500); // Example frequency of 500 Hz
  std::vector<float> cumulativeWaveform(waveformSize, 0);
  nSamples = waveformSize;

  auto pmttriggerHandle = e.getValidHandle<std::vector<sbnd::comm::pmtTrigger>>(fPMTTriggerLabel);
  auto pmtTrigger = (*pmttriggerHandle)[0];
  // Second loop: Process and accumulate the waveforms
  for (auto const &waveform : *waveforms)
  {
    std::vector<short> raw_wf = waveform.Waveform();
    std::vector<float> wf(raw_wf.size(), 0.);
    std::transform(raw_wf.begin(), raw_wf.end(), wf.begin(),
                   [](short adc)
                   { return (float)adc; });
    auto wf_corrected = subtractBaseline(wf, 2.);
    unsigned int startBin = static_cast<unsigned int>((waveform.TimeStamp() - (startTime + readoutDelay)) * 500);
    for (size_t i = 0; i < wf.size(); ++i)
    {
      if (startBin + i < cumulativeWaveform.size())
      {
        cumulativeWaveform[startBin + i] += fFinderPolarity * wf_corrected[i];
      }
    }
  }
  auto rollingSum = applyRollingSum(cumulativeWaveform);
  std::vector<std::pair<size_t, float>> peakMap;
  findPeaks(rollingSum, cumulativeWaveform, peakMap);
  std::vector<std::pair<size_t, size_t>> peakPairs;
  findPeakPairs(peakMap, false, peakPairs);
  std::vector<std::pair<float, float>> peakPairTimes;
  for (auto peakpair : peakPairs)
  {
    auto mult_win = (int)(peakpair.first / 4.);
    auto mult_win_min = std::max(0, mult_win - 30);
    auto mult_win_max = std::min((int)pmtTrigger.numPassed.size(), mult_win + 30);
    auto muon_mult = *max_element(pmtTrigger.numPassed.begin() + mult_win_min, pmtTrigger.numPassed.begin() + mult_win_max);
    if (muon_mult < fMinPMTMultiplicity)
      continue;
    float muontime = (float)peakpair.first / 500. + startTime;
    float micheltime = (float)peakpair.second / 500. + startTime;
    peakPairTimes.push_back(std::make_pair(muontime, micheltime));
    StoppingMuonTrigger stoppingmu_trigger;
    stoppingmu_trigger.MuonTime = muontime;
    stoppingmu_trigger.MichelTime = micheltime;
    stoppingmu_trigger.MuonMult = muon_mult;
    for (auto &crthit : crtHits)
    {
      if (abs(crthit.Time - muontime) * 1000. < fCRTClockSpeed)
        crthit.HasFlash = true;
      stoppingmu_trigger.CRTPlane = crthit.Plane;
      if (fUseMC)
      {
        stoppingmu_trigger.G4ID = crthit.G4ID;
        stoppingmu_trigger.G4PDG = crthit.PDG;
      }
    }
    stoppingmu_v->push_back(stoppingmu_trigger);

    if (fUseMC)
    {
      for (auto &muon : muon_tuple_vect)
      {
        if (abs(muontime - muon.Time) < 0.2)
        {
          muon.HasFlash = true;
          muon.FlashTime = muontime;
          muon.MichelFlashTime = micheltime;
        }
      }
    }
  }
  return peakPairTimes;
}

void sbnd::LateLightProducer::addPeakToMap(std::map<int, std::vector<std::pair<int, float>>> &multMap, int g4id, int opChannel, int pair_channel, float peakAmp, bool requirePositive)
{
  // Ensure the peak is valid (for Michel peaks, requirePositive enforces that the peakAmp > 0)
  if (requirePositive && peakAmp <= 0)
    return;

  // Check if the pair channel already exists in the map
  auto &mapEntry = multMap[g4id];
  auto it = std::find_if(mapEntry.begin(), mapEntry.end(),
                         [pair_channel](const std::pair<int, float> &opamp)
                         { return opamp.first == pair_channel; });

  // If paired peak does not exist, add it
  if (it == mapEntry.end())
  {
    mapEntry.push_back(std::make_pair(opChannel, peakAmp));
  }
  // If paired peak exists but current peak is lower, replace it
  else if (peakAmp < it->second)
  {
    it->first = opChannel;
    it->second = peakAmp;
  }
}

std::pair<std::vector<float>, std::vector<float>> sbnd::LateLightProducer::binAndFit(const std::vector<float> &data, size_t peakIndex, int channel)
{
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  size_t startBin = peakIndex + static_cast<size_t>(fStartTimeAfterMuon / conversionFactor);
  size_t endBin = peakIndex + static_cast<size_t>(fEndTimeAfterMuon / conversionFactor);
  size_t binWidthSamples = static_cast<size_t>(fTimeBinWidth / conversionFactor);
  size_t nBins = (fBinWaveform) ? (endBin - startBin) / binWidthSamples : endBin - startBin;

  std::vector<float> bins(nBins);
  std::vector<float> binValues(nBins);

  if (fBinWaveform)
  {
    for (size_t i = 0; i < nBins; ++i)
    {
      float binValue = std::accumulate(data.begin() + startBin + i * binWidthSamples,
                                       data.begin() + startBin + (i + 1) * binWidthSamples, 0.0) /
                       binWidthSamples;
      bins[i] = i * fTimeBinWidth;
      binValues[i] = binValue;
    }
  }
  else
  {
    for (unsigned i = startBin; i < endBin; i++)
    {
      binValues[i - startBin] = data[i];
    }
  }

  // Fit an exponential decay
  float binWidth = (fBinWaveform) ? fTimeBinWidth : conversionFactor;
  TH1D hist("hist", "hist", nBins, 0, nBins * binWidth);
  for (size_t i = 0; i < nBins; ++i)
  {
    hist.SetBinContent(i + 1, binValues[i]);
  }

  TF1 fit("fit", "[0]*exp(-x/[1])", 0, nBins * fTimeBinWidth);
  fit.SetParameters(hist.GetBinContent(0), 1000.);
  fit.SetParLimits(0, 0.5 * hist.GetBinContent(0), 2. * hist.GetBinContent(0));
  fit.SetParLimits(1, 200., 3000.);
  hist.Fit(&fit, "Q");

  float fitParam0 = fit.GetParameter(0);
  float fitParam1 = fit.GetParameter(1);

  return {binValues, {fitParam0, fitParam1}};
}

// Function to convert raw::OpDetWaveform to ROOT histogram
void sbnd::LateLightProducer::ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, int eventNumber)
{
  art::ServiceHandle<art::TFileService> tfs;
  // Retrieve necessary information from the waveform
  int channel = waveform.ChannelNumber();
  auto opType = fPDMap.pdType(fChannel);
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  int start_bin = (int)((fMuonPeakTime + 0.4 - waveform.TimeStamp()) * 1000. / conversionFactor);
  ;
  int end_bin = (int)((fMuonPeakTime + 2. - waveform.TimeStamp()) * 1000. / conversionFactor);
  // Get the channel number
  float timestamp = waveform.TimeStamp();                   // Get the timestamp
  const std::vector<short> &waveData = waveform.Waveform(); // Get the waveform data (amplitudes)
  // int numSamples = waveData.size();                         // Get the number of waveform samples
  int numSamples = end_bin - start_bin;

  // Create a unique histogram name using channel number, event number, and timestamp
  std::stringstream histName;
  histName << "waveform_ch" << channel << "_evt" << eventNumber << "_ts" << (size_t)timestamp * 1000;

  // Create the histogram with an appropriate binning (one bin per sample)
  TH1F *hist = tfs->make<TH1F>(histName.str().c_str(), histName.str().c_str(), numSamples, fMuonPeakTime + 0.4, fMuonPeakTime + 2.);

  // Fill the histogram with the waveform data
  for (int i = start_bin; i < numSamples; ++i)
  {
    hist->SetBinContent(i + 1 - start_bin, waveData[i]); // Set bin content (ROOT bins start at 1)
  }

  // Optionally, set axis labels (if needed)
  hist->GetXaxis()->SetTitle("Sample Number");
  hist->GetYaxis()->SetTitle("Amplitude");
}

DEFINE_ART_MODULE(sbnd::LateLightProducer)
