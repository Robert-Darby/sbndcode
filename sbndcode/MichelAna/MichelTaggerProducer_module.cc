////////////////////////////////////////////////////////////////////////
// Class:       MichelTaggerProducer
// Module Type: producer
// File:        MichelTaggerProducer_module.cc
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
#include "lardata/Utilities/AssociationUtil.h"

#include "canvas/Persistency/Common/Assns.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
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

namespace sbnd
{
  class MichelTaggerProducer;
}

class sbnd::MichelTaggerProducer : public art::EDProducer
{
public:
  explicit MichelTaggerProducer(fhicl::ParameterSet const &p);
  void produce(art::Event &e) override;

  // Operators
  MichelTaggerProducer(MichelTaggerProducer const &) = delete;
  MichelTaggerProducer(MichelTaggerProducer &&) = delete;
  MichelTaggerProducer &operator=(MichelTaggerProducer const &) = delete;
  MichelTaggerProducer &operator=(MichelTaggerProducer &&) = delete;

  struct MCMuon
  {
    int G4ID;
    float StartX, StartY, StartZ;
    float Time, RawAmp, SADCWAmp;
    int Mult, MichelMult;
    int MichelG4ID;
    float MichelStartX, MichelStartY, MichelStartZ;
    float MichelEnergy, MichelTime, MichelRawAmp, MichelSADCWAmp;
    long MichelDepPE;
    // float MichelTime, MichelEnergy;
    float FlashTime, MichelFlashTime;
    bool HasFlash, Stopping, EntersTPC, HasWvf;
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
  std::vector<MCMuon> muon_tuple_vect;
  std::vector<std::pair<float, float>> opMuons;
  std::vector<CRTHit> crtHits, wvfCRTHits;
  std::vector<float> muonTimes, michelTimes;
  float wvfStartTime, wvfEndTime;
  float muonThreshold, michelThreshold;
  unsigned nSamples;

  // Other variables shared between different methods.
  geo::GeometryCore const *fGeometryService;
  std::map<int, std::vector<std::pair<int, float>>> muonMultCoatMap, michelMultCoatMap, muonMultUncoatMap, michelMultUncoatMap;
  std::map<int, std::vector<std::pair<int, float>>> g4muonMultCoatMap, g4michelMultCoatMap, g4muonMultUncoatMap, g4michelMultUncoatMap;

  TTree *fTree;

  int fChannel;
  int fEventID;
  int fRun;
  int fSubRun;
  int fOpChannel, fOpChannelType;
  float fMuonPeakTime;
  float fMuonPeakAmp;
  int fMichelID;
  float fMCMichelStartX, fMCMichelStartY, fMCMichelStartZ;
  float fMichelTime, fMichelEnergy;
  float fMichelPeakAmp;
  bool fHasMichelPeak;
  float fMichelPeakTime;
  long fMichelDepPE;
  float fMuonFitParam0;
  float fMuonFitParam1;
  float fMichelFitParam0;
  float fMichelFitParam1;
  std::vector<float> fBinnedWaveformMuon;
  std::vector<float> fBinnedWaveformMichel;
  std::vector<float> fRawWaveform;
  int fMCMuonG4ID;
  float fMCMuonTime;
  float fMCMuonStartX, fMCMuonStartY, fMCMuonStartZ;
  bool fMCMuonStopping;
  bool fMCMuonEntersTPC;
  bool fMCMuonHasFlash;
  float fMuonFlashTime, fMichelFlashTime, fCRTTime;
  int fMuonMult;

  TTree *fTagTree;

  TTree *fTriggerTree;
  int fTriggerG4ID, fTriggerPDG, fTriggerPlane;
  float fMuonRawAmp, fMichelRawAmp, fMuonSADCWAmp, fMichelSADCWAmp;
  float fMCMuonRawAmp, fMCMuonSADCWAmp;
  int fMCMichelMult;
  int fMCMuonMult;
  float fMCMichelSADCWAmp, fMCMichelRawAmp;
  bool fTriggerHasFlash;
  bool fMCMuonHasWvf;
  std::vector<float> fMuonMultCoat, fMichelMultCoat, fMuonMultUncoat, fMichelMultUncoat;

  int findChannelPair(int opChannel);
  void findMCMuons(const art::Event &e);
  std::vector<float> CalcRunningAvg(std::vector<float> &wvf);
  // void GaussianSmoothing(std::vector<float> &Baseline);
  std::vector<float> subtractBaseline(const std::vector<float> &waveform, const float conversionFactor);
  std::vector<float> applyRollingSum(const std::vector<float> &waveform);
  void findPeaks(const std::vector<float> &seazrchWaveform, const std::vector<float> &waveform, std::vector<std::pair<size_t, float>> &peakIndices);
  void findPeakPairs(const std::vector<std::pair<size_t, float>> &peakIndices, const bool use_opMuons, std::vector<std::pair<size_t, size_t>> &peakPairs, int channel);
  void findCRTTimes(const art::Event &e);
  std::vector<std::pair<float, float>> findOpMuons(art::Event &e, std::unique_ptr<std::vector<sbnd::MichelTag>> &micheltag_v, std::unique_ptr<art::Assns<recob::OpFlash, sbnd::MichelTag>> &micheltag_opflash_assn_v);
  void addPeakToMap(std::map<int, std::vector<std::pair<int, float>>> &multMap, int g4id, int opChannel, int pair_channel, float peakAmp, bool requirePositive);
  std::pair<std::vector<float>, std::vector<float>> binAndFit(const std::vector<float> &data, size_t peakIndex, int channel);
  void ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, int eventNumber);
  void saveWaveformToHistogram(const std::vector<float> &waveform, int eventNumber, std::string suffix);
};

// Constructor
sbnd::MichelTaggerProducer::MichelTaggerProducer(fhicl::ParameterSet const &p)
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
  produces<std::vector<sbnd::MichelTag>>();
  produces<std::vector<raw::OpDetWaveform>>();
  produces<art::Assns<recob::OpFlash, sbnd::MichelTag>>();

  art::ServiceHandle<art::TFileService> tfs;

  fTagTree = tfs->make<TTree>("TagTree", "sbnd::MichelTagInfo");
  fTagTree->Branch("Run", &fRun);
  fTagTree->Branch("Subrun", &fSubRun);
  fTagTree->Branch("Event", &fEventID);
  fTagTree->Branch("recomuon.time", &fMuonFlashTime);
  fTagTree->Branch("recomuon.rawamp", &fMuonRawAmp);
  fTagTree->Branch("recomichel.rawamp", &fMichelRawAmp);
  fTagTree->Branch("recomuon.sadcwamp", &fMuonSADCWAmp);
  fTagTree->Branch("recomichel.sadcwamp", &fMichelSADCWAmp);
  fTagTree->Branch("recomichel.time", &fMichelFlashTime);
  fTagTree->Branch("crt.plane", &fTriggerPlane);
  fTagTree->Branch("crt.time", &fCRTTime);
  fTagTree->Branch("recomuon.mult", &fMuonMult);
  fTagTree->Branch("recomuon.multcoat", &fMuonMultCoat);
  fTagTree->Branch("recomichel.multcoat", &fMichelMultCoat);
  fTagTree->Branch("recomuon.multuncoat", &fMuonMultUncoat);
  fTagTree->Branch("recomichel.multuncoat", &fMichelMultUncoat);

  // Set branches
  fTriggerTree = tfs->make<TTree>("MCTagTree", "Waveform Producers Tree");
  fTriggerTree->Branch("Run", &fRun);
  fTriggerTree->Branch("Subrun", &fSubRun);
  fTriggerTree->Branch("Event", &fEventID);
  fTriggerTree->Branch("mcmuon.g4id", &fMCMuonG4ID);
  fTriggerTree->Branch("mcmuon.start_x", &fMCMuonStartX);
  fTriggerTree->Branch("mcmuon.start_y", &fMCMuonStartY);
  fTriggerTree->Branch("mcmuon.start_z", &fMCMuonStartZ);
  fTriggerTree->Branch("mcmichel.start_x", &fMCMichelStartX);
  fTriggerTree->Branch("mcmichel.start_y", &fMCMichelStartY);
  fTriggerTree->Branch("mcmichel.start_z", &fMCMichelStartZ);
  fTriggerTree->Branch("mcmuon.time", &fMCMuonTime);
  fTriggerTree->Branch("mcmichel.g4id", &fMichelID);
  fTriggerTree->Branch("mcmichel.energy", &fMichelEnergy);
  fTriggerTree->Branch("mcmichel.time", &fMichelTime);
  fTriggerTree->Branch("mcmichel.dep_pe", &fMichelDepPE);
  fTriggerTree->Branch("mcmuon.is_stopping", &fMCMuonStopping);
  fTriggerTree->Branch("mcmuon.enters_tpc", &fMCMuonEntersTPC);
  fTriggerTree->Branch("mcmuon.has_wvf", &fMCMuonHasWvf);
  fTriggerTree->Branch("hasTag", &fMCMuonHasFlash);
  fTriggerTree->Branch("recomuon.time", &fMuonFlashTime);
  fTriggerTree->Branch("recomuon.mult", &fMuonMult);
  fTriggerTree->Branch("recomuon.rawamp", &fMuonRawAmp);
  fTriggerTree->Branch("recomichel.rawamp", &fMichelRawAmp);
  fTriggerTree->Branch("recomuon.sadcwamp", &fMuonSADCWAmp);
  fTriggerTree->Branch("recomichel.sadcwamp", &fMichelSADCWAmp);
  fTriggerTree->Branch("mcmuon.rawamp", &fMCMuonRawAmp);
  fTriggerTree->Branch("mcmuon.sadcwamp", &fMCMuonSADCWAmp);
  fTriggerTree->Branch("mcmichel.rawamp", &fMCMichelRawAmp);
  fTriggerTree->Branch("mcmichel.rawamp", &fMCMichelSADCWAmp);
  fTriggerTree->Branch("recomichel.time", &fMichelFlashTime);
  fTriggerTree->Branch("crt.g4id", &fTriggerG4ID);
  fTriggerTree->Branch("crt.pdg", &fTriggerPDG);
  fTriggerTree->Branch("crt.plane", &fTriggerPlane);
  fTriggerTree->Branch("crt.time", &fCRTTime);
  fTriggerTree->Branch("recomuon.multcoat", &fMuonMultCoat);
  fTriggerTree->Branch("recomichel.multcoat", &fMichelMultCoat);
  fTriggerTree->Branch("recomuon.multuncoat", &fMuonMultUncoat);
  fTriggerTree->Branch("recomichel.multuncoat", &fMichelMultUncoat);
}

void sbnd::MichelTaggerProducer::produce(art::Event &e)
{
  std::unique_ptr<std::vector<sbnd::MichelTag>>
      micheltag_v(new std::vector<sbnd::MichelTag>);
  std::unique_ptr<std::vector<raw::OpDetWaveform>>
      michelwvfms_v(std::make_unique<std::vector<raw::OpDetWaveform>>());
  std::unique_ptr<art::Assns<recob::OpFlash, sbnd::MichelTag>>
      micheltag_opflash_assn_v(new art::Assns<recob::OpFlash, sbnd::MichelTag>);

  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();
  muon_tuple_vect.clear();

  muonMultCoatMap.clear();
  muonMultUncoatMap.clear();
  michelMultCoatMap.clear();
  michelMultUncoatMap.clear();
  g4muonMultCoatMap.clear();
  g4muonMultUncoatMap.clear();
  g4michelMultCoatMap.clear();
  g4michelMultUncoatMap.clear();
  fMuonMultCoat.clear();
  fMuonMultUncoat.clear();
  fMichelMultCoat.clear();
  fMichelMultUncoat.clear();

  if (fUseMC)
  {
    findMCMuons(e);
    for (const auto &mcmuon : muon_tuple_vect)
    {
      if (g4muonMultCoatMap.find(mcmuon.G4ID) == g4muonMultCoatMap.end())
      {
        std::vector<std::pair<int, float>> blank_vec;
        g4muonMultCoatMap.insert(std::make_pair(mcmuon.G4ID, blank_vec));
        g4michelMultCoatMap.insert(std::make_pair(mcmuon.G4ID, blank_vec));
        g4muonMultUncoatMap.insert(std::make_pair(mcmuon.G4ID, blank_vec));
        g4michelMultUncoatMap.insert(std::make_pair(mcmuon.G4ID, blank_vec));
      }
    }
  }
  opMuons.clear();
  opMuons = findOpMuons(e, micheltag_v, micheltag_opflash_assn_v);
  mf::LogInfo("MichelTag") << "Found " << opMuons.size() << "opmuons\n";
  for (const auto &opmuon : opMuons)
  {
    std::vector<std::pair<int, float>> blank_vec;
    muonMultCoatMap.insert(std::make_pair((int)opmuon.first * 1000, blank_vec));
    michelMultCoatMap.insert(std::make_pair((int)opmuon.first * 1000, blank_vec));
    muonMultUncoatMap.insert(std::make_pair((int)opmuon.first * 1000, blank_vec));
    michelMultUncoatMap.insert(std::make_pair((int)opmuon.first, blank_vec));
  }
  if (fVerbose)
  {
    mf::LogInfo("MichelTag") << "OpMuons: \n";
    for (auto muonmichel : opMuons)
    {
      mf::LogInfo("MichelTag") << "   " << muonmichel.first << "   " << muonmichel.second << "\n";
    }
  }
  if (fVerbose)
  {
    mf::LogInfo("MichelTag") << "MC Muons\n";
    for (auto &muon : muon_tuple_vect)
    {
      mf::LogInfo("MichelTag") << "    " << muon.G4ID << "   " << muon.Time << "   " << muon.MichelTime << "   " << muon.MichelDepPE;
    }
  }

  if (fProduceWaveforms)
  {
    auto handle = e.getHandle<std::vector<raw::OpDetWaveform>>(fInputLabel);

    if (!handle)
      mf::LogError("MichelTag") << "Invalid handle name: " << fInputLabel << "\n";

    for (const auto &waveform : *handle)
    {
      fChannel = waveform.ChannelNumber();
      int pair_channel = findChannelPair(fChannel);
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
      auto correctedWaveform = subtractBaseline(waveformData, conversionFactor);
      // Apply rolling sum to the waveform
      auto smooth_wvf = CalcRunningAvg(correctedWaveform);
      auto rollingSumWaveform = applyRollingSum(correctedWaveform);

      nSamples = waveformData.size();

      // Find if a muon and michel time is within the waveform, skip if not
      wvfEndTime = waveform.TimeStamp() + nSamples * conversionFactor / 1000.;
      wvfStartTime = waveform.TimeStamp();
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
        }

        fRawWaveform.clear();
        if (fHasMichelPeak && fMichelPeakTime - fMuonPeakTime > 2.)
        {
          int startBin = (muonPeakIndex < 50) ? 0 : muonPeakIndex - 50;
          int endBin = (muonPeakIndex + 4500 > nSamples) ? nSamples : muonPeakIndex + 4500;
          for (unsigned i = (unsigned)startBin; i < (unsigned)endBin; i++)
            fRawWaveform.push_back(correctedWaveform[i]);
        }
        for (auto micheltag : *micheltag_v)
        {
          if (abs(micheltag.MuonTime - fMuonPeakTime) > 0.05)
            continue;
          int muontime = (int)micheltag.MuonTime * 1000;
          if (opType == "pmt_coated")
          {
            addPeakToMap(muonMultCoatMap, muontime, fChannel, pair_channel, fMuonPeakAmp, false);
            addPeakToMap(michelMultCoatMap, muontime, fChannel, pair_channel, fMichelPeakAmp, true);
          }
          else if (opType == "pmt_uncoated")
          {
            addPeakToMap(muonMultUncoatMap, muontime, fChannel, pair_channel, fMuonPeakAmp, false);
            addPeakToMap(michelMultUncoatMap, muontime, fChannel, pair_channel, fMichelPeakAmp, true);
          } // If PMT uncoated
        } // Loop over sbnd::MichelTag vector

        fMCMuonG4ID = -9999;
        fMCMuonTime = std::numeric_limits<float>::max();
        fMCMuonStartX = std::numeric_limits<float>::max();
        fMCMuonStartY = std::numeric_limits<float>::max();
        fMCMuonStartZ = std::numeric_limits<float>::max();
        fMCMuonStopping = false;
        fMCMuonEntersTPC = false;
        fMCMuonHasWvf = false;
        fMCMuonHasFlash = false;
        fMuonFlashTime = -9999.;
        fMuonMult = -1;
        fMichelFlashTime = -9999.;
        fTriggerG4ID = -1;
        fMuonRawAmp = -9999.;
        fMichelRawAmp = -9999.;
        fMuonSADCWAmp = -9999.;
        fMichelSADCWAmp = -9999.;
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
                fMCMuonHasWvf = mcmuon.HasWvf;
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

        // Fill the tree
        if (fTriggerG4ID < 0)
          continue;

        if (opType == "pmt_coated")
        {
          addPeakToMap(g4muonMultCoatMap, fTriggerG4ID, fChannel, pair_channel, fMuonPeakAmp, false);
          addPeakToMap(g4michelMultCoatMap, fTriggerG4ID, fChannel, pair_channel, fMichelPeakAmp, true);
        }
        else if (opType == "pmt_uncoated")
        {
          addPeakToMap(g4muonMultUncoatMap, fTriggerG4ID, fChannel, pair_channel, fMuonPeakAmp, false);
          addPeakToMap(g4michelMultUncoatMap, fTriggerG4ID, fChannel, pair_channel, fMichelPeakAmp, true);
        } // If PMT uncoated
        if (fTriggerG4ID > 0 && fSaveHists && !peakPairs.empty() &&
            std::find(fSaveEvents.begin(), fSaveEvents.end(), fEventID) != fSaveEvents.end() &&
            fMichelPeakTime - fMuonPeakTime > 0.1)
        {
          ConvertWaveformToHistogram(waveform, fEventID);
        } // Save histograms
      } // Loop over peak pairs
    } // Loop over waveforms
  } // Loop over labels

  // Add muon trigger mults to MichelTags
  for (auto &micheltag : *micheltag_v)
  {
    int muontime = (int)micheltag.MuonTime * 1000;
    auto muoncoat_v = muonMultCoatMap[muontime];
    std::vector<float> multcoat;
    for (auto adc : muoncoat_v)
      multcoat.push_back(adc.second);
    std::sort(multcoat.begin(), multcoat.end());
    micheltag.MuonMultCoat = multcoat;

    auto muonuncoat_v = muonMultUncoatMap[muontime];
    std::vector<float> multuncoat;
    for (auto adc : muonuncoat_v)
      multuncoat.push_back(adc.second);
    std::sort(multuncoat.begin(), multuncoat.end());
    micheltag.MuonMultUncoat = multuncoat;
  }

  // Fill Trigger Tree
  if (fUseMC)
  {
    // Load sim::SimEnergyDeposits
    /* art::Handle<std::vector<sim::SimEnergyDeposit>> sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit>> sedVect;
    if (e.getByLabel(fSEDLabel, sedHandle)) // Make sure artHandle is from module
      art::fill_ptr_vector(sedVect, sedHandle);

    */
    for (auto it = g4muonMultCoatMap.begin(); it != g4muonMultCoatMap.end(); ++it)
    {
      auto g4id = it->first;
      fMuonMultCoat.clear();
      fMuonMultUncoat.clear();
      fMichelMultCoat.clear();
      fMichelMultUncoat.clear();

      if (g4muonMultCoatMap.find(g4id) != g4muonMultCoatMap.end())
      {
        for (auto muon_opamp : g4muonMultCoatMap[g4id])
          fMuonMultCoat.push_back(muon_opamp.second);
      }
      if (g4michelMultCoatMap.find(g4id) != g4michelMultCoatMap.end())
      {
        for (auto michel_opamp : g4michelMultCoatMap[g4id])
          fMichelMultCoat.push_back(michel_opamp.second);
      }
      if (g4muonMultUncoatMap.find(g4id) != g4muonMultUncoatMap.end())
      {
        for (auto muon_opamp : g4muonMultUncoatMap[g4id])
          fMuonMultUncoat.push_back(muon_opamp.second);
      }
      if (g4michelMultUncoatMap.find(g4id) != g4michelMultUncoatMap.end())
      {
        for (auto michel_opamp : g4michelMultUncoatMap[g4id])
          fMichelMultUncoat.push_back(michel_opamp.second);
      }

      if (fVerbose)
      {
        mf::LogInfo("MichelTag") << "\nFound " << fMuonMultCoat.size() << " coated muon peaks for " << g4id << "\n"
                                 << "Found " << fMichelMultCoat.size() << " coated michel peaks for " << g4id << "\n"
                                 << "Found " << fMuonMultUncoat.size() << " uncoated muon peaks for " << g4id << "\n"
                                 << "Found " << fMichelMultUncoat.size() << " uncoated michel peaks for " << g4id << "\n";
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
      fMCMuonStartX = -9999.;
      fMCMuonStartY = -9999.;
      fMCMuonStartZ = -9999.;
      fMCMichelStartX = -9999.;
      fMCMichelStartY = -9999.;
      fMCMichelStartZ = -9999.;

      fMCMuonHasFlash = false;
      fMCMuonStopping = false;
      fMCMuonEntersTPC = false;
      fMichelEnergy = -9999.;
      fMichelTime = -9999.;
      fMichelID = false;
      fMCMuonTime = -9999.;
      fMCMuonRawAmp = -9999.;
      fMCMuonSADCWAmp = -9999.;
      fMCMichelRawAmp = -9999.;
      fMCMichelSADCWAmp = -9999.;

      for (auto &muon_tuple : muon_tuple_vect)
      {
        if (muon_tuple.G4ID == g4id)
        {
          fMCMuonG4ID = g4id;
          fMCMuonTime = muon_tuple.Time;
          fMCMuonStartX = muon_tuple.StartX;
          fMCMuonStartY = muon_tuple.StartY;
          fMCMuonStartZ = muon_tuple.StartZ;

          fMichelID = muon_tuple.MichelG4ID;
          fMichelTime = muon_tuple.MichelTime;
          fMCMichelStartX = muon_tuple.MichelStartX;
          fMCMichelStartY = muon_tuple.MichelStartY;
          fMCMichelStartZ = muon_tuple.MichelStartZ;
          fMichelEnergy = muon_tuple.MichelEnergy;
          fMichelDepPE = muon_tuple.MichelDepPE;
          fMCMuonStopping = muon_tuple.Stopping;
          fMCMuonEntersTPC = muon_tuple.EntersTPC;
          fMCMuonHasFlash = muon_tuple.HasFlash;
          fMCMuonRawAmp = muon_tuple.RawAmp;
          fMCMichelRawAmp = muon_tuple.MichelRawAmp;
          fMCMuonSADCWAmp = muon_tuple.SADCWAmp;
          fMCMichelSADCWAmp = muon_tuple.MichelSADCWAmp;
        }
      }

      // Fill CRT Info
      fTriggerG4ID = -1;
      fTriggerPDG = -1;
      fTriggerPlane = -1;

      fMuonFlashTime = -9999.;
      fMichelFlashTime = -9999.;
      fMuonRawAmp = -9999.;
      fMuonSADCWAmp = -9999.;
      fMichelRawAmp = -9999.;
      fMichelSADCWAmp = -9999.;
      fMuonMult = -1;
      for (auto &crthit : crtHits)
      {
        if (crthit.G4ID == (unsigned)g4id)
        {
          fTriggerG4ID = crthit.G4ID;
          fTriggerPDG = crthit.PDG;
          fTriggerPlane = crthit.Plane;
          fCRTTime = crthit.Time;
          fTriggerHasFlash = crthit.HasFlash;
          if (crthit.HasFlash)
          {
            for (const auto &micheltag : *micheltag_v)
            {
              if (crthit.Plane == micheltag.CRTPlane &&
                  abs(crthit.Time - micheltag.MuonTime) < fCRTCoincidence / 1000.)
              {
                fMuonFlashTime = micheltag.MuonTime;
                fMichelFlashTime = micheltag.MichelTime;
                fMuonRawAmp = micheltag.MuonRawAmp;
                fMichelRawAmp = micheltag.MichelRawAmp;
                fMuonSADCWAmp = micheltag.MuonSADCWAmp;
                fMichelSADCWAmp = micheltag.MichelSADCWAmp;
                fMuonMult = micheltag.MuonMult;
              }
            }
          }
        }
      }
      fTriggerTree->Fill();
    } // Loop over g4muonMultCoatMap
  } // If UseMC

  for (const auto &micheltag : *micheltag_v)
  {
    mf::LogInfo("MichelTagger") << "Found Michel tag:" << "\n"
                                << "  Muon Time   : " << micheltag.MuonTime << "\n"
                                << "  Michel time : " << micheltag.MichelTime << "\n"
                                << "  CRT Plane   : " << micheltag.CRTPlane << "\n"
                                << "  Coat mult   : " << micheltag.MuonMultCoat.size() << "\n"
                                << "  Uncoat mult : " << micheltag.MuonMultUncoat.size() << "\n"
                                << "  Muon Raw Amp: " << micheltag.MuonRawAmp << "\n"
                                << "  Michel RawAmp: " << micheltag.MichelRawAmp << "\n"
                                << "  Muon SADCW   : " << micheltag.MuonSADCWAmp << "\n"
                                << "  Michel SADCW : " << micheltag.MichelSADCWAmp << "\n";
  }

  if (fMakeTree)
  {
    for (const auto &tag : *micheltag_v)
    {
      fMuonFlashTime = tag.MuonTime;
      fMichelFlashTime = tag.MichelTime;
      fTriggerPlane = tag.CRTPlane;
      fMuonRawAmp = tag.MuonRawAmp;
      fMuonSADCWAmp = tag.MuonSADCWAmp;
      fMichelRawAmp = tag.MichelRawAmp;
      fMichelSADCWAmp = tag.MichelSADCWAmp;
      fMuonMultCoat = tag.MuonMultCoat;
      fMuonMultUncoat = tag.MuonMultUncoat;
      fMichelMultCoat.clear();
      fMichelMultUncoat.clear();
      fTagTree->Fill();
    }
  }

  e.put(std::move(micheltag_v));
  mf::LogInfo("MichelTag") << "Saving " << michelwvfms_v->size() << " Michel waveforms";
  e.put(std::move(michelwvfms_v));
  e.put(std::move(micheltag_opflash_assn_v));
}

int sbnd::MichelTaggerProducer::findChannelPair(int opChannel)
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

void sbnd::MichelTaggerProducer::findMCMuons(const art::Event &e)
{
  art::Handle<std::vector<simb::MCParticle>> mctruthHandle;
  std::vector<art::Ptr<simb::MCParticle>> mctruthVect;
  if (e.getByLabel(fMCTruthLabel, mctruthHandle)) // Make sure artHandle is from module
    art::fill_ptr_vector(mctruthVect, mctruthHandle);

  art::Handle<std::vector<sim::SimEnergyDeposit>> sedHandle;
  std::vector<art::Ptr<sim::SimEnergyDeposit>> sedVect;
  if (e.getByLabel(fSEDLabel, sedHandle)) // Make sure artHandle is from module
    art::fill_ptr_vector(sedVect, sedHandle);

  art::Handle<std::vector<sim::SimEnergyDeposit>> sedoutHandle;
  std::vector<art::Ptr<sim::SimEnergyDeposit>> sedoutVect;
  if (e.getByLabel(fSEDOutLabel, sedoutHandle)) // Make sure artHandle is from module
    art::fill_ptr_vector(sedoutVect, sedoutHandle);

  for (const auto &mcp : mctruthVect)
  {
    if (abs(mcp->PdgCode()) != 13 || abs(mcp->T()) / 1000. > 1510.)
      continue;
    auto mu_id = mcp->TrackId();
    bool mu_decayintpc = (abs(mcp->EndX()) < 185. &&
                          abs(mcp->EndY()) < 185. &&
                          mcp->EndZ() < 485. && mcp->EndZ() > 15.);
    float mu_time = std::numeric_limits<float>::max();
    float mu_x = std::numeric_limits<float>::max();
    float mu_y = std::numeric_limits<float>::max();
    float mu_z = std::numeric_limits<float>::max();

    bool mu_enterstpc = false;
    for (unsigned i = 0; i < mcp->NumberTrajectoryPoints(); i++)
    {
      if (abs(mcp->Position(i).X()) < 200. &&
          abs(mcp->Position(i).Y()) < 200. &&
          mcp->Position(i).Z() < 500. && mcp->Position(i).Z() > 0.)
      {
        mu_time = mcp->T(i) / 1000.;
        mu_enterstpc = true;
        mu_x = mcp->Position(i).X();
        mu_y = mcp->Position(i).Y();
        mu_z = mcp->Position(i).Z();
        break;
      }
    }
    int michel_id = -1;
    float michel_time = -9999.;
    float michel_energy = -9999.;
    long michel_dep_pe = 0;
    float michel_x = std::numeric_limits<float>::max();
    float michel_y = std::numeric_limits<float>::max();
    float michel_z = std::numeric_limits<float>::max();

    if (mcp->EndProcess() == "Decay")
    {
      for (auto &mcp2 : mctruthVect)
      {
        if (mcp2->Mother() != mcp->TrackId() ||
            (mcp2->Position() - mcp->EndPosition()).Mag() > 5. ||
            abs(mcp2->PdgCode()) != 11)
          continue;
        michel_id = mcp2->TrackId();
        michel_time = mcp2->T() / 1000.;
        michel_x = mcp2->Position().X();
        michel_y = mcp2->Position().Y();
        michel_z = mcp2->Position().Z();
        michel_energy = mcp2->E() * 1000.;
        if (mu_decayintpc)
        {
          auto sed_part = std::partition(sedVect.begin(), sedVect.end(),
                                         [michel_id](art::Ptr<sim::SimEnergyDeposit> &sed)
                                         { return (sed->TrackID() == michel_id); });

          int n_sed = 0;
          for (auto sed_it = sedVect.begin(); sed_it != sed_part; sed_it++)
          {
            auto sed = *sed_it;
            michel_dep_pe += sed->NumPhotons();
            n_sed++;
          }
        }
        else
        {
          auto sed_part = std::partition(sedoutVect.begin(), sedoutVect.end(),
                                         [michel_id](art::Ptr<sim::SimEnergyDeposit> &sed)
                                         { return (sed->TrackID() == michel_id); });

          int n_sed = 0;
          for (auto sed_it = sedoutVect.begin(); sed_it != sed_part; sed_it++)
          {
            auto sed = *sed_it;
            michel_dep_pe += sed->NumPhotons();
            n_sed++;
          }
        } // if (michel_energy < 0.)
      } // Michel loop
    } // Check if muon decays
    //   continue;
    if (!mu_enterstpc)
    {
      mu_time = mcp->T() / 1000.;
      mu_x = mcp->Position().X();
      mu_y = mcp->Position().Y();
      mu_z = mcp->Position().Z();
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
    mcmuon.MichelDepPE = michel_dep_pe;
    mcmuon.MichelEnergy = michel_energy;
    mcmuon.Stopping = mu_decayintpc;
    mcmuon.EntersTPC = mu_enterstpc;
    mcmuon.StartX = mu_x;
    mcmuon.StartY = mu_y;
    mcmuon.StartZ = mu_z;
    mcmuon.MichelStartX = michel_x;
    mcmuon.MichelStartY = michel_y;
    mcmuon.MichelStartZ = michel_z;

    mcmuon.HasFlash = false;
    mcmuon.FlashTime = -9999.;
    mcmuon.MichelFlashTime = -9999.;
    mcmuon.HasWvf = false;
    muon_tuple_vect.push_back(mcmuon);
  }
  std::sort(muon_tuple_vect.begin(), muon_tuple_vect.end(),
            [](const MCMuon &muon1, const MCMuon &muon2)
            { return (muon2.Time > muon1.Time); });
}

std::vector<float> sbnd::MichelTaggerProducer::CalcRunningAvg(std::vector<float> &wvf)
{
  // int index =0;
  double Median = 0.;
  // Bury a nested loop in a thing Im not sure is any faster
  std::vector<float> Baseline(wvf.size());
  std::vector<float> smooth_wvf(wvf.size(), 0.);
  for (int i = 0; i < int(wvf.size()) - fRunningAvgSampleWidth; i++)
  {
    int EndIndex = i + fRunningAvgSampleWidth;
    double sum = 0;
    for (int j = i; j < EndIndex; j++)
    {
      if (wvf[j] < Median - 40)
        sum = sum + wvf[j]; // mask out pulses
      else
        sum = sum + Median;
    }
    Baseline[i] = sum / (EndIndex - i);
  } // CalcAvgBaseline
  // std::for_each(Baseline.begin(), Baseline.begin()+(wvf.size()-fRunningAvgSampleWidth), [wvf, &index, this](int& Val)
  //           {
  //           Val = std::accumulate(wvf.begin()+index, wvf.begin()+index+fRunningAvgSampleWidth, 0)/fRunningAvgSampleWidth;
  //           index = index+1;
  //           } );
  for (int i = int(wvf.size()) - fRunningAvgSampleWidth; i < int(wvf.size()); i++)
  {
    int EndIndex = i + fRunningAvgSampleWidth;
    if (EndIndex > int(wvf.size()))
      EndIndex = int(wvf.size());
    double sum = 0;
    for (int j = i; j < EndIndex; j++)
    {
      if (wvf[j] < Median - 40)
        sum = sum + wvf[j]; // mask out pulses
      else
        sum = sum + Median;
    }
    Baseline[i] = sum / (EndIndex - i);
  }
  for (unsigned i = 0; i < wvf.size(); i++)
    smooth_wvf[i] -= Baseline[i];
  // Baseline is all updated an can return
  return smooth_wvf;
}
/*
void sbnd::MichelTaggerProducer::GaussianSmoothing(std::vector<float> &Baseline)
{
  std::vector<int> Out(Baseline.size());
  std::vector<int> X(fGuassianConvlSize * 2 + 1);
  std::iota(X.begin(), X.end(), -fGuassianConvlSize);
  std::vector<double> Weights(fGuassianConvlSize * 2 + 1);
  double sum = 0;
  for (int i = 0; i < int(Weights.size()); i++)
  {
    Weights[i] = TMath::Exp(-TMath::Power(double(X[i]), 2.0) / (2 * TMath::Power(double(fGaussianConvlWidth), 2.0)));
    sum += Weights[i];
  }
  // Now do the convolution
  for (int i = fGuassianConvlSize + 1; i < int(Baseline.size()) - (fGuassianConvlSize + 1); i++)
  {
    double PointSum = 0;
    std::for_each(X.begin(), X.end(), [&PointSum, i, Weights, Baseline, this](int Index)
                  { PointSum += Baseline[i + Index] * Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  // Handle edges properly
  for (int i = 0; i < fGuassianConvlSize + 1; i++)
  {
    double PointSum = 0;
    std::for_each(X.begin() + fGuassianConvlSize - i, X.end(), [&PointSum, i, Weights, Baseline, this](int Index)
                  { PointSum += Baseline[i + Index] * Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  for (int i = int(Baseline.size()) - (fGuassianConvlSize + 1); i < int(Baseline.size()); i++)
  {
    double PointSum = 0;
    std::for_each(X.begin(), X.begin() + fGuassianConvlSize - (i - int(Baseline.size())), [&PointSum, i, Weights, Baseline, this](int Index)
                  { PointSum += Baseline[i + Index] * Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  // Finally copy out to baseline
  for (int i = 0; i < int(Baseline.size()); i++)
    Baseline[i] = Out[i] / sum;
} // Gaussian Smoothin
*/
std::vector<float> sbnd::MichelTaggerProducer::subtractBaseline(const std::vector<float> &waveform, const float conversionFactor)
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

std::vector<float> sbnd::MichelTaggerProducer::applyRollingSum(const std::vector<float> &waveform)
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

void sbnd::MichelTaggerProducer::findPeaks(
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

void sbnd::MichelTaggerProducer::findPeakPairs(
    const std::vector<std::pair<size_t, float>> &peakIndices,
    const bool use_opMuons,
    std::vector<std::pair<size_t, size_t>> &peakPairs,
    int channel = 6)
{
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  std::vector<std::pair<size_t, size_t>> muon_windows, michel_windows;
  size_t mu_coincidence_window, michel_coincidence_window;
  std::vector<float> muon_times, michel_times;
  if (!use_opMuons)
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
          (peakIndices[i_peak].second < fMuonMichelMaxRatio * mu_peak || use_opMuons))
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

void sbnd::MichelTaggerProducer::findCRTTimes(const art::Event &e)
{
  crtHits.clear();
  auto febData = e.getValidHandle<std::vector<sbnd::crt::FEBData>>(fCRTLabel);
  std::vector<std::pair<float, unsigned>> valid_hits;

  for (const auto &febdata : *febData)
  {
    auto mac5 = febdata.Mac5();
    auto plane = sbnd::crt::CRTCommonUtils::GetTaggerEnum(fCRTGeo.GetTaggerName(fCRTGeo.ChannelToStripName(mac5 * 32)));
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

std::vector<std::pair<float, float>> sbnd::MichelTaggerProducer::findOpMuons(
    art::Event &e,
    std::unique_ptr<std::vector<sbnd::MichelTag>> &micheltag_v,
    std::unique_ptr<art::Assns<recob::OpFlash, sbnd::MichelTag>> &micheltag_opflash_assn_v)
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
    if (std::find(fFinderOpTypes.begin(), fFinderOpTypes.end(), channelType) == fFinderOpTypes.end())
      continue;
    if (channelType == "pmt_coated" || channelType == "pmt_uncoated")
      readoutDelay = fPMTReadoutDelay;
    else
      readoutDelay = fARAReadoutDelay;

    startTime = std::min(startTime, waveformStart);
    endTime = std::max(endTime, waveformEnd);
  }

  if (fVerbose)
    mf::LogInfo("MichelTag") << "StartTime: " << startTime << "    End Time: " << endTime << "\n";
  wvfStartTime = startTime;
  wvfEndTime = endTime;
  wvfCRTHits.clear();
  for (auto crthit : crtHits)
  {
    if (crthit.Time < endTime - readoutDelay && crthit.Time > startTime - readoutDelay)
    {
      wvfCRTHits.push_back(crthit);
      if (fUseMC)
      {
        std::vector<std::pair<int, float>> blank_vec;
        if (g4muonMultCoatMap.find(crthit.G4ID) == g4muonMultCoatMap.end())
        {
          g4muonMultCoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          g4muonMultUncoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          g4michelMultCoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
          g4michelMultUncoatMap.insert(std::make_pair(crthit.G4ID, blank_vec));
        }
      }
    }
  }
  // Initialize the cumulative waveform with appropriate size based on the calculated bounds
  size_t waveformSize = static_cast<size_t>((endTime - startTime) * 500); // Example frequency of 500 Hz
  std::vector<float> cumulativeWaveform(waveformSize, 0);
  nSamples = waveformSize;

  auto pmttriggerHandle = e.getValidHandle<std::vector<sbnd::comm::pmtTrigger>>(fPMTTriggerLabel);
  if (pmttriggerHandle->size() == 0)
    mf::LogError("MichelTagger") << "Found no pmttriggers in producer: " << fPMTTriggerLabel << "\n";
  auto pmtTrigger = (*pmttriggerHandle)[0];

  // Second loop: Process and accumulate the waveforms
  for (auto const &waveform : *waveforms)
  {
    std::vector<short> raw_wf = waveform.Waveform();
    std::vector<float> wf(raw_wf.size(), 0.);
    std::transform(raw_wf.begin(), raw_wf.end(), wf.begin(),
                   [](short adc)
                   { return (float)adc; });
    unsigned int startBin = static_cast<unsigned int>((waveform.TimeStamp() - (startTime)) * 500);
    auto wf_corrected = subtractBaseline(wf, 2.);

    for (size_t i = 0; i < wf.size(); ++i)
    {
      if (startBin + i < cumulativeWaveform.size())
      {
        cumulativeWaveform[startBin + i] += fFinderPolarity * wf_corrected[i];
      }
    }
  }
  auto smooth_wvf = CalcRunningAvg(cumulativeWaveform);
  auto rollingSum = applyRollingSum(cumulativeWaveform);

  if (fUseMC)
  {
    for (auto &mcmuon : muon_tuple_vect)
    {
      if (!mcmuon.Stopping || mcmuon.MichelG4ID < 0)
        continue;
      if (mcmuon.Time < startTime || mcmuon.Time > endTime)
        continue;
      mcmuon.HasWvf = true;
      int muon_win = (int)(1000. * (mcmuon.Time - ((startTime - readoutDelay))) / (fCAENFreq * 4.));
      auto muon_win_min = std::max(0, muon_win - 30);
      auto muon_win_max = std::min((int)pmtTrigger.numPassed.size() - 1, muon_win + 30);
      mcmuon.Mult = *std::max_element(pmtTrigger.numPassed.begin() + muon_win_min, pmtTrigger.numPassed.begin() + muon_win_max);
      int michel_win = (int)(1000. * (mcmuon.MichelTime - ((startTime - readoutDelay))) / (fCAENFreq * 4.));
      auto michel_win_min = std::max(0, michel_win - 10);
      auto michel_win_max = std::min((int)pmtTrigger.numPassed.size() - 1, michel_win + 10);
      mcmuon.MichelMult = *std::max_element(pmtTrigger.numPassed.begin() + michel_win_min, pmtTrigger.numPassed.begin() + michel_win_max);

      muon_win = (int)(1000. * ((mcmuon.Time + readoutDelay) - startTime) / fCAENFreq);
      muon_win_min = std::max(muon_win - 20, 0);
      muon_win_max = std::min(muon_win + 20, (int)cumulativeWaveform.size() - 1);
      mcmuon.RawAmp = *std::max_element(cumulativeWaveform.begin() + muon_win_min, cumulativeWaveform.begin() + muon_win_max);
      mcmuon.SADCWAmp = *std::max_element(rollingSum.begin() + muon_win_min, rollingSum.begin() + muon_win_max);
      michel_win = (int)(1000. * ((mcmuon.MichelTime + readoutDelay) - startTime) / fCAENFreq);
      michel_win_min = std::max(michel_win - 20, 0);
      michel_win_max = std::min(michel_win + 20, (int)cumulativeWaveform.size() - 1);
      mcmuon.MichelRawAmp = *std::max_element(cumulativeWaveform.begin() + michel_win_min, cumulativeWaveform.begin() + michel_win_max);
      mcmuon.MichelSADCWAmp = *std::max_element(rollingSum.begin() + michel_win_min, rollingSum.begin() + michel_win_max);
    }
  } // Find MC muon, michel sadcw/raw amps

  std::vector<std::pair<size_t, float>> peakMap;
  findPeaks(rollingSum, smooth_wvf, peakMap);
  mf::LogInfo("MichelTag") << "Found " << peakMap.size() << " peaks\n";
  std::vector<std::pair<size_t, size_t>> peakPairs;
  findPeakPairs(peakMap, false, peakPairs);
  mf::LogInfo("MichelTag") << "Found " << peakPairs.size() << " peak pairs\n";
  std::vector<std::pair<float, float>> peakPairTimes;
  for (auto peakpair : peakPairs)
  {
    auto mult_win = (int)(peakpair.first / 4.);
    auto mult_win_min = std::max(0, mult_win - 30);
    auto mult_win_max = std::min((int)pmtTrigger.numPassed.size() - 1, mult_win + 30);
    auto muon_mult = *max_element(pmtTrigger.numPassed.begin() + mult_win_min, pmtTrigger.numPassed.begin() + mult_win_max);
    float muontime = (float)peakpair.first / 500. + startTime;
    float micheltime = (float)peakpair.second / 500. + startTime;
    peakPairTimes.push_back(std::make_pair(muontime, micheltime));
    MichelTag micheltag_trigger;

    micheltag_trigger.MuonTime = muontime;
    micheltag_trigger.MichelTime = micheltime;
    micheltag_trigger.MuonRawAmp = cumulativeWaveform[peakpair.first];
    micheltag_trigger.MichelRawAmp = cumulativeWaveform[peakpair.second];
    micheltag_trigger.MuonSADCWAmp = rollingSum[peakpair.first];
    micheltag_trigger.MichelSADCWAmp = rollingSum[peakpair.second];
    micheltag_trigger.MuonMult = muon_mult;
    // micheltag_trigger.MuonMult = muon_mult;
    for (auto &crthit : crtHits)
    {
      if (abs(crthit.Time - muontime) * 1000. < fCRTCoincidence)
      {
        crthit.HasFlash = true;
        micheltag_trigger.CRTPlane = crthit.Plane;
        if (fUseMC)
        {
          micheltag_trigger.G4ID = crthit.G4ID;
          micheltag_trigger.G4PDG = crthit.PDG;
        } // UseMC
      } // If crthit is time coincident with muon flash
    } // CRTHit loop
    micheltag_v->push_back(micheltag_trigger);

    // Add associations to OpFlashes within 50ns
    for (const auto &opfTag : fOpFlashLabels)
    {
      auto const &flash_h = e.getValidHandle<std::vector<recob::OpFlash>>(opfTag);
      if (!flash_h.isValid() || flash_h->empty())
        mf::LogInfo("MichelTagger") << "Don't have good flashes from producer " << opfTag << "\n";
      std::vector<art::Ptr<recob::OpFlash>> _opflash_ptr_v;
      art::fill_ptr_vector(_opflash_ptr_v, flash_h);

      for (const auto flash : _opflash_ptr_v)
      {
        if (abs(flash->Time() - (micheltag_trigger.MuonTime - readoutDelay)) > 0.05)
          continue;
        util::CreateAssn(*this, e, *micheltag_v, flash, *micheltag_opflash_assn_v);
      }
    }

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
    if (fSaveHists && std::find(fSaveEvents.begin(), fSaveEvents.end(), fEventID) != fSaveEvents.end())
    {
      unsigned hist_wvf_start = (peakpair.first > 200) ? peakpair.first - 200 : 0;
      unsigned hist_wvf_end = (peakpair.second + 2000 > nSamples - 1) ? nSamples - 1 : peakpair.second + 2000;
      std::stringstream suffix_stream;
      suffix_stream << "raw_" << peakpair.first;
      std::vector<float> wvf_slice(hist_wvf_end - hist_wvf_start);
      for (auto i = hist_wvf_start; i < hist_wvf_end; i++)
        wvf_slice[i - hist_wvf_start] = cumulativeWaveform[i];
      saveWaveformToHistogram(wvf_slice, fEventID, suffix_stream.str());
      std::stringstream suffix_stream_sadcw;
      suffix_stream_sadcw << "sadcw_" << peakpair.first;
      for (auto i = hist_wvf_start; i < hist_wvf_end; i++)
        wvf_slice[i - hist_wvf_start] = rollingSum[i];
      saveWaveformToHistogram(wvf_slice, fEventID, suffix_stream_sadcw.str());
      std::stringstream suffix_stream_smooth;
      suffix_stream_smooth << "smooth_" << peakpair.first;
      for (auto i = hist_wvf_start; i < hist_wvf_end; i++)
        wvf_slice[i - hist_wvf_start] = smooth_wvf[i];
      saveWaveformToHistogram(wvf_slice, fEventID, suffix_stream_smooth.str());
    }
  }
  return peakPairTimes;
}

void sbnd::MichelTaggerProducer::addPeakToMap(std::map<int, std::vector<std::pair<int, float>>> &multMap, int g4id, int opChannel, int pair_channel, float peakAmp, bool requirePositive)
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

std::pair<std::vector<float>, std::vector<float>> sbnd::MichelTaggerProducer::binAndFit(const std::vector<float> &data, size_t peakIndex, int channel)
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
void sbnd::MichelTaggerProducer::ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, int eventNumber)
{
  art::ServiceHandle<art::TFileService> tfs;
  // Retrieve necessary information from the waveform
  int channel = waveform.ChannelNumber();
  auto opType = fPDMap.pdType(fChannel);
  auto conversionFactor = (fPDMap.electronicsType(channel) == "daphne") ? fDaphneFreq : fCAENFreq;
  int start_bin = (fMuonPeakTime - waveform.TimeStamp() > 0.2) ? (int)(((fMuonPeakTime - 0.2) - waveform.TimeStamp()) * 1000. / conversionFactor) : 0;
  int end_bin = (int)(((fMichelPeakTime + 1.) - waveform.TimeStamp()) * 1000. / conversionFactor);
  // Get the channel number
  const std::vector<short> &waveData = waveform.Waveform(); // Get the waveform data (amplitudes)
  // int numSamples = waveData.size();
  if ((int)waveData.size() < end_bin)
    end_bin = waveData.size() - 1;
  // Get the number of waveform samples
  int numSamples = end_bin - start_bin;

  // Create a unique histogram name using channel number, event number, and timestamp
  std::stringstream histName;
  histName << "waveform_ch" << channel << "_evt" << eventNumber << "_ts_" << (size_t)fMuonPeakTime * 1000;

  // Create the histogram with an appropriate binning (one bin per sample)
  TH1F *hist = tfs->make<TH1F>(histName.str().c_str(), histName.str().c_str(), numSamples, start_bin * conversionFactor / 1000., end_bin * conversionFactor / 1000.);

  // Fill the histogram with the waveform data
  for (int i = start_bin; i < numSamples; ++i)
  {
    hist->SetBinContent(i + 1 - start_bin, waveData[i]); // Set bin content (ROOT bins start at 1)
  }

  // Optionally, set axis labels (if needed)
  hist->GetXaxis()->SetTitle("Sample Number");
  hist->GetYaxis()->SetTitle("Amplitude");
}

void sbnd::MichelTaggerProducer::saveWaveformToHistogram(const std::vector<float> &waveform, int eventNumber, std::string suffix)
{
  art::ServiceHandle<art::TFileService> tfs;
  // Retrieve necessary information from the waveform
  int start_bin = 0;
  int end_bin = waveform.size();
  // Get the channel number
  // int numSamples = waveData.size();                         // Get the number of waveform samples
  int numSamples = end_bin - start_bin;

  // Create a unique histogram name using channel number, event number, and timestamp
  std::stringstream histName;
  histName << "sumwaveform_evt_" << eventNumber << "_" << suffix;

  // Create the histogram with an appropriate binning (one bin per sample)
  TH1F *hist = tfs->make<TH1F>(histName.str().c_str(), histName.str().c_str(), numSamples, 0, numSamples);

  // Fill the histogram with the waveform data
  for (int i = start_bin; i < numSamples; ++i)
  {
    hist->SetBinContent(i + 1, waveform[i]); // Set bin content (ROOT bins start at 1)
  }

  // Optionally, set axis labels (if needed)
  hist->GetXaxis()->SetTitle("Sample Number");
  hist->GetYaxis()->SetTitle("Amplitude");
}

DEFINE_ART_MODULE(sbnd::MichelTaggerProducer)
