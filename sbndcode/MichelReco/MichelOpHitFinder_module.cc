////////////////////////////////////////////////////////////////////////
// Class:       MichelOpHitFinder
// Plugin Type: producer (Unknown Unknown)
// File:        MichelOpHitFinder_module.cc
//
// Generated at Wed Oct 18 11:51:42 2023 by Robert Darby using cetskelgen
// from  version .
//
// This module finds custom OpHits from Michel electrons taking the
// deconvolved raw::OpDetWaveforms as input and returns a vector of
// recob::OpHits associated with Michels in the event.
////////////////////////////////////////////////////////////////////////

#include "sbndcode/MichelReco/MichelOpHitFinder.h"
#include "art/Framework/Core/ModuleMacros.h"

sbnd::MichelOpHitFinder::MichelOpHitFinder(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fOpDetWaveformLabel(p.get<std::string>("OpDetWaveformLabel"))
  , fOpDetTypes(p.get<std::vector<std::string>>("OpDetTypes"))
  , fMinADC(p.get<std::vector<float>>("MinADC"))
  , fMuonDelay(p.get<float>("MuonDelay", 200.))
  , fPositivePolarity(p.get<bool>("PositivePolarity", false))
  , fMakeTree(p.get<bool>("MakeTree", false))
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSamplingFreq=clockData.OpticalClock().Frequency();//in MHz


  if(fMakeTree) initTree();
}  // sbnd::MichelOpHitFinder(fhicl::MichelOpHitFinder const& p)


void sbnd::MichelOpHitFinder::produce(art::Event& e)
{
  if(fMakeTree) {
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();
    resetTree();
  }

  //Load raw waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fOpDetWaveformLabel, wfHandle);
  if(!wfHandle.isValid()) {
    mf::LogError("MIchelOpHitFinder")
      << "Couldn't find input waveforms with label: " << fOpDetWaveformLabel << std::endl;
    throw cet::exception("MichelOpHitFinder")
      << "Couldn't find input waveforms with label: " << fOpDetWaveformLabel << std::endl;
  }

  // Store waveforms in vector
  std::vector< raw::OpDetWaveform > decoWVFVec;
  decoWVFVec.reserve(wfHandle->size());

  for(auto const& wf : *wfHandle){
    decoWVFVec.push_back(wf);
  }

  // Loop over waveforms
  for(const auto& wv : decoWVFVec) {
    size_t wfsize=wv.Waveform().size();
    std::vector<double> wave;
    wave.reserve(wfsize);
    wave.assign(wv.Waveform().begin(), wv.Waveform().end());

    // Find size, time and position of muon peak
    std::vector<double>::iterator peak_adc_pos = (fPositivePolarity) ? std::max_element(wave.begin(), wave.end()) :
                                                                       std::min_element(wave.begin(), wave.end());
    if((fPositivePolarity && *peak_adc_pos < fMinADC[0]) ||
       (!fPositivePolarity && *peak_adc_pos > fMinADC[0])) continue;
    size_t peak_adc_idx = std::distance(wave.begin(), peak_adc_pos);
    if(peak_adc_idx + int(fMuonDelay*fSamplingFreq) > wave.size()) continue;
    double muon_time = wv.TimeStamp() + (fSamplingFreq*peak_adc_idx);
  //  if(fMakeTree)
    _muon_times.push_back(muon_time);

  // ToDo: find muon peak time, skip waveforms that exceed peak time + delay
  } // Loop over waveform vector
} // void sbnd::MichelOpHitFinder::produce(art::Event& e)

void sbnd::MichelOpHitFinder::beginJob() {}
void sbnd::MichelOpHitFinder::endJob() {}

void sbnd::MichelOpHitFinder::initTree(void)
{
  art::ServiceHandle<art::TFileService> tfs;
  _mophf_tree = tfs->make<TTree>("mophftree", "Michel OpHit tree");
  _mophf_tree->Branch("run", &_run, "run/I");
  _mophf_tree->Branch("sub", &_sub, "sub/I");
  _mophf_tree->Branch("evt", &_evt, "evt/I");
  _mophf_tree->Branch("mu_time", &_muon_times, "mu_time/D");
} // sbnd::MichelOpHitFinder::initTree

void sbnd::MichelOpHitFinder::resetTree(void)
{
  _muon_times.clear();
}

DEFINE_ART_MODULE(sbnd::MichelOpHitFinder)
