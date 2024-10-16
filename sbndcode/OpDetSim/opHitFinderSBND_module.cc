////////////////////////////////////////////////////////////////////////
// Class:       opHitFinderSBND
// Module Type: producer
// File:        opHitFinderSBND_module.cc
//
// This module produces an OpHit object for light analysis
// Created by L. Paulucci, F. Marinho, and I.L. de Icaza
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
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
//#include "lardataobj/Simulation/BeamGateInfo.h"

#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
//#include "larsim/MCCheater/PhotonBackTracker.h"

#include <memory>
#include <algorithm>
#include <vector>
#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

namespace opdet {

  class opHitFinderSBND;

  class opHitFinderSBND : public art::EDProducer {
  public:
    explicit opHitFinderSBND(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    opHitFinderSBND(opHitFinderSBND const &) = delete;
    opHitFinderSBND(opHitFinderSBND &&) = delete;
    opHitFinderSBND & operator = (opHitFinderSBND const &) = delete;
    opHitFinderSBND & operator = (opHitFinderSBND &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;
    opdet::sbndPDMapAlg map; //map for photon detector types

  private:

    // Declare member data here.
    std::string fInputModuleName;
    //  art::ServiceHandle<cheat::PhotonBackTracker> pbt;
    double fSampling; //in MHz
    double fSampling_Daphne; //Daphne electronics's readout sampling frequency, in MHz
    double fBaselineSample; //in ticks
    double fPulsePolarityPMT;
    double fPulsePolarityArapuca;
    double fSaturation; //in number of p.e.
    double fArea1pePMT; //area of 1 pe in ADC*ns for PMTs
    double fArea1peSiPM; //area of 1 pe in ADC*ns for Arapucas
    bool fUseDenoising;
    int fThresholdPMT; //in ADC
    int fThresholdArapuca; //in ADC
    int fEvNumber;
    int fChNumber;
    std::string opdetType;
    std::string electronicsType;
    int threshold;
    std::vector<double> fwaveform;
    std::vector<double> outwvform;
    //int fSize;
    //int fTimePMT;         //Start time of PMT signal
    //int fTimeMax;         //Time of maximum (minimum) PMT signal
    void subtractBaseline(std::vector<double>& waveform, std::string pdtype, std::string electronicsType, double& rms);
    bool findAndSuppressPeak(std::vector<double>& waveform, size_t& timebin,
                             double& Area, double& amplitude,
                             const int& threshold, const std::string& opdetType, const std::string& electronicsType);
    void denoise(std::vector<double>& waveform, std::vector<double>& outwaveform);
    bool TV1D_denoise(std::vector<double>& waveform,
                      std::vector<double>& outwaveform,
                      const double lambda);
    void TV1D_denoise_v2(std::vector<double>& input, std::vector<double>& output,
                         unsigned int width, const double lambda);
    //std::stringstream histname;
  };

  opHitFinderSBND::opHitFinderSBND(fhicl::ParameterSet const & p)
    : EDProducer{p}
      // Initialize member data here.
  {
    fInputModuleName = p.get< std::string >("InputModule" );
    fBaselineSample  = p.get< int    >("BaselineSample"); //in ticks
    fSaturation      = p.get< double >("Saturation"   ); //in number of p.e.
    fArea1pePMT      = p.get< double >("Area1pePMT"   ); //in ADC*ns for PMTs
    fArea1peSiPM     = p.get< double >("Area1peSiPM"  ); //in ADC*ns for SiPMs
    fThresholdPMT    = p.get< double >("ThresholdPMT" ); //in ADC
    fThresholdArapuca = p.get< double>("ThresholdArapuca"); //in ADC
    fPulsePolarityPMT = p.get< int   >("PulsePolarityPMT");
    fPulsePolarityArapuca = p.get<int>("PulsePolarityArapuca");
    fUseDenoising     = p.get< bool  >("UseDenoising");

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fSampling = clockData.OpticalClock().Frequency(); // MHz
    fSampling_Daphne = p.get<double>("DaphneFrequency"); 

    // Call appropriate produces<>() functions here.
    produces<std::vector<recob::OpHit>>();
  }

  void opHitFinderSBND::produce(art::Event & e)
  {
    // Implementation of required member function here.
    fEvNumber = e.id().event();
    mf::LogInfo("opHitFinder") << "Event #" << fEvNumber;

    std::unique_ptr< std::vector< recob::OpHit > > pulseVecPtr(std::make_unique< std::vector< recob::OpHit > > ());
    fwaveform.reserve(30000); // TODO: no hardcoded value
    outwvform.reserve(30000); // TODO: no hardcoded value

    art::ServiceHandle<art::TFileService> tfs;
    art::Handle< std::vector< raw::OpDetWaveform > > wvfHandle;
    std::vector<art::Ptr<raw::OpDetWaveform>> wvfList;
    if(e.getByLabel(fInputModuleName, wvfHandle))
      art::fill_ptr_vector(wvfList, wvfHandle);

    if(!wvfHandle.isValid()) {
      mf::LogWarning("opHitFinder") << Form("Did not find any waveform");
    }

    size_t timebin = 0;
    double FWHM = 1, Area = 0, phelec, fasttotal = 3./4., rms = 0, amplitude = 0, time = 0;
    unsigned short frame = 1;
    //int histogram_number = 0;
    for(auto const& wvf_P : wvfList) {
      auto const& wvf = *wvf_P;
      if (wvf.size() == 0 ) {
        mf::LogInfo("opHitFinder") << "Empty waveform, continue.";
        continue;
      }

      fChNumber = wvf.ChannelNumber();
      opdetType = map.pdType(fChNumber);
      electronicsType = map.electronicsType(fChNumber);
      if(opdetType == "pmt_coated" || opdetType == "pmt_uncoated") {
        threshold = fThresholdPMT;
      }
      else if((opdetType == "xarapuca_vuv") || (opdetType == "xarapuca_vis")) {
        threshold = fThresholdArapuca;
      }
      else {
        mf::LogWarning("opHitFinder") << "Unexpected OpChannel: " << opdetType;
        continue;
      }

      fwaveform.resize(wvf.size());
      for(unsigned int i = 0; i < wvf.size(); i++) {
        fwaveform[i] = wvf[i];
      }

      subtractBaseline(fwaveform, opdetType, electronicsType, rms);

      if(fUseDenoising) {
        if((opdetType == "pmt_coated") || (opdetType == "pmt_uncoated")) {
        }
        else if((opdetType == "xarapuca_vuv") || (opdetType == "xarapuca_vis")) {
          denoise(fwaveform, outwvform);
        }
        else {
          mf::LogInfo("opHitFinder") << "Unexpected OpChannel: " << opdetType
                    << ", continue." << std::endl;
          std::terminate();
        }
      }

      // TODO: pass rms to this function once that's sorted. ~icaza
      while(findAndSuppressPeak(fwaveform, timebin, Area, amplitude, threshold, opdetType, electronicsType)){
        if(electronicsType == "daphne") time = wvf.TimeStamp() + (double)timebin / fSampling_Daphne;
        else time = wvf.TimeStamp() + (double)timebin / fSampling;

        if(opdetType == "pmt_coated" || opdetType == "pmt_uncoated") {
          phelec = Area / fArea1pePMT;
        }
        else if((opdetType == "xarapuca_vuv") || (opdetType == "xarapuca_vis")) {
          phelec = Area / fArea1peSiPM;
        }
        else {
          mf::LogWarning("opHitFinder")  << "Unexpected OpChannel: " << opdetType
                                         << ", continue.";
          continue;
        }

        //including hit info: OpChannel, PeakTime, PeakTimeAbs, Frame, Width, Area, PeakHeight, PE, FastToTotal
        recob::OpHit opHit(fChNumber, time, time, frame, FWHM, Area, amplitude, phelec, fasttotal);
        pulseVecPtr->emplace_back(opHit);
      } // while findAndSuppressPeak()
    } // for(auto const& wvf : (*wvfHandle)){
    e.put(std::move(pulseVecPtr));
    std::vector<double>().swap(fwaveform); // clear and release the memory of fwaveform
    std::vector<double>().swap(outwvform); // clear and release the memory of outwvform
  } // void opHitFinderSBND::produce(art::Event & e)

  DEFINE_ART_MODULE(opHitFinderSBND)

  void opHitFinderSBND::subtractBaseline(std::vector<double>& waveform,
                                         std::string pdtype, std::string electronicsType, double& rms)
  {
    double baseline = 0.0;
    rms = 0.0;
    int cnt = 0;
    double NBins=fBaselineSample;
    if (electronicsType=="daphne") NBins/=(fSampling/fSampling_Daphne);//correct the number of bins to the sampling frecuency. TODO: use a fixed time interval instead, then use the channel frequency to get the number of bins ~rodrigoa
    // TODO: this is broken it assumes that the beginning of the
    // waveform is only noise, which is not always the case. ~icaza.
    // TODO: use std::accumulate instead of this loop. ~icaza.
    for(int i = 0; i < NBins; i++) {
      baseline += waveform[i];
      rms += std::pow(waveform[i], 2);
      cnt++;
    }

    baseline = baseline / cnt;
    rms = sqrt(rms / cnt - baseline * baseline);
    rms = rms / sqrt(cnt - 1);

    if(pdtype == "pmt_coated" || pdtype == "pmt_uncoated") {
      for(unsigned int i = 0; i < waveform.size(); i++) waveform[i] = fPulsePolarityPMT * (waveform[i] - baseline);
    }
    else if((opdetType == "xarapuca_vuv") || (opdetType == "xarapuca_vis")) {
      for(unsigned int i = 0; i < waveform.size(); i++) waveform[i] = fPulsePolarityArapuca * (waveform[i] - baseline);
    }
    else {
      mf::LogWarning("opHitFinder") << "Unexpected OpChannel: " << opdetType;
      return;
    }
  }


  // TODO: pass rms to this function once that's sorted. ~icaza
  bool opHitFinderSBND::findAndSuppressPeak(std::vector<double>& waveform,
                                            size_t& timebin, double& Area,
                                            double& amplitude, const int& threshold,
                                            const std::string& opdetType,
                                            const std::string& electronicsType)
  {

    std::vector<double>::iterator max_element_it = std::max_element(waveform.begin(), waveform.end());
    amplitude = *max_element_it;
    if(amplitude < threshold) return false; // stop if there's no more peaks
    timebin = std::distance(waveform.begin(), max_element_it);

    // it_e contains the iterator to the last element in the peak
    // where waveform is above threshold
    auto it_e = std::find_if(max_element_it,
                            waveform.end(),
                            [threshold](const double& x)->bool
                              {return x < threshold;} );
    // it_s contains the iterator to the first element in the peak
    // where waveform is above threshold
    auto it_s = std::find_if(std::make_reverse_iterator(max_element_it),
                            std::make_reverse_iterator(waveform.begin()),
                            [threshold](const double& x)->bool
                              {return x < threshold;} ).base();

    // integrate the area below the peak
    // note that fSampling is in MHz and
    // we convert it to GHz here so as to
    // have an area in ADC*ns.
    Area = std::accumulate(it_s, it_e, 0.0);
    if (electronicsType == "daphne"){
    Area = Area / (fSampling_Daphne / 1000.);}
    else{
    Area = Area / (fSampling / 1000.);};

    // TODO: try to just remove this
    // TODO: better even, return iterator to last position
    std::fill(it_s, it_e, 0.0); // zeroes out that peak
    return true;
  } // bool opHitFinderSBND::findAndSuppressPeak()


  void opHitFinderSBND::denoise(std::vector<double>& waveform, std::vector<double>& outwaveform)
  {

    int wavelength = waveform.size();
    outwaveform = waveform;  // copy
    double lambda = 10.0;
    const uint retries = 5; uint try_ = 0;
    if (wavelength > 0) {
      while (try_ <= retries) {
        if (TV1D_denoise(waveform, outwaveform, lambda)) break;
        try_++;
        mf::LogInfo("opHitFinder") << try_ << "/" << retries
                                   << " Coming out of TV1D_denoise() unsuccessfully, "
                                   << "using lambda: " << lambda;
        lambda += 0.1 * lambda;
        if (try_ == retries) mf::LogWarning("opHitFinder") <<  "Couldn't denoise!";
      }
    }
    // if (wavelength > 0) TV1D_denoise_v2(waveform, outwaveform, wavelength, lambda);

    // TODO: fairly certain this for is completely redundant,
    // and if not a std::move or swap would be better. ~icaza
    for(int i = 0; i < wavelength; i++) {
      if(outwaveform[i]) waveform[i] = outwaveform[i];
    }
  } // void opHitFinderSBND::denoise()

  // TODO: this function is not robust, check if the expected input is given and put exceptions
  bool opHitFinderSBND::TV1D_denoise(std::vector<double>& waveform,
                                     std::vector<double>& outwaveform,
                                     const double lambda)
  {
    int width = waveform.size();
    int k = 0, k0 = 0; // k: current sample location, k0: beginning of current segment
    double umin = lambda, umax = -lambda; // u is the dual variable
    double vmin = waveform[0] - lambda, vmax = waveform[0] + lambda; // bounds for the segment's value
    int kplus = 0, kminus = 0; // last positions where umax=-lambda, umin=lambda, respectively
    const double twolambda = 2.0 * lambda; // auxiliary variable
    const double minlambda = -lambda; // auxiliary variable
    for (;;) { // simple loop, the exit test is inside
      while (k == width - 1) { // we use the right boundary condition
        if (umin < 0.0) { // vmin is too high -> negative jump necessary
          do outwaveform[k0++] = vmin; while (k0 <= kminus);
          umax = (vmin = waveform[kminus = k = k0]) + (umin = lambda) - vmax;
        }
        else if (umax > 0.0) { // vmax is too low -> positive jump necessary
          do outwaveform[k0++] = vmax; while (k0 <= kplus);
          umin = (vmax = waveform[kplus = k = k0]) + (umax = minlambda) - vmin;
        }
        else {
          vmin += umin / (k - k0 + 1);
          do outwaveform[k0++] = vmin; while(k0 <= k);
          return true;
        }
      } // while (k == width - 1)
      if ((umin += waveform[k + 1] - vmin) < minlambda) { // negative jump necessary
        if (k0 > width) return false;
        do outwaveform[k0++] = vmin; while (k0 <= kminus);
        vmax = (vmin = waveform[kplus = kminus = k = k0]) + twolambda;
        umin = lambda; umax = minlambda;
      }
      else if ((umax += waveform[k + 1] - vmax) > lambda) { // positive jump necessary
        if (k0 > width) return false;
        do outwaveform[k0++] = vmax; while (k0 <= kplus);
        vmin = (vmax = waveform[kplus = kminus = k = k0]) - twolambda;
        umin = lambda; umax = minlambda;
      }
      else {   //no jump necessary, we continue
        k++;
        if (k > width) return false;
        if (umin >= lambda) { // update of vmin
          vmin += (umin - lambda) / ((kminus = k) - k0 + 1);
          umin = lambda;
        }
        if (umax <= minlambda) { // update of vmax
          vmax += (umax + lambda) / ((kplus = k) - k0 + 1);
          umax = minlambda;
        }
      }
    } // for (;;)
  } // bool opHitFinderSBND::TV1D_denoise()


  void opHitFinderSBND::TV1D_denoise_v2(std::vector<double>& input, std::vector<double>& output,
                                        unsigned int width, const double lambda)
  {
    // unsigned int* indstart_low = malloc(sizeof *indstart_low * width);
    // unsigned int* indstart_up = malloc(sizeof *indstart_up * width);
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);
    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    double output_low_first = input[0] - lambda;
    double output_low_curr = output_low_first;
    double output_up_first = input[0] + lambda;
    double output_up_curr = output_up_first;
    const double twolambda = 2.0 * lambda;
    if (width == 1) {
      output[0] = input[0];
      return;
    }
    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; i++) {
      if (input[i] >= output_low_curr) {
        if (input[i] <= output_up_curr) {
          output_up_curr += (input[i] - output_up_curr) / (i - indstart_up[j_up] + 1);
          output[indjseg] = output_up_first;
          while ((j_up > jseg) && (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
            output_up_curr += (output[ind] - output_up_curr) *
                              ((double)(indstart_up[j_up--] - ind) / (i - ind + 1));
          if (j_up == jseg) {
            while ((output_up_curr <= output_low_first) && (jseg < j_low)) {
              indjseg2 = indstart_low[++jseg];
              output_up_curr += (output_up_curr - output_low_first) *
                                ((double)(indjseg2 - indjseg) / (i - indjseg2 + 1));
              while (indjseg < indjseg2) output[indjseg++] = output_low_first;
              output_low_first = output[indjseg];
            }
            output_up_first = output_up_curr;
            indstart_up[j_up = jseg] = indjseg;
          }
          else output[indstart_up[j_up]] = output_up_curr;
        }
        else
          output_up_curr = output[i] = input[indstart_up[++j_up] = i];
        output_low_curr += (input[i] - output_low_curr) / (i - indstart_low[j_low] + 1);
        output[indjseg] = output_low_first;
        while ((j_low > jseg) && (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
          output_low_curr += (output[ind] - output_low_curr) *
                             ((double)(indstart_low[j_low--] - ind) / (i - ind + 1));
        if (j_low == jseg) {
          while ((output_low_curr >= output_up_first) && (jseg < j_up)) {
            indjseg2 = indstart_up[++jseg];
            output_low_curr += (output_low_curr - output_up_first) *
                               ((double)(indjseg2 - indjseg) / (i - indjseg2 + 1));
            while (indjseg < indjseg2) output[indjseg++] = output_up_first;
            output_up_first = output[indjseg];
          }
          if ((indstart_low[j_low = jseg] = indjseg) == i) output_low_first = output_up_first - twolambda;
          else output_low_first = output_low_curr;
        }
        else output[indstart_low[j_low]] = output_low_curr;
      }
      else {
        output_up_curr += ((output_low_curr = output[i] = input[indstart_low[++j_low] = i]) -
                           output_up_curr) / (i - indstart_up[j_up] + 1);
        output[indjseg] = output_up_first;
        while ((j_up > jseg) && (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
          output_up_curr += (output[ind] - output_up_curr) *
                            ((double)(indstart_up[j_up--] - ind) / (i - ind + 1));
        if (j_up == jseg) {
          while ((output_up_curr <= output_low_first) && (jseg < j_low)) {
            indjseg2 = indstart_low[++jseg];
            output_up_curr += (output_up_curr - output_low_first) *
                              ((double)(indjseg2 - indjseg) / (i - indjseg2 + 1));
            while (indjseg < indjseg2) output[indjseg++] = output_low_first;
            output_low_first = output[indjseg];
          }
          if ((indstart_up[j_up = jseg] = indjseg) == i) output_up_first = output_low_first + twolambda;
          else output_up_first = output_up_curr;
        }
        else output[indstart_up[j_up]] = output_up_curr;
      }
    }
    /* here i==width (with value the actual width minus one) */
    if (input[i] + lambda <= output_low_curr) {
      while (jseg < j_low) {
        indjseg2 = indstart_low[++jseg];
        while (indjseg < indjseg2) output[indjseg++] = output_low_first;
        output_low_first = output[indjseg];
      }
      while (indjseg < i) output[indjseg++] = output_low_first;
      output[indjseg] = input[i] + lambda;
    }
    else if (input[i] - lambda >= output_up_curr) {
      while (jseg < j_up) {
        indjseg2 = indstart_up[++jseg];
        while (indjseg < indjseg2) output[indjseg++] = output_up_first;
        output_up_first = output[indjseg];
      }
      while (indjseg < i) output[indjseg++] = output_up_first;
      output[indjseg] = input[i] - lambda;
    }
    else {
      output_low_curr += (input[i] + lambda - output_low_curr) / (i - indstart_low[j_low] + 1);
      output[indjseg] = output_low_first;
      while ((j_low > jseg) && (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
        output_low_curr += (output[ind] - output_low_curr) *
                           ((double)(indstart_low[j_low--] - ind) / (i - ind + 1));
      if (j_low == jseg) {
        if (output_up_first >= output_low_curr)
          while (indjseg <= i) output[indjseg++] = output_low_curr;
        else {
          output_up_curr += (input[i] - lambda - output_up_curr) / (i - indstart_up[j_up] + 1);
          output[indjseg] = output_up_first;
          while ((j_up > jseg) && (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
            output_up_curr += (output[ind] - output_up_curr) *
                              ((double)(indstart_up[j_up--] - ind) / (i - ind + 1));
          while (jseg < j_up) {
            indjseg2 = indstart_up[++jseg];
            while (indjseg < indjseg2) output[indjseg++] = output_up_first;
            output_up_first = output[indjseg];
          }
          indjseg = indstart_up[j_up];
          while (indjseg <= i) output[indjseg++] = output_up_curr;
        }
      }
      else {
        while (jseg < j_low) {
          indjseg2 = indstart_low[++jseg];
          while (indjseg < indjseg2) output[indjseg++] = output_low_first;
          output_low_first = output[indjseg];
        }
        indjseg = indstart_low[j_low];
        while (indjseg <= i) output[indjseg++] = output_low_curr;
      }
    }
    // free(indstart_low);
    // free(indstart_up);
  }// void opHitFinderSBND::TV1D_denoise_v2()

} // namespace opdet
