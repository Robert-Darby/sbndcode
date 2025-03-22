#ifndef MICHELRECOUTILS_H
#define MICHELRECOUTILS_H

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include <algorithm>
#include <numeric>
#include <map>
#include <vector>
#include <cmath>
#include <limits>

namespace MichelRecoUtil
{
    int findChannelPair(const std::vector<int>& pair1, const std::vector<int>& pair2, const std::vector<int>& unpaired, int opChannel);
    std::vector<float> CalcRunningAvg(const int fRunngingAvgSampleWidth, std::vector<float> &wvf);
    std::vector<float> subtractBaseline(const std::vector<float> &waveform, const float conversionFactor);
    std::vector<float> applyRollingSum(const std::vector<float> &waveform);
    void findPeaks(const std::vector<float> &seazrchWaveform, const std::vector<float> &waveform, const float michelThreshold, std::vector<std::pair<size_t, float>> &peakIndices);
    void ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, const float fMuonPeakTime, const float fMichelPeakTime, const float conversionFactor, int eventNumber);

}

#endif
