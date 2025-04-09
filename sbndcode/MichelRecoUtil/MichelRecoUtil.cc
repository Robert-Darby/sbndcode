#include "sbndcode/MichelRecoUtil/MichelRecoUtil.hh"

int MichelRecoUtil::findChannelPair(const std::vector<int> &pair1, const std::vector<int> &pair2, const std::vector<int> &unpaired, const int opChannel)
{
    int pair_channel = -1;
    int pair_idx;
    auto pair_it = std::find(pair1.begin(), pair1.end(), opChannel); // Fixing typo fOpChannel to opChannel
    if (std::find(unpaired.begin(), unpaired.end(), opChannel) != unpaired.end())
        return opChannel;
    if (pair_it != pair1.end())
    {
        pair_idx = std::distance(pair1.begin(), pair_it);
        pair_channel = pair2[pair_idx];
    }
    else
    {
        pair_it = std::find(pair2.begin(), pair2.end(), opChannel);
        pair_idx = std::distance(pair2.begin(), pair_it);
        pair_channel = pair1[pair_idx];
    }
    return pair_channel; // Moved outside the else block to fix error
}

std::vector<float> MichelRecoUtil::CalcRunningAvg(const int fRunningAvgSampleWidth, std::vector<float> &wvf)
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
        smooth_wvf[i] = wvf[i] + Baseline[i];
    // Baseline is all updated an can return
    return Baseline;
}

std::vector<float> MichelRecoUtil::subtractBaseline(const std::vector<float> &waveform, const float conversionFactor)
{
    size_t baselineSampleCount = static_cast<size_t>(100 / conversionFactor);
    float baseline = std::accumulate(waveform.begin(), waveform.begin() + baselineSampleCount, 0.0f) / baselineSampleCount;

    // Subtract baseline from waveform
    std::vector<float> correctedWaveform(waveform.size());
    std::transform(waveform.begin(), waveform.end(), correctedWaveform.begin(),
                   [baseline](float val)
                   { return val - baseline; });

    return correctedWaveform;
}

std::vector<float> MichelRecoUtil::applyRollingSum(const std::vector<float> &waveform)
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

void MichelRecoUtil::findPeaks(
    const std::vector<float> &searchWaveform,
    const std::vector<float> &waveform,
    const float michelThreshold,
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

void MichelRecoUtil::ConvertWaveformToHistogram(const raw::OpDetWaveform &waveform, const float fMuonPeakTime, const float fMichelPeakTime, const float conversionFactor, int eventNumber)
{
    art::ServiceHandle<art::TFileService> tfs;
    // Retrieve necessary information from the waveform
    int channel = waveform.ChannelNumber();
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
