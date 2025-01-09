#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>

void muRatioFinder(const char* inputFileList) {
    // Define the histograms
    TH1F* lifetimes = new TH1F("lifetimes", "Muon Lifetimes", 50, 0, 10); // 0-10 with 0.2 bin width
    TH1F* alpha_hist = new TH1F("alpha", "Alpha values", 50, 0, 10);

    // Open the text file containing the list of ROOT files
    std::ifstream fileList(inputFileList);
    std::string rootFileName;
    std::vector<std::string> rootFiles;

    while (fileList >> rootFileName) {
        rootFiles.push_back(rootFileName);
    }

    // Variables to hold tree branches
    bool Is_stopping;
    bool HasFlash;
    double MuonFlashTime;
    double MichelFlashTime;
    double timestamp;
    double michel_time;

    // Loop over the root files
    for (const auto& fileName : rootFiles) {
        TFile* file = TFile::Open(fileName.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            continue;
        }

        // Get the TTree from the directory "latelight"
        TDirectory* dir = file->GetDirectory("latelight");
        if (!dir) {
            std::cerr << "Directory 'latelight' not found in file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        TTree* tree = (TTree*)dir->Get("TriggerTree");
        if (!tree) {
            std::cerr << "TTree 'TriggerTree' not found in directory 'latelight' of file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        // Set the branch addresses
        tree->SetBranchAddress("Is_stopping", &Is_stopping);
        tree->SetBranchAddress("HasFlash", &HasFlash);
        tree->SetBranchAddress("MuonFlashTime", &MuonFlashTime);
        tree->SetBranchAddress("MichelFlashTime", &MichelFlashTime);
        tree->SetBranchAddress("timestamp", &timestamp);
        tree->SetBranchAddress("michel_time", &michel_time);

        // Loop over the events
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);

            // Apply the criteria
            if (Is_stopping && HasFlash &&
                TMath::Abs(MuonFlashTime - timestamp) < 0.05 &&
                TMath::Abs(MichelFlashTime - michel_time) < 0.05) {

                // Calculate lifetime and fill the histogram
                double lifetime = MichelFlashTime - timestamp;
                lifetimes->Fill(lifetime);
            }
        }

        file->Close();
    }

    // Now we perform the fitting for each bin
    double w = lifetimes->GetBinWidth(1);
    for (int bin = 1; bin <= lifetimes->GetNbinsX(); ++bin) {
        double t_i = lifetimes->GetBinCenter(bin);
        double t_min = t_i - 0.5 * w;
        double t_max = t_i + 0.5 * w;

        // Define the fitting function for each bin
        TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-x/[1]) + [2]*TMath::Exp(-x/[3])", t_min, t_max);
        fitFunc->SetParameters(1.0, 2.2, 1.0, 2.0);  // Initial guesses for C, tau+, alpha, tau-

        // Perform the fit on the data within this bin
        lifetimes->Fit(fitFunc, "RQ", "", t_min, t_max);

        // Extract the alpha parameter and fill the alpha histogram
        double alpha = fitFunc->GetParameter(2);
        alpha_hist->SetBinContent(bin, alpha);

        // Clean up the fit function
        delete fitFunc;
    }

    // Save the histograms to a new ROOT file
    TFile* outputFile = new TFile("muRatio.root", "RECREATE");
    lifetimes->Write();
    alpha_hist->Write();
    outputFile->Close();

    // Clean up
    delete lifetimes;
    delete alpha_hist;

    std::cout << "Histograms saved to muRatio.root" << std::endl;
}

