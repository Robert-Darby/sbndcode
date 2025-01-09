#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TTree.h"

namespace sbnd
{

    class MichelAnalyser;
}
class sbnd::MichelAnalyser : public art::EDAnalyzer
{
public:
    explicit MichelAnalyser(fhicl::ParameterSet const &p);

    void analyze(art::Event const &e) override;
    void beginJob() override;

    MichelAnalyser(MichelAnalyser const &) = delete;
    MichelAnalyser(MichelAnalyser &&) = delete;
    MichelAnalyser &operator=(MichelAnalyser const &) = delete;
    MichelAnalyser &operator=(MichelAnalyser &&) = delete;

private:
    TTree *fTree;

    // Variables for TTree
    int fPdgCode;
    float fStartX, fStartY, fStartZ, fEndX, fEndY, fEndZ;
    float fEnergy, fTime;
    float fDirX, fDirY, fDirZ;
    float fEnterX, fEnterY, fEnterZ, fEnterTime, fEnterEnergy;
    int fMichelPDG;
    float fMichelStartX, fMichelStartY, fMichelStartZ;
    float fMichelEndX, fMichelEndY, fMichelEndZ;
    float fMichelEnergy, fMichelTime;
    float fMichelDirX, fMichelDirY, fMichelDirZ;

    // Configuration
    const art::InputTag fMCParticleLabel;
    const float fDetectorXMin, fDetectorXMax, fDetectorYMin, fDetectorYMax, fDetectorZMin, fDetectorZMax;

    // Helper function to check if particle enters the detector
    bool EntersDetector(float x, float y, float z);
};

sbnd::MichelAnalyser::MichelAnalyser(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      fMCParticleLabel(p.get<art::InputTag>("MCParticleLabel", "largeant")),
      fDetectorXMin(-200), fDetectorXMax(200),
      fDetectorYMin(-200), fDetectorYMax(200),
      fDetectorZMin(0), fDetectorZMax(500) {}

void sbnd::MichelAnalyser::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("MCParticleTree", "MCParticle Analysis Tree");

    fTree->Branch("PdgCode", &fPdgCode, "PdgCode/I");
    fTree->Branch("StartX", &fStartX, "StartX/F");
    fTree->Branch("StartY", &fStartY, "StartY/F");
    fTree->Branch("StartZ", &fStartZ, "StartZ/F");
    fTree->Branch("EndX", &fEndX, "EndX/F");
    fTree->Branch("EndY", &fEndY, "EndY/F");
    fTree->Branch("EndZ", &fEndZ, "EndZ/F");
    fTree->Branch("Energy", &fEnergy, "Energy/F");
    fTree->Branch("Time", &fTime, "Time/F");
    fTree->Branch("DirX", &fDirX, "DirX/F");
    fTree->Branch("DirY", &fDirY, "DirY/F");
    fTree->Branch("DirZ", &fDirZ, "DirZ/F");
    fTree->Branch("EnterX", &fEnterX, "EnterX/F");
    fTree->Branch("EnterY", &fEnterY, "EnterY/F");
    fTree->Branch("EnterZ", &fEnterZ, "EnterZ/F");
    fTree->Branch("EnterTime", &fEnterTime, "EnterTime/F");
    fTree->Branch("EnterEnergy", &fEnterEnergy, "EnterEnergy/F");

    fTree->Branch("MichelPDG", &fMichelPDG, "MichelPDG/I");
    fTree->Branch("MichelStartX", &fMichelStartX, "MichelStartX/F");
    fTree->Branch("MichelStartY", &fMichelStartY, "MichelStartY/F");
    fTree->Branch("MichelStartZ", &fMichelStartZ, "MichelStartZ/F");
    fTree->Branch("MichelEndX", &fMichelEndX, "MichelEndX/F");
    fTree->Branch("MichelEndY", &fMichelEndY, "MichelEndY/F");
    fTree->Branch("MichelEndZ", &fMichelEndZ, "MichelEndZ/F");
    fTree->Branch("MichelEnergy", &fMichelEnergy, "MichelEnergy/F");
    fTree->Branch("MichelTime", &fMichelTime, "MichelTime/F");
    fTree->Branch("MichelDirX", &fMichelDirX, "MichelDirX/F");
    fTree->Branch("MichelDirY", &fMichelDirY, "MichelDirY/F");
    fTree->Branch("MichelDirZ", &fMichelDirZ, "MichelDirZ/F");
}

void sbnd::MichelAnalyser::analyze(art::Event const &e)
{
    auto const &mcParticles = *e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);

    for (auto const &particle : mcParticles)
    {
        fPdgCode = particle.PdgCode();

        // Select only muons (PDG code 13) and protons (PDG code 2212)
        if (abs(fPdgCode) == 13 || fPdgCode == 2212)
        {
            fStartX = particle.Vx();
            fStartY = particle.Vy();
            fStartZ = particle.Vz();
            fEndX = particle.EndX();
            fEndY = particle.EndY();
            fEndZ = particle.EndZ();
            fEnergy = particle.E();
            fTime = particle.T();
            fDirX = particle.Px() / particle.P();
            fDirY = particle.Py() / particle.P();
            fDirZ = particle.Pz() / particle.P();

            // Check if particle enters detector
            fEnterX = fEnterY = fEnterZ = fEnterTime = fEnterEnergy = -999;
            for (unsigned i = 0; i < particle.NumberTrajectoryPoints(); ++i)
            {
                if (EntersDetector(particle.Vx(i), particle.Vy(i), particle.Vz(i)))
                {
                    fEnterX = particle.Vx(i);
                    fEnterY = particle.Vy(i);
                    fEnterZ = particle.Vz(i);
                    fEnterTime = particle.T(i);
                    fEnterEnergy = particle.E(i);
                    break;
                }
            }

            // If particle is a muon, look for Michel electron (PDG code 11 or -11)
            if (abs(fPdgCode) == 13)
            {
                fMichelPDG = fMichelStartX = fMichelStartY = fMichelStartZ = fMichelEndX = fMichelEndY = fMichelEndZ = -999;
                fMichelEnergy = fMichelTime = fMichelDirX = fMichelDirY = fMichelDirZ = -999;

                for (auto const &daughter : mcParticles)
                {
                    if ((abs(daughter.PdgCode()) == 11 || abs(daughter.PdgCode()) == -11) && daughter.Mother() == particle.TrackId())
                    {
                        fMichelPDG = daughter.PdgCode();
                        fMichelStartX = daughter.Vx();
                        fMichelStartY = daughter.Vy();
                        fMichelStartZ = daughter.Vz();
                        fMichelEndX = daughter.EndX();
                        fMichelEndY = daughter.EndY();
                        fMichelEndZ = daughter.EndZ();
                        fMichelEnergy = daughter.E();
                        fMichelTime = daughter.T();
                        fMichelDirX = daughter.Px() / daughter.P();
                        fMichelDirY = daughter.Py() / daughter.P();
                        fMichelDirZ = daughter.Pz() / daughter.P();
                        break;
                    }
                }
            }

            // Fill the TTree with the current particle data
            fTree->Fill();
        }
    }
}

bool sbnd::MichelAnalyser::EntersDetector(float x, float y, float z)
{
    return (fDetectorXMin < x && x < fDetectorXMax &&
            fDetectorYMin < y && y < fDetectorYMax &&
            fDetectorZMin < z && z < fDetectorZMax);
}

DEFINE_ART_MODULE(sbnd::MichelAnalyser)
