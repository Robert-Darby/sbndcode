////////////////////////////////////////////////////////////////////////
// Class:       MichelVisMap
// Plugin Type: analyzer (Unknown Unknown)
// File:        MichelVisMap_module.cc
//
// Generated at Wed Feb 26 07:34:45 2025 by Robert Darby using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// Obj includes
#include "canvas/Persistency/Common/Assns.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "sbnobj/SBND/Trigger/MichelTag.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// ROOT Includes
#include "TTree.h"

namespace sbnd
{
  class MichelVisMap;
}

class sbnd::MichelVisMap : public art::EDAnalyzer
{
public:
  explicit MichelVisMap(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelVisMap(MichelVisMap const &) = delete;
  MichelVisMap(MichelVisMap &&) = delete;
  MichelVisMap &operator=(MichelVisMap const &) = delete;
  MichelVisMap &operator=(MichelVisMap &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  struct MCMuon
  {
    int G4ID;
    float StartX, StartY, StartZ;
    float Time, RawAmp, SADCWAmp;
    int Mult, MichelMult;
    int MichelG4ID;
    float MichelStartX, MichelStartY, MichelStartZ;
    float MichelEnergy, MichelLength, MichelTime, MichelRawAmp, MichelSADCWAmp;
    long MichelDepPE;
    // float MichelTime, MichelEnergy;
    float FlashTime, MichelFlashTime;
    bool HasFlash, Stopping, EntersTPC;
  };

  struct PDInfo
  {
    int Channel;
    int Type;
    float X, Y, Z;

    float NMCPE;
    float NRecoPE;
  };

private:
  // Declare member data here.
  // Globals
  std::vector<MCMuon>
      muon_tuple_vect;

  // Tree variables
  TTree *fTree;
  int evt_, sub_, run_;
  int mcmuong4id_;
  float mcmuonstartx_, mcmuonstarty_, mcmuonstartz_;
  float mcmuonendx_, mcmuonendy_, mcmuonendz_;
  float mcmuonstarttime_, mcmuonendtime_;

  int mcmichelg4id_;
  float mcmichelstartx_, mcmichelstarty_, mcmichelstartz_;
  float mcmichelendx_, mcmichelendy_, mcmichelendz_;
  float mcmichelstarttime_, mcmichelendtime_;
  float mcmichelenergy_, mcmichellength_;
  float mcmicheltotalpe_, mcmichelvispe_;

  float micheltagmuontime_, micheltagmicheltime_;

  std::vector<int> pdchannel_;
  std::vector<int> pdtype_;
  std::vector<float> pdx_, pdy_, pdz_;
  std::vector<float> pdmcphotons_;
  std::vector<float> pdrecophotons_;

  const std::string fMCTruthLabel;
  const std::string fSEDLabel;
  const std::string fSEDOutLabel;

  const bool fUseMC;

  // Functions
  void resetVars();
  void findMCMuons(const art::Event &e);
  void fillMCVars(const MCMuon &mcmuon);
};

sbnd::MichelVisMap::MichelVisMap(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
      fMCTruthLabel(p.get<std::string>("MCTruthLabel")),
      fSEDLabel(p.get<std::string>("SEDLabel")),
      fSEDOutLabel(p.get<std::string>("SEDOutLabel")),

      fUseMC(p.get<bool>("UseMC", false))

// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::MichelVisMap::analyze(art::Event const &e)
{
  // Implementation of required member function here.
  evt_ = e.id().event();
  run_ = e.run();
  sub_ = e.subRun();

  if (fUseMC)
    findMCMuons(e);

  for (const auto &mcmuon : muon_tuple_vect)
  {
    resetVars();
    fillMCVars(mcmuon);
    fTree->Fill();
  }
}

void sbnd::MichelVisMap::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("vismap_tree", "Michel Visibility Map Metrics");

  fTree->Branch("mcmuong4id", &mcmuong4id_, "mcmuong4id/I");
  fTree->Branch("mcmuonstartx", &mcmuonstartx_, "mcmuonstartx/F");
  fTree->Branch("mcmuonstarty", &mcmuonstarty_, "mcmuonstarty/F");
  fTree->Branch("mcmuonstartz", &mcmuonstartz_, "mcmuonstartz/F");
  fTree->Branch("mcmuonendx", &mcmuonendx_, "mcmuonendx/F");
  fTree->Branch("mcmuonendy", &mcmuonendy_, "mcmuonendy/F");
  fTree->Branch("mcmuonendz", &mcmuonendz_, "mcmuonendz/F");
  fTree->Branch("mcmuonstarttime", &mcmuonstarttime_, "mcmuonstarttime/F");
  fTree->Branch("mcmuonendtime", &mcmuonendtime_, "mcmuonendtime/F");

  fTree->Branch("mcmichelg4id", &mcmichelg4id_, "mcmichelg4id/I");
  fTree->Branch("mcmichelstartx", &mcmichelstartx_, "mcmichelstartx/F");
  fTree->Branch("mcmichelstarty", &mcmichelstarty_, "mcmichelstarty/F");
  fTree->Branch("mcmichelstartz", &mcmichelstartz_, "mcmichelstartz/F");
  fTree->Branch("mcmichelendx", &mcmichelendx_, "mcmichelendx/F");
  fTree->Branch("mcmichelendy", &mcmichelendy_, "mcmichelendy/F");
  fTree->Branch("mcmichelendz", &mcmichelendz_, "mcmichelendz/F");
  fTree->Branch("mcmichelstarttime", &mcmichelstarttime_, "mcmichelstarttime/F");
  fTree->Branch("mcmichelendtime", &mcmichelendtime_, "mcmichelendtime/F");
  fTree->Branch("mcmichelenergy", &mcmichelenergy_, "mcmichelenergy/F");
  fTree->Branch("mcmichellength", &mcmichellength_, "mcmichellength/F");
  fTree->Branch("mcmicheltotalpe", &mcmicheltotalpe_, "mcmicheltotalpe/F");
  fTree->Branch("mcmichelvispe", &mcmichelvispe_, "mcmichelvispe/F");

  fTree->Branch("micheltagmuontime", &micheltagmuontime_, "micheltagmuontime/F");
  fTree->Branch("micheltagmicheltime", &micheltagmicheltime_, "micheltagmicheltime/F");

  fTree->Branch("pdchannel", &pdchannel_);
  fTree->Branch("pdtype", &pdtype_);
  fTree->Branch("pdx", &pdx_);
  fTree->Branch("pdy", &pdy_);
  fTree->Branch("pdz", &pdz_);
  fTree->Branch("pdmcphotons", &pdmcphotons_);
  fTree->Branch("pdrecophotons", &pdrecophotons_);
}

void sbnd::MichelVisMap::endJob()
{
  // Implementation of optional member function here.
}

void sbnd::MichelVisMap::resetVars()
{
  mcmuong4id_ = -9999;
  mcmuonstartx_ = mcmuonstarty_ = mcmuonstartz_ = -9999.;
  mcmuonendx_ = mcmuonendy_ = mcmuonendz_ = -9999.;
  mcmuonstarttime_ = mcmuonendtime_ = -9999.;

  mcmichelg4id_ = -9999;
  mcmichelstartx_ = mcmichelstarty_ = mcmichelstartz_ = -9999.;
  mcmichelendx_ = mcmichelendy_ = mcmichelendz_ = -9999.;
  mcmichelstarttime_ = mcmichelendtime_ = -9999.;
  mcmichelenergy_ = mcmichellength_ = -9999.;
  mcmicheltotalpe_ = mcmichelvispe_ = -9999.;

  micheltagmuontime_ = micheltagmicheltime_ = -9999.;

  pdchannel_.clear();
  pdtype_.clear();
  pdx_.clear();
  pdy_.clear();
  pdz_.clear();
  pdmcphotons_.clear();
  pdrecophotons_.clear();
}

void sbnd::MichelVisMap::findMCMuons(const art::Event &e)
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
    float michel_length = -9999.;
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
        michel_length = mcp2->Trajectory().TotalLength();
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

    MCMuon mcmuon;
    mcmuon.G4ID = mu_id;
    mcmuon.Time = mu_time;
    mcmuon.MichelG4ID = michel_id;
    mcmuon.MichelTime = michel_time;
    mcmuon.MichelDepPE = michel_dep_pe;
    mcmuon.MichelEnergy = michel_energy;
    mcmuon.MichelLength = michel_length;
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
    muon_tuple_vect.push_back(mcmuon);
  }
  std::sort(muon_tuple_vect.begin(), muon_tuple_vect.end(),
            [](const MCMuon &muon1, const MCMuon &muon2)
            { return (muon2.Time > muon1.Time); });
}

void sbnd::MichelVisMap::fillMCVars(const MCMuon &mcmuon)
{
  mcmuong4id_ = mcmuon.G4ID;
  mcmuonstartx_ = mcmuon.StartX;
  mcmuonstarty_ = mcmuon.StartY;
  mcmuonstartz_ = mcmuon.StartZ;
  mcmuonstarttime_ = mcmuon.Time;

  mcmichelg4id_ = mcmuon.MichelG4ID;
  mcmichelstartx_ = mcmuon.MichelStartX;
  mcmichelstarty_ = mcmuon.MichelStartY;
  mcmichelstartz_ = mcmuon.MichelStartZ;
  mcmichelstarttime_ = mcmuon.MichelTime;
  mcmichelenergy_ = mcmuon.MichelEnergy;
  mcmichellength_ = mcmuon.MichelLength;
  mcmicheltotalpe_ = mcmuon.MichelDepPE;
  mcmichelvispe_ = 0.;
}

DEFINE_ART_MODULE(sbnd::MichelVisMap)
