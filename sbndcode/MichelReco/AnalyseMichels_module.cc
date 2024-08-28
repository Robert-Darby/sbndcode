//////////////////////////////////////////////////////////////////
//
// 			TO DO 					// 
//////////////////////////////////////////////////////////////////

// GET START POS FOR MC/RECO MICHEL AND RECO MUON
// GET TRAJECTORY ANGLE FOR MC/RECO MUON/MICHEL
// CREATE BOOLEAN FOR IF VERTEX IS AT WRONG END OF MUON
// N Clusters for reco particles
// Find neutrino interactrion vertex

///////////////////////////////////////////////////////////////////////
// Class:       AnalyseMichels
// Plugin Type: Analyser (Unknown Unknown)
// File:        AnalyseMichels_module.cc
//
// Generated at Wed Oct 13 08:28:13 2021 by Edward Tyley using cetskelgen
// from  version .
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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"   // Find associations as pointers
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include <stdio.h>
#include <stdlib.h>
#include "lardataobj/RecoBase/Slice.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbndaq-artdaq-core/Obj/SBND/Coincidence.hh"

#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/SimPhotons.h"

// ROOT includes
#include <TH1F.h>
#include <TTree.h>
#include<THStack.h>
#include<TMath.h>

// STL includes
#include <string>
#include <vector>
#include <iostream>

namespace sbnd {
class AnalyseMichels;
}



class sbnd::AnalyseMichels : public art::EDAnalyzer {
  public:
  explicit AnalyseMichels(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseMichels(AnalyseMichels const&) = delete;
  AnalyseMichels(AnalyseMichels&&) = delete;
  AnalyseMichels& operator=(AnalyseMichels const&) = delete;
  AnalyseMichels& operator=(AnalyseMichels&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  private:

  void rotateVector(TVector3 *vector);  // Rotates components of 3-vector clockwise e.g. XTZ -> ZXY

  void ResetVars();

  void FillMC(const art::Ptr<simb::MCParticle>& mcp, const art::Event& e, std::vector<art::Ptr<simb::MCParticle>>& mctruthVect);
  void FindRecoMichelShower(const art::Event& e);
  void FindRecoMichelTrack(const art::Event& e);
  void FindRecoMichel(const art::Event& e);
  void FillTrigTree(const art::Event& e);
  void FillPulseTree(const art::Event& e);
  void MatchMCRecoLight(const art::Event& e);

  // Create out output tree
  TTree* fTree;

  // Create our output histograms

  TH1D *fMCMichelEnergyHist;
  TH1D *fMCMuonEnergyHist;
  TH1D *fRecoMuonEnergyHist;
  TH1D *fRecoMichelEnergyHist;
  TH1D *fMCMichelThetaHist;
  TH1D *fMCMuonThetaHist;
  TH1D *fRecoMichelThetaHist;
  TH1D *fRecoMuonThetaHist;
  TH1D *fMCMichelPhiHist;
  TH1D *fMCMuonPhiHist;
  TH1D *fRecoMichelPhiHist;
  TH1D *fRecoMuonPhiHist;
  TH1D *fMCMichelRelThetaHist;
  TH1D *fRecoMichelRelThetaHist;
  TH1D *fMCMichelRelPhiHist;
  TH1D *fRecoMichelRelPhiHist;
  TH1D *fMCMichelLengthHist;


  // Event Variables
  int fRun, fSubRun;
  int fEventID;
  int fNPFParticles;
  std::vector<int> fPFParticlePDG;
  int fEventNHits;
  int fEventNClusterHits;		// No hits in event associated with a cluster
  int fNRecoElectrons;
  std::vector<int> fMCPDG;
  std::vector<int> fMCTrackID;
  std::vector<int> fNHitsInPFP;
  std::vector<int> fNTrueHitsInPFP;
  std::vector<int> fNClusters;
  std::vector<int> fPFPIDs;
  bool fWrongMuonEnd;
  std::vector<int> fClusterID;
  std::vector<std::string> fClusterPlane;
  int fNMuons;
  float fRecoMuonMichelDist;
  std::vector<bool> fPFPIsPrimary;
  std::vector<int> fPFPMother;
  int fNMuonsWithTrueHits;
  float fMaxSADCX;
  float fMaxSADCY;
  float fMaxSADCZ;
  float fMaxSADCDist;  // Distance from hit with highest SADC integral to michel/muon vertex
  float fMaxAmpX;
  float fMaxAmpY;
  float fMaxAmpZ;
  float fMaxAmpDist;   // Distnace from hit with max peak amplitude to michel/muon vertex
  float fAmpMean;
  float fAmpSigma;
  float fSADCMean;
  float fSADCSigma;
  int fNSlices;
  int fNMichelNonClust;	// No. of non-clustered Michel hits in event

  // MC Michel
  int fMCMichelID;
  float fMCMichelEnergy;
  float fMCMichelTheta;
  float fMCMichelPhi;
  float fMCMichelRelTheta;
  float fMCMichelRelPhi;
  TVector3 *fMCMichelVect;
  float fMCMichelLength;
  float fMCMichelStartX;
  float fMCMichelStartY;
  float fMCMichelStartZ;
  float fMCMichelStartT;
  float fMCMichelRelAngle;
  int fMCMichelNPoints;
  float fMCMichelEnergyFrac;		// Fraction of end MC muon energy taken by michel
  std::string fMCMichelStartProcess;
  float fMCMichelMamophtTime;
  int fMCMichelMamophtADC;

  float fPurity;
  float fCompleteness;
  int fNHitsInRecoMichel;
  int fNTotalMichelHits;
  int fNShowers;
  int fNTracks;
  bool fRecoMichel;
  float fShowerPurity;
  float fTrackPurity;
  int fShowerBestPlane;
  bool fIsTrack;
  bool fIsShower;
  int fNTrueHitsInRecoMichel;
  float fRecoMichelEnergy;
  float fRecoMichelTheta;
  float fRecoMichelPhi;
  float fRecoMichelRelTheta;
  float fRecoMichelRelPhi;
  float fRecoMichelRelAngle;
  TVector3 *fRecoMichelVect;
  float fRecoMichelStartX;
  float fRecoMichelStartY;
  float fRecoMichelStartZ;
  int fRecoMichelID;
  int fRecoMichelMother;
  int fNRecoMichelSpacePoints;
  float fRecoMichelLength;
  float fRecoMichelAzi;
  float fRecoMichelZen;
  float fRecoMichelCloseProximity;		// Distance between closest two hits
  int fRecoMichelNClusters;
  float fRecoMichelEnergySigma;
  std::vector<double> fRecoMicheldEdX;
  std::vector<float> fRecoMichelHitIntegral;
  std::vector<float> fRecoMichelHitPeakTime;
  std::vector<float> fRecoMichelHitNElectrons;
  std::vector<int> fRecoMichelHitPlane;
  float fRecoMichelEnergyU;
  float fRecoMichelEnergyV;
  float fRecoMichelEnergyW;
  float fRecoMichelEnergyDiffU;
  float fRecoMichelEnergyDiffV;
  float fRecoMichelEnergyDiffW;
  float fRecoMichelIntegralMaxU;
  float fRecoMichelIntegralMeanU;
  float fRecoMichelIntegralSigmaU;
  float fRecoMichelIntegralMaxV;
  float fRecoMichelIntegralMeanV;
  float fRecoMichelIntegralSigmaV;
  float fRecoMichelIntegralMaxW;
  float fRecoMichelIntegralMeanW;
  float fRecoMichelIntegralSigmaW;
  float fRecoMichelMultiplicityMaxU;
  float fRecoMichelMultiplicityMeanU;
  float fRecoMichelMultiplcitySigmaU;
  float fRecoMichelGOFMaxU;
  float fRecoMichelGOFSigmaU;
  float fRecoMichelGOFMeanU;
  float fRecoMichelMultiplicityMaxV;
  float fRecoMichelMultiplicityMeanV;
  float fRecoMichelMultiplcitySigmaV;
  float fRecoMichelGOFMaxV;
  float fRecoMichelGOFSigmaV;
  float fRecoMichelGOFMeanV;
  float fRecoMichelMultiplicityMaxW;
  float fRecoMichelMultiplicityMeanW;
  float fRecoMichelMultiplcitySigmaW;
  float fRecoMichelGOFMaxW;
  float fRecoMichelGOFSigmaW;
  float fRecoMichelGOFMeanW;
  std::vector<float> fRecoMichelhitMultiplicity;
  std::vector<float> fRecoMichelHitGOF;
  float fRecoMichelEDiffMin;
  float fRecoMichelEDiffMinGOF;
  float fRecoMichelEDiffMinMult;
  float fRecoMichelPurity;
  float fRecoMichelPurityU;
  float fRecoMichelPurityV;
  float fRecoMichelPurityW;
  float fRecoMichelCompleteness;
  float fRecoMichelCompletenessW;
  float fRecoMichelEPurityU;
  float fRecoMichelEPurityV;
  float fRecoMichelEPurityW;
  float fRecoMichelEPurity;
  std::vector<int> fRecoMichelPlaneIndex;
  std::vector<float> fRecoMichelEnergyDiffVect;
  std::vector<float> fRecoMichelEPurityVect;
  std::vector<float> fRecoMichelEnergyVect;
  std::vector<float> fRecoMichelPurityVect;
  std::vector<float> fRecoMichelcompletenessVect;
  std::vector<float> fRecoMichelGOFMeanVect;
  std::vector<float> fRecoMichelGOFSigmaVect;
  std::vector<float> fRecoMichelGOFMaxVect;
    std::vector<float> fRecoMichelMultMeanVect;
  std::vector<float> fRecoMichelMultSigmaVect;
  std::vector<float> fRecoMichelMultMaxVect;
  std::vector<float> fRecoMichelIntegralMeanVect;
  std::vector<float> fRecoMichelIntegralSigmaVect;
  std::vector<float> fRecoMichelIntegralMaxVect;
  std::vector<int> fRecoMichelNHitsPlane;
  TVector3 *fRecoMichelStartVect;

  // MC Muon
  int fMCMuonPDG;
  int fMCMuonG4ID;
  int fNDeltas;
  int fNPoints;
  int fNElectronsWithTrueHits;
  float fMCMuonEnergy;
  float fMCMuonTheta;
  float fMCMuonPhi;
  TVector3 *fMCMuonVect;
  float fMCMuonEndX;
  float fMCMuonEndY;
  float fMCMuonEndZ;
  float fMCMuonEndT;
  bool fMCMuonStopping;
  float fMCMuonEndEnergy;
  float fMCMuonEndPx;
  float fMCMuonEndPy;
  float fMCMuonEndPz;
  float fMCMuonGenX, fMCMuonGenY, fMCMuonGenZ, fMCMuonGenE;
  float fMCMuonTime;
  float fMCMuonStartX;
  float fMCMuonStartY;
  float fMCMuonStartZ;
  float fMCMuonStartT;
  float fMCMuonBendiness;
  bool fMCMuonTrigger;
  unsigned fMCMuonTriggerID;
  float fMCMuonMinDist;
  bool fMCMuonEntersTPC;
  std::string fMCMuonEndProcess;
  float fMCMuonMamophtTime;
  int fMCMuonMamophtADC;

  // Reco Muon
  int fNRecoMuonHits;
  float fRecoMuonEnergy;
  float fRecoMuonTheta;
  float fRecoMuonPhi;
  TVector3 *fRecoMuonVect;
  int fRecoMuonID;
 float fRecoMuonStartX;
  float fRecoMuonStartY;
  float fRecoMuonStartZ;
  float fRecoMuonEndX;
  float fRecoMuonEndY;
  float fRecoMuonEndZ;
  int fRecoMuonCheckHits;
  int fNRecoMuonSpacePoints;
  bool fRecoMuonIsPrimary;
  int fRecoMuonMother;
  int fNTrueHitsInRecoMuon;
  float fRecoMuonAzi;
  float fRecoMuonZen;
  int fRecoMuonNClusters;
  std::vector<int> fRecoMuonNHitsPlane;
  std::vector<float> fRecoMuonEndBendinessVect;
  TVector3 *fRecoMuonEndVect;

  // NonClustered Hits
  std::vector<float> fNonClustX;
  std::vector<float> fNonClustY;
  std::vector<float> fNonClustZ;
  std::vector<float> fNonClustDist;		// Distance of non-clustered hit to MC Michel start
  std::vector<float> fNonClustMaxAmp;
  std::vector<float> fNonClustSADC;
  std::vector<int> fNonClustPlane;
  std::vector<bool> fNonClustIsMichel;

  calo::CalorimetryAlg fCalorimetryAlg;


  // SoftwareTrigger Waveforms
  TTree* fTrigTree;
  int fTrigID;
  float fTrigTime;
  std::vector<int> fTrigMultVec;
  std::vector<float> fTrigADCVec, fTrigADCRiseVec;
  unsigned fCoincID;
  int fCoincPlanes;
  float fCoincCRTTime, fCoincPMTTime;

  // Pulses
  TTree* fPulseTree;
  int fPulseTriggerID;
  std::vector<int> fPulseCh;
  std::vector<float> fPulseChBaseline, fPulseChBaselineSigma;
  std::vector<float> fPulsePeak, fPulseArea, fPulseTStart, fPulseTPeak, fPulseTEnd, fPulsePE;


  // Light Matching Tree
  TTree* fOpTree;
  std::vector<int> fOpCh;
  std::vector<float> fRecoMichelPhotons, fRecoMuonPhotons;
  std::vector<float> fMCMichelPhotons, fMCMuonPhotons;
  std::vector<unsigned> fRecoMichelNOpHits, fRecoMuonNOpHits;

  // Temporary variables
  int countpfps;
  std::string lastchar;

  // Define input labels
  const std::string fPFParticleLabel;
  const std::string fHitLabel;
  const std::string fTrackLabel;
  const std::string fShowerLabel;
  const std::string fCaloLabel;
  const std::string fMCTruthLabel;
  const std::string fClusterLabel;
  const std::string fHitClusterLabel;
  const std::string fHitShowerLabel;
  const std::string fHitTrackLabel;
  const std::string fHitSpacePointLabel;
  const std::string fTrackCaloLabel;
  const int sUseWPlaneOnly;
  const std::string fSliceLabel;
  const float fRecombinationFactor;
  const bool fUseTrack;
  const bool fUseShower;

  const std::string fSoftwareTriggerLabel;
  const std::string fCoincidenceLabel;

  const std::vector<std::string> fOpHitLabels;
  const std::vector<float> fOpHitDelays;
  const std::vector<std::string> fSimPhotonsLabel;
  const std::string fMaMOpHTLabel;
  // Declare member data here.
    // BackTrackerService
//    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
};

sbnd::AnalyseMichels::AnalyseMichels(fhicl::ParameterSet const& p)
    : EDAnalyzer { p }, fMCMichelVect(nullptr), fRecoMichelVect(nullptr), fRecoMichelStartVect(nullptr), fMCMuonVect(nullptr), fRecoMuonVect(nullptr), fRecoMuonEndVect(nullptr)
    // Initialise out input labels by reading the fhicl parameters
    , fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
    , fPFParticleLabel(p.get<std::string>("PFParticleLabel"))
   , fHitLabel(p.get<std::string>("HitLabel")) 
   , fTrackLabel(p.get<std::string>("TrackLabel"))
   , fShowerLabel(p.get<std::string>("ShowerLabel"))
    , fCaloLabel(p.get<std::string>("CalorimetryLabel"))
   , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
   , fClusterLabel(p.get<std::string>("ClusterLabel"))
   , fHitClusterLabel(p.get<std::string>("HitClusterLabel"))
  , fHitShowerLabel(p.get<std::string>("HitShowerLabel"))
  , fHitTrackLabel(p.get<std::string>("HitTrackLabel"))
  , fHitSpacePointLabel(p.get<std::string>("HitSpacePointLabel"))
  , fTrackCaloLabel(p.get<std::string>("TrackCaloLabel"))
  , sUseWPlaneOnly(p.get<int>("UseWPlaneOnly"))
  , fSliceLabel(p.get<std::string>("SliceLabel"))
  , fRecombinationFactor(p.get<float>("RecombinationFactor"))
  , fUseTrack(p.get<bool>("UseTrack"))
  , fUseShower(p.get<bool>("UseShower"))
  , fSoftwareTriggerLabel(p.get<std::string>("SoftwareTriggerLabel"))
  , fCoincidenceLabel(p.get<std::string>("CoincidenceLabel"))
  , fOpHitLabels(p.get<std::vector<std::string>>("OpHitLabels"))
  , fOpHitDelays(p.get<std::vector<float>>("OpHitDelays"))
  , fSimPhotonsLabel(p.get<std::vector<std::string>>("SimPhotonsLabel"))
  , fMaMOpHTLabel(p.get<std::string>("MaMOpHTLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::AnalyseMichels::rotateVector(TVector3 *vector) {
  float temp = vector->Y();
  vector->SetY(vector->Z());
  vector->SetZ(temp);
  temp = vector->X();
  vector->SetX(vector->Y());
  fMCMuonVect->SetY(temp);
}



void sbnd::AnalyseMichels::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fRun = e.run(); fSubRun = e.subRun();
  fEventID = e.id().event();


  // Load the PFParticles from pandora
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  // If there are no PFParticles then give up and skip the event
  if (pfpVec.empty())
    return;

   // Accessing hits
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr <recob::Hit> > hitVect;
  if(e.getByLabel(fHitLabel, hitHandle))
    art::fill_ptr_vector(hitVect, hitHandle);

  fEventNHits = hitVect.size();

  // Access SpacePoints associated with each hit
  art::FindManyP<recob::SpacePoint> spacepointAssoc(hitVect, e, fHitSpacePointLabel);

  // Accessing slices
  art::Handle< std::vector<recob::Slice> > sliceHandle;
  std::vector<art::Ptr <recob::Slice> > sliceVect;
  if(e.getByLabel(fSliceLabel, sliceHandle))
    art::fill_ptr_vector(sliceVect, sliceHandle);

  fNSlices = sliceVect.size();

  // Selecting reco michel
  std::vector< art::Ptr<recob::Track> > michelTrack;
  std::vector< art::Ptr<recob::Shower> > michelShower;
  std::vector< art::Ptr<recob::Hit> > michelHitsNotRecod;

  // Select reco muon
  std::vector <art::Ptr<recob::Track> > muonTrack;

  // If no hits then skip event
  if (hitVect.empty())
    return;

  fNPFParticles = pfpVec.size();

  // Accessing MCParticles
  art::Handle< std::vector<simb::MCParticle> > mctruthHandle;
  std::vector< art::Ptr<simb::MCParticle> > mctruthVect;
  if(e.getByLabel(fMCTruthLabel, mctruthHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(mctruthVect, mctruthHandle);

  // Count other electrons produced by parent muon

  for (auto const &mcp: mctruthVect) {
    fMCPDG.push_back(mcp->PdgCode());
    fMCTrackID.push_back(mcp->TrackId());
  }
   
  // Get G4ID of stopping muons
  std::vector<int> muon_ids;
  for(auto const mcp : mctruthVect) {
    if(abs(mcp->PdgCode()) != 13) continue;
    muon_ids.push_back(mcp->TrackId());
  }
  for(auto const &mcp: mctruthVect) {
    if(std::find(muon_ids.begin(), muon_ids.end(), mcp->TrackId()) == muon_ids.end()) continue;
    ResetVars();
    FillMC(mcp, e, mctruthVect);
    // if(fMCMichelID == 0) continue;
    FindRecoMichel(e);
    if(fIsShower) FindRecoMichelShower(e);
    else if(fIsTrack) FindRecoMichelTrack(e);
    MatchMCRecoLight(e);
    fTree->Fill();
  }
  FillTrigTree(e);
  FillPulseTree(e);
}

void sbnd::AnalyseMichels::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Add branches to the TTree
  // Event
  fTree->Branch("run",				&fRun);
  fTree->Branch("sub,				&fSubRun");
  fTree->Branch("event.ID", 			&fEventID);
  fTree->Branch("event.NPFParticles", 		&fNPFParticles);
  fTree->Branch("event.PFParticlePDG", 		&fPFParticlePDG);
  fTree->Branch("event.NHits", 			&fEventNHits);
  fTree->Branch("event.NClusterHits",		&fEventNClusterHits);
  fTree->Branch("event.MCPDG", 			&fMCPDG);
  fTree->Branch("event.MCTrackID", 		&fMCTrackID);
  fTree->Branch("event.PFPNHits", 		&fNHitsInPFP);
  fTree->Branch("event.PFPNTrueHits", 		&fNTrueHitsInPFP);
  fTree->Branch("event.PFPNCluster", 		&fNClusters);
  fTree->Branch("event.PFPID",			&fPFPIDs);
  fTree->Branch("event.NRecoElectrons",		&fNRecoElectrons);
  fTree->Branch("event.WrongRecoMuonEnd",	&fWrongMuonEnd);
  fTree->Branch("event.ClusterID",		&fClusterID);
  fTree->Branch("event.ClusterPlane",		&fClusterPlane);
  fTree->Branch("event.NRecoMuons",		&fNMuons);
  fTree->Branch("event.RecoMuonMichelDist",	&fRecoMuonMichelDist);
  fTree->Branch("event.PFPIsPrimary",		&fPFPIsPrimary);
  fTree->Branch("event.PFPMother",		&fPFPMother);
  fTree->Branch("event.NElectronsWithTrueHits",	&fNElectronsWithTrueHits);
  fTree->Branch("event.NMuonsWithTrueHits",	&fNMuonsWithTrueHits);
  fTree->Branch("event.MaxSADCX",		&fMaxSADCX);
  fTree->Branch("event.MaxSADCY",               &fMaxSADCY);
  fTree->Branch("event.MaxSADCZ",		&fMaxSADCZ);
  fTree->Branch("event.MaxSADCDist",               &fMaxSADCDist);
  fTree->Branch("event.MaxAmpX",               &fMaxAmpX);
  fTree->Branch("event.MaxAmpY",               &fMaxAmpY);
  fTree->Branch("event.MaxAmpZ",               &fMaxAmpZ);
  fTree->Branch("event.MaxAmpDist",            	&fMaxAmpDist);
  fTree->Branch("event.AmpMean",		&fAmpMean);
  fTree->Branch("event.AmpSigma",		&fAmpSigma);
  fTree->Branch("event.SADCMean",		&fSADCMean);
  fTree->Branch("event.SADCSigma",		&fSADCSigma);
  fTree->Branch("event.NSlices",		&fNSlices);
  fTree->Branch("event.NNonclustMichelHits",	&fNMichelNonClust);

  // MC Michel
  fTree->Branch("mcMichel.ID", 			&fMCMichelID);
  fTree->Branch("mcMichel.Energy", 		&fMCMichelEnergy);
  fTree->Branch("mcMichel.Theta", 		&fMCMichelTheta);
  fTree->Branch("mcMichel.Phi", 		&fMCMichelPhi);
  fTree->Branch("mcMichel.RelTheta", 		&fMCMichelRelTheta);
  fTree->Branch("mcMichel.RelPhi", 		&fMCMichelRelPhi);
  fTree->Branch("mcMichel.Vect", 		&fMCMichelVect);
  fTree->Branch("mcMichel.Length", 		&fMCMichelLength);
  fTree->Branch("mcMichel.StartT",    &fMCMichelStartT);
  fTree->Branch("mcMichel.StartX",		&fMCMichelStartX);
  fTree->Branch("mcMichel.StartY",		&fMCMichelStartY);
  fTree->Branch("mcMichel.StartZ",		&fMCMichelStartZ);
  fTree->Branch("mcMichel.RelAngle",		&fMCMichelRelAngle);
  fTree->Branch("mcMichel.NPoints",		&fMCMichelNPoints);
  fTree->Branch("mcMichel.EnergyFrac",		&fMCMichelEnergyFrac);
  fTree->Branch("mcMichel.Process",		&fMCMichelStartProcess);
  fTree->Branch("mcMichel.MamophtTime",		&fMCMichelMamophtTime);
  fTree->Branch("mcMichel.MamophtADC",          &fMCMichelMamophtADC);
 
 // Reco Michel
  fTree->Branch("recoMichel.NHits", 		&fNHitsInRecoMichel);
  fTree->Branch("recoMichel.TotalHits", 	&fNTotalMichelHits);
  fTree->Branch("recoMichel.NShowers", 		&fNShowers);
  fTree->Branch("recoMichel.NTracks", 		&fNTracks);
  fTree->Branch("recoMichel.Exists", 		&fRecoMichel);
  fTree->Branch("recoMichel.ShowerPurity", 	&fShowerPurity);
  fTree->Branch("recoMichel.TrackPurity", 	&fTrackPurity);
  fTree->Branch("recoMichel.ShowerBestPlane", 	&fShowerBestPlane);
  fTree->Branch("recoMichel.IsTrack", 		&fIsTrack);
  fTree->Branch("recoMichel.IsShower", 		&fIsShower);
  fTree->Branch("recoMichel.NTrueHits", 	&fNTrueHitsInRecoMichel);
  fTree->Branch("recoMichel.Energy", 		&fRecoMichelEnergy);
  fTree->Branch("recoMichel.Theta", 		&fRecoMichelTheta);
  fTree->Branch("recoMichel.Phi", 		&fRecoMichelPhi);
  fTree->Branch("recoMichel.RelTheta", 		&fRecoMichelRelTheta);
  fTree->Branch("recoMichel.RelPhi", 		&fRecoMichelRelPhi);
  fTree->Branch("recoMichel.RelAngle",		&fRecoMichelRelAngle);
  fTree->Branch("recoMichel.Vect", 		&fRecoMichelVect);
  fTree->Branch("recoMichel.StartX",		&fRecoMichelStartX);
  fTree->Branch("recoMichel.StartY",		&fRecoMichelStartY);
  fTree->Branch("recoMichel.StartZ",		&fRecoMichelStartZ);
  fTree->Branch("recoMichel.ID", 		&fRecoMichelID);
  fTree->Branch("recoMichel.Mother",		&fRecoMichelMother);
  fTree->Branch("recoMichel.NSpacePoints",	&fNRecoMichelSpacePoints);
  fTree->Branch("recoMichel.Length",		&fRecoMichelLength);
  fTree->Branch("recoMichel.Azimuth",		&fRecoMichelAzi);
  fTree->Branch("recoMichel.Zenith",		&fRecoMichelZen);
  fTree->Branch("recoMichel.ClosestHitDist",	&fRecoMichelCloseProximity);
  fTree->Branch("recoMichel.NClusters",		&fRecoMichelNClusters);
  fTree->Branch("recoMichel.EnergySigma",	&fRecoMichelEnergySigma);
  fTree->Branch("recoMichel.dEdX",		&fRecoMicheldEdX);
  fTree->Branch("recoMichel.HitIntegral",	&fRecoMichelHitIntegral);
  fTree->Branch("recoMichel.HitPeakTime",	&fRecoMichelHitPeakTime);
  fTree->Branch("recoMichel.HitNElectrons",	  &fRecoMichelHitNElectrons);
  fTree->Branch("recoMichel.HitPlane",		&fRecoMichelHitPlane);
  fTree->Branch("recoMIchel.EnergyU",		&fRecoMichelEnergyU);
  fTree->Branch("recoMIchel.EnergyV", 		 &fRecoMichelEnergyV);
  fTree->Branch("recoMIchel.EnergyW",		  &fRecoMichelEnergyW);
  fTree->Branch("recoMichel.EnergyDiffU",	&fRecoMichelEnergyDiffU);
  fTree->Branch("recoMichel.EnergyDiffV",       &fRecoMichelEnergyDiffV);
  fTree->Branch("recoMichel.EnergyDiffW",       &fRecoMichelEnergyDiffW);
  fTree->Branch("recoMichel.IntegralMaxU",	&fRecoMichelIntegralMaxU);
  fTree->Branch("recoMichel.IntegralMeanU",	&fRecoMichelIntegralMeanU);
  fTree->Branch("recoMichel.IntegralSigmaU",	&fRecoMichelIntegralSigmaU);
  fTree->Branch("recoMichel.IntegralMaxV",      &fRecoMichelIntegralMaxV);
  fTree->Branch("recoMichel.IntegralMeanV",     &fRecoMichelIntegralMeanV);
  fTree->Branch("recoMichel.IntegralSigmaV",    &fRecoMichelIntegralSigmaV);
  fTree->Branch("recoMichel.IntegralMaxW",      &fRecoMichelIntegralMaxW);
  fTree->Branch("recoMichel.IntegralMeanW",     &fRecoMichelIntegralMeanW);
  fTree->Branch("recoMichel.IntegralSigmaW",    &fRecoMichelIntegralSigmaW);
  fTree->Branch("recoMichel.MultiplicityMaxV",   &fRecoMichelMultiplicityMaxV);
  fTree->Branch("recoMichel.MultiplicityMeanV",   &fRecoMichelMultiplicityMeanV);
  fTree->Branch("recoMichel.MultiplicitySigmaV",   &fRecoMichelMultiplcitySigmaV);
  fTree->Branch("recoMichel.GOFMaxV",		   &fRecoMichelGOFMaxV);
  fTree->Branch("recoMichel.GOFMeanV",		   &fRecoMichelGOFSigmaV);
  fTree->Branch("recoMichel.GOFSigmaV",		   &fRecoMichelGOFMeanV);
  fTree->Branch("recoMichel.MultiplicityMaxU",	   &fRecoMichelMultiplicityMaxU);
  fTree->Branch("recoMichel.MultiplicityMeanU",	   &fRecoMichelMultiplicityMeanU);
  fTree->Branch("recoMichel.MultiplicitysigmaU",   &fRecoMichelMultiplcitySigmaU);
  fTree->Branch("recoMichel.GOFMaxU",		   &fRecoMichelGOFMaxU);
  fTree->Branch("recoMichel.GOFMeanU",		   &fRecoMichelGOFSigmaU);
  fTree->Branch("recoMichel.GOFSigmaU",		   &fRecoMichelGOFMeanU);
  fTree->Branch("recoMichel.MultiplicityMaxW",	   &fRecoMichelMultiplicityMaxW);
  fTree->Branch("recoMichel.MultiplicityMeanW",	   &fRecoMichelMultiplicityMeanW);
  fTree->Branch("recoMichel.MultiplicitySigmaW",	  &fRecoMichelMultiplcitySigmaW);
  fTree->Branch("recoMichel.GOFMaxW",		   &fRecoMichelGOFMaxW);
  fTree->Branch("recoMichel.GOFMeanW",		   &fRecoMichelGOFSigmaW);
  fTree->Branch("recoMichel.GOFSigmaW",		   &fRecoMichelGOFMeanW);
  fTree->Branch("recoMichel.HitMultiplicity",	&fRecoMichelhitMultiplicity);
  fTree->Branch("recoMichel.HitGOF",		&fRecoMichelHitGOF);
  fTree->Branch("recoMichel.EDiffMin",		&fRecoMichelEDiffMin);
  fTree->Branch("recoMichel.EDiffMinMult",	&fRecoMichelEDiffMinMult);
  fTree->Branch("recoMichel.EDiffMinGOF",	&fRecoMichelEDiffMinGOF);
  fTree->Branch("recoMichel.PurityU",		&fRecoMichelPurityU);
  fTree->Branch("recoMichel.PurityV",           &fRecoMichelPurityV);
  fTree->Branch("recoMichel.PurityW",           &fRecoMichelPurityW);
  fTree->Branch("recoMichel.Purity",		&fRecoMichelPurity);
  fTree->Branch("recoMichel.Completeness",	&fRecoMichelCompleteness);
  fTree->Branch("recoMichel.CompletenessW",      &fRecoMichelCompletenessW);
  fTree->Branch("recoMichel.EPurityU",		&fRecoMichelEPurityU);
  fTree->Branch("recoMichel.EPurityV",          &fRecoMichelEPurityV);
  fTree->Branch("recoMichel.EPurityW",          &fRecoMichelEPurityW);
  fTree->Branch("recoMichel.PlaneIndex",	&fRecoMichelPlaneIndex);
  fTree->Branch("recoMichel.EnergyVect",		&fRecoMichelEnergyVect);
  fTree->Branch("recoMichel.EnergyDiffVect",		&fRecoMichelEnergyDiffVect);
  fTree->Branch("recoMichel.EPurityVect",		&fRecoMichelEPurityVect);
  fTree->Branch("recoMichel.EPurity",			&fRecoMichelEPurity);
  fTree->Branch("recoMichel.PurityVect",          	&fRecoMichelPurityVect);
  fTree->Branch("recoMichel.CompletenessVect",          &fRecoMichelcompletenessVect);
  fTree->Branch("recoMichel.GOFMeanVect",          	&fRecoMichelGOFMeanVect);
  fTree->Branch("recoMichel.GOFSigmaVect",          	&fRecoMichelGOFSigmaVect);
  fTree->Branch("recoMichel.GOFMaxVect",          	&fRecoMichelGOFMaxVect);
    fTree->Branch("recoMichel.MultMeanVect",          	&fRecoMichelMultMeanVect);
  fTree->Branch("recoMichel.MultSigmaVecy",          	&fRecoMichelMultSigmaVect);
  fTree->Branch("recoMichel.MultMaxVect",          	&fRecoMichelMultMaxVect);
  fTree->Branch("recoMichel.IntegralMeanVect",          	&fRecoMichelIntegralMeanVect);
  fTree->Branch("recoMichel.IntegralSigmaVect",          	&fRecoMichelIntegralSigmaVect);
  fTree->Branch("recoMichel.IntegralMaxVect",		&fRecoMichelIntegralMaxVect);
  fTree->Branch("recoMichel.NHitsVect",			&fRecoMichelNHitsPlane);
  fTree->Branch("recoMichel.StartDir",			&fRecoMichelStartVect);

  // MC Muon
  fTree->Branch("mcMuon.PDG",			&fMCMuonPDG);
  fTree->Branch("mcMuon.G4ID",			&fMCMuonG4ID);
  fTree->Branch("mcMuon.NDeltas", 		&fNDeltas);
  fTree->Branch("mcMuon.NPoints", 		&fNPoints);
  fTree->Branch("mcMuon.Energy", 		&fMCMuonEnergy);
  fTree->Branch("mcMuon.Theta", 		&fMCMuonTheta);
  fTree->Branch("mcMuon.Phi", 			&fMCMuonPhi);
  fTree->Branch("mcMuon.Vect", 			&fMCMuonVect);
  fTree->Branch("mcMuon.EndX",			&fMCMuonEndX);
  fTree->Branch("mcMuon.EndY",			&fMCMuonEndY);
  fTree->Branch("mcMuon.EndZ",			&fMCMuonEndZ);
  fTree->Branch("mcMuon.EndT",      &fMCMuonEndT);
  fTree->Branch("mcMuon.Stopping",   &fMCMuonStopping);
  fTree->Branch("mcMuon.EndEnergy",		&fMCMuonEndEnergy);
  fTree->Branch("mcMuon.EndPx",			&fMCMuonEndPx);
  fTree->Branch("mcMuon.EndPy",			&fMCMuonEndPy);
  fTree->Branch("mcMuon.EndPz",			&fMCMuonEndPz);
  fTree->Branch("mcMuon.T",		&fMCMuonTime);
  fTree->Branch("mcMuon.GenX",		&fMCMuonGenX);
  fTree->Branch("mcMuon.GenY",		&fMCMuonGenY);
  fTree->Branch("mcMuon.GenZ",		&fMCMuonGenZ);
  fTree->Branch("mcMuon.GenE",    &fMCMuonGenE);
  fTree->Branch("mcMuon.StartX",		&fMCMuonStartX);
  fTree->Branch("mcMuon.StartY",		&fMCMuonStartY);
  fTree->Branch("mcMuon.StartZ",		&fMCMuonStartZ);
  fTree->Branch("mcMuon.Bendiness",		&fMCMuonBendiness);
  fTree->Branch("mcMuon.Trigger",		&fMCMuonTrigger);
  fTree->Branch("mcMuon.TriggerID",		&fMCMuonTriggerID);
  fTree->Branch("mcMuon.EndProcess",		&fMCMuonEndProcess);
  fTree->Branch("mcMuon.EntersTPC",		&fMCMuonEntersTPC);
  fTree->Branch("mcMuon.MinDist",		&fMCMuonMinDist);
  fTree->Branch("mcMuonMamophtTime",          &fMCMuonMamophtTime);
  fTree->Branch("mcMuonMamophtADC",          &fMCMuonMamophtADC);

  // Reco Muon
  fTree->Branch("recoMuon,NHits",		&fNRecoMuonHits);
  fTree->Branch("recoMuon.Energy", 		&fRecoMuonEnergy);
  fTree->Branch("recoMuon.Theta", 		&fRecoMuonTheta);
  fTree->Branch("recoMuon.Phi", 		&fRecoMuonPhi);
  fTree->Branch("recoMuon.Vect", 		&fRecoMuonVect);
  fTree->Branch("recoMuon.ID", 			&fRecoMuonID);
  fTree->Branch("recoMuon.StartX", 		&fRecoMuonStartX);
  fTree->Branch("recoMuon.StartY",		&fRecoMuonStartY);
  fTree->Branch("recoMuon.StartZ",		&fRecoMuonStartZ);
  fTree->Branch("recoMuon.EndX",		&fRecoMuonEndX);
  fTree->Branch("recoMuon.EndY",		&fRecoMuonEndY);
  fTree->Branch("recoMuon.EndZ",		&fRecoMuonEndZ);
  fTree->Branch("recoMuon.CheckHits",		&fRecoMuonCheckHits);
  fTree->Branch("recoMuon.IsPrimary",		&fRecoMuonIsPrimary);
  fTree->Branch("recoMuon.Mother",		&fRecoMuonMother);
  fTree->Branch("recoMuon.NTrueHits",		&fNTrueHitsInRecoMuon);
  fTree->Branch("recoMuon.Azimuth",		&fRecoMuonAzi);
  fTree->Branch("recoMuon.Zenith",		&fRecoMuonZen);
  fTree->Branch("recoMuon.EndX",		&fRecoMuonEndX);
  fTree->Branch("recoMuon.EndY",		&fRecoMuonEndY);
  fTree->Branch("recoMuon.EndZ",		&fRecoMuonEndZ);
  fTree->Branch("recoMuon.NClusters",		&fRecoMuonNClusters);
  fTree->Branch("recoMuon.NHitsVect",		&fRecoMuonNHitsPlane);
  fTree->Branch("recoMuon.EndBendinessVect",	&fRecoMuonEndBendinessVect);
  fTree->Branch("recoMuon.EndDir",		&fRecoMuonEndVect);

  fTree->Branch("nonClust.X",			&fNonClustX);
  fTree->Branch("nonClust.Y",                   &fNonClustY);
  fTree->Branch("nonClust.Z",                   &fNonClustZ);
 fTree->Branch("nonClust.Dist",                   &fNonClustDist);
  fTree->Branch("nonClust.MaxAmp",                &fNonClustMaxAmp);
  fTree->Branch("nonClust.SADC",                   &fNonClustSADC);
  fTree->Branch("nonClust.Plane",                   &fNonClustPlane);
  fTree->Branch("nonClust.IsMichel",		&fNonClustIsMichel);

  fTree->Branch("countpfpf",				&countpfps);
  fTree->Branch("lastchar",				&lastchar);

  fTrigTree = tfs->make<TTree>("trig_tree", "Output Tree");
  fTrigTree->Branch("run",		&fRun);
  fTrigTree->Branch("sub",		&fSubRun);
  fTrigTree->Branch("evt",		&fEventID);
  fTrigTree->Branch("id",		&fTrigID);
  fTrigTree->Branch("trig_time",	&fTrigTime);
  fTrigTree->Branch("mult_vec",		&fTrigMultVec);
  fTrigTree->Branch("adc_vec",		&fTrigADCVec);
  fTrigTree->Branch("adc_rise_vec",	&fTrigADCRiseVec);
  fTrigTree->Branch("coinc_id",		&fCoincID);
  fTrigTree->Branch("coinc_planes",     &fCoincPlanes);
  fTrigTree->Branch("coinc_crttimes",   &fCoincCRTTime);
  fTrigTree->Branch("coinc_pmttimes",   &fCoincPMTTime);

  fPulseTree = tfs->make<TTree>("pulse_tree", "Output Tree");
  fPulseTree->Branch("run",              &fRun);
  fPulseTree->Branch("sub",              &fSubRun);
  fPulseTree->Branch("evt",              &fEventID);
  fPulseTree->Branch("id",               &fPulseTriggerID);
  fPulseTree->Branch("ch",		&fPulseCh);
  fPulseTree->Branch("ch_base",		&fPulseChBaseline);
  fPulseTree->Branch("ch_sigma",	&fPulseChBaselineSigma);
  fPulseTree->Branch("pulse_peak",	&fPulsePeak);
  fPulseTree->Branch("pulse_area",	&fPulseArea);
  fPulseTree->Branch("pulse_pe",	&fPulsePE);
  fPulseTree->Branch("pulse_tstart",	&fPulseTStart);
  fPulseTree->Branch("pulse_tend",	&fPulseTEnd);
  fPulseTree->Branch("pulse_tpeak",	&fPulseTPeak);

  fOpTree = tfs->make<TTree>("op_tree", "Output Tree");
  fOpTree->Branch("run",              &fRun);
  fOpTree->Branch("sub",              &fSubRun);
  fOpTree->Branch("evt",              &fEventID);
  fOpTree->Branch("muon_id",		&fMCMuonG4ID);
  fOpTree->Branch("michel_id",		&fMCMichelID);
  fOpTree->Branch("ch",			&fOpCh);
  fOpTree->Branch("michel_mc_phot",		&fMCMichelPhotons);
  fOpTree->Branch("michel_reco_phot",		&fRecoMichelPhotons);
  fOpTree->Branch("muon_mc_phot",             &fMCMuonPhotons);
  fOpTree->Branch("muon_reco_phot",           &fRecoMuonPhotons);
  fOpTree->Branch("muon_nophits",		&fRecoMuonNOpHits);
  fOpTree->Branch("michel_nophits",               &fRecoMichelNOpHits);

  fMCMichelEnergyHist = tfs->make<TH1D>("mcMichelEnergyHist", "Energy of MC Michels; Energy; Events", 40, 2, 1);
  fMCMuonEnergyHist = tfs->make<TH1D>("mcMuonEnergyHist", "Energy of MC muons; Energy; Events", 40, 2, 1);
  fRecoMuonEnergyHist = tfs->make<TH1D>("recoMichelEnergyHist", "Energy of reconstructed Michels; Energy; Events", 40, 2, 1);
  fRecoMichelEnergyHist = tfs->make<TH1D>("recoMuonEnergyHist", "Energy of reconstructed muons; Energy; Events", 40, 2, 1);
  fMCMichelThetaHist = tfs->make<TH1D>("mcMichelThetahist", "Theta of MC Michels; Theta; Events", 40, 2, 1);
  fMCMuonThetaHist = tfs->make<TH1D>("mcMuonThetaHist", "Theta of MC muons; Theta; Events", 40, 2, 1);
  fRecoMichelThetaHist = tfs->make<TH1D>("recoMichelThetaHist", "Theta of reconstructed Michels; Theta; Events", 40, 2, 1);
  fRecoMuonThetaHist = tfs->make<TH1D>("recoMuonThetaHist", "Theta of reconstructed muons; Theta; Events", 40, 2, 1);
  fMCMichelPhiHist = tfs->make<TH1D>("mcMichelPhiHist", "Phi of MC Michels; Phi; Events", 40, 2, 1);
  fMCMuonPhiHist = tfs->make<TH1D>("mcMuonPhiHist", "Phi of MC muons; Phi; Events", 40, 2, 1);
  fRecoMichelPhiHist = tfs->make<TH1D>("recoMichelPhiHist", "Phi of reconstructed Michels; Phi; Events", 40, 2, 1);
  fRecoMuonPhiHist = tfs->make<TH1D>("recoMuonHist", "Phi of reconstructed muons; Phi; Events", 40, 2, 1);
  fMCMichelRelThetaHist = tfs->make<TH1D>("mcMichelRelThetaHist", "Theta difference between MC Michel and muon; Theta; Events", 40, 2, 1);
  fRecoMichelRelThetaHist = tfs->make<TH1D>("recoMichelRelThetaHist", "Theta difference between reconstructed Michel and muon; Theta; Events", 40, 2, 1);
  fMCMichelRelPhiHist = tfs->make<TH1D>("mcMichelRelPhiHist", "Phi difference between MC Michel and muon; Phi; Events", 40, 2, 1);
  fRecoMichelRelPhiHist = tfs->make<TH1D>("recoMichelRelHist", "Phi difference between reconstructed Michel and muon; Phi; Events", 40, 2, 1);
  fMCMichelLengthHist = tfs->make<TH1D>("mcMichelLengthHist", "Length of MC Michels; Length; Events", 40, 4, 3);
}

void sbnd::AnalyseMichels::endJob()
{}

void sbnd::AnalyseMichels::ResetVars()
{
  // Reset all of our variables to 0 or empty vectors
  // This ensures things are not kept from the previous event
  fNPFParticles = 0;
  fNDeltas = 0;
  fMCMichelID = -1;
  fNTotalMichelHits= 0;
  fNHitsInRecoMichel = 0;
  fEventNHits = 0;
  fEventNClusterHits = 0;
  fPurity = 0;
  fCompleteness = 0;
  fNTrueHitsInRecoMichel = 0;
  fNPoints = 0;
  fNShowers = 0;
  fNTracks = 0;
  fRecoMichel = false;
  fIsTrack = false;
  fIsShower = false;
  fShowerPurity = 0;
  fTrackPurity = 0;
  fShowerBestPlane = 0;
  fRecoMichelID = 0;
  fRecoMuonID = -1;
  fNRecoMuonHits = 0;
  fNRecoElectrons = 0;
  fWrongMuonEnd = false;
  fMCMichelStartX = 0;
  fMCMichelStartY = 0;
  fMCMichelStartZ = 0;
  fMCMichelStartT = 0.;
  fRecoMichelStartX = 0;
  fRecoMichelStartY = 0;
  fRecoMichelStartZ = 0;
  fMCMuonPDG = -9999;
  fMCMuonG4ID = -9999;
  fMCMuonEndX = 0;
  fMCMuonEndY = 0;
  fMCMuonEndZ = 0;
  fMCMuonEndT = 0.;
  fMCMuonStopping = false;
  fRecoMuonStartX = 0;
  fRecoMuonStartY = 0;
  fRecoMuonStartZ = 0;
  fRecoMuonEndX = 0;
  fRecoMuonEndY = 0;
  fRecoMuonEndZ = 0;
  fRecoMichelMother = -1;
  fRecoMuonCheckHits = 0;
  fNRecoMichelSpacePoints = 0;
  fRecoMichelNHitsPlane = {-1, -1, -1};
  fNMuons = 0;
  fRecoMuonIsPrimary = false;
  fRecoMuonMother = -1;
  fRecoMuonMichelDist = 0;
  fNTrueHitsInRecoMuon = 0;
  fNMuonsWithTrueHits = 0;
  fNElectronsWithTrueHits = 0;
  fRecoMichelLength = 0;
  fRecoMichelAzi = 0;
  fRecoMichelZen = 0;
  fRecoMuonTheta = 0;
  fRecoMuonPhi = 0;
  fRecoMuonAzi = 0;
  fRecoMuonZen = 0;
  fMaxSADCX = 0;
  fMaxSADCY = 0;
  fMaxSADCZ = 0;
  fMaxSADCDist = 0;
  fMaxAmpX = 0;
  fMaxAmpY = 0;
  fMaxAmpZ = 0;
  fMaxAmpDist = 0;
  fMCMichelRelAngle = 0;
  fMCMichelNPoints = 0;
  fMCMuonEndEnergy = 0;
  fMCMichelEnergyFrac = 0;
  fMCMuonEndPx = 0;
  fMCMuonEndPy = 0;
  fMCMuonEndPz = 0;
  fMCMuonGenX = 0.; fMCMuonGenY = 0.; fMCMuonGenZ= 0.; fMCMuonGenE = -9999.;
  fMCMuonTime = 10000.;
  fMCMuonStartX = 1000;
  fMCMuonStartY = 1000;
  fMCMuonStartZ = 1000;
  fRecoMichelCloseProximity = -1;
  fRecoMichelNClusters = 0;
  fRecoMuonNClusters = 0;
  fMCMuonBendiness = 0;
  fRecoMichelEnergySigma = 0;
  fRecoMicheldEdX.clear();
  fAmpMean = 0;
  fAmpSigma = 0;
  fSADCMean = 0;
  fSADCSigma = 0;
  fNMichelNonClust = 0;
   fRecoMichelMultiplicityMaxV = 0;
   fRecoMichelMultiplicityMeanV = 0;
   fRecoMichelMultiplcitySigmaV = 0;
   fRecoMichelGOFMaxV = 0;
   fRecoMichelGOFSigmaV = 0;
   fRecoMichelGOFMeanV = 0;
   fRecoMichelMultiplicityMaxU = 0;
   fRecoMichelMultiplicityMeanU = 0;
   fRecoMichelMultiplcitySigmaU = 0;
   fRecoMichelGOFMaxU = 0;
   fRecoMichelGOFSigmaU = 0;
   fRecoMichelGOFMeanU = 0;
   fRecoMichelMultiplicityMaxW = 0;
   fRecoMichelMultiplicityMeanW = 0;
   fRecoMichelMultiplcitySigmaW = 0;
   fRecoMichelGOFMaxW = 0;
   fRecoMichelGOFSigmaW = 0;
   fRecoMichelGOFMeanW = 0;
  countpfps = 0;
  fNSlices = 0;
  fRecoMuonEndBendinessVect = {-1, -1, -1};
  fRecoMichelhitMultiplicity.clear();
  fRecoMichelHitGOF.clear();
  fRecoMichelPurityU = 0;
  fRecoMichelPurityV = 0;
  fRecoMichelPurityW = 0;
  fRecoMichelPurity = -1.;
  fRecoMichelCompleteness = -1.;
  fRecoMichelCompletenessW = -1.;
  fIsTrack = false;
  fIsShower = false;
  fRecoMichelEPurityU = 0;
  fRecoMichelEPurityV = 0;
  fRecoMichelEPurityW = 0;
  fMCMuonEntersTPC = false;
  fMCMuonMinDist = 99999.;
  fMCMuonMamophtTime = -99999.;
  fMCMichelMamophtTime = -99999.;
  fMCMichelMamophtADC = -9999;
  fMCMuonMamophtADC = -9999;
  fRecoMichelPlaneIndex = {0, 1, 2};
  fRecoMichelEPurityVect = {-1, -1, -1};
  fRecoMichelEPurity = 0;
  fRecoMichelEnergyVect = {-1 , -1 , -1};
  fRecoMichelEnergyDiffVect = {-1 , -1 , -1};
  fRecoMichelPurityVect = {-1, -1, -1};
  fRecoMichelcompletenessVect = {-1, -1, -1};
  fRecoMichelGOFMeanVect = {-1, -1, -1};
  fRecoMichelGOFSigmaVect = {-1, -1, -1};
  fRecoMichelGOFMaxVect = {-1, -1, -1};
  fRecoMichelMultMeanVect = {-1, -1, -1};
  fRecoMichelMultSigmaVect = {-1, -1, -1};
  fRecoMichelMultMaxVect = {-1, -1, -1};
  fRecoMichelIntegralMeanVect = {-1, -1, -1};
  fRecoMichelIntegralSigmaVect = {-1, -1, -1};
  fRecoMichelIntegralMaxVect = {-1, -1, -1};

  fRecoMichelEnergyDiffU = 0;
  fRecoMichelEnergyDiffV = 0;
  fRecoMichelEnergyDiffW = 0;
  fRecoMichelIntegralMaxU = 0;
  fRecoMichelIntegralMeanU = 0;
  fRecoMichelIntegralSigmaU = 0;
  fRecoMichelIntegralMaxV = 0;
  fRecoMichelIntegralMeanV = 0;
  fRecoMichelIntegralSigmaV = 0;
  fRecoMichelIntegralMaxW = 0;
  fRecoMichelIntegralMeanW = 0;
  fRecoMichelIntegralSigmaW = 0;
  fRecoMichelEDiffMin = 0;
  fRecoMichelEDiffMinGOF = 0;
  fRecoMichelEDiffMinMult = 0;

  fRecoMichelHitIntegral.clear();
  fRecoMichelHitPeakTime.clear();
  fRecoMichelHitNElectrons.clear();
  fRecoMichelHitPlane.clear();
  fRecoMichelEnergyU = 0;
  fRecoMichelEnergyV = 0;
  fRecoMichelEnergyW = 0;

  fPFParticlePDG.clear();
  fMCPDG.clear();
  fMCTrackID.clear();
  fNHitsInPFP.clear();
  fNTrueHitsInPFP.clear();
  fNClusters.clear();
  fPFPIDs.clear();
  fClusterID.clear();
  fClusterPlane.clear();

  fMCMichelEnergy = 0;
  fMCMuonEnergy = 0;
  fRecoMuonEnergy = 0;
  fRecoMichelEnergy = 0;
  fMCMichelTheta = 0;
  fMCMuonTheta = 0;
  fRecoMichelTheta = 0;
  fMCMuonPhi = 0;
  fRecoMichelPhi = 0;
  fRecoMuonPhi = 0;
  fMCMichelRelTheta = 0;
  fRecoMichelRelTheta = 0;
  fMCMichelRelPhi = 0;
  fRecoMichelRelPhi = 0;
  fRecoMichelRelAngle = 0.0;
  fMCMichelLength = 0;
   lastchar = "a";
  fPFPIsPrimary.clear();
  fPFPMother.clear();

  fRecoMuonEndX = 0;
  fRecoMuonEndY = 0;
  fRecoMuonEndZ = 0;

  fNonClustX.clear();
  fNonClustY.clear();
  fNonClustZ.clear();
  fNonClustDist.clear();             // Distance of non-clustered hit to MC Michel start
  fNonClustMaxAmp.clear();
  fNonClustSADC.clear();
  fNonClustPlane.clear();
  fNonClustIsMichel.clear();
}


void sbnd::AnalyseMichels::FillMC(
  const art::Ptr<simb::MCParticle>& mcp,
  const art::Event& e, 
  std::vector<art::Ptr<simb::MCParticle>>& mctruthVect)
{
  fMCMuonPDG = mcp->PdgCode();
  fMCMuonG4ID = mcp->TrackId();
  fNPoints = mcp->NumberTrajectoryPoints();
  unsigned int i =  fNPoints/2;
  unsigned int j = fNPoints-1;
  unsigned int firstquart = fNPoints/4;
  unsigned int lastquart = 3*fNPoints/4;
  TVector3 tempvect = mcp->Position(0).Vect();
  TVector3 mcmuonend = mcp->EndPosition().Vect();
  fMCMuonGenX = mcp->Position().X();
  fMCMuonGenY = mcp->Position().Y();
  fMCMuonGenZ = mcp->Position().Z();
  fMCMuonGenE = mcp->E(0);
  fMCMuonEndProcess = mcp->EndProcess();

  float mu_trig_time = mcp->T() / 1000.;
  for(int pos=0;pos<fNPoints;pos++) {
    if((abs(mcp->Position(pos).X()) < 200.) &&
       (abs(mcp->Position(pos).Y()) < 200.) &&
       (mcp->Position(pos).Z() > 0.) &&
       (mcp->Position(pos).Z() < 500.) && !fMCMuonEntersTPC) {
      fMCMuonTime = mcp->T(pos) / 1000.; // us relative to beam spill
      fMCMuonStartX = mcp->Position(pos).X();
      fMCMuonStartY = mcp->Position(pos).Y();
      fMCMuonStartZ = mcp->Position(pos).Z();
      fMCMuonEnergy = mcp->E(pos);
      fMCMuonEntersTPC = true;
      mu_trig_time = mcp->T(pos) / 1000.;
      break;
    }
  }
  TVector3 vect1 = (mcp->Position(firstquart).Vect() - mcp->Position(0).Vect()).Unit();
  TVector3 vect2 = (mcp->Position(i).Vect() - mcp->Position(firstquart).Vect()).Unit();
  TVector3 vect3 = (mcp->Position(lastquart).Vect() - mcp->Position(i).Vect()).Unit();
  TVector3 vect4 = (mcp->Position(j).Vect() - mcp->Position(lastquart).Vect()).Unit();
  fMCMuonBendiness = (vect1.Dot(vect2) + vect2.Dot(vect3) + vect3.Dot(vect4))/3;

  *fMCMuonVect = mcp->Position(i).Vect() - mcp->Position(j).Vect();
  TVector3 fMCMuonEndVect = mcp->Position(j).Vect();
  fMCMuonEndX = mcp->EndX();
  fMCMuonEndY = mcp->EndY();
  fMCMuonEndZ = mcp->EndZ();
  fMCMuonEndT = mcp->T(mcp->NumberTrajectoryPoints()-1) / 1000.;
  fMCMuonStopping = ((fMCMuonEndX) > 200 || abs(fMCMuonEndY) > 200 || fMCMuonEndZ < 0. || fMCMuonEndZ > 500) ? false : true;
  TVector3 fMCMuonStartVect = mcp->Position(0).Vect();
  rotateVector(fMCMuonVect);
  fMCMuonTheta = (fMCMuonVect->Theta()) * 180 / M_PI -90; 		// Minus sign is there as muons are coming down;
  rotateVector(fMCMuonVect);
  fMCMuonPhi = (fMCMuonVect->Theta()) * 180 / M_PI -90;
  rotateVector(fMCMuonVect);
  fMCMuonEnergy = mcp->E() * 1000;				// Scale energy to MeV
  fMCMuonEndEnergy = mcp->EndE() * 1000;
  fMCMuonEnergyHist->Fill(fMCMuonEnergy);
  fMCMuonThetaHist->Fill(fMCMuonTheta);
  fMCMuonPhiHist->Fill(fMCMuonPhi);
  fMCMuonEndPx = mcp->EndPx();
  fMCMuonEndPy = mcp->EndPy();
  fMCMuonEndPz = mcp->EndPz();
  
  for(auto& mcp2 : mctruthVect) {
    if(abs(mcp2->PdgCode()) != 11 || mcp2->Mother() != fMCMuonG4ID ||
       abs(mcp->EndX() - mcp2->Position().X()) > 5. ||
       abs(mcp->EndY() - mcp2->Position().Y()) > 5. ||
       abs(mcp->EndZ() - mcp2->Position().Z()) > 5. ) continue;
    fMCMichelStartProcess = mcp2->Process();
    fMCMichelNPoints = mcp2->NumberTrajectoryPoints();  
    unsigned int i  = mcp2->NumberTrajectoryPoints() -1;
    unsigned int j = mcp2->NumberTrajectoryPoints() / 2;
    fMCMichelID = mcp2->TrackId();
    fMCMichelEnergy = mcp2->E() * 1000;
    fMCMichelEnergyHist->Fill(fMCMichelEnergy);
    fMCMichelStartT = mcp2->T() / 1000.;
    TVector3 fMCMichelStartVect = mcp2->Position(0).Vect();
    fMCMichelStartX = fMCMichelStartVect.X();
    fMCMichelStartY = fMCMichelStartVect.Y();
      fMCMichelStartZ = fMCMichelStartVect.Z();
    if(fMCMichelNPoints < 10) {
     *fMCMichelVect = mcp2->Position(0).Vect() - mcp2->Position(i).Vect();
    } else {*fMCMichelVect = mcp2->Position(0).Vect() - mcp2->Position(j).Vect();
    } 
    rotateVector(fMCMichelVect);
    fMCMichelTheta = (fMCMichelVect->Theta()) * 180 / M_PI - 90;
    rotateVector(fMCMichelVect);
    fMCMichelPhi = (fMCMichelVect->Theta()) * 180 / M_PI - 90;
    rotateVector(fMCMichelVect);
    fMCMichelLength = (mcp2->Position(0).Vect() - mcp2->Position(i).Vect()).Mag();
    fMCMichelLengthHist->Fill(fMCMichelLength);
    fMCMichelRelTheta = fMCMichelTheta - fMCMuonTheta;
    fMCMichelRelPhi = fMCMichelPhi - fMCMuonPhi;
    fMCMichelRelThetaHist->Fill(fMCMichelRelTheta);
    fMCMichelRelPhiHist->Fill(fMCMichelRelPhi);
    fMCMichelEnergyFrac = fMCMichelEnergy / fMCMuonEndEnergy;
  }

  // Find if there is a coincident software trigger
  // Accessing MCParticles
  art::Handle< std::vector<sbnd::trigger::pmtSoftwareTrigger> > trigHandle;
  std::vector< art::Ptr<sbnd::trigger::pmtSoftwareTrigger> > trigVect;
  if(e.getByLabel(fSoftwareTriggerLabel, trigHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(trigVect, trigHandle);

  fMCMuonTrigger = false;
  fMCMuonTriggerID = -1;
  for(unsigned i = 0; i < trigVect.size(); i++) {
    auto trig  = trigVect[i];
    if((abs((double)trig->trig_ts / 1000. - 1510. - mu_trig_time) < 1.)) {
      fMCMuonTrigger = true;
      fMCMuonTriggerID = i;
    }
  }
  // Accessing MCParticles
  art::Handle< std::vector<std::pair<int, short> > > mamophtHandle;
  std::vector< art::Ptr<std::pair<int, short>> > mamophtVect;
  if(e.getByLabel(fMaMOpHTLabel, mamophtHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(mamophtVect, mamophtHandle);

  fMCMuonTrigger = false;
  for(unsigned i = 0; i < mamophtVect.size(); i+=2) {
    auto trig  = mamophtVect[i];
    auto mich_trig = mamophtVect[i+1];
    std::cout << ((float)trig->first / 500) - 1510. << "\n";
    if((abs(((float)trig->first / 500.) - 1510. -  mu_trig_time) < 1.)) {
      fMCMuonMamophtTime = ((float)trig->first / 500.) - 1510.;
      fMCMuonMamophtADC = trig->second;
      fMCMichelMamophtTime = ((float)mich_trig->first / 500.) - 1510.;
      fMCMichelMamophtADC = mich_trig->second;
    }
  }
}

void sbnd::AnalyseMichels::FindRecoMichelShower(
  const art::Event& e)
{
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  art::FindManyP<recob::Shower> pfpShowerAssoc(pfpVec, e, fShowerLabel);

  for(auto& pfp : pfpVec) {
    fPFPIDs.push_back(pfp->Self());
    if(pfp->Self() != unsigned(fRecoMichelID)) continue;
    std::vector< art::Ptr<recob::Shower> > pfpShowers = pfpShowerAssoc.at(pfp.key());
    fNShowers = pfpShowers.size();
    if(pfpShowers.size() == 0) continue;
    auto& shw = pfpShowers[0];
    TVector3 shwstart = shw->ShowerStart();
    fRecoMichelStartX = shwstart.X();
    fRecoMichelStartY = shwstart.Y();
    fRecoMichelStartZ = shwstart.Z();
    fRecoMichelLength = shw->Length();
    fShowerBestPlane = shw->best_plane();
    std::vector<float> recomichelstartvect = {1.0, 1.0, 1.0};
//    recomichelstartvect[0] = shw->Direction().X();
//    recomichelstartvect[1] = shw->Direction().Y();
//    recomichelstartvect[2] = shw->Direction().Z();
    fRecoMichelStartVect->SetXYZ(recomichelstartvect[0], recomichelstartvect[1], recomichelstartvect[2]);
//    fRecoMichelStartVect.Unit();
//    fRecoMichelRelAngle = fRecoMuonEndVect->Angle(&fRecoMichelStartVect) * 180.0 / M_PI;
    int ShowerBestPlane = fShowerBestPlane;
    if(sUseWPlaneOnly==1) ShowerBestPlane = 2;
    fRecoMichelEnergy = shw->Energy().at(ShowerBestPlane);
    TVector3 showerdir = shw->Direction();
    rotateVector(&showerdir);
    fRecoMichelTheta = showerdir.Theta() * 180 / M_PI + 90.0;
    rotateVector(&showerdir);
    fRecoMichelPhi = showerdir.Theta() * 180 / M_PI + 90.0;
    fRecoMichelEnergySigma = shw->EnergyErr().at(ShowerBestPlane);
    fRecoMicheldEdX = shw->dEdx();
    fRecoMichelEnergyVect[0] = shw->Energy().at(0);
    fRecoMichelEnergyVect[1] = shw->Energy().at(1);
    fRecoMichelEnergyVect[2] = shw->Energy().at(2);
    if(fRecoMichelEnergyVect[0] > 0) fRecoMichelEnergyDiffVect[0] = fRecoMichelEnergyVect[0] - fMCMichelEnergy;
    if(fRecoMichelEnergyVect[1] > 0) fRecoMichelEnergyDiffVect[1] = fRecoMichelEnergyVect[1] - fMCMichelEnergy;
    if(fRecoMichelEnergyVect[2] > 0) fRecoMichelEnergyDiffVect[2] = fRecoMichelEnergyVect[2] - fMCMichelEnergy;
  }
}

void sbnd::AnalyseMichels::FindRecoMichelTrack(
  const art::Event& e)
{
  art::Handle<std::vector<recob::Track>> trackHandle;
  std::vector<art::Ptr<recob::Track>> trackVec;
  if (e.getByLabel(fTrackLabel, trackHandle))
    art::fill_ptr_vector(trackVec, trackHandle);

//  art::FindManyP<recob::Hit> hitShowerAssoc(showerVec, e, fHitShowerLabel);

  for(auto& trk : trackVec) {
    if(trk->ID() != fRecoMichelID) continue;
    fRecoMichelStartX = trk->Start().X();
    fRecoMichelStartY = trk->Start().Y();
    fRecoMichelStartZ = trk->Start().Z();      
  }
}

void sbnd::AnalyseMichels::FindRecoMichel(
  const art::Event& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVec;
  if (e.getByLabel(fHitLabel, hitHandle))
    art::fill_ptr_vector(hitVec, hitHandle);

  float event_nmichelhits = 0., event_nmuonhits = 0.,
        event_nmichelhits_w = 0;
  for(auto& hit : hitVec) {
    int hit_id = RecoUtils::TrueParticleID(clockData, hit);
    if(hit_id == fMCMichelID) {
      event_nmichelhits++;
      if(hit->View()==2) event_nmichelhits_w++;
    }
    else if(hit_id == fMCMuonG4ID) event_nmuonhits++;
  }

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVec;
  if (e.getByLabel(fClusterLabel, clusterHandle))
    art::fill_ptr_vector(clusterVec, clusterHandle);

  art::FindManyP<recob::Cluster> pfpClusterAssoc(pfpVec, e, fClusterLabel);
  art::FindManyP<recob::Hit> hitClusterAssoc(clusterVec, e, fHitClusterLabel);

  float min_purity = -0., min_muonhits = 0.;
  for(auto& pfp : pfpVec) {
    std::vector< art::Ptr<recob::Cluster> > pfpClusters = pfpClusterAssoc.at(pfp.key());
    fNClusters.push_back(pfpClusters.size());
    if(pfpClusters.empty()) continue;
    std::vector<art::Ptr <recob::Hit> > pfpHits;
    for(const auto& clust : pfpClusters) {     // For each cluster
      std::vector< art::Ptr<recob::Hit> > clusterHits = hitClusterAssoc.at(clust.key());
      pfpHits.insert(pfpHits.end(), clusterHits.begin(), clusterHits.end());
    }
    float pfp_hits = 0., pfp_nmichelhits = 0., 
          pfp_nhits_w = 0., pfp_nmichelhits_w = 0.,
          pfp_nmuonhits=0.;
    for (const art::Ptr<recob::Hit> &hit : pfpHits) {           // Loop over individual hits in each PFP
      pfp_hits++;
//      art::FindManyP<recob::SpacePoint> spacepointHitAssoc(pfpHits, e, fHitSpacePointLabel);
//      std::vector< art::Ptr<recob::SpacePoint> > hitPoints = spacepointHitAssoc.at(hit.key());
//      if(hitPoints.empty()) continue;
      if(hit->View()==2) pfp_nhits_w++;
      int hitid = RecoUtils::TrueParticleID(clockData, hit);
      if(hitid==fMCMichelID) {
        pfp_nmichelhits++;
        if(hit->View()==2) pfp_nmichelhits_w++;
      }
      else if(hitid==fMCMuonG4ID) pfp_nmuonhits++;
    }
    float pfp_pur = pfp_nmichelhits / pfp_hits;
    float pfp_cmp = pfp_nmichelhits / event_nmichelhits;
    if(pfp_pur > min_purity && pfp_pur > 0.5) {
      if((!fUseTrack && abs(pfp->PdgCode()) == 13) ||
         (!fUseShower && abs(pfp->PdgCode()) == 11)) continue;
      fRecoMichelPurity = pfp_pur;
      fRecoMichelCompleteness = pfp_cmp;
      fRecoMichelPurityW = pfp_nmichelhits_w / pfp_nhits_w;
      fRecoMichelCompletenessW = pfp_nmichelhits_w / event_nmichelhits_w;
      fNTrueHitsInPFP.push_back(pfp_nmichelhits);
      fNHitsInPFP.push_back(pfp_hits);
      fRecoMichelID = pfp->Self();
      if(abs(pfp->PdgCode()) == 11) fIsShower = true;
      else if(abs(pfp->PdgCode()) == 13) fIsTrack = true;
    }
    if(pfp_nmuonhits > min_muonhits) {
      min_muonhits = pfp_nmuonhits;
      fNRecoMuonHits = event_nmuonhits;
      fRecoMuonNClusters = pfpClusters.size();
      fNTrueHitsInRecoMuon = pfp_nmichelhits;
      fRecoMuonID = pfp->Self();
    }
  }
}
void sbnd::AnalyseMichels::FillTrigTree(
  const art::Event& e)

{
  // Accessing MCParticles
  art::Handle< std::vector<sbnd::trigger::pmtSoftwareTrigger> > trigHandle;
  std::vector< art::Ptr<sbnd::trigger::pmtSoftwareTrigger> > trigVect;
  if(e.getByLabel(fSoftwareTriggerLabel, trigHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(trigVect, trigHandle);

  art::Handle< std::vector<sbnd::trigger::Coincidence> > coincHandle;
  std::vector< art::Ptr<sbnd::trigger::Coincidence> > coincVect;
  if(e.getByLabel(fCoincidenceLabel, coincHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(coincVect, coincHandle);


  fMCMuonTrigger = false;
  fMCMuonTriggerID = -1;
  for(unsigned i = 0; i < trigVect.size(); i++) {
    // Reset vectors
    fTrigMultVec.clear(); fTrigMultVec.resize(5120);
    fTrigADCVec.clear(); fTrigADCVec.resize(5120);
    fTrigADCRiseVec.clear(); fTrigADCRiseVec.resize(5120);

    auto coinc = coincVect[i];
    fCoincID = coinc->ID;
    fCoincPlanes = coinc->Planes;
    fCoincCRTTime = coinc->CRTTime;
    fCoincPMTTime = coinc->PMTTime;
    auto trig  = trigVect[i];
    fTrigID = i;
    fTrigTime = (float)trig->trig_ts / 1000. - 1510.;
    int sum  = 0;
    for(unsigned i = 0; i < trig->wvfVec.size() ; i++) {
      fTrigMultVec[i] = trig->multVec[i];
      fTrigADCVec[i] = trig->wvfVec[i];
      if(i==0) {fTrigADCRiseVec[0] = 0;}
      else {
       if(trig->wvfVec[i]>trig->wvfVec[i-1]) {
         sum += trig->wvfVec[i] - trig->wvfVec[i-1];
         fTrigADCRiseVec[i] = sum;
       }
       else {
         sum = 0;
         fTrigADCRiseVec[i] = 0;
       }
     }
    }
    fTrigTree->Fill();
  } 
}

void sbnd::AnalyseMichels::FillPulseTree(
  const art::Event& e)
{
  // Accessing MCParticles
  art::Handle< std::vector<sbnd::trigger::pmtSoftwareTrigger> > trigHandle;
  std::vector< art::Ptr<sbnd::trigger::pmtSoftwareTrigger> > trigVect;                                            if(e.getByLabel(fSoftwareTriggerLabel, trigHandle))     // Make sure artHandle is from mo$
    art::fill_ptr_vector(trigVect, trigHandle);

//  std::cout << "N Triggers: " << trigVect.size() << "\n";
  for(unsigned i=0; i<trigVect.size();i++) {
    auto trig = trigVect[i];
    fPulseTriggerID = (int)i;
    fPulseCh.clear(); fPulseChBaseline.clear(); fPulseChBaselineSigma.clear();
    fPulsePeak.clear(); fPulsePE.clear(); fPulseArea.clear();
    fPulseTStart.clear(); fPulseTEnd.clear(); fPulseTPeak.clear();
//    std::cout << "  N PMTs: " << trig->pmtInfoVec.size() << "\n";
    for(auto ch : trig->pmtInfoVec) {
      for(auto& pulse : ch.pulseVec) {
        if(pulse.pe < 25.) continue;
        fPulseCh.push_back(ch.channel);
        fPulseChBaseline.push_back(ch.baseline);
        fPulseChBaselineSigma.push_back(ch.baselineSigma);
        fPulsePeak.push_back(pulse.peak);
        fPulseArea.push_back(pulse.area);
        fPulsePE.push_back(pulse.pe);
        fPulseTStart.push_back(pulse.t_start);
        fPulseTEnd.push_back(pulse.t_end);
        fPulseTPeak.push_back(pulse.t_peak);
      }
    }
    if(fPulsePeak.empty()) continue;
    fPulseTree->Fill();
  }
}

void sbnd::AnalyseMichels::MatchMCRecoLight(const art::Event& e)
{
  constexpr int kNumChannels = 312;

  // Reset photon vectors
  fOpCh.clear();
  for(unsigned ch=0; ch<312; ch++) fOpCh.push_back(ch);
  fRecoMichelPhotons.clear(); fRecoMichelPhotons.resize(kNumChannels);
  fMCMichelPhotons.clear(); fMCMichelPhotons.resize(kNumChannels);
  fRecoMuonPhotons.clear(); fRecoMuonPhotons.resize(kNumChannels);
  fMCMuonPhotons.clear(); fMCMuonPhotons.resize(kNumChannels);
  fRecoMichelNOpHits.clear(); fRecoMichelNOpHits.resize(kNumChannels);
  fRecoMuonNOpHits.clear(); fRecoMuonNOpHits.resize(kNumChannels);


  // Process each OpHit input tag
  for(unsigned i=0; i<fOpHitLabels.size(); i++) {
    art::InputTag inputTag(fOpHitLabels[i]);

    if (auto opHitsHandle = e.getValidHandle<std::vector<recob::OpHit>>(inputTag)) {
      for (const auto& hit : *opHitsHandle) {
//        if(abs(hit.PeakTime() - fMCMuonTime) < 10.) std::cout << fMCMuonTime << "	" << hit.PeakTime() << "	" << hit.PE() << "\n";
        if (std::abs(hit.PeakTime() - fOpHitDelays[i] - fMCMichelStartT) <= 0.1) {
          int opChannel = hit.OpChannel();
          if (opChannel >= 0 && opChannel < kNumChannels) {
            fRecoMichelPhotons[opChannel] += hit.PE();
            fRecoMichelNOpHits[opChannel]++;
          }
        }
        else if (std::abs(hit.PeakTime() - fOpHitDelays[i] - fMCMuonTime) <= 0.1) {
          int opChannel = hit.OpChannel();
          if (opChannel >= 0 && opChannel < kNumChannels) {
            fRecoMuonPhotons[opChannel] += hit.PE();
            fRecoMuonNOpHits[opChannel]++;
          }
        }
      }
    } else {
      std::cerr << "Error: No OpHits found for input tag: " << inputTag << std::endl;
    }
  }

  // Retrieve and process the SimPhotons
  std::vector<art::Handle<std::vector<sim::SimPhotons>>> fPhotonHandles;
  fPhotonHandles = e.getMany<std::vector<sim::SimPhotons>>();
  for(auto& handle : fPhotonHandles) {
    if (handle.isValid()) {
      for (const auto& photons : *handle) {
        int opChannel = photons.OpChannel();
        if (opChannel >= 0 && opChannel < kNumChannels) {
          for (const auto& photon : photons) {
            if (photon.MotherTrackID == fMCMichelID) {
              fMCMichelPhotons[opChannel]++;
            }
            else if (photon.MotherTrackID == fMCMuonG4ID) {
              fMCMuonPhotons[opChannel]++;
            }
          }
        }
      }
    } else {
      std::cerr << "Error: No SimPhotons found for input tag: " << fSimPhotonsLabel[0] << std::endl;
    }
  }

  // Vector to store the channel information
  std::vector<std::pair<int, std::pair<double, int>>> channelData;

  // Collect data for each channel with recorded values
  for (size_t channel = 0; channel < kNumChannels; ++channel) {
    channelData.emplace_back(channel, std::make_pair(fRecoMichelPhotons[channel], fMCMichelPhotons[channel]));
  }

  // Output the recorded information
//  for (const auto& data : channelData) {
//    int channel = data.first;
//    double totalPEs = data.second.first;
//    int onePhotonCount = data.second.second;

//    std::cout << "Optical Channel: " << channel
//              << ", Total PEs: " << totalPEs
//              << ", OnePhoton Count: " << onePhotonCount
//              << std::endl;
//  }
  fOpTree->Fill();
}

DEFINE_ART_MODULE(sbnd::AnalyseMichels)
