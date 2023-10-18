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

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"   // Find associations as pointers
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "sbncode/OpDet/PDMapAlg.h"
#include "lardataobj/RecoBase/OpHit.h"
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

////////////////////////////////////////////////////////////////////////
// Class:       pmtSoftwareTriggerProducer
// Plugin Type: producer (Unknown Unknown)
// File:        pmtSoftwareTriggerProducer_module.cc
//
// Generated at Thu Feb 17 13:22:51 2022 by Patrick Green using cetskelgen
// from  version .

// Module to implement software trigger metrics to the PMT Trigger simulation
// Input: artdaq fragment output from the pmtArtdaqFragmentProducer.cc module
// Calculates various PMT metrics for every event (that passes the hardware trigger)
// Output: sbnd::trigger::pmtSoftwareTrigger data product 

// More information can be found at:
// https://sbnsoftware.github.io/sbndcode_wiki/SBND_Trigger
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
#include "art_root_io/TFileService.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

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

  void FillMC(const art::Event& e, const art::Ptr<simb::MCParticle>& mcp, std::vector<art::Ptr<simb::MCParticle>>& mctruthVect);
  void FindRecoMichelShower(const art::Event& e);
  void FindRecoMichelTrack(const art::Event& e);
  void FindRecoMichel(const art::Event& e);
  void FindOpticalData(const art::Event& e);
  void FindFragmentData(const art::Event& e);

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
  int fEventID, fEventRun, fEventSub;
  int fNPFParticles;
  int fEventNOpHits;
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

  double fMichelOpHits;
  std::vector<double> fOpHitTimes;
  std::vector<double> fOpHitWidths;
  std::vector<double> fOpHitPEs;
  std::vector<double> fWVFChannel;
  std::vector<double> fWVFBaseline;
  std::vector<double> fWVFSigma;
  std::vector<std::string> fWVFDet;

  // MC Michel
  int fMCMichelID;
  double fMCMichelTime;
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
  float fMCMichelRelAngle;
  int fMCMichelNPoints;
  float fMCMichelEnergyFrac;		// Fraction of end MC muon energy taken by michel
  float fMCMichelTotalADC;

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
  float fRecoMichel2DVertexDist;
  float fRecoMichel3DVertexDist;
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
  double fMCMuonTime;
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
  float fMCMuonEndEnergy;
  float fMCMuonEndPx;
  float fMCMuonEndPy;
  float fMCMuonEndPz;
  float fMCMuonStartX;
  float fMCMuonStartY;
  float fMCMuonStartZ;
  float fMCMuonBendiness;

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

  // PMT Fragments
  std::vector<int> fFragmentID;
  std::vector<float> fFragmentTimeStamp;
  std::vector<int> fWvfmsStartBin;
  std::vector<std::vector<uint16_t>> fWvfmsVec;

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
  const std::string fOpHitLabel;
  const std::string fPMTFragmentLabel;
  const unsigned fWvfmLength;
  const int sUseWPlaneOnly;
  const std::string fSliceLabel;
  const std::unique_ptr<opdet::PDMapAlg> fPDMapAlgPtr;
  const float fRecombinationFactor;
  const bool fUseTrack;
  const bool fUseShower;
  const double fTimeStart;
  const double fTimeEnd;
  const bool fMakeHists;

  // Declare member data here.
  void FillMC(const art::Event& e, art::Ptr<simb::MCParticle>& mcp, std::vector<art::Ptr<simb::MCParticle>>& mctruthVect);
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
  , fOpHitLabel(p.get<std::string>("OpHitLabel"))
  , fPMTFragmentLabel(p.get<std::string>("PMTFragmentLabel", "fragmentProducer"))
  , fWvfmLength(p.get<unsigned>("FragmentLength", 5120))
  , sUseWPlaneOnly(p.get<int>("UseWPlaneOnly"))
  , fSliceLabel(p.get<std::string>("SliceLabel"))
  , fPDMapAlgPtr(art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg")))
  , fRecombinationFactor(p.get<float>("RecombinationFactor"))
  , fUseTrack(p.get<bool>("UseTrack"))
  , fUseShower(p.get<bool>("UseShower"))
  , fTimeStart(p.get<double>("TimeStart", 0.))
  , fTimeEnd(p.get<double>("TimeEnd", 10.))
  , fMakeHists(p.get<bool>("MakeHists", true))
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
  fEventID = e.id().event();
  fEventSub = e.subRun();
  fEventRun = e.run();

// std::cout << __FILE__ << "::" << __func__ << "():[" << __LINE__ << "]\t\n";

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
//std::cout << __FILE__ << "::" << __func__ << "():[" << __LINE__ << "]\t\n";
   
  // Get ID of Michel - always the last electron

  for(auto const &mcp: mctruthVect) {
    if(abs(mcp->PdgCode()) != 13 || abs(mcp->EndX()) > 200. || abs(mcp->EndY()) > 200. || mcp->EndZ() < 0. || mcp->EndZ() > 500.) continue;
    ResetVars();
    FillMC(e, mcp, mctruthVect);
    if(fMCMichelID == 0) continue;
    FindRecoMichel(e);
    if(fIsShower) FindRecoMichelShower(e);
    else if(fIsTrack) FindRecoMichelTrack(e);
    FindOpticalData(e);
    FindFragmentData(e);
    fTree->Fill();
  }
}

void sbnd::AnalyseMichels::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Add branches to the TTree
  // Event
  fTree->Branch("event.Run",                     &fEventRun);
  fTree->Branch("event.Sub",                     &fEventSub);
  fTree->Branch("event.ID",                     &fEventID);
  fTree->Branch("event.ID", 			&fEventID);
  fTree->Branch("event.NPFParticles", 		&fNPFParticles);
  fTree->Branch("event.NOpHits",		&fEventNOpHits);
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
  fTree->Branch("event.OpHitTimes",		&fOpHitTimes);
  fTree->Branch("event.OpHitWidths",		&fOpHitWidths);
  fTree->Branch("event.OpHitPEs",		&fOpHitPEs);
  fTree->Branch("event.MichelPEs",		&fMichelOpHits);
  fTree->Branch("event.WVFChannel",		&fWVFChannel);
  fTree->Branch("event.WVFBaseline",            &fWVFBaseline);
  fTree->Branch("event.WVFSigma",               &fWVFSigma);
  fTree->Branch("event.WVFDet",             	&fWVFDet);

  // MC Michel
  fTree->Branch("mcMichel.ID", 			&fMCMichelID);
  fTree->Branch("mcMichel.Time",		&fMCMichelTime);
  fTree->Branch("mcMichel.Energy", 		&fMCMichelEnergy);
  fTree->Branch("mcMichel.Theta", 		&fMCMichelTheta);
  fTree->Branch("mcMichel.Phi", 		&fMCMichelPhi);
  fTree->Branch("mcMichel.RelTheta", 		&fMCMichelRelTheta);
  fTree->Branch("mcMichel.RelPhi", 		&fMCMichelRelPhi);
  fTree->Branch("mcMichel.Vect", 		&fMCMichelVect);
  fTree->Branch("mcMichel.Length", 		&fMCMichelLength);
  fTree->Branch("mcMichel.StartX",		&fMCMichelStartX);
  fTree->Branch("mcMichel.StartY",		&fMCMichelStartY);
  fTree->Branch("mcMichel.StartZ",		&fMCMichelStartZ);
  fTree->Branch("mcMichel.RelAngle",		&fMCMichelRelAngle);
  fTree->Branch("mcMichel.NPoints",		&fMCMichelNPoints);
  fTree->Branch("mcMichel.EnergyFrac",		&fMCMichelEnergyFrac);
  fTree->Branch("mcMichel.TotalADC",		&fMCMichelTotalADC);

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
  fTree->Branch("recoMichel.2DVertexDist",	&fRecoMichel2DVertexDist);
  fTree->Branch("recoMichel.3DVertexDist",      &fRecoMichel3DVertexDist);
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
  fTree->Branch("mcMuon.Time",			&fMCMuonTime);
  fTree->Branch("mcMuon.NDeltas", 		&fNDeltas);
  fTree->Branch("mcMuon.NPoints", 		&fNPoints);
  fTree->Branch("mcMuon.Energy", 		&fMCMuonEnergy);
  fTree->Branch("mcMuon.Theta", 		&fMCMuonTheta);
  fTree->Branch("mcMuon.Phi", 			&fMCMuonPhi);
  fTree->Branch("mcMuon.Vect", 			&fMCMuonVect);
  fTree->Branch("mcMuon.EndX",			&fMCMuonEndX);
  fTree->Branch("mcMuon.EndY",			&fMCMuonEndY);
  fTree->Branch("mcMuon.EndZ",			&fMCMuonEndZ);
  fTree->Branch("mcMuon.EndEnergy",		&fMCMuonEndEnergy);
  fTree->Branch("mcMuon.EndPx",			&fMCMuonEndPx);
  fTree->Branch("mcMuon.EndPy",			&fMCMuonEndPy);
  fTree->Branch("mcMuon.EndPz",			&fMCMuonEndPz);
  fTree->Branch("mcMuon.StartX",		&fMCMuonStartX);
  fTree->Branch("mcMuon.StartY",		&fMCMuonStartY);
  fTree->Branch("mcMuon.StartZ",		&fMCMuonStartZ);
  fTree->Branch("mcMuon.Bendiness",		&fMCMuonBendiness);

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
  fEventNOpHits = -9999;
  fNDeltas = 0;
  fMCMichelID = 0;
  fMCMichelTime = -9999.;
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
  fOpHitTimes.clear();
  fOpHitWidths.clear();
  fOpHitPEs.clear();
  fMichelOpHits = 0.;
  fWVFChannel.clear();
  fWVFBaseline.clear();
  fWVFSigma.clear();
  fWVFDet.clear();
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
  fMCMichelTotalADC = 0.;
  fRecoMichel2DVertexDist = -9999.;
  fRecoMichel3DVertexDist = -9999.;
  fRecoMichelStartX = 0;
  fRecoMichelStartY = 0;
  fRecoMichelStartZ = 0;
  fMCMuonPDG = -9999;
  fMCMuonTime = -9999.;
  fMCMuonG4ID = -9999;
  fMCMuonEndX = 0;
  fMCMuonEndY = 0;
  fMCMuonEndZ = 0;
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
  fMCMuonStartX = 0;
  fMCMuonStartY = 1000;
  fMCMuonStartZ = 0;
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

  fFragmentID.clear();
  fFragmentTimeStamp.clear();
  fWvfmsStartBin.clear();
  fWvfmsVec.clear();;
}


void sbnd::AnalyseMichels::FillMC(
  const art::Event& e,
  const art::Ptr<simb::MCParticle>& mcp,
  std::vector<art::Ptr<simb::MCParticle>>& mctruthVect)
{
  fMCMuonPDG = mcp->PdgCode();
  fMCMuonG4ID = mcp->TrackId();
  fNPoints = mcp->NumberTrajectoryPoints();
  for(int pos=0;pos<fNPoints;pos++) {
    if((abs(mcp->Position(pos).X()) < 200.) &&
       (abs(mcp->Position(pos).Y()) < 200.) &&
       (mcp->Position(pos).Z() > 0.) &&
       (mcp->Position(pos).Z() < 500.)) {
      fMCMuonTime = mcp->T(pos); 
      fMCMuonStartX = mcp->Position(pos).X();
      fMCMuonStartY = mcp->Position(pos).Y();
      fMCMuonStartZ = mcp->Position(pos).Z();
      break;
    }
  }
  unsigned int i =  fNPoints/2;
  unsigned int j = fNPoints-1;
  unsigned int firstquart = fNPoints/4;
  unsigned int lastquart = 3*fNPoints/4;
  TVector3 tempvect = mcp->Position(0).Vect();
  TVector3 mcmuonend = mcp->EndPosition().Vect();

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
    if(abs(mcp2->PdgCode()) != 11 || mcp2->Mother() != fMCMuonG4ID) continue;
    fMCMichelTime = mcp2->T();
    fMCMichelNPoints = mcp2->NumberTrajectoryPoints();  
    unsigned int i  = mcp2->NumberTrajectoryPoints() -1;
    unsigned int j = mcp2->NumberTrajectoryPoints() / 2;
    fMCMichelID = mcp2->TrackId();
    fMCMichelEnergy = mcp2->E() * 1000;
    fMCMichelEnergyHist->Fill(fMCMichelEnergy);
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
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVec;
  if (e.getByLabel(fHitLabel, hitHandle))
    art::fill_ptr_vector(hitVec, hitHandle);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  for(auto& hit : hitVec) {
    int hit_id = RecoUtils::TrueParticleID(clockData, hit);
    if(hit_id == fMCMichelID &&
       hit->View()==2) fMCMichelTotalADC += hit->SummedADC();
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
    fRecoMichel2DVertexDist = sqrt((fRecoMichelStartY - fMCMichelStartY)*(fRecoMichelStartY - fMCMichelStartY) +
                                   (fRecoMichelStartZ - fMCMichelStartZ)*(fRecoMichelStartZ - fMCMichelStartZ));
    fRecoMichel3DVertexDist = sqrt((fRecoMichel2DVertexDist* fRecoMichel2DVertexDist) +
                                   (fRecoMichelStartX - fMCMichelStartX)*(fRecoMichelStartX - fMCMichelStartX));
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
    fRecoMichel2DVertexDist = sqrt((fRecoMichelStartY - fMCMichelStartY)*(fRecoMichelStartY - fMCMichelStartY) +
                                   (fRecoMichelStartZ - fMCMichelStartZ)*(fRecoMichelStartZ - fMCMichelStartZ));
    fRecoMichel3DVertexDist = sqrt((fRecoMichel2DVertexDist* fRecoMichel2DVertexDist) +
                                   (fRecoMichelStartX - fMCMichelStartX)*(fRecoMichelStartX - fMCMichelStartX));
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

void sbnd::AnalyseMichels::FindOpticalData(
  const art::Event& e)
{
  art::Handle<std::vector<recob::OpHit>> ophitHandle;
  std::vector<art::Ptr<recob::OpHit>> ophitVec;
  if (e.getByLabel(fOpHitLabel, ophitHandle))
    art::fill_ptr_vector(ophitVec, ophitHandle);

  std::cout << "NOpHits: " << ophitVec.size() << std::endl;
  fEventNOpHits = ophitVec.size();

  std::vector<recob::OpHit> ophits;
  for(unsigned i=0; i<ophitVec.size();i++) {
    auto& oph = *(ophitVec)[i];
    if(oph.StartTime() > fTimeEnd ||
       oph.StartTime() < fTimeStart) continue;
    ophits.emplace_back(oph);
    fOpHitTimes.push_back(oph.StartTime());
    fOpHitWidths.push_back(oph.Width());
    fOpHitPEs.push_back(oph.PE());
    if((oph.StartTime()*1000. > fMCMichelTime) && (oph.StartTime()*1000 < fMCMichelTime + 200)) fMichelOpHits += oph.PE();
  }
  std::sort(ophits.begin(), ophits.end(),
    [this] (const recob::OpHit oph1, recob::OpHit oph2)
    { return oph1.StartTime() < oph2.StartTime(); });

  double start_time = ophits[0].StartTime();
  double end_time = ophits[ophits.size()-1].StartTime();
  std::cout << "StartTime: " << start_time << std::endl;
  std::cout << "EndTime: " << end_time << std::endl;
  auto clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  double tickPeriod = clockData.OpticalClock().TickPeriod();
  unsigned bins = unsigned(1./tickPeriod * (fTimeEnd - fTimeStart));

  if(fMakeHists) {
    std::stringstream histname;
    histname.str(std::string());
    histname << "ophits_hist_" << fEventID;
  
    art::ServiceHandle<art::TFileService> tfs;
    TH1D* ophits_hist = tfs->make<TH1D>(histname.str().c_str(), "ophit time profile", bins, start_time, end_time);
    for(auto& oph : ophits) ophits_hist->Fill(oph.StartTime(), oph.PE());
  }
}

void sbnd::AnalyseMichels::FindFragmentData(
  const art::Event& e)
{
  art::Handle<std::vector<artdaq::Fragment>> fragHandle;
  std::vector<art::Ptr<artdaq::Fragment>> fragVec;
  if (e.getByLabel(fPMTFragmentLabel, fragHandle))
    art::fill_ptr_vector(fragVec, fragHandle);

  std::vector<std::vector<uint16_t>> fWvfmsVec;
  fWvfmsVec.resize(15*8); // 15 pmt channels per fragment, 8 fragments per trigger
  for(unsigned fragmentIdx = 0; fragmentIdx<fragVec.size(); fragmentIdx++) {
    const auto& frag = fragVec[fragmentIdx];
    int fragId = static_cast<int>(frag->fragmentID()); 
    fFragmentID.push_back(fragId);
    sbndaq::CAENV1730Fragment bb(*frag);
    auto const* md = bb.Metadata();
    // access timestamp
    uint32_t timestamp = md->timeStampNSec;
    fFragmentTimeStamp.push_back(timestamp);
    double fFragTimeStamp = timestamp - fTimeStart;
    int startbin = (fFragTimeStamp >= 1000) ? 0 : int(500 - (fFragTimeStamp)/2); // units of bins
    fWvfmsStartBin.push_back(startbin);

    // access waveforms in fragment and save
    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag->dataBeginBytes() 
                 + sizeof(sbndaq::CAENV1730EventHeader));
    const uint16_t* value_ptr =  data_begin;
    uint16_t value = 0;

    // channel offset
    size_t nChannels = 15; // 15 pmts per fragment
    size_t ch_offset = 0;
    for (size_t i_ch = 0; i_ch < nChannels; ++i_ch) {
      fWvfmsVec[i_ch + nChannels*fragId].resize(fWvfmLength);
      ch_offset = (size_t)(i_ch * fWvfmLength);
      //--loop over waveform samples
      for(size_t i_t = 0; i_t < fWvfmLength; ++i_t){ 
        value_ptr = data_begin + ch_offset + i_t; // pointer arithmetic
        value = *(value_ptr);
        fWvfmsVec[i_ch + nChannels*fragId][i_t] = value;
      } //--end loop samples
    } //--end loop channels
  } // Waveform handle loop
}

DEFINE_ART_MODULE(sbnd::AnalyseMichels)
