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
// Class:       AnalyseWaveforms
// Plugin Type: Analyser (Unknown Unknown)
// File:        AnalyseWaveforms_module.cc
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
#include "lardataobj/RawData/OpDetWaveform.h"
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
class AnalyseWaveforms;
}



class sbnd::AnalyseWaveforms : public art::EDAnalyzer {
  public:
  explicit AnalyseWaveforms(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyseWaveforms(AnalyseWaveforms const&) = delete;
  AnalyseWaveforms(AnalyseWaveforms&&) = delete;
  AnalyseWaveforms& operator=(AnalyseWaveforms const&) = delete;
  AnalyseWaveforms& operator=(AnalyseWaveforms&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  private:

  void rotateVector(TVector3 *vector);  // Rotates components of 3-vector clockwise e.g. XTZ -> ZXY

  void ResetVars();

  void FillMC(const art::Ptr<simb::MCParticle>& mcp, std::vector<art::Ptr<simb::MCParticle>>& mctruthVect);
  void FindWVFData(const art::Event& e);

  // Create out output tree
  TTree* fTree;

  // Event Variables
  int fRun;
  int fSubrun;
  int fEventID;
  double fEventTriggerTime;
  std::vector<int> fMCPDG;
  std::vector<int> fMCTrackID;

  std::vector<double> fWVFChannel;
  std::vector<double> fWVFBaseline;
  std::vector<double> fWVFSigma;
  std::vector<std::string> fWVFDet;
  std::vector<double> fWVFTime;
  std::vector<short> fWVFADCMax;
  std::vector<float> fWVFADCMaxTime;

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
  float fMCMichelEnergyFrac;            // Fraction of end MC muon energy taken by michel


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

  // Temporary variables
  int countpfps;
  std::string lastchar;

  // Define input labels
  const std::string fMCTruthLabel;
  const std::string fWVFLabel;
  const std::unique_ptr<opdet::PDMapAlg> fPDMapAlgPtr;
  const double fStartTime;
  const double fEndTime;

  // Declare member data here.
  void FillMC(art::Ptr<simb::MCParticle>& mcp, std::vector<art::Ptr<simb::MCParticle>>& mctruthVect);
    // BackTrackerService
//    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
};

sbnd::AnalyseWaveforms::AnalyseWaveforms(fhicl::ParameterSet const& p)
    : EDAnalyzer { p }, fMCMichelVect(nullptr), fMCMuonVect(nullptr)
    //Initialise out input labels by reading the fhicl parameters
  , fMCTruthLabel(p.get<std::string>("MCTruthLabel"))
  , fWVFLabel(p.get<std::string>("WVFLabel"))
  , fPDMapAlgPtr(art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg")))
  , fStartTime(p.get<double>("StartTime", -1.))
  , fEndTime(p.get<double>("EndTime", 10.))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::AnalyseWaveforms::rotateVector(TVector3 *vector) {
  float temp = vector->Y();
  vector->SetY(vector->Z());
  vector->SetZ(temp);
  temp = vector->X();
  vector->SetX(vector->Y());
  fMCMuonVect->SetY(temp);
}



void sbnd::AnalyseWaveforms::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEventID = e.id().event();

// std::cout << __FILE__ << "::" << __func__ << "():[" << __LINE__ << "]\t\n";

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

//  for(auto const &mcp: mctruthVect) {
//    if(abs(mcp->PdgCode()) != 13 || abs(mcp->EndX()) > 200. || abs(mcp->EndY()) > 200. || mcp->EndZ() < 0. || mcp->EndZ() > 500.) continue;
//    ResetVars();
//    FillMC(mcp, mctruthVect);
//    if(fMCMichelID == 0) continue;
    FindWVFData(e);
    fTree->Fill();
//  }
}

void sbnd::AnalyseWaveforms::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  // Get the TFileService to create out output tree for us
  fTree = tfs->make<TTree>("tree", "Output Tree");

  // Add branches to the TTree
  // Event
  fTree->Branch("event.ID", 			&fEventID);
  fTree->Branch("event.TriggerTime",		&fEventTriggerTime);
  fTree->Branch("event.MCPDG", 			&fMCPDG);
  fTree->Branch("event.MCTrackID", 		&fMCTrackID);
  fTree->Branch("event.WVFChannel",		&fWVFChannel);
  fTree->Branch("event.WVFBaseline",            &fWVFBaseline);
  fTree->Branch("event.WVFSigma",               &fWVFSigma);
  fTree->Branch("event.WVFDet",             	&fWVFDet);
  fTree->Branch("event.WVFTime",		&fWVFTime);
  fTree->Branch("event.WVFADCMax",		&fWVFADCMax);
  fTree->Branch("event.WVFADCMaxTime",		&fWVFADCMaxTime);

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

  fTree->Branch("countpfpf",				&countpfps);
  fTree->Branch("lastchar",				&lastchar);
}

void sbnd::AnalyseWaveforms::endJob()
{}

void sbnd::AnalyseWaveforms::ResetVars()
{
  // Reset all of our variables to 0 or empty vectors
  // This ensures things are not kept from the previous event
  // Event Variables
  fEventID = -9999;
  fEventTriggerTime = -9999.;

  fMCPDG.clear();
  fMCTrackID.clear();

  fWVFChannel.clear();
  fWVFBaseline.clear();
  fWVFSigma.clear();
  fWVFDet.clear();
  fWVFADCMax.clear();
  fWVFADCMaxTime.clear();

  // MC Michel
  fMCMichelID = -9999;
  fMCMichelTime = -9999.;
  fMCMichelEnergy = -9999.;
  fMCMichelTheta = -9999.;
  fMCMichelPhi = -9999.;
  fMCMichelRelTheta = -9999.;
  fMCMichelRelPhi = -9999.;
  fMCMichelLength = -9999.;
  fMCMichelStartX = -9999.;
  fMCMichelStartY = -9999.;
  fMCMichelStartZ = -9999.;
  fMCMichelRelAngle = -999.;
  fMCMichelNPoints = -9999;
  fMCMichelEnergyFrac = -9999.;            // Fraction of end MC muon energy taken by michel


  // MC Muon
  fMCMuonPDG = -9999.;
  fMCMuonTime = -9999.;
  fMCMuonG4ID = -9999.;
  fNDeltas = -9999.;
  fNPoints = -9999.;
  fNElectronsWithTrueHits = -9999.;
  fMCMuonEnergy = -9999.;
  fMCMuonTheta = -9999.;
  fMCMuonPhi = -9999.;
  fMCMuonEndX = -9999.;
  fMCMuonEndY = -9999.;
  fMCMuonEndZ = -9999.;
  fMCMuonEndEnergy = -9999.;
  fMCMuonEndPx = -9999.;
  fMCMuonEndPy = -9999.;
  fMCMuonEndPz = -9999.;
  fMCMuonStartX = -9999.;
  fMCMuonStartY = -9999.;
  fMCMuonStartZ = -9999.;
  fMCMuonBendiness = -9999.;
}


void sbnd::AnalyseWaveforms::FillMC(const art::Ptr<simb::MCParticle>& mcp,
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
    fMCMichelRelTheta = fMCMichelTheta - fMCMuonTheta;
    fMCMichelRelPhi = fMCMichelPhi - fMCMuonPhi;
    fMCMichelEnergyFrac = fMCMichelEnergy / fMCMuonEndEnergy;
  }
}

void sbnd::AnalyseWaveforms::FindWVFData(
  const art::Event& e)
{
  art::ServiceHandle<art::TFileService> tfs;

  art::Handle<std::vector<raw::OpDetWaveform>> wvfHandle;
  e.getByLabel(fWVFLabel, wvfHandle);

  std::vector<raw::OpDetWaveform> wvfms;
  wvfms.reserve(wvfHandle->size());
  for(const auto& wv : *wvfHandle) wvfms.push_back(wv);

  int n_coated = 0, n_uncoated = 0;
  for(const auto& wv: wvfms) {
    unsigned ch = wv.ChannelNumber();
    fWVFChannel.push_back(ch);
    if(fPDMapAlgPtr->pdType(ch) == "pmt_coated") n_coated++;
    if(fPDMapAlgPtr->pdType(ch) == "pmt_uncoated") n_uncoated++;
  }

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  double sampling_rate = double(1./clockData.OpticalClock().Frequency());
  std::cout << "EVent Trigger Time: " << clockData.TriggerTime() << std::endl;
  fEventTriggerTime =  clockData.TriggerTime();

  for(const auto& wv: wvfms) {
    unsigned ch = wv.ChannelNumber();
    fWVFChannel.push_back(ch);
    fWVFDet.push_back(fPDMapAlgPtr->pdType(ch));
    fWVFTime.push_back(wv.TimeStamp());
    auto& wvf = wv.Waveform();
    double sum = 0.;
    auto min_adc = std::min_element(wvf.begin(), wvf.end());
    fWVFADCMax.push_back(*min_adc);
    fWVFADCMaxTime.push_back(wv.TimeStamp() + sampling_rate*std::distance(wvf.begin(), min_adc));
    for(unsigned ix = 0 ; ix<wvf.size();ix++) {
      sum += wvf[ix];
    }
    fWVFBaseline.push_back(sum/wvf.size());
  }
}

DEFINE_ART_MODULE(sbnd::AnalyseWaveforms)
