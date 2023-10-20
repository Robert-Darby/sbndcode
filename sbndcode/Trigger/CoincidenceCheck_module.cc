////////////////////////////////////////////////////////////////////////
// Class:       CoincidenceCheck
// Plugin Type: filter (Unknown Unknown)
// File:        CoincidenceCheck_module.cc
//
// Generated at Thu Oct 19 15:06:07 2023 by Robert Darby using cetskelgen
// from  version .
//
// Thise module checks for coincidence between CRT/PMT trigger times
// within a given window. Can be used for CRT and PDS or CRT only.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"


#include <memory>

class CoincidenceCheck;


class CoincidenceCheck : public art::EDFilter {
public:
  explicit CoincidenceCheck(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CoincidenceCheck(CoincidenceCheck const&) = delete;
  CoincidenceCheck(CoincidenceCheck&&) = delete;
  CoincidenceCheck& operator=(CoincidenceCheck const&) = delete;
  CoincidenceCheck& operator=(CoincidenceCheck&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::map<size_t, std::vector<float>> crtMapMaker(
    artdaq::Fragment& frag,
    std::map<size_t, std::vector<float>>& crtTimeMap);
  // Declare member data here.

  // FHICL Params
  const float fCoincidenceWindow;
  const std::vector<size_t> fCRTPlanes;
  const std::string fCRTLogic;
  const bool fCRTOnly;
  const bool fVerbose;
};


CoincidenceCheck::CoincidenceCheck(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  , fCoincidenceWindow(p.get<float>("CoincidenceWindow", 0.24))
  , fCRTPlanes(p.get<std::vector<size_t>>("CRTPlanes"))
  , fCRTLogic(p.get<std::string>("CRTLogic", "AND"))
  , fCRTOnly(p.get<bool>("CRTOnly", false))
  , fVerbose(p.get<bool>("Verbose", false))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool CoincidenceCheck::filter(art::Event& e)
{
  std::map<size_t, std::vector<float>> crtTimeMap;
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  // loop over fragment handles
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0) continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      // container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() == sbndaq::detail::FragmentType::BERNCRTV2){
          if (fVerbose)     std::cout << "    Found " << contf.block_count() << " CRT Fragments in container " << std::endl;
          for (size_t ii = 0; ii < contf.block_count(); ++ii) crtMapMaker(*contf[ii].get(), crtTimeMap);
        }
      } // Fragment loop
    } // if CRT Fragment
    else {
      // normal fragment
      for (auto frag : *handle){
        if (frag.type()==sbndaq::detail::FragmentType::BERNCRTV2) {
          crtMapMaker(frag, crtTimeMap);
        }
      }
    } // Fragment loop
  } // Fragment handle loop 
  return false;
}

void CoincidenceCheck::beginJob()
{
  // Implementation of optional member function here.
}

void CoincidenceCheck::endJob()
{
  // Implementation of optional member function here.
}

std::map<size_t, std::vector<float>> CoincidenceCheck::crtMapMaker(
  artdaq::Fragment& frag,
  std::map<size_t, std::vector<float>>& crtTimeMap)
{
  // use  fragment ID to get plane information
  sbndaq::BernCRTFragmentV2 bern_fragment(frag);
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
  auto thisone = frag.fragmentID();  uint plane = (thisone & 0x0700) >> 8;
  if (plane>7) {std::cout << "bad plane value " << plane << std::endl; plane=0;}
  // Check if this fragment is considered in the coincidence check
  if(fCRTPlanes[plane-1] == 2) return crtTimeMap;
  for(unsigned int iHit = 0; iHit < md->hits_in_fragment(); iHit++) {
    sbndaq::BernCRTHitV2 const* bevt = bern_fragment.eventdata(iHit);
    // require that this is data and not clock reset (0xC), and that the ts1 time is valid (0x2)
    auto thisflag = bevt->flags;
    if (thisflag & 0x2 && !(thisflag & 0xC) ) {
      // check ts1 for beam window
      auto thistime=bevt->ts1;
      crtTimeMap[thisone-1].push_back(thistime);
    }
  }
  return crtTimeMap;
} // CoincidenceCheck::crtMapMaker

DEFINE_ART_MODULE(CoincidenceCheck)
