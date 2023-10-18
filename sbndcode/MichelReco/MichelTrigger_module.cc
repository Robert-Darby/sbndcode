#include <iostream>
#include <algorithm>

#include "TGeoManager.h"
#include "TVector3.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
#include <string>

namespace filt{

  class MichelTriggerFilter : public art::EDFilter {
    public:
      explicit MichelTriggerFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:

   std::string fPMTTriggerLabel;
   short fMultiplicity;
  };


  MichelTriggerFilter::MichelTriggerFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
  }

  void MichelTriggerFilter::reconfigure(fhicl::ParameterSet const& pset){
    fPMTTriggerLabel = pset.get<std::string>("PMTTriggerLabel", "pmttriggerproducer");
    fMultiplicity = pset.get<short>("Multiplicity", 5);
  }


  bool MichelTriggerFilter::filter(art::Event & e)
  {
    art::Handle< std::vector< sbnd::comm::pmtTrigger > > triggerHandle;
    e.getByLabel(fPMTTriggerLabel, triggerHandle);
    std::cout << "N triggers: " << triggerHandle->size() << std::endl;
    if(!triggerHandle.isValid()) return false;

    for(auto const& trigger : (*triggerHandle)) {
      for (size_t idx = 0; idx < trigger.numPassed.size(); idx++) {
        if (trigger.numPassed[idx] >= fMultiplicity) {
          return true;
        }
      } // trigger handle loop
    }
    return false;
  }

  void MichelTriggerFilter::beginJob() {
    std::cout << "Running MichelFilter" << std::endl;
  }

  DEFINE_ART_MODULE(MichelTriggerFilter)

}

