
set(      MODULE_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
			                larreco::Calorimetry
                           larreco::RecoAlg
                           lardata::RecoObjects
                           larpandora::LArPandoraInterface
                           sbndcode_CRTUtils
                           sbndcode_CRT
                           sbnobj::Common_CRT
                           nusimdata::SimulationBase
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support
                           art_root_io::TFileService_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
        	            artdaq_core::artdaq-core_Data
	                    artdaq_core::artdaq-core_Utilities
                	    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	                    sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        	            sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
 
                           messagefacility::MF_MessageLogger
                           ROOT::Core
                           ROOT::Tree
                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           sbndcode_RecoUtils
                           sbndcode_OpDetSim)

# cet_build_plugin(AnalyzeMichels art::Module SOURCE AnalyzeMichels_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(AnalyseMichels art::Module SOURCE AnalyseMichels_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(MaMOpHT art::Module SOURCE MaMOpHT_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(LateLightAnalyser art::Module SOURCE LateLightAnalyser_module.cc LIBRARIES ${MODULE_LIBRARIES} )
# install_headers()
install_fhicl()
install_source()

