install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirectory(recomb_variations)
add_subdirectory(legacy)
add_subdirectory(mu_lifetime_systematics)
