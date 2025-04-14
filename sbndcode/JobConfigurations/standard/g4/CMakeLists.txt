install_fhicl()

FILE(GLOB fcl_files *.fcl)
install_source(EXTRAS  ${fcl_files} )

add_subdirecttry(mu_lifetime_systematics)
add_subdirectory(recomb_variations)
add_subdirectory(legacy)
