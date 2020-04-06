#set(CMAKE_PREFIX_PATH /Users/nevay/physics/reps/bdsim-develop-install/lib/cmake/bdsim)
set(CMAKE_PREFIX_PATH /Users/nevay/physics/reps/bdsim-tc-sixtrack-link-install/lib/cmake/bdsim)
find_package(BDSIM REQUIRED)
include_directories(${BDSIM_INCLUDE_DIR}/bdsim)
set(G4COLLIMATION_SOURCES BDSSixtrackInterface.cc)
add_library(g4collimation STATIC ${G4COLLIMATION_SOURCES})
target_link_libraries(g4collimation PRIVATE gmadShared bdsim)

