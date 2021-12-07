# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

#sets the project name
PROJECT(myproject)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(INCLUDES
  AddRooDataSet.h
  MLAsymmKing.h  
  MLException.h
  )

set(SOURCES
  src/AddRooDataSet.cxx
  src/MLAsymmKing.cxx
  )

#################
# Configuration #
#################


## ============================================================================
##  ROOT
## ============================================================================

# Look for ROOT
message("Looking for ROOT > 6.0...")
find_package(ROOT 6 REQUIRED COMPONENTS
  RooFitCore
  TreePlayer
  MathCore
  Gui
  pthread
  MathMore
  )
list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})
message(STATUS "Found ROOT!")
include(${ROOT_USE_FILE})
message(STATUS "ROOT include directories: ${ROOT_INCLUDE_DIRS}")
include_directories("${ROOT_INCLUDE_DIRS}")
message(STATUS "Using ROOT libraries ${ROOT_LIBRARIES}")


## ============================================================================
##  EIGEN
## ============================================================================

include_directories(/work/LAS/amandajw-lab/veritas/eigen-3.3.9)


## ============================================================================
##  PKG-CONFIG
## ============================================================================

##Find PkgConfig which is required for VEGAS to work
find_package(PkgConfig REQUIRED)
if(PKG_CONFIG_FOUND)
  message(STATUS "PkgConfig found.")
else(NOT PKG_CONFIG_FOUND)
  message(FATAL_ERROR "PkgConfig not found")
endif(PKG_CONFIG_FOUND)


include_directories(/work/LAS/amandajw-lab/users/achrmy/macros/RooAddDataSetPlotCombineNoise/include)


add_definitions(
  -D__STDC_LIMIT_MACROS
  -D__STDC_CONSTANT_MACROS
  -DNOHDF5
  -DENABLE_ML_ROOTFILE_IO
  -Wno-deprecated
  )


add_executable(AddRooDataSet src/main.cpp ${SOURCES} ${ROOT_INCLUDE_DIRS})
target_link_libraries(AddRooDataSet ${ROOT_LIBRARIES})
 # ${SOURCES}
 # ${INCLUDES}
 # ${ROOT_LIBRARIES}
 # )

#install(TARGETS MLRooDataStore DESTINATION ${CMAKE_INSTALL_BINDIR})
#install(MLRooDataStore ${CMAKE_INSTALL_BINDIR})

#requests that the executable be built
#add_executable(MLRooDataStore ${SOURCES})
#cmake ./ -DVEGAS_BASE_DIR=$VEGAS
#cmake --build .
#rm -r CMakeCache.txt CMakeFiles cmake_install.cmake example.txt
