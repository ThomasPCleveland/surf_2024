cmake_minimum_required(VERSION 3.29)
project(KID_DMD)

# ----------------------------------------------------------------------------
# Find G4CMP package -- NOTE: Assumes user has run g4cmp_env.*sh script!
#
set(CMAKE_PREFIX_PATH $ENV{G4CMPINSTALL}/../..)
find_package(G4CMP REQUIRED)
include(${G4CMP_USE_FILE})

# ----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
include($ENV{G4CMPINSTALL}/FindGeant4.cmake)

# ----------------------------------------------------------------------------
# Setup include directories and compile definitions
# NOTE: Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${G4CMP_INCLUDE_DIRS})
include(${G4CMP_USE_FILE})
include(${Geant4_USE_FILE})

# ----------------------------------------------------------------------------
# RPATH stuff
#
# Automatically add paths of linked libraries to RPATH *unless* they are in the
# build directory
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set the RPATH of the installed objects *unless* they are installing to a
# system directory, in which case they shouldn't need an explicitly set RPATH.
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/lib" is_SYSTEM_DIR)

if("${is_SYSTEM_DIR}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# ----------------------------------------------------------------------------
# Build stuff
#
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Action_Initialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Config_Manager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Config_Messenger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Detector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Electrode.cc
    # ${CMAKE_CURRENT_SOURCE_DIR}/src/Event_Action.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kaplan_QP.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Phonon.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Custom_Physics.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Physics_List.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Primary_Action.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Quasiparticle.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Stepping_Action.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQ_Sensitivity.cc
)

if(USE_GEANT4_STATIC_LIBS)
    add_library(KID_DMD_lib STATIC ${SOURCES})
else()
    add_library(KID_DMD_lib SHARED ${SOURCES})
endif()

set_target_properties(KID_DMD_lib PROPERTIES OUTPUT_NAME g4cmpKID_DMD)

message("G4CMP Libraries: ")
message(${G4CMP_LIBRARIES})

target_link_libraries(KID_DMD_lib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(KID_DMD main.cc)
target_link_libraries(KID_DMD KID_DMD_lib)

install(TARGETS KID_DMD_lib DESTINATION lib)
install(TARGETS KID_DMD DESTINATION bin)
