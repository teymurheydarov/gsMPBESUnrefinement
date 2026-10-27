######################################################################
## CMakeLists.txt --- gsMPBESUnrefinement
## This file is part of the G+Smo library.
##
## Unrefinement of planar multipatch THB-spline parameterisations in
## the MPBES basis, with the fitting stages (least squares, linear and
## non-linear optimisation) formulated with CasADi
######################################################################

## gsMPBESUnrefinement module
project(gsMPBESUnrefinement)

# Apply G+Smo config
include(gsConfig)

## The MPBES basis (src/gsMPBES.*) is part of this module; a separate
## module defining the same class must not be enabled at the same time
if(TARGET gsMPBES)
  message(FATAL_ERROR "${PROJECT_NAME} contains the MPBES basis (gsMPBES); "
    "remove gsMPBES from GISMO_OPTIONAL")
endif()

## CasADi (https://web.casadi.org): the optimisation stages use its
## matrix expressions and automatic differentiation; the binary
## distributions (pip wheel, conda package, release archive) ship a
## casadi-config.cmake, point casadi_DIR to the directory holding it.
find_package(casadi REQUIRED)
get_target_property(${PROJECT_NAME}_CASADI_INCLUDE casadi::casadi INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(${PROJECT_NAME}_CASADI_CONFIGS casadi::casadi IMPORTED_CONFIGURATIONS)
set(${PROJECT_NAME}_CASADI_LIBRARY "")
foreach(cfg ${${PROJECT_NAME}_CASADI_CONFIGS})
  # on Windows the linker needs the import library, not the DLL
  if(WIN32 AND NOT ${PROJECT_NAME}_CASADI_LIBRARY)
    get_target_property(${PROJECT_NAME}_CASADI_LIBRARY casadi::casadi IMPORTED_IMPLIB_${cfg})
  endif()
  if(NOT ${PROJECT_NAME}_CASADI_LIBRARY)
    get_target_property(${PROJECT_NAME}_CASADI_LIBRARY casadi::casadi IMPORTED_LOCATION_${cfg})
  endif()
endforeach()
if(NOT ${PROJECT_NAME}_CASADI_LIBRARY)
  get_target_property(${PROJECT_NAME}_CASADI_LIBRARY casadi::casadi IMPORTED_LOCATION)
endif()
if(NOT ${PROJECT_NAME}_CASADI_LIBRARY)
  message(FATAL_ERROR "${PROJECT_NAME}: the location of the CasADi library was not found in the casadi::casadi target")
endif()
message(STATUS "${PROJECT_NAME}: CasADi ${casadi_VERSION} (${${PROJECT_NAME}_CASADI_LIBRARY})")

# the library is linked into gismo, the headers are visible to the module sources
set(gismo_LINKER ${gismo_LINKER} ${${PROJECT_NAME}_CASADI_LIBRARY}
  CACHE INTERNAL "gismo extra linker objects")
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${${PROJECT_NAME}_CASADI_INCLUDE}
  CACHE INTERNAL "gismo include directories")
include_directories(${${PROJECT_NAME}_CASADI_INCLUDE})

## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules" )

target_include_directories(${PROJECT_NAME} PRIVATE ${${PROJECT_NAME}_CASADI_INCLUDE})
# CasADi < 3.8 defines operator<<(ostream, std::pair) in namespace std, as gsSortedVector.h
# does; G+Smo skips its own definition when TORCH_VERSION is defined
if(casadi_VERSION VERSION_LESS 3.8)
  target_compile_definitions(${PROJECT_NAME} PRIVATE TORCH_VERSION=0)
endif()

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

# Symlink include dir (headers are in /src)
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" )

# add example files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# add unittests
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")
