cmake_minimum_required(VERSION 3.1)

project(Position-Based-Dynamics VERSION 0.0.1)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(SofaMisc REQUIRED)
find_package(SofaOpenglVisual REQUIRED)

FILE(GLOB_RECURSE SOURCES_FILES ${CMAKE_SOURCE_DIR}/*.cpp)
FILE(GLOB_RECURSE HEADER_FILES ${CMAKE_SOURCE_DIR}/*.hpp)
FILE(GLOB_RECURSE INLINE_FILES ${CMAKE_SOURCE_DIR}/*.inl)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${INLINE_FILES} ${SOURCES_FILES})
target_link_libraries(${PROJECT_NAME} SofaCore SofaSimulationCore SofaSimulationCommon SofaBaseMechanics SofaTopologyMapping SofaRigid SofaOpenglVisual)
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_MYPLUGIN")
