# ------------------------------------------------------------------------------
# NestdaqTFSlicer Library
# ------------------------------------------------------------------------------
add_library(NestdaqTFSlicer SHARED
    TFSlicerUnpacker.cxx
)

target_compile_features(NestdaqTFSlicer PUBLIC cxx_std_17)

target_include_directories(NestdaqTFSlicer
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/unpacker
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(NestdaqTFSlicer
    FairMQ::FairMQ
    FairLogger::FairLogger
)

# KTimer is header-only (included in source), so we don't compile it here to avoid duplication.
# It is available to other targets via include_directories in top CMake.
