# ------------------------------------------------------------------------------
# OnlineAnalysisSlicer
# ------------------------------------------------------------------------------
add_executable(OnlineAnalysisSlicer
    OnlineAnalysisSlicer.cxx
)

target_include_directories(OnlineAnalysisSlicer
    PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/unpacker
    ${CMAKE_SOURCE_DIR}/sub/NestdaqAnalyzer/NestdaqTFSlicer
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(OnlineAnalysisSlicer
    PRIVATE
    NestdaqTFSlicer
    FairMQ::FairMQ
#    FairMQStateMachine
    FairLogger::FairLogger
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

# ------------------------------------------------------------------------------
# OnlineAnalysisTFBuilder
# ------------------------------------------------------------------------------
# NestdaqUnpacker does not seem to export a library.
# We compile its sources directly here.
set(NESTDAQ_UNPACKER_SRC_DIR ${CMAKE_SOURCE_DIR}/unpacker)

add_executable(OnlineAnalysisTFBuilder
    OnlineAnalysisTFBuilder.cxx
    ${NESTDAQ_UNPACKER_SRC_DIR}/NodeUnpacker.cc
    ${NESTDAQ_UNPACKER_SRC_DIR}/LeafProcessor.cc
    ${NESTDAQ_UNPACKER_SRC_DIR}/LeafDecoder.cc
    ${NESTDAQ_UNPACKER_SRC_DIR}/NodeHeaderDecoder.cc
)

target_include_directories(OnlineAnalysisTFBuilder
    PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/unpacker
    ${CMAKE_SOURCE_DIR}/sub/NestdaqAnalyzer/NestdaqTFSlicer
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(OnlineAnalysisTFBuilder
    PRIVATE
    NestdaqTFSlicer
    FairMQ::FairMQ
#    FairMQStateMachine
    FairLogger::FairLogger
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

install(TARGETS OnlineAnalysisSlicer OnlineAnalysisTFBuilder RUNTIME DESTINATION bin)
