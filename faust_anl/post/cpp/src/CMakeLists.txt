#add_executable(real_w_mixed real_w_mixed.cpp)
#target_link_libraries(real_w_mixed PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

add_executable(real_w_mixed real_w_mixed_hist.cpp formatting.cpp)
target_link_libraries(real_w_mixed PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

add_executable(rescale_similar_e_rel rescale_similar_e_rel.cpp)
target_link_libraries(rescale_similar_e_rel PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)


set(SRC_LIST 3h_correlation_correction_to_3a.cpp formatting.cpp)
add_executable(3h_correlation_correction_to_3a ${SRC_LIST})
target_link_libraries(3h_correlation_correction_to_3a PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

set(SRC_LIST 3a_bkgnd.cpp formatting.cpp)
add_executable(3a_bkgnd ${SRC_LIST})
target_link_libraries(3a_bkgnd PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

set(SRC_LIST 3a_kde.cpp formatting.cpp)
add_executable(3a_kde ${SRC_LIST})
target_link_libraries(3a_kde PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

set(SRC_LIST source_momentum_v_e_rel.cpp formatting.cpp)
add_executable(source_momentum_v_e_rel ${SRC_LIST})
target_link_libraries(source_momentum_v_e_rel PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

set(SRC_LIST pir_2023.cpp formatting.cpp)
add_executable(pir_2023 ${SRC_LIST})
target_link_libraries(pir_2023 PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)

set(SRC_LIST vector_diagram.cpp formatting.cpp)
add_executable(vector_diagram ${SRC_LIST})
target_link_libraries(vector_diagram PRIVATE ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree)
