function(add_particle_tracker_test target source)
  add_executable(${target} ${source})
  target_link_libraries(${target}
    PRIVATE particle_transport_framework::particle_tracker)
  add_test(NAME ${target} COMMAND ${target})
  set_tests_properties(${target} PROPERTIES TIMEOUT 30)
endfunction()

add_particle_tracker_test(test_geography test_geography.cpp)
add_particle_tracker_test(test_datetime test_datetime.cpp)
add_particle_tracker_test(test_fixed_point_record test_fixed_point_record.cpp)
add_particle_tracker_test(test_fixed_point_gzip test_fixed_point_gzip.cpp)
add_particle_tracker_test(test_pkd2_roundtrip test_pkd2_roundtrip.cpp)
add_particle_tracker_test(test_interpolation test_interpolation.cpp)
