cmake_minimum_required(VERSION 3.0.2)
project(ik_constraint2_body_contact_sample)

add_compile_options(-std=c++14)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS
  choreonoid_viewer
  ik_constraint2_body_contact
  prioritized_inverse_kinematics_solver2
  roslib
  )

find_package(choreonoid REQUIRED)

catkin_package(
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${CHOREONOID_INCLUDE_DIRS}
)

add_library(CnoidIKConstraintBodyContactSamplePlugin SHARED
  src/CnoidIKConstraintBodyContactSamplePlugin.cpp
  src/sample0.cpp
  src/sample1.cpp
  )

target_link_libraries(CnoidIKConstraintBodyContactSamplePlugin
  ${catkin_LIBRARIES}
  ${CHOREONOID_BODY_PLUGIN_LIBRARIES}
  )
set_target_properties(CnoidIKConstraintBodyContactSamplePlugin PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CHOREONOID_PLUGIN_DIR})
