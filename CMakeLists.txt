cmake_minimum_required(VERSION 3.16)
add_definitions(-w)

project(quad)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS Off) 
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build/bin)

add_compile_definitions(EIGEN_STACK_ALLOCATION_LIMIT=0)

set(casadi_INCLUDE_DIRS /usr/local/include/casadi)
set(mujoco_INCLUDE_DIRS /usr/local/include/mujoco)
set(eigen_INCLUDE_DIRS /usr/local/include/eigen3)
set(pinocchio_INCLUDE_DIRS /usr/local/include/pinocchio)

find_package(casadi REQUIRED)
find_package(mujoco REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(pinocchio REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Boost REQUIRED)

include_directories(
    ./
    ${casadi_INCLUDE_DIRS}
    ${mujoco_INCLUDE_DIRS}
    ${pinocchio_INCLUDE_DIRS}
    ${eigen_INCLUDE_DIRS}
    application
    simulate
    estimate
    utils
    control
    control/wbc 
    control/wbc/TaskContact 
    control/nmpc/acados_cpp_demo
)

include_directories($ENV{ACADOS_SOURCE_DIR}/include)
include_directories($ENV{ACADOS_SOURCE_DIR}/include/blasfeo/include)
include_directories($ENV{ACADOS_SOURCE_DIR}/include/hpipm/include)
include_directories($ENV{ACADOS_SOURCE_DIR}/include/acados)
include_directories(control/nmpc/acados_python_demo/c_generated_code)

link_directories($ENV{ACADOS_SOURCE_DIR}/lib)

aux_source_directory(application APP_SRC)
aux_source_directory(simulate SIM_SRC)
aux_source_directory(estimate ESTM_SRC)
aux_source_directory(utils UTILS_SRC)
aux_source_directory(control CONTROL_SRC)
aux_source_directory(control/wbc WBC_SRC)
aux_source_directory(control/wbc/TaskContact WBC_TASK_CONTACT_SRC)
aux_source_directory(control/nmpc/acados_cpp_demo NMPC_SRC)

add_executable(main main.cpp
                    ${APP_SRC}
                    ${SIM_SRC}
                    ${ESTM_SRC}
                    ${UTILS_SRC}
                    ${WBC_SRC}
                    ${WBC_TASK_CONTACT_SRC}
                    ${CONTROL_SRC}
                    ${NMPC_SRC}
)

target_link_libraries(main
    casadi
    mujoco::mujoco
    glfw
    OsqpEigen::OsqpEigen
    pinocchio
    yaml-cpp::yaml-cpp
    acados
    hpipm 
    blasfeo
    ${CMAKE_SOURCE_DIR}/control/nmpc/acados_python_demo/c_generated_code/libacados_ocp_solver_srbd.so
    ${CMAKE_SOURCE_DIR}/control/nmpc/acados_python_demo/c_generated_code/libacados_sim_solver_srbd.so
)

