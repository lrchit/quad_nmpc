
#pragma once

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

extern "C" {
#include "acados/utils/math.h"
#include "acados_c/sim_interface.h"

#include "acados_sim_solver_srbd.h"
}

using namespace std;
using namespace Eigen;

class QuadNmpcInterpolator {
public:
  QuadNmpcInterpolator();
  ~QuadNmpcInterpolator();

  void getSolution(const VectorXd &currentState,
                   const VectorXd &currentInputCommand,
                   VectorXd &nextStateCommand);

  void copyArray2Eigen(VectorXd &target, const double *source, int len,
                       int startIndexEigen, int startIndexArray);
  void copyEigen2Array(double *target, const VectorXd &source, int len,
                       int startIndexArray, int startIndexEigen);

private:
  int status; // acados operation state
  int horizons;
  int nx;
  int nu;

  srbd_sim_solver_capsule *srbd_sim_solver_capsule_;
  sim_config *sim_config_;
  void *sim_dims_;
  sim_in *sim_in_;
  sim_out *sim_out_;
  sim_solver *sim_solver_;

  VectorXd lastOptimalState;
};
