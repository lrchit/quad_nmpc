
#pragma once

#include <chrono>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

extern "C" {
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"

#include "acados_solver_srbd.h"
}

using namespace std;
using namespace Eigen;

class QuadNmpc {
public:
  QuadNmpc();
  ~QuadNmpc();

  void nmpcUpdate(const VectorXd &currentState, const VectorXd &gaitTable,
                  const VectorXd &desiredComState,
                  const VectorXd &desiredFootZPos,
                  const VectorXd &desiredFootZVel);

  void getSolution(VectorXd &optimalState, VectorXd &optimalControl);

private:
  void copyArray2Eigen(VectorXd &target, const double *source, int len,
                       int startIndexEigen, int startIndexArray);
  void copyEigen2Array(double *target, const VectorXd &source, int len,
                       int startIndexArray, int startIndexEigen);

  Vector2d getFootProjectReference(const double yaw, Vector2d pos);

  int status; // acados operation state
  int horizons;
  int nx;
  int nu;

  double reference[48];
  int referenceIndex[48];

  double mass;
  double g;

  bool useTimeStepOpti;

  Vector4d Xshoulder;
  Vector4d Yshoulder;

  srbd_solver_capsule *acados_ocp_capsule;
  ocp_nlp_config *nlp_config;
  ocp_nlp_dims *nlp_dims;
  ocp_nlp_in *nlp_in;
  ocp_nlp_out *nlp_out;

  VectorXd lastOptimalState;
  VectorXd lastOptimalControl;
};
