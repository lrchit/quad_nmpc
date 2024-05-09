
#pragma once

#include "kinematics.h"
#include "orientation_tools.h"
#include "state.h"
#include <yaml-cpp/yaml.h>

using namespace Eigen;

class QuadEstm {
public:
  QuadEstm(std::vector<LegKinematics> legKin);

  void callStateEstimator(QuadState &state, double *sensordata);

private:
  void linearKFComputeState(QuadState &state, double *sensordata);
  void cheaterComputeState(QuadState &state, double *sensordata);

  bool cheaterMode;

  Matrix<double, 18, 1> xhat_;
  Matrix<double, 12, 1> ps_;
  Matrix<double, 12, 1> vs_;
  Matrix<double, 18, 18> A_;
  Matrix<double, 18, 18> Q0_;
  Matrix<double, 18, 18> P_;
  Matrix<double, 28, 28> R0_;
  Matrix<double, 18, 3> B_;
  Matrix<double, 28, 18> C_;

  double processNoisePimu;
  double processNoiseVimu;
  double processNoisePfoot;
  double sensorNoisePimuRelFoot;
  double sensorNoiseVimuRelFoot;
  double sensorNoiseZFoot;

  double footRadius_ = 0.005;

  std::vector<LegKinematics> legKin_;
};