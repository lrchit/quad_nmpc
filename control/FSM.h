#ifndef QUADFSM_H
#define QUADFSM_H

#pragma once
#include "Gait.h"

#include "body_optimization.h"
#include "estimate.h"
#include "kinematics.h"

#include "nmpc_caller.h"
#include "nmpc_integrator.h"

#include "orientation_tools.h"
#include "state.h"
#include <yaml-cpp/yaml.h>

#include "WbcCtrl.hpp"

#define pi 3.1415926535

class QuadFSM {
public:
  QuadFSM();
  ~QuadFSM(){};

  void mainProgram();
  void computeNmpc();
  void computeWbc();
  void bodyOptimize();
  void stateEstimate(double *sensordata);
  Matrix<double, 12, 1> getJointTorque();

  bool useWbc;
  bool bodyOptiNeeded;
  bool nmpcUpdateNeeded;
  bool wbcUpdateNeeded;

  QuadState stateCur;

private:
  void computeJointTorque();
  void updateMpcData();
  void updateNmpc();
  void updateWbcData();

  // stateEstimate，mpc，wbc，body_opti
  QuadEstm *estimater;
  QuadWbc *wbcCaller;
  QuadBodyOpti *bodyOptimizer;
  QuadNmpc *nmpcCaller;
  QuadNmpcInterpolator *nmpcIntegrator;

  std::vector<LegKinematics> legKin;
  WbcData wbcData;

  Matrix3d kp_, kd_;
  Matrix<double, 3, 4> footForceFromKin;

  int counter;

  double dt;
  double dtmpc;
  int horizon;
  int iterationBetweenMpc;
  int iterationBetweenBodyOpti;

  // for body_opti
  Matrix<double, 12, 1> footHoldDes;
  double weightX, weightY, weightZ;

  // gait
  Gait *gait, *currentGait;
  OffsetDurationGait *trotting, *standing, *flying, *bounding, *trot_climbing;
  int iterationCounter;
  int gaitIterationCounter;

  Vector<double, 6> trajIntegrate;

  // nmpc
  VectorXd currentState;
  VectorXd desiredComState;
  VectorXd desiredFootZPos;
  VectorXd desiredFootZVel;
  VectorXd gaitTable;
  VectorXd currentfootVel;
  VectorXd optimalControl;
  VectorXd optimalState;
  MatrixXd phaseTable;
};
#endif
