#pragma once
#include <Eigen/Dense>

struct QuadState {
  QuadState() {
    rpy.setZero();
    pos.setZero();
    omegaWorld.setZero();
    linVel.setZero();
    linAcc.setZero();
    rotMat.setZero();
    footPos.setZero();
    footVel.setZero();
    footPosWorld.setZero();
    footVelWorld.setZero();
    legQPos.setZero();
    legQVel.setZero();
    contactPhase.setZero();
    grfRef.setZero();
    jointTorque.setZero();
  }

  Eigen::Vector3d rpy;
  Eigen::Vector3d pos;
  // this 3 are in world frame
  Eigen::Vector3d omegaWorld;
  Eigen::Vector3d linVel;
  Eigen::Vector3d linAcc;

  Eigen::Matrix3d rotMat;
  Eigen::Matrix<double, 3, 4> footPos;
  Eigen::Matrix<double, 3, 4> footVel;
  Eigen::Matrix<double, 3, 4> footPosWorld;
  Eigen::Matrix<double, 3, 4> footVelWorld;
  Eigen::Matrix<double, 3, 4> legQPos;
  Eigen::Matrix<double, 3, 4> legQVel;

  Eigen::Matrix<double, 4, 1> contactPhase;

  Eigen::Matrix<double, 12, 1> grfRef;
  Eigen::Matrix<double, 12, 1> jointTorque;
};
