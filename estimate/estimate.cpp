#include "estimate.h"
#include <chrono>
#include <iostream>

QuadEstm::QuadEstm(std::vector<LegKinematics> legKin) {
  double dt = 0.002;
  xhat_.setZero();
  ps_.setZero();
  vs_.setZero();
  A_.setZero();
  A_.block(0, 0, 3, 3) = Matrix3d::Identity();
  A_.block(0, 3, 3, 3) = dt * Matrix3d::Identity();
  A_.block(3, 3, 3, 3) = Matrix3d::Identity();
  A_.block(6, 6, 12, 12) = Matrix<double, 12, 12>::Identity();
  B_.setZero();
  B_.block(3, 0, 3, 3) = dt * Matrix3d::Identity();
  Matrix<double, Dynamic, Dynamic> C1(3, 6);
  C1 << Matrix3d::Identity(), Matrix3d::Zero();
  Matrix<double, Dynamic, Dynamic> C2(3, 6);
  C2 << Matrix3d::Zero(), Matrix3d::Identity();
  C_.setZero();
  C_.block(0, 0, 3, 6) = C1;
  C_.block(3, 0, 3, 6) = C1;
  C_.block(6, 0, 3, 6) = C1;
  C_.block(9, 0, 3, 6) = C1;
  C_.block(0, 6, 12, 12) = -1 * Matrix<double, 12, 12>::Identity();
  C_.block(12, 0, 3, 6) = C2;
  C_.block(15, 0, 3, 6) = C2;
  C_.block(18, 0, 3, 6) = C2;
  C_.block(21, 0, 3, 6) = C2;
  C_(27, 17) = 1;
  C_(26, 14) = 1;
  C_(25, 11) = 1;
  C_(24, 8) = 1;
  P_.setIdentity();
  P_ = 100 * P_;
  Q0_.setIdentity();
  Q0_.block(0, 0, 3, 3) = (dt / 20) * Matrix3d::Identity();
  Q0_.block(3, 3, 3, 3) = (dt * 9.8 / 20) * Matrix3d::Identity();
  Q0_.block(6, 6, 12, 12) = dt * Matrix<double, 12, 12>::Identity();
  R0_.setIdentity();

  legKin_ = legKin;

  YAML::Node config = YAML::LoadFile("../estimate/esti_config.yaml");
  processNoisePimu = config["imu_process_noise_position"].as<double>();
  processNoiseVimu = config["imu_process_noise_velocity"].as<double>();
  processNoisePfoot = config["foot_process_noise_position"].as<double>();
  sensorNoisePimuRelFoot = config["foot_sensor_noise_position"].as<double>();
  sensorNoiseVimuRelFoot = config["foot_sensor_noise_velocity"].as<double>();
  sensorNoiseZFoot = config["foot_height_sensor_noise"].as<double>();

  cheaterMode = config["cheaterMode"].as<bool>();
}

// 调用状态估计器
void QuadEstm::callStateEstimator(QuadState &state, double *sensordata) {
  if (cheaterMode)
    cheaterComputeState(state, sensordata);
  else
    linearKFComputeState(state, sensordata);
}

// 线性卡尔曼滤波
void QuadEstm::linearKFComputeState(QuadState &state, double *sensordata) {

  // std::cout << "************* LinearKF *************" << std::endl;

  // --- get rpy ---
  int quat_sensor_adr = 3;
  Vector<double, 4> q(
      *(sensordata + quat_sensor_adr + 0), *(sensordata + quat_sensor_adr + 1),
      *(sensordata + quat_sensor_adr + 2), *(sensordata + quat_sensor_adr + 3));
  state.rpy = ori::quatToRPY(q);
  state.rotMat = ori::quaternionToRotationMatrix(q);

  // --- get angvel ---
  int angvel_sensor_adr = 10;
  for (int i = 0; i < 3; ++i) {
    state.omegaWorld(i) = sensordata[angvel_sensor_adr + i];
  }
  state.omegaWorld = state.rotMat.transpose() * state.omegaWorld;

  // --- get linacc ---
  int linacc_sensor_adr = 13;
  for (int i = 0; i < 3; ++i) {
    state.linAcc(i) = sensordata[linacc_sensor_adr + i];
  }
  state.linAcc = state.rotMat.transpose() * state.linAcc;

  // --- get legQPos ---
  int qpos_sensor_adr = 16;
  for (int i = 0; i < 4; ++i) {
    Vector3d one_legQPos;
    one_legQPos << *(sensordata + qpos_sensor_adr + 3 * i),
        *(sensordata + qpos_sensor_adr + 3 * i + 1),
        *(sensordata + qpos_sensor_adr + 3 * i + 2);
    state.legQPos.col(i) = one_legQPos;
  }

  // --- get legQVel ---
  int qvel_sensor_adr = 28;
  for (int i = 0; i < 4; ++i) {
    Vector3d one_legQVel;
    one_legQVel << *(sensordata + qvel_sensor_adr + 3 * i),
        *(sensordata + qvel_sensor_adr + 3 * i + 1),
        *(sensordata + qvel_sensor_adr + 3 * i + 2);
    state.legQVel.col(i) = one_legQVel;
  }

  // --- contact ---
  int contact_sensor_adr = 40;
  bool isContact[4];
  for (int i = 0; i < 4; ++i) {
    if (*(sensordata + contact_sensor_adr + i) > 0) {
      isContact[i] = true;
    } else {
      isContact[i] = false;
    }
  }

  Matrix<double, 18, 18> Q = Matrix<double, 18, 18>::Identity();
  Q.block(0, 0, 3, 3) = Q0_.block(0, 0, 3, 3) * processNoisePimu;
  Q.block(3, 3, 3, 3) = Q0_.block(3, 3, 3, 3) * processNoiseVimu;
  Q.block(6, 6, 12, 12) = Q0_.block(6, 6, 12, 12) * processNoisePfoot;

  Matrix<double, 28, 28> R = Matrix<double, 28, 28>::Identity();
  R.block(0, 0, 12, 12) = R0_.block(0, 0, 12, 12) * sensorNoisePimuRelFoot;
  R.block(12, 12, 12, 12) = R0_.block(12, 12, 12, 12) * sensorNoiseVimuRelFoot;
  R.block(24, 24, 4, 4) = R0_.block(24, 24, 4, 4) * sensorNoiseZFoot;

  int qindex = 0;
  int rindex1 = 0;
  int rindex2 = 0;
  int rindex3 = 0;

  Vector3d g(0, 0, -9.81);
  Matrix3d Rbod = state.rotMat.transpose();
  Vector3d a = state.linAcc + g;
  Vector<double, 4> pzs = Vector<double, 4>::Zero();
  Vector<double, 4> trusts = Vector<double, 4>::Zero();
  Vector3d p0, v0;
  p0 << xhat_[0], xhat_[1], xhat_[2];
  v0 << xhat_[3], xhat_[4], xhat_[5];

  for (int i = 0; i < 4; ++i) {
    int i1 = 3 * i;
    Vector3d p_rel, dp_rel, ddp_rel;
    legKin_[i].forwardKinFrame(state.legQPos.col(i), state.legQVel.col(i),
                               Vector3d::Zero(), p_rel, dp_rel, ddp_rel, 14);

    state.footPos.col(i) = p_rel;
    state.footVel.col(i) = dp_rel;

    Vector3d p_f = Rbod * p_rel;
    Vector3d dp_f =
        Rbod * ((state.rotMat * state.omegaWorld).cross(p_rel) + dp_rel);

    qindex = 6 + i1;
    rindex1 = i1;
    rindex2 = 12 + i1;
    rindex3 = 24 + i;

    double high_suspect_number(100);

    Q.block(qindex, qindex, 3, 3) = (isContact[i] ? 1. : high_suspect_number) *
                                    Q.block(qindex, qindex, 3, 3);
    R.block(rindex1, rindex1, 3, 3) = 1 * R.block(rindex1, rindex1, 3, 3);
    R.block(rindex2, rindex2, 3, 3) =
        (isContact[i] ? 1. : high_suspect_number) *
        R.block(rindex2, rindex2, 3, 3);
    R(rindex3, rindex3) =
        (isContact[i] ? 1. : high_suspect_number) * R(rindex3, rindex3);

    ps_.segment(i1, 3) = -p_f;
    vs_.segment(i1, 3) = -dp_f;

    pzs[i] = isContact[i] ? footRadius_ : xhat_[2] - p_f[2];
  }

  Matrix<double, 28, 1> y;
  y << ps_, vs_, pzs;
  xhat_ = A_ * xhat_ + B_ * a;
  Matrix<double, 18, 18> At = A_.transpose();
  Matrix<double, 18, 18> Pm = A_ * P_ * At + Q;
  Matrix<double, 18, 28> Ct = C_.transpose();
  Matrix<double, 28, 1> yModel = C_ * xhat_;
  Matrix<double, 28, 1> ey = y - yModel;
  Matrix<double, 28, 28> S = C_ * Pm * Ct + R;

  // todo compute LU only once
  Matrix<double, 28, 1> S_ey = S.lu().solve(ey);
  xhat_ += Pm * Ct * S_ey;

  Matrix<double, 28, 18> S_C = S.lu().solve(C_);
  P_ = (Matrix<double, 18, 18>::Identity() - Pm * Ct * S_C) * Pm;

  Matrix<double, 18, 18> Pt = P_.transpose();
  P_ = (P_ + Pt) / 2;

  if (P_.block(0, 0, 2, 2).determinant() > 0.000001) {
    P_.block(0, 2, 2, 16).setZero();
    P_.block(2, 0, 16, 2).setZero();
    P_.block(0, 0, 2, 2) /= 10;
  }

  // update pos, linVel, footPos
  state.pos = xhat_.block(0, 0, 3, 1);
  state.linVel = xhat_.block(3, 0, 3, 1);
  for (int i = 0; i < 4; ++i) {
    state.footPosWorld.col(i) =
        state.rotMat.transpose() * state.footPos.col(i) + state.pos;
    state.footVelWorld.col(i) =
        state.rotMat.transpose() * state.footVel.col(i) + state.linVel;

    if (isContact[i]) {
      state.footPosWorld(2, i) = 0.005;
    }
  }

  // // 与真实值误差
  // Vector3d pos_real =
  //     Vector3d({*sensordata, *(sensordata + 1), *(sensordata + 2)});
  // std::cout << "pos_err = \n" << state.pos - pos_real << std::endl;

  // std::cout << "pz_esti = " << state.footPosWorld(2, 0)
  //           << state.footPosWorld(2, 1) << state.footPosWorld(2, 2)
  //           << state.footPosWorld(2, 3) << std::endl;
}

// 直接读的传感器
void QuadEstm::cheaterComputeState(QuadState &state, double *sensordata) {

  // std::cout << "************* cheaterMode *************" << std::endl;

  // --- get pos ---
  int pos_sensor_adr = 0;
  state.pos = Vector3d({*sensordata, *(sensordata + 1), *(sensordata + 2)});

  // --- get rpy ---
  int quat_sensor_adr = 3;
  Vector<double, 4> q(
      *(sensordata + quat_sensor_adr + 0), *(sensordata + quat_sensor_adr + 1),
      *(sensordata + quat_sensor_adr + 2), *(sensordata + quat_sensor_adr + 3));
  state.rpy = ori::quatToRPY(q);
  state.rotMat = ori::quaternionToRotationMatrix(q);

  // --- get linvel ---
  int linvel_sensor_adr = 7;
  for (int i = 0; i < 3; ++i) {
    state.linVel(i) = sensordata[linvel_sensor_adr + i];
  }

  // --- get angvel ---
  int angvel_sensor_adr = 10;
  for (int i = 0; i < 3; ++i) {
    state.omegaWorld(i) = sensordata[angvel_sensor_adr + i];
  }
  state.omegaWorld = state.omegaWorld;

  // --- get legQPos ---
  int qpos_sensor_adr = 16;
  for (int i = 0; i < 4; ++i) {
    Vector3d one_legQPos;
    one_legQPos << *(sensordata + qpos_sensor_adr + 3 * i),
        *(sensordata + qpos_sensor_adr + 3 * i + 1),
        *(sensordata + qpos_sensor_adr + 3 * i + 2);
    state.legQPos.col(i) = one_legQPos;
  }

  // --- get legQVel ---
  int qvel_sensor_adr = 28;
  for (int i = 0; i < 4; ++i) {
    Vector3d one_legQVel;
    one_legQVel << *(sensordata + qvel_sensor_adr + 3 * i),
        *(sensordata + qvel_sensor_adr + 3 * i + 1),
        *(sensordata + qvel_sensor_adr + 3 * i + 2);
    state.legQVel.col(i) = one_legQVel;
  }

  // --- get footPos, footVel ---
  for (int i = 0; i < 4; ++i) {
    Vector3d p_rel, dp_rel, ddp_rel;
    legKin_[i].forwardKinFrame(state.legQPos.col(i), state.legQVel.col(i),
                               Vector3d::Zero(), p_rel, dp_rel, ddp_rel, 14);
    state.footPos.col(i) = p_rel;
    state.footPosWorld.col(i) =
        state.rotMat.transpose() * state.footPos.col(i) + state.pos;
    state.footVel.col(i) = dp_rel;
    state.footVelWorld.col(i) =
        state.rotMat.transpose() * state.footVel.col(i) + state.linVel;
  }
}