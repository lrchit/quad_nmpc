#include "BodyPosTask.hpp"
// (X, Y, Z)
#include "orientation_tools.h"
BodyPosTask::BodyPosTask(const FBDynModel *robot) : Task(3), _robot_sys(robot) {
  TK::Jt_ = DMat::Zero(TK::dim_task_, 18); // configure space dim
  // TK::Jt_.block(0, 3, 3, 3).setIdentity();
  TK::Jt_.block(0, 0, 3, 3).setIdentity();
  TK::JtDotQdot_ = DVec::Zero(TK::dim_task_);

  _Kp_kin = DVec::Constant(TK::dim_task_, 1.);
  _Kp = DVec::Constant(TK::dim_task_, 50.);
  _Kd = DVec::Constant(TK::dim_task_, 1.0);
}

BodyPosTask::~BodyPosTask() {}

bool BodyPosTask::_UpdateCommand(const void *pos_des, const DVec &vel_des,
                                 const DVec &acc_des) {
  Vec3 *pos_cmd = (Vec3 *)pos_des;
  Vec3 link_pos = _robot_sys->_state.bodyPosition;
  // std::cout << (*pos_cmd).transpose() <<std::endl;

  Quat quat = _robot_sys->_state.bodyOrientation;
  Mat3 Rot = ori::quaternionToRotationMatrix(quat);

  SVec curr_vel = _robot_sys->_state.bodyVelocity;
  curr_vel.tail(3) = Rot.transpose() * curr_vel.tail(3);

  // X, Y, Z
  for (int i(0); i < 3; ++i) {
    TK::pos_err_[i] = _Kp_kin[i] * ((*pos_cmd)[i] - link_pos[i]);
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * ((*pos_cmd)[i] - link_pos[i]) +
                     _Kd[i] * (TK::vel_des_[i] - curr_vel[i + 3]) +
                     TK::acc_des_[i];
  }

  // std::cout << TK::op_cmd_.transpose() <<std::endl;
  // Quat quat = _robot_sys->_state.bodyOrientation;
  // Mat3 Rot = ori::quaternionToRotationMatrix(quat);
  // TK::pos_err_ = Rot * TK::pos_err_;
  // TK::vel_des_ = Rot * TK::vel_des_;
  // TK::acc_des_ = Rot * TK::acc_des_;
  return true;
}

bool BodyPosTask::_UpdateTaskJacobian() {
  Quat quat = _robot_sys->_state.bodyOrientation;
  Mat3 Rot = ori::quaternionToRotationMatrix(quat);
  // TK::Jt_.block(0, 3, 3, 3) = Rot.transpose();
  TK::Jt_.block(0, 0, 3, 3) = Rot.transpose();
  // TK::Jt_.block(0,3, 3,3) = Rot;
  // TK::Jt_.block(0,3, 3,3) = Rot*TK::Jt_.block(0,3,3,3);
  return true;
}

bool BodyPosTask::_UpdateTaskJDotQdot() { return true; }