#include "BodyOriTask.hpp"
// (Rx, Ry, Rz)
#include <orientation_tools.h>

BodyOriTask::BodyOriTask(const FBDynModel *robot) : Task(3), _robot_sys(robot) {
  TK::Jt_ = DMat::Zero(TK::dim_task_, 18); // 18 is configure space size
  // TK::Jt_.block(0, 0, 3, 3).setIdentity();
  TK::Jt_.block(0, 3, 3, 3).setIdentity();
  TK::JtDotQdot_ = DVec::Zero(TK::dim_task_);

  _Kp_kin = DVec::Constant(TK::dim_task_, 1.); // not very useful
  _Kp = DVec::Constant(TK::dim_task_, 50.);    // task space kp
  _Kd = DVec::Constant(TK::dim_task_, 1.);     // task space kd
}

BodyOriTask::~BodyOriTask() {}

// oritation is presented in quaterion
// angle vel, acc is presented in Vec3
bool BodyOriTask::_UpdateCommand(const void *pos_des, const DVec &vel_des,
                                 const DVec &acc_des) {
  Quat *ori_cmd = (Quat *)pos_des;
  Quat link_ori = (_robot_sys->_state.bodyOrientation);

  Quat link_ori_inv;
  link_ori_inv[0] = link_ori[0];
  link_ori_inv[1] = -link_ori[1];
  link_ori_inv[2] = -link_ori[2];
  link_ori_inv[3] = -link_ori[3];
  // link_ori_inv /= link_ori.norm();

  // Explicit because operational space is in global frame
  Quat ori_err = ori::quatProduct(*ori_cmd, link_ori_inv);

  // std::cout <<"[RPY]"<< ori::quatToRPY(link_ori).transpose() << std::endl;
  // std::cout <<"[RPY]"<< ori::quatToRPY(*ori_cmd).transpose() << std::endl;

  if (ori_err[0] < 0.) {
    ori_err *= (-1.);
  }
  Vec3 ori_err_so3;
  ori::quaternionToso3(ori_err, ori_err_so3);

  Mat3 Rot = ori::quaternionToRotationMatrix(link_ori);
  SVec curr_vel = _robot_sys->_state.bodyVelocity;
  curr_vel.head(3) =
      Rot.transpose() *
      curr_vel.head(3); //_state.bodyVelocity是body系的需要转为global

  // Configuration space: Local
  // Operational Space: Global

  // Rx, Ry, Rz
  for (int i(0); i < 3; ++i) {
    TK::pos_err_[i] = _Kp_kin[i] * ori_err_so3[i];
    TK::vel_des_[i] = vel_des[i];
    TK::acc_des_[i] = acc_des[i];

    TK::op_cmd_[i] = _Kp[i] * ori_err_so3[i] +
                     _Kd[i] * (TK::vel_des_[i] - curr_vel[i]) + TK::acc_des_[i];
  }
  return true;
}

bool BodyOriTask::_UpdateTaskJacobian() {
  Quat quat = _robot_sys->_state.bodyOrientation;
  Mat3 Rot = ori::quaternionToRotationMatrix(quat);
  // TK::Jt_.block(0, 0, 3, 3) = Rot.transpose();
  TK::Jt_.block(0, 3, 3, 3) = Rot.transpose();
  // pretty_print(Rot, std::cout, "Rot mat");
  return true;
}

bool BodyOriTask::_UpdateTaskJDotQdot() { return true; }