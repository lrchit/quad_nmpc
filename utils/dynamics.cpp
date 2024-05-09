#include "dynamics.hpp"

#include <orientation_tools.h>

void FBDynModel::_updateMBkinmatics() {
  // update kinmatics update Jc Jcdqd update pFoot update vFoot
  pinocchio::Data::Matrix6x J(6, model.nv);

  Vec3 xdd;
  Vec3 xd;
  Vec3 x;

  J.setZero();

  xdd.setZero();
  xd.setZero();
  x.setZero();

  // Vec3 xdd;
  // Vec18 qddZero;
  // xdd.setZero();
  // qddZero.setZero();

  pinocchio::forwardKinematics(model, data, q, v);
  pinocchio::computeJointJacobians(model, data); // or model, data, q
  pinocchio::computeJointJacobiansTimeVariation(model, data, q, v);

  pinocchio::updateFramePlacements(model, data);

  for (int i = 0; i < _foots.size(); i++) {
    int frame_id = model.getFrameId(_foots[i]);

    J.setZero();
    pinocchio::getFrameJacobian(model, data, frame_id,
                                pinocchio::LOCAL_WORLD_ALIGNED,
                                J); // LOCAL_WORLD_ALIGNED
    _Jc[i] = J.block(0, 0, 3, model.nv);

    // this is also correct
    // xdd = pinocchio::getFrameAcceleration(model, data, frame_id,
    //                                       pinocchio::LOCAL_WORLD_ALIGNED)
    //           .linear();
    // _Jcdqd[i] = xdd;
    J.setZero();
    pinocchio::getFrameJacobianTimeVariation(model, data, frame_id,
                                             pinocchio::LOCAL_WORLD_ALIGNED, J);
    // this is not the member qd here, actually the qd is v
    _Jcdqd[i] = J.block(0, 0, 3, model.nv) * v;

    xd = pinocchio::getFrameVelocity(model, data, frame_id,
                                     pinocchio::LOCAL_WORLD_ALIGNED)
             .linear();
    _vGC[i] = xd;

    x = data.oMf[frame_id].translation();
    _pGC[i] = x;

    // std::cout << _vGC[i].transpose() << std::endl;
  }
  // std::cout << _Jc.size() << std::endl;

  // std::cout << "end" << std::endl;

  // Vec3 pG = data.oMf[model.getFrameId("trunk")].translation();
  // std::cout << pG.transpose() <<std::endl;

  // Vec3 rpytheta
  // =pinocchio::rpy::matrixToRpy(data.oMf[model.getFrameId("trunk")].rotation());
  // std::cout << rpytheta.transpose() <<std::endl;
}

void FBDynModel::_updateMBdynamics() {
  _A = pinocchio::crba(model, data, q);
  _A.triangularView<Eigen::StrictlyLower>() =
      _A.transpose().triangularView<Eigen::StrictlyLower>();
  // std::cout <<"[A]"<< std::endl<< _A.block(0,0,6,6) << std::endl;
  _Ainv = _A.inverse();
  _coriolis = pinocchio::computeCoriolisMatrix(model, data, q, v) * v;
  _grav = pinocchio::computeGeneralizedGravity(model, data, q);
}

void FBDynModel::updateModel(const FBModelState &state) {
  _state = state;

  _full_config.segment(0, 3) = _state.bodyPosition;
  Vec3 rpy = ori::quatToRPY(_state.bodyOrientation);
  _full_config.segment(3, 3) = rpy;
  _full_config.segment(6, 12) = _state.q;

  // std::cout <<"[RPY]"<< rpy.transpose() <<std::endl;
  // std::cout << _full_config.transpose() <<std::endl;

  q.segment(0, 3) = _state.bodyPosition;

  Quat bodyQuat; // in pinocchio quat is [xyzw] but eigen [wxyz]
  bodyQuat.segment(0, 3) = _state.bodyOrientation.segment(1, 3);
  bodyQuat(3) = _state.bodyOrientation(0);

  q.segment(3, 4) = bodyQuat;
  q.segment(7, 12) = _state.q;

  v.segment(0, 3) = _state.bodyVelocity.segment(3, 3); // local frame linear vel
  v.segment(3, 3) =
      _state.bodyVelocity.segment(0, 3); // local frame angular vel
  v.segment(6, 12) = _state.qd;

  // Mat3 RotMat = ori::quaternionToRotationMatrix(_state.bodyOrientation);
  // qd = v;
  // qd.segment(0, 3) = RotMat.transpose() * v.segment(0, 3);

  _updateMBkinmatics();
  _updateMBdynamics();
}