#include "kinematics.h"

void LegKinematics::forwardKinFrame(Eigen::Vector3d q, Eigen::Vector3d dq,
                                    Eigen::Vector3d ddq, Eigen::Vector3d &x,
                                    Eigen::Vector3d &dx, Eigen::Vector3d &ddx,
                                    pinocchio::FrameIndex id) {
  pinocchio::forwardKinematics(model, data, q, dq, ddq);
  pinocchio::updateFramePlacements(model, data);
  // pinocchio::framesForwardKinematics(model, data, q);
  x = data.oMf[id].translation();
  dx = pinocchio::getFrameVelocity(model, data, id,
                                   pinocchio::LOCAL_WORLD_ALIGNED)
           .linear();

  ddx = pinocchio::getFrameAcceleration(model, data, id,
                                        pinocchio::LOCAL_WORLD_ALIGNED)
            .linear();

  // std::cout<<'\n'<<"x: "<<x.transpose()<< '\n'<<"dx "<< \
    // dx.transpose()<<'\n'<<"ddx"<<ddx.transpose()<<std::endl;
}

void LegKinematics::inverseKinFrame(Eigen::Vector3d &q, Eigen::Vector3d &dq,
                                    Eigen::Vector3d &ddq, Eigen::Vector3d x,
                                    Eigen::Vector3d dx, Eigen::Vector3d ddx,
                                    pinocchio::FrameIndex id,
                                    Eigen::Vector3d q_init) {
  Eigen::Vector3d temp_q = q_init;
  const double eps = 1e-3;
  const int IT_MAX = 1000;
  const double DT = 0.1;

  pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), x);
  bool success = false;

  Eigen::Vector3d err;
  pinocchio::Data::Matrix6x J(6, model.nv);
  pinocchio::Data::Matrix6x dJ(6, model.nv);
  // Eigen::Matrix<double,6,6> dJ;
  J.setZero();
  dJ.setZero();
  int i = 0;
  for (i = 0; i < IT_MAX; i++) {
    pinocchio::forwardKinematics(model, data, temp_q);
    pinocchio::computeJointJacobians(
        model, data,
        temp_q); // if we use this ,we can get frame or joint jacobian

    pinocchio::updateFramePlacements(model, data);
    pinocchio::getFrameJacobian(model, data, id, pinocchio::LOCAL_WORLD_ALIGNED,
                                J); // cal this, you need forwardKinematics

    pinocchio::SE3 T_temp = data.oMf[id];
    err = oMdes.translation() -
          T_temp.translation(); // must compute the forward kinematics

    if (err.norm() < eps) {
      success = true;
      break;
    }
    if (i >= IT_MAX) {
      success = false;
      break;
    }

    Eigen::Vector3d delta_q = J.block(0, 0, 3, 3).inverse() * err;
    temp_q = pinocchio::integrate(model, temp_q, delta_q * DT);

    //  if(!(i%10))
    //  std::cout << i << ": error = " << err.transpose() << std::endl;
  }

  // if (success) {
  //   std::cout << "Convergence achieved!" << std::endl;
  // } else {
  //   std::cout << "\nWarning: the iterative algorithm has not reached "
  //                "convergence to the desired precision"
  //             << std::endl;
  // }
  q = temp_q;

  dq = J.block(0, 0, 3, 3).inverse() * dx;
  pinocchio::computeJointJacobiansTimeVariation(
      model, data, q, dq); // if we use this .we can get frame or joint
                           // derviation of jocabians
  pinocchio::getFrameJacobianTimeVariation(
      model, data, id, pinocchio::LOCAL_WORLD_ALIGNED,
      dJ); // this can repace getframejacobian

  ddq = J.block(0, 0, 3, 3).inverse() * (ddx - (dJ * dq).block(0, 0, 3, 3));
  // std::cout<<"dj"<<dJ<<std::endl;
  // std::cout << "\nresult: " << q.transpose() << dq.transpose() << "ddq"
  //           << ddq.transpose() << std::endl;
  // std::cout << "    i =   " << i << std::endl;
  // std::cout << "\nfinal error: " << err.transpose() << std::endl;
}

Eigen::Matrix3d LegKinematics::getJacobian(Eigen::Vector3d q,
                                           pinocchio::FrameIndex id) {
  pinocchio::Data::Matrix6x J(6, model.nv);
  // pinocchio::Data::Matrix6x dJ(6,model.nv);
  // Eigen::Matrix<double,6,6> dJ;
  J.setZero();
  pinocchio::forwardKinematics(model, data, q);

  pinocchio::computeJointJacobians(
      model, data, q); // if we use this ,we can get frame or joint jocobians
  pinocchio::updateFramePlacements(model, data);

  pinocchio::getFrameJacobian(model, data, id, pinocchio::LOCAL_WORLD_ALIGNED,
                              J); // cal this, you need forwardKinematics

  return J.block(0, 0, 3, 3);
}
