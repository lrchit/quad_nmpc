/* related header files */
#pragma once
#include <pinocchio/fwd.hpp> // always include it before any other header

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/spatial/explog.hpp"

#include "pseudoInverse.h"

using namespace Eigen;

class LegKinematics {
public:
  LegKinematics(const std::string urdf_filename) {
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data temp_data(model);

    data = temp_data;
  }
  ~LegKinematics(){};

  void forwardKinFrame(Eigen::Vector3d q, Eigen::Vector3d dq,
                       Eigen::Vector3d ddq, Eigen::Vector3d &x,
                       Eigen::Vector3d &dx, Eigen::Vector3d &ddx,
                       pinocchio::FrameIndex id);

  void
  inverseKinFrame(Eigen::Vector3d &q, Eigen::Vector3d &dq, Eigen::Vector3d &ddq,
                  Eigen::Vector3d x, Eigen::Vector3d dx, Eigen::Vector3d ddx,
                  pinocchio::FrameIndex id,
                  Eigen::Vector3d q_init = Eigen::Vector3d(0, 0.92, -1.84));

  Eigen::Matrix3d getJacobian(Eigen::Vector3d q, pinocchio::FrameIndex id);

private:
  pinocchio::Model model;
  pinocchio::Data data;
};
