#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#include "cppTypes.h"
#include <string>
#include <vector>

#include <pinocchio/fwd.hpp> // always include it before any other header

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/spatial/explog.hpp"

#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/multibody/model.hpp"

#include "pinocchio/math/rpy.hpp"

struct FBModelState {
  Quat bodyOrientation; // world frame [wxyz]
  Vec3 bodyPosition;    // world frame
  SVec bodyVelocity;    // local frame [ang_vel, linVel]
  Vec12 q;  // q here is 12 dof and the order should be "FL", "FR", "RL", "RR"
  Vec12 qd; // qd here is 12 dof
};

class FBDynModel {
public:
  // vector<D6Mat<T>> _J;
  // vector<SVec<T>> _Jdqd;
  FBDynModel(
      const std::string urdf_filename = "../unitree_a1_desc/urdf/a1.urdf") {

    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data temp_data(model);
    data = temp_data;

    _Jc.resize(4);
    _Jcdqd.resize(4);
    _pGC.resize(4);
    _vGC.resize(4);

    // std::cout<<
    // model.getFrameId("FL_hip")<<"," <<
    // model.getFrameId("FL_thigh")<<"," <<
    // model.getFrameId("FL_calf")<<"," <<
    // model.getFrameId("FL_foot")<<"," <<
    // model.getFrameId("FR_foot")<<"," <<
    // model.getFrameId("RL_foot")<<"," <<
    // model.getFrameId("RR_foot")<<"," <<
    // std::endl;
  }

  ~FBDynModel() {}

  std::vector<D3Mat> _Jc; // contact jacobian list (foot point jacobian)
  D3Mat _Jcd;             // contact jacobian time derivation
  std::vector<Vec3> _Jcdqd;
  FBModelState _state;
  std::vector<Vec3> _pGC; // foot point position (global)
  std::vector<Vec3> _vGC; // foot point velocity (global)

  pinocchio::Model model;
  pinocchio::Data data;

  DMat _A;
  DMat _Ainv;
  DVec _grav;
  DVec _coriolis;
  Vec18 _full_config;

  void updateModel(const FBModelState &state);

protected:
  Vec19
      q; // first 7 [global_base_position, global_base_quaternion] in pinocchio
  Vec18 v; // first 6 [local_base_velocity_linear, local_base_velocity_angular]
           // in pinocchio
  Vec18 qd;

  std::vector<std::string> _foots = {"FL_foot", "FR_foot", "RL_foot",
                                     "RR_foot"};

  void _updateMBkinmatics();
  void _updateMBdynamics();
};

#endif