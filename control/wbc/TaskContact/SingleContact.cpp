#include "SingleContact.hpp"

// [ Fx, Fy, Fz ]
SingleContact::SingleContact(const FBDynModel *robot, int pt)
    : ContactSpec(3), _max_Fz(500.), _contact_pt(pt), _dim_U(6) {
  Contact::idx_Fz_ = 2;
  robot_sys_ = robot;
  Contact::Jc_ = DMat(Contact::dim_contact_, 18); // configure space dim
  Contact::JcDotQdot_ = DVec::Zero(Contact::dim_contact_);
  Contact::Uf_ = DMat::Zero(_dim_U, Contact::dim_contact_);

  double mu(0.4);

  Contact::Uf_(0, 2) = 1.;

  Contact::Uf_(1, 0) = 1.;
  Contact::Uf_(1, 2) = mu;
  Contact::Uf_(2, 0) = -1.;
  Contact::Uf_(2, 2) = mu;

  Contact::Uf_(3, 1) = 1.;
  Contact::Uf_(3, 2) = mu;
  Contact::Uf_(4, 1) = -1.;
  Contact::Uf_(4, 2) = mu;

  // Upper bound of normal force
  Contact::Uf_(5, 2) = -1.;
}

SingleContact::~SingleContact() {}

bool SingleContact::_UpdateJc() {
  Contact::Jc_ = robot_sys_->_Jc[_contact_pt];

  // Quat quat = robot_sys_->_state.bodyOrientation;
  // Mat3 Rot = ori::quaternionToRotationMatrix(quat);
  // Contact::Jc_.block(0,3, 3,3) = Rot*Contact::Jc_.block(0,3,3,3);

  // Contact::Jc_.block(0,0, 3,3) = Rot.transpose()*Contact::Jc_.block(0,0,3,3);

  return true;
}

bool SingleContact::_UpdateJcDotQdot() {
  Contact::JcDotQdot_ = robot_sys_->_Jcdqd[_contact_pt];
  // std::cout << Contact::JcDotQdot_.transpose() << std::endl;
  return true;
}

bool SingleContact::_UpdateUf() { return true; }

bool SingleContact::_UpdateInequalityVector() {
  Contact::ieq_vec_ = DVec::Zero(_dim_U);
  Contact::ieq_vec_[5] = -_max_Fz;
  return true;
}
