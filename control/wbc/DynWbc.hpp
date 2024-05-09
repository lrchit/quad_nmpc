#ifndef DYN_WBC
#define DYN_WBC

#include "ContactSpec.hpp"
#include "Task.hpp"
#include "pseudoInverse.h"
#include <OsqpEigen/OsqpEigen.h>
#include <cppTypes.h>
#include <vector>

class DynWbcExtraData {
public:
  // Output
  DVec _opt_result; // optimal deltafr deltaf
  DVec _qddot;      // final qddot
  DVec _Fr;         // final wbc reaction force

  // Input
  DVec _W_floating; // floating base acceration offset weight matrix
  DVec _W_rf;       // reaction force offset weight matrix

  DynWbcExtraData() {}
  ~DynWbcExtraData() {}
};

class DynWbc {
public:
  DynWbc(int num_qdot, const std::vector<ContactSpec *> *contact_list,
         const std::vector<Task *> *task_list)
      : num_act_joint_(num_qdot - 6), num_qdot_(num_qdot), _dim_floating(6) {
    Sa_ = DMat::Zero(num_act_joint_, num_qdot_);
    Sv_ = DMat::Zero(6, num_qdot_);

    Sa_.block(0, 6, num_act_joint_, num_act_joint_).setIdentity();
    Sv_.block(0, 0, 6, 6).setIdentity();

    _contact_list = contact_list;
    _task_list = task_list;

    _eye = DMat::Identity(num_qdot_, num_qdot_); // q size eye
    _eye_floating = DMat::Identity(_dim_floating, _dim_floating);
  }

  ~DynWbc() {}

  virtual void UpdateSetting(const DMat &A, const DMat &Ainv, const DVec &cori,
                             const DVec &grav, void *extra_setting = NULL);

  virtual void MakeTorque(DVec &cmd, void *extra_input);

protected:
  // full rank fat matrix only (inverse matrix with weight)
  void _WeightedInverse(const DMat &J, const DMat &Winv, DMat &Jinv,
                        double threshold = 0.0001) {
    DMat lambda(J * Winv * J.transpose());
    DMat lambda_inv;
    pseudoInverse(lambda, threshold, lambda_inv);
    Jinv = Winv * J.transpose() * lambda_inv;
  }

  int num_act_joint_; // Actuated joint num
  int num_qdot_;      // Virtual joint num

  DMat Sa_; // Actuated joint select matrix
  DMat Sv_; // Virtual joint select matrix

  DMat A_;
  DMat Ainv_;
  DVec cori_;
  DVec grav_;

  bool b_updatesetting_;

private:
  const std::vector<ContactSpec *> *_contact_list;
  const std::vector<Task *> *_task_list;

  void _SetEqualityConstraint(const DVec &qddot);
  void _SetInEqualityConstraint();
  void _ContactBuilding();

  void _SolveQP(); // solve QP with _Hessian, _gradient, _dyn_CE, _dyn_ce0,
                   // _dyn_CI, _dyn_ci0, _z

  void _GetSolution(const DVec &qddot, DVec &cmd);
  void _SetCost();
  void _SetOptimizationSize();

  int _dim_opt;     // Contact pt delta, delta_f,delta_fr
  int _dim_eq_cstr; // equality constraints

  int _dim_rf; // all reaction force dim
  int _dim_Uf; // all reaction force constraints matrix row

  int _dim_floating; // floating base joints dim

  DynWbcExtraData *_data;

  DMat _Hessian;
  DVec _gradient;
  DVec _z; // optimal deltaf deltafr

  DMat _dyn_CE;
  DVec _dyn_ce0;
  DMat _dyn_CI;
  DVec _dyn_ci0;

  DMat _eye;
  DMat _eye_floating;

  DMat _S_delta;
  DMat _Uf;
  DVec _Uf_ieq_vec;

  DMat _Jc;
  DVec _JcDotQdot;
  DVec _Fr_des;

  DMat _B;
  DVec _c;
  DVec task_cmd_;

  OsqpEigen::Solver qp_solver; //[qpsolver]
};
#endif