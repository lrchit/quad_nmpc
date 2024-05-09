#include "DynWbc.hpp"
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>

void DynWbc::UpdateSetting(const DMat &A, const DMat &Ainv, const DVec &cori,
                           const DVec &grav, void *extra_setting) {
  A_ = A;
  Ainv_ = Ainv;
  cori_ = cori;
  grav_ = grav;
  b_updatesetting_ = true;

  (void)extra_setting;
}

void DynWbc::_SetOptimizationSize() {
  // Dimension
  _dim_rf = 0;
  _dim_Uf = 0;
  for (int i = 0; i < (*_contact_list).size(); ++i) {
    _dim_rf += (*_contact_list)[i]->getDim();
    _dim_Uf += (*_contact_list)[i]->getDimRFConstraint();
  }

  _dim_opt = _dim_floating + _dim_rf; // dim of optimal output
  _dim_eq_cstr = _dim_floating;       // dim of equal contriant

  _Hessian = DMat::Zero(_dim_opt, _dim_opt);
  _gradient = DVec::Zero(_dim_opt);

  // Eigen Matrix Setting
  _dyn_CE = DMat::Zero(_dim_eq_cstr, _dim_opt);
  _dyn_ce0 = DVec(_dim_eq_cstr);
  if (_dim_rf > 0) {
    _dyn_CI = DMat::Zero(_dim_Uf, _dim_opt);
    _dyn_ci0 = DVec(_dim_Uf);

    _Jc = DMat(_dim_rf, num_qdot_);
    _JcDotQdot = DVec(_dim_rf);
    _Fr_des = DVec(_dim_rf);

    _Uf = DMat::Zero(_dim_Uf, _dim_rf);
    _Uf_ieq_vec = DVec(_dim_Uf);
  } else {
    _dyn_CI = DMat::Zero(1, _dim_opt);
    _dyn_ci0 = DVec::Zero(1); // may be error
  }
}

void DynWbc::_SetCost() {
  // Set Cost
  int idx_offset = 0;
  for (int i = 0; i < _dim_floating; ++i) {
    _Hessian(i + idx_offset, i + idx_offset) = _data->_W_floating[i];
  }
  idx_offset += _dim_floating;
  for (int i = 0; i < _dim_rf; ++i) {
    _Hessian(i + idx_offset, i + idx_offset) = _data->_W_rf[i];
  }
}

void DynWbc::_ContactBuilding() {
  DMat Uf;
  DVec Uf_ieq_vec;
  // Initial
  DMat Jc;
  DVec JcDotQdot;
  int dim_accumul_rf, dim_accumul_uf;
  (*_contact_list)[0]->getContactJacobian(Jc);
  (*_contact_list)[0]->getJcDotQdot(JcDotQdot);
  (*_contact_list)[0]->getRFConstraintMtx(Uf);
  (*_contact_list)[0]->getRFConstraintVec(Uf_ieq_vec);

  dim_accumul_rf =
      (*_contact_list)[0]
          ->getDim(); // acculmulate reaction force dim for for loop
  dim_accumul_uf =
      (*_contact_list)[0]
          ->getDimRFConstraint(); // acculmulate reaction force constraint
                                  // matrix dim for for loop

  _Jc.block(0, 0, dim_accumul_rf, num_qdot_) = Jc;
  _JcDotQdot.head(dim_accumul_rf) = JcDotQdot;
  _Uf.block(0, 0, dim_accumul_uf, dim_accumul_rf) = Uf;
  _Uf_ieq_vec.head(dim_accumul_uf) = Uf_ieq_vec;
  _Fr_des.head(dim_accumul_rf) = (*_contact_list)[0]->getRFDesired();

  int dim_new_rf, dim_new_uf;

  for (int i = 1; i < (*_contact_list).size(); ++i) {
    (*_contact_list)[i]->getContactJacobian(Jc);
    (*_contact_list)[i]->getJcDotQdot(JcDotQdot);

    dim_new_rf = (*_contact_list)[i]->getDim();
    dim_new_uf = (*_contact_list)[i]->getDimRFConstraint();

    // Jc append
    _Jc.block(dim_accumul_rf, 0, dim_new_rf, num_qdot_) = Jc;

    // JcDotQdot append
    _JcDotQdot.segment(dim_accumul_rf, dim_new_rf) = JcDotQdot;

    // Uf
    (*_contact_list)[i]->getRFConstraintMtx(Uf);
    _Uf.block(dim_accumul_uf, dim_accumul_rf, dim_new_uf, dim_new_rf) = Uf;

    // Uf inequality vector
    (*_contact_list)[i]->getRFConstraintVec(Uf_ieq_vec);
    _Uf_ieq_vec.segment(dim_accumul_uf, dim_new_uf) = Uf_ieq_vec;

    // Fr desired
    _Fr_des.segment(dim_accumul_rf, dim_new_rf) =
        (*_contact_list)[i]->getRFDesired();
    dim_accumul_rf += dim_new_rf;
    dim_accumul_uf += dim_new_uf;
  }
}

void DynWbc::_SetInEqualityConstraint() {
  _dyn_CI.block(0, _dim_floating, _dim_Uf, _dim_rf) = _Uf;
  _dyn_ci0 = _Uf_ieq_vec - _Uf * _Fr_des;
}

void DynWbc::_SetEqualityConstraint(const DVec &qddot) {
  if (_dim_rf > 0) {
    _dyn_CE.block(0, 0, _dim_eq_cstr, _dim_floating) =
        A_.block(0, 0, _dim_floating, _dim_floating);
    _dyn_CE.block(0, _dim_floating, _dim_eq_cstr, _dim_rf) =
        -Sv_ * _Jc.transpose();
    _dyn_ce0 = -Sv_ * (A_ * qddot + cori_ + grav_ - _Jc.transpose() * _Fr_des);

    // std::cout << "[FrMPC]"<<_Fr_des.transpose() <<std::endl;
    // std::cout << "[-Sv*JcT*_Fr_des]"<< (-Sv_*_Jc.transpose() *
    // _Fr_des).transpose() << std::endl;
  } else {
    _dyn_CE.block(0, 0, _dim_eq_cstr, _dim_floating) =
        A_.block(0, 0, _dim_floating, _dim_floating);
    _dyn_ce0 = -Sv_ * (A_ * qddot + cori_ + grav_);
  }
}

void DynWbc::_GetSolution(const DVec &qddot, DVec &cmd) {
  DVec tot_tau;
  if (_dim_rf > 0) {
    _data->_Fr = DVec(_dim_rf);
    // get Reaction forces
    for (int i = 0; i < _dim_rf; ++i)
      _data->_Fr[i] = _z[i + _dim_floating] + _Fr_des[i];

    // DVec qddrand = qddot;
    // qddrand(0) += 10;
    // qddrand(1) += 10;
    // qddrand(2) += 10;
    // qddrand(3) += 10;
    // qddrand(4) += 10;
    // qddrand(5) += 10;
    // DVec temp = Sv_ * (A_ * qddot + cori_ + grav_ );
    // std::cout << temp.transpose() << std::endl;

    // std::cout <<"[qddot]"<< qddot.transpose().segment(0 ,6) << std::endl;
    // std::cout << "[fr]" << _data->_Fr.transpose() << std::endl;
    // std::cout << "[frdes]" << _Fr_des.transpose() << std::endl;
    tot_tau = A_ * qddot + cori_ + grav_ - _Jc.transpose() * _data->_Fr;
    // A_ * qddot + cori_ + grav_ - _Jc.transpose() * _Fr_des;

  } else {
    tot_tau = A_ * qddot + cori_ + grav_;
  }
  _data->_qddot = qddot;
  cmd = tot_tau.tail(num_act_joint_);
}

void DynWbc::MakeTorque(DVec &cmd, void *extra_input) {
  if (!b_updatesetting_) {
    printf("[Wanning] WBIC setting is not done\n");
  }
  if (extra_input)
    _data = static_cast<DynWbcExtraData *>(extra_input); //_dyn_wbc_data

  // resize G, g0, CE, ce0, CI, ci0
  _SetOptimizationSize();
  _SetCost();

  DVec qddot_pre;
  DMat JcBar;
  DMat Npre;

  if (_dim_rf > 0) {
    // Contact Setting
    _ContactBuilding();

    // Set inequality constraints
    _SetInEqualityConstraint();
    _WeightedInverse(_Jc, Ainv_, JcBar);
    qddot_pre = JcBar * (-_JcDotQdot);
    // std::cout <<"[qddot_pre]"<< qddot_pre.transpose() <<std::endl;
    // std::cout <<"[_JcDotQdot]"<< _JcDotQdot.transpose() <<std::endl;
    Npre = _eye - JcBar * _Jc;
    // std::cout <<"[JcBar]"<< JcBar <<std::endl;
    // pretty_print(JcBar, std::cout, "JcBar");
    // pretty_print(_JcDotQdot, std::cout, "JcDotQdot");
    // pretty_print(qddot_pre, std::cout, "qddot 1");
  } else {
    qddot_pre = DVec::Zero(num_qdot_);
    Npre = _eye;
  }

  // Task
  Task *task;
  DMat Jt, JtBar, JtPre;
  DVec JtDotQdot, xddot;

  for (int i = 0; i < (*_task_list).size(); ++i) {
    // std::cout << i << std ::endl;
    task = (*_task_list)[i];

    task->getTaskJacobian(Jt);
    task->getTaskJacobianDotQdot(JtDotQdot);
    task->getCommand(xddot);
    // if (i == 0)
    // std::cout << "[xddotori]"<<xddot.transpose()<<std::endl;
    // if (i == 1) {
    //     xddot(0) = 120;
    //     xddot(1) = 10;
    //     xddot(2) = 100;
    // }

    // if (i == 1)
    // std::cout << "[xddotpos]"<<xddot.transpose()<<std::endl;

    // std::cout << task -> getPosError().transpose()<<std::endl;

    JtPre = Jt * Npre;
    _WeightedInverse(JtPre, Ainv_, JtBar);

    qddot_pre += JtBar * (xddot - JtDotQdot - Jt * qddot_pre);
    Npre = Npre * (_eye - JtBar * JtPre);
  }

  // std::cout <<"[qddcmd]"<< qddot_pre.segment(6,12).transpose() << std::endl;
  // std::cout <<"[qddcmd]"<< qddot_pre.segment(0,6).transpose() << std::endl;
  // Set equality constraints
  _SetEqualityConstraint(qddot_pre);

  // printf("G:\n");
  // std::cout<<G<<std::endl;
  // printf("g0:\n");
  // std::cout<<g0<<std::endl;

  // Optimization
  // Timer timer;
  // T f = solve_quadprog(G, g0, CE, ce0, CI, ci0, z);
  // std::cout<<"\n wbic old time: "<<timer.getMs()<<std::endl;
  // (void)f;
  _SolveQP();
  // pretty_print(qddot_pre, std::cout, "qddot_cmd");
  for (int i = 0; i < _dim_floating; ++i)
    qddot_pre[i] += _z[i];
  _GetSolution(qddot_pre, cmd);
  // std::cout << "cmd" << cmd.transpose()  << std::endl;

  _data->_opt_result = DVec(_dim_opt);
  for (int i = 0; i < _dim_opt; ++i) {
    _data->_opt_result[i] = _z[i];
  }
}

void DynWbc::_SolveQP() {
  /*
  QP problem:
  min 1 / 2 * _z.T * _Hessian * _z + _gradient.transpos() * _z
  s.t.
  _dyn_CE * z = _dyn_ce0
  _dyn_CI * z >= _dyn_ci0
  */

  DMat linear_constraints(_dyn_CE.rows() + _dyn_CI.rows(), _dyn_CE.cols());
  linear_constraints << _dyn_CE, _dyn_CI;

  DVec lower_bound(_dyn_ce0.rows() + _dyn_ci0.rows());
  DVec upper_bound(_dyn_ce0.rows() + _dyn_ci0.rows());
  lower_bound.segment(0, _dyn_ce0.rows()) = _dyn_ce0;
  upper_bound.segment(0, _dyn_ce0.rows()) = _dyn_ce0;

  lower_bound.segment(_dyn_ce0.rows(), _dyn_ci0.rows()) = _dyn_ci0;
  upper_bound.segment(_dyn_ce0.rows(), _dyn_ci0.rows()) =
      Eigen::MatrixXd::Constant(_dyn_ci0.rows(), 1, OsqpEigen::INFTY);
  Eigen::SparseMatrix<double> linear_constraints_sparse;
  linear_constraints_sparse = linear_constraints.sparseView();
  Eigen::SparseMatrix<double> Hessian_sparse;
  Hessian_sparse = _Hessian.sparseView();

  if (!qp_solver.isInitialized()) {
    qp_solver.settings()->setVerbosity(false);
    qp_solver.settings()->setWarmStart(false);
    qp_solver.data()->setNumberOfVariables(_dyn_CE.cols());
    qp_solver.data()->setNumberOfConstraints(_dyn_CE.rows() + _dyn_CI.rows());
    qp_solver.data()->setLinearConstraintsMatrix(linear_constraints_sparse);
    qp_solver.data()->setHessianMatrix(Hessian_sparse);
    qp_solver.data()->setGradient(_gradient);
    qp_solver.data()->setLowerBound(lower_bound);
    qp_solver.data()->setUpperBound(upper_bound);
    qp_solver.initSolver();
  } else {
    // qp_solver.clearSolver();
    qp_solver.settings()->setVerbosity(false);
    qp_solver.settings()->setWarmStart(true);

    qp_solver.data()->setNumberOfVariables(_dyn_CE.cols());
    qp_solver.data()->setNumberOfConstraints(_dyn_CE.rows() + _dyn_CI.rows());

    qp_solver.data()->clearHessianMatrix();
    qp_solver.data()->clearLinearConstraintsMatrix();

    qp_solver.data()->setLinearConstraintsMatrix(linear_constraints_sparse);
    qp_solver.data()->setHessianMatrix(Hessian_sparse);

    qp_solver.data()->setLowerBound(lower_bound);
    qp_solver.data()->setUpperBound(upper_bound);

    // qp_solver.initSolver();
  }
  qp_solver.solveProblem();

  // Vec18 temp1;
  // temp1.setZero();
  // temp1(6) = 1;
  // DVec temp = Sv_ * (A_ * qdd + cori_ + grav_);
  // std::cout << temp.transpose() << std::endl;

  // std::cout << "[QP]" <<_z.transpose() * Hessian_sparse * _z << std::endl;
  // std::cout << "[lower_bound]"<<lower_bound.transpose() << std::endl;
  // std::cout << "[middle]"<< (linear_constraints_sparse * _z).transpose() <<
  // std::endl; std::cout << "[upper_bound]"<<upper_bound.transpose() <<
  // std::endl;

  // std::cout << upper_bound.transpose() <<std::endl;
  _z = qp_solver.getSolution();

  if (!qp_solver.solve()) {
    std::cout << "[DynWbc] fail to solve qp!" << std::endl;
  }
  // std::cout << _z.transpose() << std::endl;
}