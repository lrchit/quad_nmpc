#include "KinWbc.hpp"
#include "pseudoInverse.h"

KinWBC::KinWBC(int num_qdot)
    : threshold_(0.001), num_qdot_(num_qdot), num_act_joint_(num_qdot - 6) {
  I_mtx = DMat::Identity(num_qdot_, num_qdot_);
}

/*curr_config is current joint position*/
bool KinWBC::FindConfiguration(const DVec &curr_config,
                               const std::vector<Task *> &task_list,
                               const std::vector<ContactSpec *> &contact_list,
                               DVec &jpos_cmd, DVec &jvel_cmd) {

  // Contact Jacobian Setup
  DMat Nc(num_qdot_, num_qdot_);
  Nc.setIdentity();
  if (contact_list.size() > 0) {
    DMat Jc, Jc_i;
    contact_list[0]->getContactJacobian(Jc);
    int num_rows = Jc.rows();

    // stack Jc of every contact task into a complete Jc
    for (int i = 1; i < contact_list.size(); ++i) {
      contact_list[i]->getContactJacobian(Jc_i);
      int num_new_rows = Jc_i.rows();
      Jc.conservativeResize(num_rows + num_new_rows, num_qdot_);
      Jc.block(num_rows, 0, num_new_rows, num_qdot_) = Jc_i;
      num_rows += num_new_rows;
    }

    // Projection Matrix
    _BuildProjectionMatrix(Jc, Nc);
  }

  // First Task (not contact force task)
  DVec delta_q,
      qdot; // after contact force task delta_q = qdotcmd = 0,qddotcmd != 0
  DMat Jt, JtPre, JtPre_pinv, N_nx, N_pre; // N_nx: Ni|i - 1|N_pre: Ni-1

  Task *task = task_list[0];
  task->getTaskJacobian(Jt);
  JtPre = Jt * Nc;
  _PseudoInverse(JtPre, JtPre_pinv);

  delta_q = JtPre_pinv * (task->getPosError());
  qdot = JtPre_pinv * (task->getDesVel());

  DVec prev_delta_q = delta_q;
  DVec prev_qdot = qdot;

  _BuildProjectionMatrix(JtPre, N_nx);
  N_pre = Nc * N_nx;

  for (int i = 1; i < task_list.size(); ++i) {
    task = task_list[i];

    task->getTaskJacobian(Jt);
    JtPre = Jt * N_pre;

    _PseudoInverse(JtPre, JtPre_pinv);
    delta_q =
        prev_delta_q + JtPre_pinv * (task->getPosError() - Jt * prev_delta_q);
    qdot = prev_qdot + JtPre_pinv * (task->getDesVel() - Jt * prev_qdot);

    // For the next task
    _BuildProjectionMatrix(JtPre, N_nx);
    N_pre *= N_nx;
    prev_delta_q = delta_q;
    prev_qdot = qdot;
  }
  for (int i = 0; i < num_act_joint_; ++i) {
    jpos_cmd[i] = curr_config[i + 6] + delta_q[i + 6];
    jvel_cmd[i] = qdot[i + 6];
  }
  return true;
}

void KinWBC::_BuildProjectionMatrix(const DMat &J, DMat &N) {
  DMat J_pinv;
  _PseudoInverse(J, J_pinv);
  N = I_mtx - J_pinv * J;
}

void KinWBC::_PseudoInverse(const DMat J, DMat &Jinv) {
  pseudoInverse(J, threshold_, Jinv);
}
