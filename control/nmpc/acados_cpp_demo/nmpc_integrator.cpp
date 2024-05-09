
#include "nmpc_integrator.h"

void QuadNmpcInterpolator::copyArray2Eigen(VectorXd &target,
                                           const double *source, int len,
                                           int startIndexEigen,
                                           int startIndexArray) {
  for (int i = 0; i < len; i++) {
    target(i + startIndexEigen) = source[i + startIndexArray];
  }
}

void QuadNmpcInterpolator::copyEigen2Array(double *target,
                                           const VectorXd &source, int len,
                                           int startIndexArray,
                                           int startIndexEigen) {
  for (int i = 0; i < len; i++) {
    target[i + startIndexArray] = source(i + startIndexEigen);
  }
}

// 构造函数对求解器进行初始化，基本没啥可改的
QuadNmpcInterpolator::QuadNmpcInterpolator() {
  // create a capsule according to the pre-defined model
  srbd_sim_solver_capsule_ = srbd_acados_sim_solver_create_capsule();

  // optimizer
  status = srbd_acados_sim_create(srbd_sim_solver_capsule_);

  // 判断一下是否初始化成功，失败则立即退出
  if (status) {
    printf("srbd_acados_sim_create() optimater returned status %d. Exiting.\n",
           status);
    exit(1);
  }

  // some important structure of ocp
  sim_config_ = srbd_acados_get_sim_config(srbd_sim_solver_capsule_);
  sim_dims_ = srbd_acados_get_sim_dims(srbd_sim_solver_capsule_);
  sim_in_ = srbd_acados_get_sim_in(srbd_sim_solver_capsule_);
  sim_out_ = srbd_acados_get_sim_out(srbd_sim_solver_capsule_);
  sim_solver_ = srbd_acados_get_sim_solver(srbd_sim_solver_capsule_);

  horizons = 51;
  nx = 24;
  nu = 24;

  lastOptimalState.setZero(nx);
}

QuadNmpcInterpolator::~QuadNmpcInterpolator() {}

// 求解非线性优化问题并返回最优状态和输入轨迹
void QuadNmpcInterpolator::getSolution(const VectorXd &currentState,
                                       const VectorXd &currentInputCommand,
                                       VectorXd &nextStateCommand) {

  double current_state[nx], current_input_command[nu];
  copyEigen2Array(current_state, currentState, nx, 0, 0);
  copyEigen2Array(current_input_command, currentInputCommand, nu, 0, 0);

  sim_in_set(sim_config_, sim_dims_, sim_in_, "x", current_state);
  sim_in_set(sim_config_, sim_dims_, sim_in_, "u", current_input_command);

  status = srbd_acados_sim_solve(srbd_sim_solver_capsule_);

  // 判断一下是否有解，无解的时候要保留上次的结果，直至有解的时候才能更新
  if (status != ACADOS_SUCCESS) {
    // printf("acados_solve() failed with status %d.\n", status);
    nextStateCommand = lastOptimalState;
  } else {
    // get solution
    double optimalStateArray[nx];
    sim_out_get(sim_config_, sim_dims_, sim_out_, "x", optimalStateArray);
    copyArray2Eigen(lastOptimalState, optimalStateArray, nx, 0, 0);

    nextStateCommand = lastOptimalState;
  }
}