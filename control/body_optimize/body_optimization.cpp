
#include "body_optimization.h"

QuadBodyOpti::QuadBodyOpti() {

  side_sign_x.push_back(1);
  side_sign_x.push_back(1);
  side_sign_x.push_back(-1);
  side_sign_x.push_back(-1);

  side_sign_y.push_back(-1);
  side_sign_y.push_back(1);
  side_sign_y.push_back(-1);
  side_sign_y.push_back(1);
}

QuadBodyOpti::~QuadBodyOpti() {}

// rpy转为旋转矩阵
casadi::SX QuadBodyOpti::rpy2rot(casadi::SX rpy) {
  casadi::SX Rx = casadi::SX(3, 3);
  casadi::SX Ry = casadi::SX(3, 3);
  casadi::SX Rz = casadi::SX(3, 3);
  casadi::SX R = casadi::SX(3, 3);

  Rz(0, 0) = cos(rpy(2));
  Rz(0, 1) = -sin(rpy(2));
  Rz(1, 0) = sin(rpy(2));
  Rz(1, 1) = cos(rpy(2));
  Rz(2, 2) = 1;

  Ry(0, 0) = cos(rpy(1));
  Ry(0, 2) = sin(rpy(1));
  Ry(1, 1) = 1;
  Ry(2, 0) = -sin(rpy(1));
  Ry(2, 2) = cos(rpy(1));

  Rx(0, 0) = 1;
  Rx(1, 1) = cos(rpy(0));
  Rx(1, 2) = -sin(rpy(0));
  Rx(2, 1) = sin(rpy(0));
  Rx(2, 2) = cos(rpy(0));

  R = casadi::SX::mtimes({Rz, Ry, Rx});

  return R;
}

//创建求解器
void QuadBodyOpti::setBodyOptiWeight(double weightX, double weightY,
                                     double weightZ) {

  // 读取参数
  casadi::SX weight = casadi::SX(3, 3);
  weight(0, 0) = weightX;
  weight(1, 1) = weightY;
  weight(2, 2) = weightZ;

  // Variables，状态和控制输入为变量
  casadi::SX obj = casadi::SX::sym("obj", 5);
  // Parameters，含有足端位置,yaw
  casadi::SX p = casadi::SX::sym("p", 12 + 1);

  casadi::SX rotation_mat =
      rpy2rot(casadi::SX::vertcat({obj(casadi::Slice(0, 2)), p(12)}));

  casadi::SX Foot_default = casadi::SX({dx, dy, dz});

  // 代价函数
  casadi::SX f = casadi::SX::sym("cost_fun");
  f = 0;
  for (int i = 0; i < 4; ++i) {
    f += casadi::SX::mtimes(
        {(Foot_default - p(casadi::Slice(i * 3, i * 3 + 3)) +
          obj(casadi::Slice(2, 5)) +
          casadi::SX::mtimes(
              {rotation_mat, casadi::SX({side_sign_x[i] * length,
                                         side_sign_y[i] * width, 0})}))
             .T(),
         weight,
         (Foot_default - p(casadi::Slice(i * 3, i * 3 + 3)) +
          obj(casadi::Slice(2, 5)) +
          casadi::SX::mtimes(
              {rotation_mat, casadi::SX({side_sign_x[i] * length,
                                         side_sign_y[i] * width, 0})}))});
  }

  //构建求解器
  casadi::SXDict nlp_prob = {{"x", obj}, {"p", p}, {"f", f}};
  std::string solver_name = "ipopt";
  casadi::Dict nlp_opts;
  nlp_opts["expand"] = true;
  nlp_opts["ipopt.max_iter"] = 100;
  nlp_opts["ipopt.linear_solver"] = "ma27";
  nlp_opts["ipopt.print_level"] = 0;
  nlp_opts["print_time"] = 0;
  nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-4;
  nlp_opts["ipopt.acceptable_tol"] = 1e-4;
  nlp_opts["ipopt.warm_start_entire_iterate"] = "yes";
  solver = nlpsol("nlpsol", solver_name, nlp_prob, nlp_opts);

  // // Generate C code for the NLP functions
  // solver.generate_dependencies("nlp");
  // // Just-in-time compilation?
  // bool jit = false;
  // if (jit) {
  //   // Create a new NLP solver instance using just-in-time compilation
  //   solver = casadi::nlpsol("nlpsol", solver_name, "nlp.c");
  // } else {
  //   // Compile the c-code
  //   int flag = system("gcc -fPIC -shared -O3 nlp.c -o nlp.so");
  //   casadi_assert(flag == 0, "Compilation failed");

  //   // Create a new NLP solver instance from the compiled code
  //   solver = casadi::nlpsol("nlpsol", solver_name, "nlp.so");
  // }

  // casadi::Dict nlp_opts;
  // nlp_opts["expand"] = true;
  // // nlp_opts["max_iter"] = 10)
  // // nlp_opts["verbose"] = true;
  // // nlp_opts["linear_solver"] = "ma57";
  // nlp_opts["hessian_approximation"] = "exact";
  // // nlp_opts["derivative_test"] = "second-order";
  // // Specify QP solver
  // nlp_opts["qpsol"] = "nlpsol";
  // nlp_opts["qpsol_options.nlpsol"] = "ipopt";
  // nlp_opts["qpsol_options.error_on_fail"] = false;
  // nlp_opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0;
  // nlp_opts["qpsol_options.nlpsol_options.print_time"] = 0;
  // // Allocate NLP solver and buffers
  // solver = nlpsol("nlpsol", "sqpmethod", nlp_prob, nlp_opts);
}

void QuadBodyOpti::optiSolution(Eigen::Matrix<double, 6, 1> current_states,
                                Eigen::Matrix<double, 12, 1> p_feet) {

  //求解参数设置
  // Initial guess and bounds for the optimization variablese
  std::vector<double> x0;
  x0.push_back(current_states(0));
  x0.push_back(current_states(1));
  for (int i = 0; i < 3; ++i) {
    x0.push_back(current_states(3 + i));
  }

  // Nonlinear bounds
  std::vector<double> lbx;
  std::vector<double> ubx;
  double pFoot_x_min =
      fmin(fmin(fmin(fmin(fmin(p_feet(0) + p_feet(3), p_feet(0) + p_feet(6)),
                          p_feet(0) + p_feet(9)),
                     p_feet(3) + p_feet(6)),
                p_feet(3) + p_feet(9)),
           p_feet(6) + p_feet(9));
  double pFoot_x_max =
      fmax(fmax(fmax(fmax(fmax(p_feet(0) + p_feet(3), p_feet(0) + p_feet(6)),
                          p_feet(0) + p_feet(9)),
                     p_feet(3) + p_feet(6)),
                p_feet(3) + p_feet(9)),
           p_feet(6) + p_feet(9));
  double pFoot_y_min =
      fmin(fmin(fmin(fmin(fmin(p_feet(1) + p_feet(4), p_feet(1) + p_feet(7)),
                          p_feet(1) + p_feet(10)),
                     p_feet(4) + p_feet(7)),
                p_feet(4) + p_feet(10)),
           p_feet(7) + p_feet(10));
  double pFoot_y_max =
      fmax(fmax(fmax(fmax(fmax(p_feet(1) + p_feet(4), p_feet(1) + p_feet(7)),
                          p_feet(1) + p_feet(10)),
                     p_feet(4) + p_feet(7)),
                p_feet(4) + p_feet(10)),
           p_feet(7) + p_feet(10));

  lbx.push_back(-0.8);
  lbx.push_back(-0.8);
  lbx.push_back(pFoot_x_min);
  lbx.push_back(pFoot_y_min);
  lbx.push_back(0.2);
  ubx.push_back(0.8);
  ubx.push_back(0.8);
  ubx.push_back(pFoot_x_max);
  ubx.push_back(pFoot_y_max);
  ubx.push_back(0.36);

  // Original parameter values
  std::vector<double> p0;
  for (int leg = 0; leg < 4; leg++) {
    for (int i = 0; i < 3; i++) {
      p0.push_back(p_feet(i + 3 * leg));
    }
  }
  p0.push_back(current_states(2));

  arg["lbx"] = lbx;
  arg["ubx"] = ubx;
  arg["x0"] = x0;
  arg["p"] = p0;
  res = solver(arg);

  // std::cout << "Objective: " << res.at("f") << std::endl;

  std::vector<double> res_all(res.at("x"));
  for (int i = 0; i < 5; i++) {
    opt_obj(i) = res_all.at(i);
  }
}
