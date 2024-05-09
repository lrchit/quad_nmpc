
#pragma once

#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

class QuadBodyOpti {
public:
  QuadBodyOpti();
  ~QuadBodyOpti();
  //定义求解器
  void setBodyOptiWeight(double weightX, double weightY, double weightZ);
  //优化求解
  void optiSolution(Eigen::Matrix<double, 6, 1> current_states,
                    Eigen::Matrix<double, 12, 1> p_feet);
  //获取解
  Eigen::Matrix<double, 5, 1> get_opt_obj() { return opt_obj; }

private:
  double dx = 0;
  double dy = 0;
  double dz = -0.28;
  double length = 0.183;
  double width = 0.13205;
  std::vector<double> side_sign_x;
  std::vector<double> side_sign_y;

  Eigen::Matrix<double, 5, 1> opt_obj;

  casadi::Function solver;               //求解器
  std::map<std::string, casadi::DM> res; //求解结果
  std::map<std::string, casadi::DM> arg; //求解参数

  casadi::SX rpy2rot(casadi::SX rpy);
};