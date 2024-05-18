
#include "nmpc_caller.h"

#include <yaml-cpp/yaml.h>

QuadNmpc::QuadNmpc() {
  // create a capsule according to the pre-defined model
  acados_ocp_capsule = srbd_acados_create_capsule();

  // optimizer
  status = srbd_acados_create(acados_ocp_capsule);
  if (status) {
    printf("srbd_acados_create() optimater returned status %d. Exiting.\n",
           status);
    exit(1);
  }

  // some important structure of ocp
  nlp_config = srbd_acados_get_nlp_config(acados_ocp_capsule);
  nlp_dims = srbd_acados_get_nlp_dims(acados_ocp_capsule);
  nlp_in = srbd_acados_get_nlp_in(acados_ocp_capsule);
  nlp_out = srbd_acados_get_nlp_out(acados_ocp_capsule);

  horizons = nlp_dims->N;
  nx = *nlp_dims->nx;
  nu = *nlp_dims->nu;

  lastOptimalState.setZero(24 * (horizons + 1));
  lastOptimalControl.setZero(24 * horizons);

  // params
  YAML::Node config = YAML::LoadFile("../control/nmpc_config.yaml");

  mass = config["mass"].as<double>();
  g = config["g"].as<double>();
  Xshoulder << 0.183, 0.183, -0.183, -0.183;
  Yshoulder << -0.11205, 0.11205, -0.11205, 0.11205;

  useTimeStepOpti = config["useTimeStepOpti"].as<bool>();

  // com task
  double parameters[70];
  int weightParamNum = 22;
  parameters[0] = config["q1"].as<double>();
  parameters[1] = config["q2"].as<double>();
  parameters[2] = config["q3"].as<double>();
  parameters[3] = config["q4"].as<double>();
  parameters[4] = config["q5"].as<double>();
  parameters[5] = config["q6"].as<double>();
  parameters[6] = config["q7"].as<double>();
  parameters[7] = config["q8"].as<double>();
  parameters[8] = config["q9"].as<double>();
  parameters[9] = config["q10"].as<double>();
  parameters[10] = config["q11"].as<double>();
  parameters[11] = config["q12"].as<double>();
  // foot pos task
  parameters[12] = config["q13"].as<double>();
  parameters[13] = config["q14"].as<double>();
  parameters[14] = config["q15"].as<double>();
  // regularize grf
  parameters[15] = config["r1"].as<double>();
  parameters[16] = config["r2"].as<double>();
  parameters[17] = config["r3"].as<double>();
  // foot vel task
  parameters[18] = config["r4"].as<double>();
  parameters[19] = config["r5"].as<double>();
  parameters[20] = config["r6"].as<double>();
  // coefficient of terminal node
  parameters[21] = config["alpha"].as<double>();
  for (int i = 0; i < horizons + 1; ++i) {
    srbd_acados_update_params(acados_ocp_capsule, i, parameters, 70);
  }

  // using this to only change desired
  for (int i = 0; i < 48; ++i) {
    referenceIndex[i] = weightParamNum + i;
  }
}

QuadNmpc::~QuadNmpc() {}

void QuadNmpc::copyArray2Eigen(VectorXd &target, const double *source, int len,
                               int startIndexEigen, int startIndexArray) {
  for (int i = 0; i < len; i++) {
    target(i + startIndexEigen) = source[i + startIndexArray];
  }
}

void QuadNmpc::copyEigen2Array(double *target, const VectorXd &source, int len,
                               int startIndexArray, int startIndexEigen) {
  for (int i = 0; i < len; i++) {
    target[i + startIndexArray] = source(i + startIndexEigen);
  }
}

// default foot pos
Vector2d QuadNmpc::getFootProjectReference(const double yaw, Vector2d pos) {
  Matrix2d Rz = Matrix2d::Zero();
  Rz(0, 0) = cos(yaw);
  Rz(0, 1) = -sin(yaw);
  Rz(1, 0) = sin(yaw);
  Rz(1, 1) = cos(yaw);

  return Rz * pos;
}

// update for desired and current state
void QuadNmpc::nmpcUpdate(const VectorXd &currentState,
                          const VectorXd &gaitTable,
                          const VectorXd &desiredComState,
                          const VectorXd &desiredFootZPos,
                          const VectorXd &desiredFootZVel) {
  // std::cout << "currentState\n" << currentState.transpose() << std::endl;
  // std::cout << "desiredComState\n" << desiredComState.transpose() <<
  // std::endl; std::cout << "desiredFootZPos\n" << desiredFootZPos.transpose()
  // << std::endl; std::cout << "desiredFootZVel\n" <<
  // desiredFootZVel.transpose()
  // << std::endl;

  // --- first node should be current state ---
  double currentStateArray[24];
  copyEigen2Array(currentStateArray, currentState, 24, 0, 0);
  ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx",
                                currentStateArray);
  ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx",
                                currentStateArray);

  // ---constraints ---
  for (int i = 0; i < horizons + 1; i++) {
    if (i == 0) { // constraint_0
      double lh[29], uh[29];
      for (int j = 0; j < 4; j++) {
        // friction cone
        lh[7 * j + 0] = -10000;
        lh[7 * j + 1] = -10000;
        lh[7 * j + 2] = 0;
        lh[7 * j + 3] = 0;

        uh[7 * j + 0] = 0;
        uh[7 * j + 1] = 0;
        uh[7 * j + 2] = 10000;
        uh[7 * j + 3] = 10000;

        // z-axis force limit
        if (gaitTable[i * 4 + j] == 0) {
          lh[7 * j + 4] = 0;
          uh[7 * j + 4] = 0;
        } else {
          lh[7 * j + 4] = 0;
          uh[7 * j + 4] = 200;
        }
        // complementary
        lh[7 * j + 5] = -10000;
        uh[7 * j + 5] = 0;
        lh[7 * j + 6] = 0;
        uh[7 * j + 6] = 10000;
      }
      // slack variable
      lh[28] = 0;
      uh[28] = 10000;

      ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
      ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    } else if (i < horizons + 1) { // constraint
      double lh[33], uh[33];
      for (int j = 0; j < 4; j++) {
        // z-axis pos of foot
        lh[8 * j + 0] = 0.005;
        uh[8 * j + 0] = 0.3;
        // friction cone
        lh[8 * j + 1] = -10000;
        lh[8 * j + 2] = -10000;
        lh[8 * j + 3] = 0;
        lh[8 * j + 4] = 0;

        uh[8 * j + 1] = 0;
        uh[8 * j + 2] = 0;
        uh[8 * j + 3] = 10000;
        uh[8 * j + 4] = 10000;

        // z-axis force limit
        if (gaitTable[i * 4 + j] == 0) {
          lh[8 * j + 5] = 0;
          uh[8 * j + 5] = 0;
        } else {
          lh[8 * j + 5] = 0;
          uh[8 * j + 5] = 200;
        }
        // complementary
        lh[8 * j + 6] = -10000;
        uh[8 * j + 6] = 0;
        lh[8 * j + 7] = 0;
        uh[8 * j + 7] = 10000;
      }
      // slack varible
      lh[32] = 0;
      uh[32] = 10000;

      ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
      ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
  }

  // ---parameter ---
  Vector<double, 8> footProjectReference;
  for (int i = 0; i < 4; ++i) {
    footProjectReference.segment(2 * i, 2) = getFootProjectReference(
        currentState[8], Vector2d(Xshoulder[i], Yshoulder[i]));
  }
  for (int j = 0; j < horizons + 1; j++) {
    // state ref
    for (int k = 0; k < 12; k++) {
      reference[k] = desiredComState[k];
    }
    for (int k = 0; k < 4; k++) {
      reference[12 + 3 * k + 0] = footProjectReference[2 * k + 0];
      reference[12 + 3 * k + 1] = footProjectReference[2 * k + 1];
      reference[12 + 3 * k + 2] = desiredFootZPos[4 * j + k];
    }
    // control ref
    double grfAverage = mass * g /
                        (gaitTable[0 + 4 * j] + gaitTable[1 + 4 * j] +
                         gaitTable[2 + 4 * j] + gaitTable[3 + 4 * j]);
    for (int k = 0; k < 4; k++) {
      reference[24 + 3 * k + 0] = 0;
      reference[24 + 3 * k + 1] = 0;
      // reference[24 + 3 * k + 2] = gaitTable[k + 4 * j] * grfAverage;
      reference[24 + 3 * k + 2] = 0;
    }
    for (int k = 0; k < 4; ++k) {
      reference[36 + 3 * k + 0] = 0;
      reference[36 + 3 * k + 1] = 0;
      reference[36 + 3 * k + 2] = desiredFootZVel[4 * j + k];
    }

    srbd_acados_update_params_sparse(acados_ocp_capsule, j, referenceIndex,
                                     reference, 48);
  }
}

// solve and get optimal trajectory
void QuadNmpc::getSolution(VectorXd &optimalState, VectorXd &optimalControl) {
  status = srbd_acados_solve(acados_ocp_capsule);

  // return last solution if solver failed
  if (status != ACADOS_SUCCESS) {
    // printf("acados_solve() failed with status %d.\n", status);
    optimalState = lastOptimalState;
    optimalControl = lastOptimalControl;
  } else {
    // get the optimized solution
    double optimalStateArray[24];
    double optimalControlArray[25];

    for (int i = 0; i < horizons; i++) {
      ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, i, "x", optimalStateArray);
      ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, i, "u",
                      optimalControlArray);

      copyArray2Eigen(lastOptimalState, optimalStateArray, 24, i * 24, 0);
      copyArray2Eigen(lastOptimalControl, optimalControlArray, 25, i * 25, 0);
    }
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, horizons, "x",
                    optimalStateArray);
    copyArray2Eigen(lastOptimalState, optimalStateArray, 24, horizons * 24, 0);

    optimalState = lastOptimalState;
    optimalControl = lastOptimalControl;
  }
  // std::cout << "control = \n"
  //           << optimalControl.head(25).transpose() << std::endl;
  // std::cout << "state = \n"
  //           << optimalState.segment(12, 12).transpose() << std::endl;
}