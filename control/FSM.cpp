#include "FSM.h"

void rpySafetyCheck(double &angle_des, double angle_act) {
  if (angle_des - angle_act > pi * 0.5) {
    angle_des -= 2 * pi;
  } else if (angle_des - angle_act < -pi * 0.5) {
    angle_des += 2 * pi;
  }
}

/**
 * @brief initialize
 */
QuadFSM::QuadFSM() {

  YAML::Node config = YAML::LoadFile("../control/nmpc_config.yaml");

  horizon = config["horizon"].as<int>();
  iterationBetweenMpc = config["iterationBetweenMpc"].as<int>();
  iterationBetweenBodyOpti = config["iterationBetweenBodyOpti"].as<int>();

  useWbc = config["useWbc"].as<bool>();

  weightX = config["body_opti_x"].as<double>();
  weightY = config["body_opti_y"].as<double>();
  weightZ = config["body_opti_z"].as<double>();

  // legKin.resize(4);
  legKin.push_back(LegKinematics("../unitree_a1_desc/urdf/FR.urdf"));
  legKin.push_back(LegKinematics("../unitree_a1_desc/urdf/FL.urdf"));
  legKin.push_back(LegKinematics("../unitree_a1_desc/urdf/RR.urdf"));
  legKin.push_back(LegKinematics("../unitree_a1_desc/urdf/RL.urdf"));

  kp_.setZero();
  kd_.setZero();
  kp_(0, 0) = config["kp_cartesian_x"].as<double>();
  kp_(1, 1) = config["kp_cartesian_y"].as<double>();
  kp_(2, 2) = config["kp_cartesian_z"].as<double>();
  kd_(0, 0) = config["kd_cartesian_x"].as<double>();
  kd_(1, 1) = config["kd_cartesian_y"].as<double>();
  kd_(2, 2) = config["kd_cartesian_z"].as<double>();

  dt = 0.002;
  dtmpc = dt * iterationBetweenMpc;

  footForceFromKin.setZero();

  trotting = new OffsetDurationGait(
      (horizon + 1) / 2,
      Vector<int, 4>(0, 0.5 * (horizon + 1) / 2, 0.5 * (horizon + 1) / 2, 0),
      Vector<int, 4>(0.5 * (horizon + 1) / 2, 0.5 * (horizon + 1) / 2,
                     0.5 * (horizon + 1) / 2, 0.5 * (horizon + 1) / 2),
      "Trotting");
  standing = new OffsetDurationGait(
      (horizon + 1) / 2, Vector<int, 4>(0, 0, 0, 0),
      Vector<int, 4>((horizon + 1) / 2, (horizon + 1) / 2, (horizon + 1) / 2,
                     (horizon + 1) / 2),
      "Standing");

  currentGait = standing;

  counter = 0;
  iterationCounter = 0;
  gaitIterationCounter = 0;
  nmpcUpdateNeeded = false;
  wbcUpdateNeeded = false;
  bodyOptiNeeded = false;

  estimater = new QuadEstm(legKin);
  wbcCaller = new QuadWbc();
  bodyOptimizer = new QuadBodyOpti();
  nmpcCaller = new QuadNmpc();
  nmpcIntegrator = new QuadNmpcInterpolator();

  bodyOptimizer->setBodyOptiWeight(weightX, weightY, weightZ);

  currentState.resize(24);
  desiredComState.resize(12);
  desiredFootZPos.resize(4 * (horizon + 1));
  desiredFootZVel.resize(4 * (horizon + 1));
  gaitTable.resize(4 * (horizon + 1));
  phaseTable.resize(4, horizon);
  optimalState.resize(24 * (horizon + 1));
  optimalControl.resize(24 * horizon);
  currentfootVel.resize(12);

  trajIntegrate.setZero();
}

/**
 * @brief main loop
 */
void QuadFSM::mainProgram() {

  counter++;

  if (counter >= 300 && counter < 2500000000) {
    iterationCounter++;
    // contact schedule
    currentGait->setIterations(iterationBetweenMpc, iterationCounter);
    if (currentGait == standing) {
      stateCur.contactPhase << 0.5, 0.5, 0.5, 0.5;
    } else {
      stateCur.contactPhase = currentGait->getContactState();
    }

    // call mpc if needed
    if ((iterationCounter % iterationBetweenMpc) == 0) {
      updateNmpc();
      nmpcUpdateNeeded = true;
    }
  } else {
    stateCur.contactPhase << 0.5, 0.5, 0.5, 0.5;
  }
}

void QuadFSM::stateEstimate(double *sensordata) {
  estimater->callStateEstimator(stateCur, sensordata);
}

// update nmpc
void QuadFSM::updateNmpc() {
  // --- gait switch --
  if (counter < 1500) {
    gait = standing;
  } else if (counter < 20000) {
    gait = trotting;
  } else {
    gait = standing;
  }
  gait->setIterations(iterationBetweenMpc, iterationCounter);
  Vector4d nextGaitPhase = gait->getSwingState();
  if ((nextGaitPhase[0] == 0 || nextGaitPhase[0] == 1) &&
      (nextGaitPhase[1] == 0 || nextGaitPhase[1] == 1) &&
      (nextGaitPhase[2] == 0 || nextGaitPhase[2] == 1) &&
      (nextGaitPhase[3] == 0 || nextGaitPhase[3] == 1)) {
    currentGait = gait;
  }

  // --- current state ---
  currentState.segment(0, 3) = stateCur.pos;
  currentState.segment(3, 3) = stateCur.linVel;
  currentState.segment(6, 3) = stateCur.rpy;
  currentState.segment(9, 3) = stateCur.omegaWorld;
  for (int i = 0; i < 4; ++i) {
    currentState.segment(12 + 3 * i, 3) = stateCur.footPosWorld.col(i);
  }

  // --- current foot vel ---
  for (int i = 0; i < 4; ++i) {
    currentfootVel.segment(3 * i, 3) = stateCur.footVelWorld.col(i);
  }

  // --- gait table --
  int *mpcTable;
  currentGait->setIterations(iterationBetweenMpc, iterationCounter);
  mpcTable = currentGait->getMpcTable();
  for (int i = 0; i < (horizon + 1) / 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      gaitTable[4 * i + j] = mpcTable[4 * i + j];
    }
  }
  currentGait->setIterations(iterationBetweenMpc,
                             iterationCounter +
                                 iterationBetweenMpc * (horizon + 1) / 2);
  mpcTable = currentGait->getMpcTable();
  for (int i = 0; i < (horizon + 1) / 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      gaitTable[4 * (horizon + 1) / 2 + 4 * i + j] = mpcTable[4 * i + j];
    }
  }

  // --- phase table ---
  for (int i = 0; i < horizon; ++i) {
    currentGait->setIterations(
        iterationBetweenMpc, iterationCounter + (i + 1) * iterationBetweenMpc);
    phaseTable.col(i) = currentGait->getSwingState();
  }

  // --- desired foot pos and vel ---
  double zPosDefault = 0.005;
  double swingHeight = 0.05;
  double zVelLift = 0.0;
  double zVelDown = 0.0;
  auto cubicSpine = [zPosDefault, swingHeight, zVelLift,
                     zVelDown](double phase) {
    double a, b, c, d, x0, x1, v0, v1, T, t;
    T = 1;
    if (phase < 0.5) {
      x0 = zPosDefault;
      x1 = zPosDefault + swingHeight;
      v0 = zVelLift;
      v1 = 0;
      t = 2 * phase;
    } else {
      x0 = zPosDefault + swingHeight;
      x1 = zPosDefault;
      v0 = 0;
      v1 = zVelDown;
      t = 2 * phase - 1;
    }

    a = (2 * x0 - 2 * x1 + v0 * T + v1 * T) / T / T / T;
    b = (-3 * x0 + 3 * x1 - 2 * v0 * T - v1 * T) / T / T;
    c = v0;
    d = x0;
    return a * t * t * t + b * t * t + c * t + d;
  };
  auto cubicSpineFirstDerivate = [zPosDefault, swingHeight, zVelLift,
                                  zVelDown](double phase) {
    double a, b, c, d, x0, x1, v0, v1, T, t;
    T = 1;
    if (phase < 0.5) {
      x0 = zPosDefault;
      x1 = zPosDefault + swingHeight;
      v0 = zVelLift;
      v1 = 0;
      t = 2 * phase;
    } else {
      x0 = zPosDefault + swingHeight;
      x1 = zPosDefault;
      v0 = 0;
      v1 = zVelDown;
      t = 2 * phase - 1;
    }

    a = (2 * x0 - 2 * x1 + v0 * T + v1 * T) / T / T / T;
    b = (-3 * x0 + 3 * x1 - 2 * v0 * T - v1 * T) / T / T;
    c = v0;
    return 3 * a * t * t + 2 * b * t + c;
  };
  auto cubicBezier = [zPosDefault, swingHeight, zVelLift,
                      zVelDown](double phase) {
    double ZPos;
    if (phase < 0.5) {
      phase = 2 * phase;
      ZPos = zPosDefault + pow(phase, 2) * (3 - 2 * phase) * swingHeight;
    } else {
      phase = 2 * phase - 1;
      ZPos = swingHeight + zPosDefault +
             pow(phase, 2) * (3 - 2 * phase) * (-swingHeight);
    }
    return ZPos;
  };
  auto cubicBezierFirstDerivate = [zPosDefault, swingHeight, zVelLift,
                                   zVelDown](double phase) {
    double ZVel;
    if (phase < 0.5) {
      phase = 2 * phase;
      ZVel = 6 * phase * (1 - phase) * swingHeight;
    } else {
      phase = 2 * phase - 1;
      ZVel = 6 * phase * (1 - phase) * (-swingHeight);
    }
    return ZVel;
  };

  for (int i = 0; i < horizon + 1; ++i) {
    // --- desired foot pos ---
    if (i == 0) {
      desiredFootZPos.segment(4 * i, 4) << stateCur.footPosWorld(2, 0),
          stateCur.footPosWorld(2, 1), stateCur.footPosWorld(2, 2),
          stateCur.footPosWorld(2, 3);
    } else {
      for (int j = 0; j < 4; ++j) {
        if (phaseTable(j, i - 1) > 0 && phaseTable(j, i - 1) <= 1) { // swing
          desiredFootZPos[4 * i + j] = cubicSpine(phaseTable(j, i - 1));
        } else { // stance
          desiredFootZPos[4 * i + j] = desiredFootZPos[4 * (i - 1) + j];
        }
      }
    }
    // --- desired foot vel ---
    for (int j = 0; j < 4; ++j) {
      if (phaseTable(j, i) > 0 && phaseTable(j, i) <= 1) { // swing
        desiredFootZVel[4 * i + j] = cubicSpineFirstDerivate(phaseTable(j, i));
      } else { // stance
        desiredFootZVel[4 * i + j] = 0;
      }
    }
  }

  // --- desired com state ---
  if (currentGait == trotting) { // trot
    desiredComState.segment(0, 12) << 0, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    // // yaw motion
    // desiredComState[8] = currentState[8];
    // double yawVel = 0.2;
    // desiredComState[8] += yawVel * horizon * iterationBetweenMpc * dt;
    // desiredComState[11] = yawVel;
    rpySafetyCheck(desiredComState[8], currentState[8]);
  } else if (currentGait == standing) { // stand
    desiredComState.segment(0, 12) << 0, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  }

  // --- update nmpc ---
  nmpcCaller->nmpcUpdate(currentState, gaitTable, desiredComState,
                         desiredFootZPos, desiredFootZVel);
}

void QuadFSM::computeNmpc() {
  nmpcUpdateNeeded = false;

  // --- compute torque ---
  nmpcCaller->getSolution(optimalState, optimalControl);
  stateCur.grfRef = optimalControl.head(12);

  if (!useWbc) {
    for (int i = 0; i < 4; ++i) {
      if (gaitTable[i] == 0) {
        footForceFromKin.col(i) =
            kp_ * (optimalState.segment(12 + 24 + 3 * i, 3) -
                   currentState.segment(12 + 3 * i, 3)) +
            kd_ * (optimalControl.segment(12 + 3 * i, 3) -
                   currentfootVel.segment(0 + 3 * i, 3));
      } else {
        footForceFromKin.col(i) << 0, 0, 0;
      }
    }
    computeJointTorque();
  } else {
    wbcUpdateNeeded = true;
  }
}

void QuadFSM::updateWbcData() {
  wbcData.contact_state.resize(4);
  wbcData.pFoot_des.resize(4);
  wbcData.vFoot_des.resize(4);
  wbcData.aFoot_des.resize(4);
  wbcData.Fr_des.resize(4);

  // update for wbc
  wbcData.state.bodyOrientation = ori::rpyToQuat(stateCur.rpy);
  wbcData.state.bodyPosition = stateCur.pos;

  Matrix3d RotMat = stateCur.rotMat;

  wbcData.state.bodyVelocity.segment(0, 3) = RotMat * stateCur.omegaWorld;
  wbcData.state.bodyVelocity.segment(3, 3) = RotMat * stateCur.linVel;

  std::vector<int> transLeg = {1, 0, 3, 2}; //[todo]
  for (int i = 0; i < 4; ++i) {
    wbcData.state.q.segment(transLeg[i] * 3, 3) = stateCur.legQPos.col(i);
    wbcData.state.qd.segment(transLeg[i] * 3, 3) = stateCur.legQVel.col(i);
  }

  // --- interpolate for wbc reference ---
  VectorXd currentInputCommand(24), nextStateCommand(24);
  currentInputCommand = optimalControl.head(24);
  // [todo] which performs better?
  nextStateCommand = optimalState.segment(24, 24);
  // nmpcIntegrator->getSolution(optimalState.head(24), currentInputCommand,
  //                             nextStateCommand);

  wbcData.pBody_des = nextStateCommand.segment(0, 3);
  wbcData.vBody_des = nextStateCommand.segment(3, 3);
  wbcData.aBody_des = stateCur.linAcc;
  wbcData.pBodyOri_des = ori::rpyToQuat(nextStateCommand.segment(6, 3));
  wbcData.vBodyOri_des = nextStateCommand.segment(9, 3);

  currentGait->setIterations(iterationBetweenMpc, iterationCounter);

  for (int i = 0; i < 4; ++i) {
    wbcData.pFoot_des[transLeg[i]] = nextStateCommand.segment(12 + 3 * i, 3);

    wbcData.vFoot_des[transLeg[i]] = optimalControl.segment(12 + 3 * i, 3);
    wbcData.aFoot_des[transLeg[i]] = Vector3d::Zero();

    if (currentGait->getContactState()[i] > 0 &&
        currentGait->getContactState()[i] <= 1) {
      wbcData.contact_state[transLeg[i]] = 1;
    } else {
      wbcData.contact_state[transLeg[i]] = 0;
    }

    wbcData.Fr_des[transLeg[i]] = optimalControl.segment(3 * i, 3);
  }
}

void QuadFSM::computeWbc() {
  // update
  updateWbcData();

  // --- solve the qp ---
  wbcCaller->run(wbcData, stateCur.jointTorque);
}

// compute joint torque using only mpc
void QuadFSM::computeJointTorque() {
  for (int i = 0; i < 4; ++i) {
    Matrix3d jacobian = legKin[i].getJacobian(stateCur.legQPos.col(i), 14);
    stateCur.jointTorque.segment(3 * i, 3) =
        jacobian.transpose() *
        (footForceFromKin.col(i) -
         stateCur.rotMat * stateCur.grfRef.segment(3 * i, 3));
  }
}

void QuadFSM::bodyOptimize() {
  bodyOptiNeeded = false;

  // body_optimization
  Matrix<double, 6, 1> current_state;
  current_state << stateCur.rpy(0), stateCur.rpy(1), stateCur.rpy(2),
      stateCur.pos(0), stateCur.pos(1), stateCur.pos(2);
  bodyOptimizer->optiSolution(current_state, footHoldDes);

  Eigen::Matrix<double, 5, 1> obj;
  obj = bodyOptimizer->get_opt_obj();
  // std::cout << "obj = " << obj.transpose() << std::endl;
  double filter = 0.1;
  trajIntegrate(0) = trajIntegrate(0) * (1 - filter) + obj(0) * filter;
  trajIntegrate(1) = trajIntegrate(1) * (1 - filter) + obj(1) * filter;
  trajIntegrate(5) = trajIntegrate(5) * (1 - filter) + obj(4) * filter;
}

Matrix<double, 12, 1> QuadFSM::getJointTorque() { return stateCur.jointTorque; }
