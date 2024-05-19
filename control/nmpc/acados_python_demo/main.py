from srbd_model import SRBD_Model
from srbd_ocp_setting import SRBD_Optimizer

import numpy as np
import yaml


def readParam(path):
    with open(path, 'r', encoding='utf-8') as f:
        data = f.read()
        result = yaml.load(data, Loader=yaml.FullLoader)

    parameters = np.zeros(22)
    parameters[0] = result["q1"]
    parameters[1] = result["q2"]
    parameters[2] = result["q3"]
    parameters[3] = result["q4"]
    parameters[4] = result["q5"]
    parameters[5] = result["q6"]
    parameters[6] = result["q7"]
    parameters[7] = result["q8"]
    parameters[8] = result["q9"]
    parameters[9] = result["q10"]
    parameters[10] = result["q11"]
    parameters[11] = result["q12"]
    parameters[12] = result["q13"]
    parameters[13] = result["q14"]
    parameters[14] = result["q15"]
    parameters[15] = result["r1"]
    parameters[16] = result["r2"]
    parameters[17] = result["r3"]
    parameters[18] = result["r4"]
    parameters[19] = result["r5"]
    parameters[20] = result["r6"]
    parameters[21] = result["alpha"]
    return parameters


def getFootProjectReference(yaw, pos):

    Rz = np.eye(2)
    Rz[0, 0] = np.cos(yaw)
    Rz[0, 1] = -np.sin(yaw)
    Rz[1, 0] = np.sin(yaw)
    Rz[1, 1] = np.cos(yaw)

    return np.hstack((np.dot(Rz, pos), [0.005]))


def solverSet(opt):
    infValue = 1e+9
    # parameters
    parametersArray = readParam("../../nmpc_config.yaml")
    weightIndex = np.zeros(22)
    for i in range(22):
        weightIndex[i] = i
    for i in range(horizon):
        opt.solver.set_params_sparse(i, weightIndex, parametersArray)

    # constraint bounds
    lh_0 = np.array([-infValue, -infValue, 0, 0, 0, -infValue, 0] * 4 + [0])
    uh_0 = np.array([0, 0, infValue, infValue, 200, 0, infValue] * 4 + [1])
    lh = np.array([0.005, -infValue, -infValue, 0, 0, 0, -infValue, 0] * 4 +
                  [0])
    uh = np.array([0.3, 0, 0, infValue, infValue, 200, 0, infValue] * 4 + [1])
    for i in range(horizon):
        if i == 0:
            opt.solver.constraints_set(i, "lh", lh_0)
            opt.solver.constraints_set(i, "uh", uh_0)
        else:
            opt.solver.constraints_set(i, "lh", lh)
            opt.solver.constraints_set(i, "uh", uh)

    # reference
    Xshoulder = np.array([0.183, 0.183, -0.183, -0.183])
    Yshoulder = np.array([-0.11205, 0.11205, -0.11205, 0.11205])
    currentState = np.array([
        0, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.183, -0.11205, 0.005, 0.183,
        0.11205, 0.005, -0.183, -0.11205, 0.005, -0.183, 0.11205, 0.005
    ])
    desiredComState = np.array([0, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    footPosReference = []
    for i in range(4):
        footPosDefault = np.array([Xshoulder[i], Yshoulder[i]])
        footPosReference = np.hstack(
            (getFootProjectReference(currentState[8],
                                     footPosDefault), footPosReference))
    referenceIndex = np.zeros(48)
    weightParamNum = 22
    for i in range(48):
        referenceIndex[i] = weightParamNum + i
    referenceArray = np.hstack(
        (desiredComState, footPosReference, [0, 0, 30] * 4, [0, 0, 0] * 4))
    for i in range(horizon + 1):
        opt.solver.set_params_sparse(i, referenceIndex, referenceArray)


if __name__ == "__main__":
    horizon = 51
    dtmpc = 0.02
    dt = 0.002

    srbd_model_ = SRBD_Model()
    t_horizon_ = dtmpc * horizon
    n_nodes_ = horizon
    generateCode_ = True

    opt = SRBD_Optimizer(model=srbd_model_,
                         t_horizon=t_horizon_,
                         n_nodes=n_nodes_,
                         solverMode="explicit",
                         generateCode=generateCode_)  # explicit/implicit

    if generateCode_ != True:
        solverSet(opt)

        opt.solver.solve()
        for i in range(horizon):
            print(opt.solver.get(i, "u"))
        for i in range(horizon + 1):
            print(opt.solver.get(i, "x"))
