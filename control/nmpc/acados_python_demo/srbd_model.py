import numpy as np
import casadi as ca
from acados_template import AcadosModel

# stateVars: [pos, vel(in world frame), rpy, omega(in world frame), footPos(in world frame)]
# controlVars: [grf(in world frame), footVel(in world frame)]


def rpy2RotMat(rpy):
    Rx = ca.SX.eye(3)
    Ry = ca.SX.eye(3)
    Rz = ca.SX.eye(3)
    R = ca.SX.eye(3)

    Rz[0, 0] = ca.cos(rpy[2])
    Rz[0, 1] = -ca.sin(rpy[2])
    Rz[1, 0] = ca.sin(rpy[2])
    Rz[1, 1] = ca.cos(rpy[2])

    Ry[0, 0] = ca.cos(rpy[1])
    Ry[0, 2] = ca.sin(rpy[1])
    Ry[2, 0] = -ca.sin(rpy[1])
    Ry[2, 2] = ca.cos(rpy[1])

    Rx[1, 1] = ca.cos(rpy[0])
    Rx[1, 2] = -ca.sin(rpy[0])
    Rx[2, 1] = ca.sin(rpy[0])
    Rx[2, 2] = ca.cos(rpy[0])

    R = Rz @ Ry @ Rx

    return R


def omega2Rpydot(rpy, omegaWorld):
    rpyDot = ca.SX(3, 1)
    transMat = ca.SX.eye(3)

    # transMat[0, 0] = 1
    # transMat[0, 1] = ca.sin(rpy[0]) * ca.tan(rpy[1])
    # transMat[0, 2] = ca.cos(rpy[0]) * ca.tan(rpy[1])
    # transMat[1, 0] = 0
    # transMat[1, 1] = ca.cos(rpy[0])
    # transMat[1, 2] = ca.sin(rpy[0])
    # transMat[2, 0] = 0
    # transMat[2, 1] = ca.sin(rpy[0]) / ca.cos(rpy[1])
    # transMat[2, 2] = ca.cos(rpy[0]) / ca.cos(rpy[1])

    transMat[0, 0] = ca.cos(rpy[2]) / ca.cos(rpy[1])
    transMat[0, 1] = ca.sin(rpy[2]) / ca.cos(rpy[1])
    transMat[0, 2] = 0
    transMat[1, 0] = -ca.sin(rpy[2])
    transMat[1, 1] = ca.cos(rpy[2])
    transMat[1, 2] = 0
    transMat[2, 0] = ca.cos(rpy[2]) * ca.tan(rpy[1])
    transMat[2, 1] = ca.sin(rpy[2]) * ca.tan(rpy[1])
    transMat[2, 2] = 1

    rpyDot = transMat @ omegaWorld

    return rpyDot


def vec2Mat(vector):
    mat = ca.SX.zeros(3, 3)

    mat[0, 1] = -vector[2]
    mat[0, 2] = vector[1]
    mat[1, 0] = vector[2]
    mat[1, 2] = -vector[0]
    mat[2, 0] = -vector[1]
    mat[2, 1] = vector[0]

    return mat


# srbd model
class SRBD_Model(object):

    def __init__(self, ):
        model = AcadosModel()
        model.name = "srbd"

        stateVars = ca.SX.sym("stateVars", 24)
        stateVarsDot = ca.SX.sym("stateVarsDot", 24)
        controlVars = ca.SX.sym("controlVars", 24)

        weight = ca.SX.sym("weight", 22)
        stateReference = ca.SX.sym("stateReference", 24)
        controlReference = ca.SX.sym("controlReference", 24)

        self.stateVars = stateVars
        self.rotMat = rpy2RotMat(self.stateVars[6:9])
        self.controlVars = controlVars
        self.weight = weight
        self.stateReference = stateReference
        self.controlReference = controlReference

        model.f_impl_expr = stateVarsDot - self.dynamics(
            stateVars, controlVars)
        model.f_expl_expr = self.dynamics(stateVars, controlVars)
        model.disc_dyn_expr = self.rk4()
        model.x = stateVars
        model.xdot = stateVarsDot
        model.u = controlVars
        model.p = ca.vertcat(weight, stateReference, controlReference)

        model.con_h_expr_0 = self.constraints_0()
        model.con_h_expr = self.constraints()
        # model.con_h_expr_e = self.constraints_e()

        self.model = model

        constraintBounds = ca.types.SimpleNamespace()
        constraintBounds.lower_bound = np.zeros(36)
        constraintBounds.upper_bound = np.zeros(36)
        constraintBounds_0 = ca.types.SimpleNamespace()
        constraintBounds_0.lower_bound = np.zeros(32)
        constraintBounds_0.upper_bound = np.zeros(32)
        constraintBounds_e = ca.types.SimpleNamespace()
        constraintBounds_e.lower_bound = np.zeros(24)
        constraintBounds_e.upper_bound = np.zeros(24)
        self.constraintBounds = constraintBounds
        self.constraintBounds_0 = constraintBounds_0
        self.constraintBounds_e = constraintBounds_e

    def getCostMatrix(self, ):
        Q_p = np.array([self.weight[0], self.weight[1], self.weight[2]])
        Q_v = np.array([self.weight[3], self.weight[4], self.weight[5]])
        Q_q = np.array([self.weight[6], self.weight[7], self.weight[8]])
        Q_w = np.array([self.weight[9], self.weight[10], self.weight[11]])
        Q_f = np.array([self.weight[12], self.weight[13], self.weight[14]])
        Q = np.concatenate((Q_p, Q_v, Q_q, Q_w, Q_f, Q_f, Q_f, Q_f))

        R = np.array([
            self.weight[15], self.weight[16], self.weight[17], self.weight[15],
            self.weight[16], self.weight[17], self.weight[15], self.weight[16],
            self.weight[17], self.weight[15], self.weight[16], self.weight[17],
            self.weight[18], self.weight[19], self.weight[20], self.weight[18],
            self.weight[19], self.weight[20], self.weight[18], self.weight[19],
            self.weight[20], self.weight[18], self.weight[19], self.weight[20]
        ])

        return np.diag(np.concatenate((Q, R)))

    def getTerminalCostMatrix(self, ):
        Q_p_e = np.array([self.weight[0], self.weight[1], self.weight[2]])
        Q_v_e = np.array([self.weight[3], self.weight[4], self.weight[5]])
        Q_q_e = np.array([self.weight[6], self.weight[7], self.weight[8]])
        Q_w_e = np.array([self.weight[9], self.weight[10], self.weight[11]])
        Q_f_e = np.array([self.weight[12], self.weight[13], self.weight[14]])

        return np.diag((np.concatenate(
            (Q_p_e, Q_v_e, Q_q_e, Q_w_e, Q_f_e, Q_f_e, Q_f_e,
             Q_f_e)))) * self.weight[21]

    def cost(self, ):
        W = self.getCostMatrix()
        err = (ca.vertcat(
            self.stateVars[0:12], self.stateVars[12:14] - self.stateVars[0:2],
            self.stateVars[14], self.stateVars[15:17] - self.stateVars[0:2],
            self.stateVars[17], self.stateVars[18:20] - self.stateVars[0:2],
            self.stateVars[20], self.stateVars[21:23] - self.stateVars[0:2],
            self.stateVars[23], self.controlVars) -
               ca.vertcat(self.stateReference, self.controlReference))

        return err.T @ W @ err

    def cost_e(self, ):
        Q_e = self.getTerminalCostMatrix()
        err = self.stateVars - self.stateReference

        return err.T @ Q_e @ err

    def constraints(self, ):
        mu = 0.4
        constraints = ca.SX.sym("constraints", 36)
        for i in range(4):
            # foot ZPos
            constraints[i * 9 + 0] = self.stateVars[i * 3 + 12 + 2]
            # friction
            constraints[i * 9 + 1] = (self.controlVars[i * 3 + 0] -
                                      mu * self.controlVars[i * 3 + 2])
            constraints[i * 9 + 2] = (self.controlVars[i * 3 + 1] -
                                      mu * self.controlVars[i * 3 + 2])
            constraints[i * 9 + 3] = (self.controlVars[i * 3 + 0] +
                                      mu * self.controlVars[i * 3 + 2])
            constraints[i * 9 + 4] = (self.controlVars[i * 3 + 1] +
                                      mu * self.controlVars[i * 3 + 2])
            constraints[i * 9 + 5] = self.controlVars[i * 3 + 2]
            # foot vel
            constraints[i * 9 + 6] = self.controlVars[i * 3 + 12 + 0]
            constraints[i * 9 + 7] = self.controlVars[i * 3 + 12 + 1]
            constraints[i * 9 + 8] = self.controlVars[i * 3 + 12 + 2]

        return constraints

    def constraints_0(self, ):
        mu = 0.4
        constraints_0 = ca.SX.sym("constraints_0", 32)
        for i in range(4):
            constraints_0[i * 8 + 0] = (self.controlVars[i * 3 + 0] -
                                        mu * self.controlVars[i * 3 + 2])
            constraints_0[i * 8 + 1] = (self.controlVars[i * 3 + 1] -
                                        mu * self.controlVars[i * 3 + 2])
            constraints_0[i * 8 + 2] = (self.controlVars[i * 3 + 0] +
                                        mu * self.controlVars[i * 3 + 2])
            constraints_0[i * 8 + 3] = (self.controlVars[i * 3 + 1] +
                                        mu * self.controlVars[i * 3 + 2])
            constraints_0[i * 8 + 4] = self.controlVars[i * 3 + 2]
            # foot vel
            constraints_0[i * 8 + 5] = self.controlVars[i * 3 + 12 + 0]
            constraints_0[i * 8 + 6] = self.controlVars[i * 3 + 12 + 1]
            constraints_0[i * 8 + 7] = self.controlVars[i * 3 + 12 + 2]

        return constraints_0

    def constraints_e(self, ):
        constraints_e = ca.SX.sym("constraints_e", 24)

        for i in range(24):
            constraints_e[i] = self.stateVars[i]

        return constraints_e

    def dynamics(self, stateVars, controlVars):

        inertia = ca.SX.zeros(3, 3)
        inertia[0, 0] = 0.0158533
        inertia[1, 1] = 0.0377999
        inertia[2, 2] = 0.0456542
        mass = 12
        g = [0, 0, 9.81]

        inertia_inv = ca.inv(inertia)

        # switch to global inertia
        inertia = self.rotMat @ inertia @ self.rotMat.T
        inertia_inv = self.rotMat.T @ inertia_inv @ self.rotMat

        tau = (ca.cross(ca.vertcat(stateVars[12:15] - stateVars[0:3]),
                        controlVars[0:3], -1) +
               ca.cross(ca.vertcat(stateVars[15:18] - stateVars[0:3]),
                        controlVars[3:6], -1) +
               ca.cross(ca.vertcat(stateVars[18:21] - stateVars[0:3]),
                        controlVars[6:9], -1) +
               ca.cross(ca.vertcat(stateVars[21:24] - stateVars[0:3]),
                        controlVars[9:12], -1))
        f = (controlVars[0:3] + controlVars[3:6] + controlVars[6:9] +
             controlVars[9:12])

        state_space = ca.vertcat(
            stateVars[3:6], f / mass - g,
            omega2Rpydot(stateVars[6:9], stateVars[9:12]),
            inertia_inv @ (tau - vec2Mat(self.stateVars[9:12])
                           @ (inertia @ self.stateVars[9:12])),
            controlVars[12:15], controlVars[15:18], controlVars[18:21],
            controlVars[21:24])

        return state_space

    def rk4(self, ):
        k1 = 0.02 * self.dynamics(self.stateVars, self.controlVars)
        k2 = 0.02 * self.dynamics(self.stateVars + k1 / 2, self.controlVars)
        k3 = 0.02 * self.dynamics(self.stateVars + k2 / 2, self.controlVars)
        k4 = 0.02 * self.dynamics(self.stateVars + k3, self.controlVars)

        stateNext = self.stateVars + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        return stateNext
