import os
import sys
import shutil
import errno
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSim, AcadosSimSolver


def safe_mkdir_recursive(directory, overwrite=False):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                raise
    else:
        if overwrite:
            try:
                shutil.rmtree(directory)
            except:
                print("Error while removing directory {}".format(directory))


class SRBD_Optimizer(object):

    def __init__(self, model, t_horizon, n_nodes, solverMode="explicit"):

        self.T = t_horizon
        self.N = n_nodes

        # dims
        nx = model.model.x.size()[0]
        self.nx = nx
        nu = model.model.u.size()[0]
        self.nu = nu
        ny = nx + nu
        n_params = model.model.p.size()[0]

        acados_source_path = os.environ["ACADOS_SOURCE_DIR"]
        sys.path.insert(0, acados_source_path)

        # OCP solver
        ocp = AcadosOcp()
        ocp.acados_include_path = acados_source_path + "/include"
        ocp.acados_lib_path = acados_source_path + "/lib"
        # simulator
        sim = AcadosSim()
        sim.acados_include_path = acados_source_path + '/include'
        sim.acados_lib_path = acados_source_path + '/lib'

        ocp.solver_options.tf = self.T
        sim.solver_options.T = self.T / self.N
        if (solverMode == "explicit"):
            ocp.solver_options.integrator_type = "ERK"
            sim.solver_options.integrator_type = "ERK"
        elif (solverMode == "implicit"):
            ocp.solver_options.integrator_type = "IRK"
            ocp.solver_options.collocation_type = "GAUSS_LEGENDRE"
            sim.solver_options.integrator_type = "IRK"
            sim.solver_options.collocation_type = "GAUSS_LEGENDRE"
        else:
            print("ERROR: 'solverMode' should be 'explicit' or 'implicit'")
            exit()

        ocp.model = model.model

        ocp.dims.N = self.N
        ocp.dims.nx = self.nx
        ocp.dims.nu = self.nu

        # param initialize
        ocp.dims.np = n_params
        ocp.parameter_values = np.zeros(n_params)

        # # cost using LQ function
        # ocp.cost.cost_type = "LINEAR_LS"
        # ocp.cost.cost_type_e = "LINEAR_LS"
        # ocp.cost.W = model.getCostMatrix()
        # ocp.cost.W_e = model.getTerminalCostMatrix()
        # ocp.cost.Vx = np.zeros((ny, nx))
        # ocp.cost.Vx[:nx, :nx] = np.eye(nx)
        # ocp.cost.Vu = np.zeros((ny, nu))
        # ocp.cost.Vu[-nu:, -nu:] = np.eye(nu)
        # ocp.cost.Vx_e = np.eye(nx)
        # x_ref = np.zeros(nx)
        # u_ref = np.zeros(nu)
        # ocp.cost.yref = np.concatenate((x_ref, u_ref))
        # ocp.cost.yref_e = x_ref

        # use this for external cost
        x_ref = np.zeros(nx)
        ocp.cost.cost_type = "EXTERNAL"
        ocp.cost.cost_type_e = "EXTERNAL"
        ocp.model.cost_expr_ext_cost = model.cost()
        ocp.model.cost_expr_ext_cost_e = model.cost_e()

        # constraint bounds
        ocp.constraints.lh_0 = model.constraintBounds_0.lower_bound
        ocp.constraints.uh_0 = model.constraintBounds_0.upper_bound
        ocp.constraints.lh = model.constraintBounds.lower_bound
        ocp.constraints.uh = model.constraintBounds.upper_bound
        # ocp.constraints.lh_e = model.constraintBounds_e.lower_bound
        # ocp.constraints.uh_e = model.constraintBounds_e.upper_bound
        ocp.constraints.x0 = x_ref

        # solver options，using HPIPM solver，using Gauss-Newton for Hessian
        ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
        ocp.solver_options.ext_fun_compile_flags = "-O3"
        ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
        ocp.solver_options.qp_solver_warm_start = 2  # this does't make a big difference
        # ocp.solver_options.qp_solver_ric_alg = 0
        # ocp.solver_options.qp_solver_cond_ric_alg = 0
        # ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
        # ocp.solver_options.levenberg_marquardt = 0.1
        # ocp.solver_options.regularize_method = "MIRROR"
        ocp.solver_options.sim_method_num_stages = 4
        ocp.solver_options.print_level = 0
        ocp.solver_options.hpipm_mode = "SPEED_ABS"
        ocp.solver_options.nlp_solver_type = "SQP_RTI"
        ocp.solver_options.nlp_solver_max_iter = 20
        ocp.solver_options.nlp_solver_tol_stat = 1e-6
        ocp.solver_options.nlp_solver_tol_eq = 1e-6
        ocp.solver_options.nlp_solver_tol_ineq = 1e-6
        ocp.solver_options.nlp_solver_tol_comp = 1e-6

        # # use ddp solver
        # ocp.translate_to_feasibility_problem(True, False)
        # ocp.solver_options.nlp_solver_type = "DDP"

        # generate ocp solver
        ocp.code_export_directory = "c_generated_code"
        json_file = os.path.join("./" + model.model.name + "_acados_ocp.json")
        self.solver = AcadosOcpSolver(
            ocp,
            json_file=json_file)  # build=False, generate=False, verbose=False)

        # set model
        sim.model = model.model

        sim.dims.nx = nx
        sim.dims.nu = nu

        # interrupts(unnecessary)
        # To be added later

        # param initialize
        sim.dims.np = n_params
        sim.parameter_values = np.zeros(n_params)

        # solver options
        sim.solver_options.sim_method_num_stages = 4
        sim.solver_options.ext_fun_compile_flags = "-O3"

        # generate simulator
        sim.code_export_directory = "c_generated_code"
        json_file = os.path.join("./" + model.model.name + "_acados_sim.json")
        self.integrator = AcadosSimSolver(
            sim,
            json_file=json_file)  # generate=False, build=False, verbose=False)
