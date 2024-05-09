from srbd_model import SRBD_Model
from srbd_ocp_setting import SRBD_Optimizer

if __name__ == "__main__":
    horizon = 51
    dtmpc = 0.02
    dt = 0.002

    srbd_model_ = SRBD_Model()
    t_horizon_ = dtmpc * horizon
    n_nodes_ = horizon

    opt = SRBD_Optimizer(
        model=srbd_model_,
        t_horizon=t_horizon_,
        n_nodes=n_nodes_,
        solverMode="multiple shooting")  # multiple shooting/collocation
