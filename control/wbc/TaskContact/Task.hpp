#ifndef WBC_TASK
#define WBC_TASK

#include <cppTypes.h>

#define TK Task
// 运动任务抽象类
class Task {
public:
  Task(int dim)
      : b_set_task_(false), dim_task_(dim), op_cmd_(dim), pos_err_(dim),
        vel_des_(dim), acc_des_(dim) {}

  virtual ~Task() {}

  void getCommand(DVec &op_cmd) { op_cmd = op_cmd_; } // get ddot{x}_{i}^{des}
  void getTaskJacobian(DMat &Jt) { Jt = Jt_; }
  void getTaskJacobianDotQdot(DVec &JtDotQdot) { JtDotQdot = JtDotQdot_; }

  bool UpdateTask(const void *pos_des, const DVec &vel_des,
                  const DVec &acc_des) {
    _UpdateTaskJacobian();
    _UpdateTaskJDotQdot();
    _UpdateCommand(pos_des, vel_des, acc_des);
    _AdditionalUpdate();
    b_set_task_ = true;
    return true;
  }

  bool IsTaskSet() { return b_set_task_; }
  int getDim() { return dim_task_; }
  void UnsetTask() { b_set_task_ = false; }

  const DVec &getPosError() { return pos_err_; }
  const DVec &getDesVel() { return vel_des_; }
  const DVec &getDesAcc() { return acc_des_; }

protected:
  // Update op_cmd_ : ddot{x}_{i}^{des}
  virtual bool _UpdateCommand(const void *pos_des, const DVec &vel_des,
                              const DVec &acc_des) = 0;
  // Update Jt_
  virtual bool _UpdateTaskJacobian() = 0;
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot() = 0;
  // Additional Update (defined in child classes)
  virtual bool _AdditionalUpdate() = 0;

  bool b_set_task_;
  size_t dim_task_;

  DVec op_cmd_;
  DVec JtDotQdot_;
  DMat Jt_;

  DVec pos_err_;
  DVec vel_des_;
  DVec acc_des_;
};

#endif