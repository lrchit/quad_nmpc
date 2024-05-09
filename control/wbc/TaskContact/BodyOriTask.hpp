#ifndef BODY_ORIENTATION_TASK
#define BODY_ORIENTATION_TASK
// (Rx, Ry, Rz)
#include "cppTypes.h"
#include "dynamics.hpp"
#include <Task.hpp>

class BodyOriTask : public Task {
public:
  BodyOriTask(const FBDynModel *);
  virtual ~BodyOriTask();

  DVec _Kp_kin;
  DVec _Kp, _Kd;

protected:
  // Update op_cmd_
  virtual bool _UpdateCommand(const void *pos_des, const DVec &vel_des,
                              const DVec &acc_des);
  // Update Jt_
  virtual bool _UpdateTaskJacobian();
  // Update JtDotQdot_
  virtual bool _UpdateTaskJDotQdot();
  virtual bool _AdditionalUpdate() { return true; }

  int link_idx_;
  bool virtual_depend_;
  const FBDynModel *_robot_sys;
};

#endif
