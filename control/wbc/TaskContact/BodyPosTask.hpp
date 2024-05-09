#ifndef BODY_POS_TASK
#define BODY_POS_TASK

// (X, Y, Z)
#include "dynamics.hpp"
#include <Task.hpp>

class BodyPosTask : public Task {
public:
  BodyPosTask(const FBDynModel *robot);
  virtual ~BodyPosTask();

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

  const FBDynModel *_robot_sys;
};

#endif
