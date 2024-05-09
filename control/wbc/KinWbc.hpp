#include "ContactSpec.hpp"
#include "Task.hpp"
#include "cppTypes.h"
#include <vector>

class KinWBC {
public:
  KinWBC(int num_qdot);
  ~KinWBC() {}

  bool FindConfiguration(const DVec &curr_config,
                         const std::vector<Task *> &task_list,
                         const std::vector<ContactSpec *> &contact_list,
                         DVec &jpos_cmd, DVec &jvel_cmd);

  DMat Ainv_;

private:
  void _PseudoInverse(const DMat J, DMat &Jinv);
  void _BuildProjectionMatrix(const DMat &J, DMat &N);

  double threshold_;
  int num_qdot_;
  int num_act_joint_;
  DMat I_mtx;
};