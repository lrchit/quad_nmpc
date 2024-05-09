#pragma once
#include "FSM.h"
/**
 * @class QuadApp
 * @brief 单例设计
 *
 */
class QuadApp {
public:
  static QuadApp *instance() {
    if (inst_ptr == nullptr) {
      inst_ptr = new QuadApp;
    }
    return inst_ptr;
  }
  void Start();

  QuadFSM FSM;

  bool flag_control_update = false;

private:
  static QuadApp *inst_ptr;
};
