#include "WbcCtrl.hpp"
#include "cppTypes.h"
#include "orientation_tools.h"
#include <yaml-cpp/yaml.h>

QuadWbc::QuadWbc()
    : _full_config(18), // [error]?plus 7?
      _tau_ff(12), _des_jpos(12), _des_jvel(12)
// _model(urdf_filename)
{
  _iter = 0;
  _full_config.setZero();

  _kin_wbc = new KinWBC(18);

  _dyn_wbc = new DynWbc(18, &(_contact_list), &(_task_list));
  _dyn_wbc_data = new DynWbcExtraData();

  _body_ori_task = new BodyOriTask(&(_model));
  _body_pos_task = new BodyPosTask(&(_model));

  _foot_contact[0] = new SingleContact(&(_model), 0);
  _foot_contact[1] = new SingleContact(&(_model), 1);
  _foot_contact[2] = new SingleContact(&(_model), 2);
  _foot_contact[3] = new SingleContact(&(_model), 3);

  _foot_task[0] = new LinkPosTask(&(_model), 0);
  _foot_task[1] = new LinkPosTask(&(_model), 1);
  _foot_task[2] = new LinkPosTask(&(_model), 2);
  _foot_task[3] = new LinkPosTask(&(_model), 3);

  _ParameterSetup();
}

QuadWbc::~QuadWbc() {
  delete _kin_wbc;
  delete _dyn_wbc;
  delete _dyn_wbc_data;

  typename std::vector<Task *>::iterator iter = _task_list.begin();
  while (iter < _task_list.end()) {
    delete (*iter);
    ++iter;
  }
  _task_list.clear();

  typename std::vector<ContactSpec *>::iterator iter2 = _contact_list.begin();
  while (iter2 < _contact_list.end()) {
    delete (*iter2);
    ++iter2;
  }
  _contact_list.clear();
}

void QuadWbc::run(const WbcData &input_data, Vec12 &joint_tau) {
  ++_iter;

  _UpdateModel(input_data);

  _ContactTaskUpdate(input_data);

  _ComputeWBC();

  Vec12 joint_torques_temp;

  _UpdateLegCMD(joint_torques_temp, input_data);

  std::vector<int> transLeg = {1, 0, 3, 2};

  for (int i = 0; i < 4; i++) {
    joint_tau.segment(3 * i, 3) =
        joint_torques_temp.segment(3 * transLeg[i], 3);
  }
}

void QuadWbc::_UpdateModel(const WbcData &input_data) {
  _model.updateModel(input_data.state);
  _full_config = _model._full_config;
  _A = _model._A;
  _Ainv = _model._Ainv;
  _coriolis = _model._coriolis;
  _grav = _model._grav;
  // std::cout << "[A]" << _A.block(3, 3, 3, 3) << std::endl;
  // std::cout << _grav.transpose() <<std::endl;
}

void QuadWbc::_UpdateLegCMD(Vec12 &joint_tau, const WbcData &input_data) {
  //[todo] the order of q qv
  //更新要发给腿部的指令 that is tau_ff + pd
  // now the joint tau is listed as the order of q in Mqdd+Cq+g = St + JF
  double err_jpos = 0;
  double err_jvel = 0;
  // std::cout<< std::endl <<"errjpos";
  for (int leg = 0; leg < 4; leg++) {
    for (int jt = 0; jt < 3; jt++) {
      err_jpos = _des_jpos(3 * leg + jt) - input_data.state.q(3 * leg + jt);
      err_jvel = _des_jvel(3 * leg + jt) - input_data.state.qd(3 * leg + jt);
      joint_tau(3 * leg + jt) = _tau_ff(3 * leg + jt) +
                                _Kp_joint[jt] * err_jpos +
                                _Kd_joint[jt] * err_jvel;
    }
  }

  // std::cout << _tau_ff.transpose() << std::endl;
  // std::cout << joint_tau.segment(0,3).transpose() << std::endl;
  // std::cout <<"[desjpos]" <<_des_jpos.segment(0,3).transpose() <<
  // std::endl;//input_data.state.q.transpose() - std::cout <<"[realjpos]" <<
  // input_data.state.q.segment(0,3).transpose() << std::endl;//-
}

void QuadWbc::_ComputeWBC() {
  _kin_wbc->FindConfiguration(_full_config, _task_list, _contact_list,
                              _des_jpos, _des_jvel);

  // std::cout <<"[_full_config]"<< _full_config.transpose() <<std::endl;
  // std::cout <<"[_des_jpos]"<< _des_jpos.segment(6,3).transpose() <<std::endl;
  // std::cout <<"[_des_jvel]"<< _des_jvel.segment(6,3).transpose() <<std::endl;

  _dyn_wbc->UpdateSetting(_A, _Ainv, _coriolis, _grav);

  // std::cout << "_A:" << _A << std::endl;
  // std::cout << "grav" << _grav.transpose() <<std::endl;

  _dyn_wbc->MakeTorque(_tau_ff, _dyn_wbc_data);
  // std::cout <<"[TAU_FF]"<< _tau_ff.segment(0,3).transpose() <<std::endl;
}

void QuadWbc::_ContactTaskUpdate(const WbcData &input_data) {
  // update contact specialize and kinmatics task

  _CleanTaskContact();

  Vec3 zero_vec3;
  zero_vec3.setZero();

  _body_ori_task->UpdateTask(&(input_data.pBodyOri_des),
                             input_data.vBodyOri_des, zero_vec3); // error

  // std::cout <<"[RPY]"<< ori::quatToRPY(input_data.pBodyOri_des).transpose()
  // << std::endl;

  _body_pos_task->UpdateTask(&(input_data.pBody_des), input_data.vBody_des,
                             input_data.aBody_des);

  _task_list.push_back(_body_ori_task);
  _task_list.push_back(_body_pos_task);

  for (int leg = 0; leg < 4; leg++) {
    if (input_data.contact_state[leg] > 0) { // contact
      // std::cout << "[des_fr_contact]" <<input_data.Fr_des[leg].transpose()<<
      // std::endl;
      _foot_contact[leg]->setRFDesired(input_data.Fr_des[leg]);
      _foot_contact[leg]->UpdateContactSpec();
      _contact_list.push_back(_foot_contact[leg]);
    } else { // no contact
      _foot_task[leg]->UpdateTask(&(input_data.pFoot_des[leg]),
                                  input_data.vFoot_des[leg],
                                  input_data.aFoot_des[leg]);
      _task_list.push_back(_foot_task[leg]);
    }
  }
}

void QuadWbc::_CleanTaskContact() {
  _contact_list.clear();
  _task_list.clear();
}

void QuadWbc::_ParameterSetup() {
  YAML::Node config = YAML::LoadFile("../control/wbc_config.yaml");

  _Kp_joint[0] = config["kp_joint_abad"].as<double>();
  _Kp_joint[1] = config["kp_joint_hip"].as<double>();
  _Kp_joint[2] = config["kp_joint_knee"].as<double>();
  _Kd_joint[0] = config["kd_joint_abad"].as<double>();
  _Kd_joint[1] = config["kd_joint_hip"].as<double>();
  _Kd_joint[2] = config["kd_joint_knee"].as<double>();

  ((BodyOriTask *)_body_ori_task)->_Kp[0] =
      config["kp_body_ori_x"].as<double>();
  ((BodyOriTask *)_body_ori_task)->_Kp[1] =
      config["kp_body_ori_y"].as<double>();
  ((BodyOriTask *)_body_ori_task)->_Kp[2] =
      config["kp_body_ori_z"].as<double>();
  ((BodyOriTask *)_body_ori_task)->_Kd[0] =
      config["kd_body_ori_x"].as<double>();
  ((BodyOriTask *)_body_ori_task)->_Kd[1] =
      config["kd_body_ori_y"].as<double>();
  ((BodyOriTask *)_body_ori_task)->_Kd[2] =
      config["kd_body_ori_z"].as<double>();

  ((BodyPosTask *)_body_pos_task)->_Kp[0] =
      config["kp_body_pos_x"].as<double>();
  ((BodyPosTask *)_body_pos_task)->_Kp[1] =
      config["kp_body_pos_y"].as<double>();
  ((BodyPosTask *)_body_pos_task)->_Kp[2] =
      config["kp_body_pos_z"].as<double>();
  ((BodyPosTask *)_body_pos_task)->_Kd[0] =
      config["kd_body_pos_x"].as<double>();
  ((BodyPosTask *)_body_pos_task)->_Kd[1] =
      config["kd_body_pos_y"].as<double>();
  ((BodyPosTask *)_body_pos_task)->_Kd[2] =
      config["kd_body_pos_z"].as<double>();

  for (int i = 0; i < 4; i++) {
    ((LinkPosTask *)_foot_task[i])->_Kp[0] = config["kp_foot_x"].as<double>();
    ((LinkPosTask *)_foot_task[i])->_Kp[1] = config["kp_foot_y"].as<double>();
    ((LinkPosTask *)_foot_task[i])->_Kp[2] = config["kp_foot_z"].as<double>();
    ((LinkPosTask *)_foot_task[i])->_Kd[0] = config["kd_foot_x"].as<double>();
    ((LinkPosTask *)_foot_task[i])->_Kd[1] = config["kd_foot_y"].as<double>();
    ((LinkPosTask *)_foot_task[i])->_Kd[2] = config["kd_foot_z"].as<double>();
  }

  _dyn_wbc_data->_W_floating =
      DVec::Constant(6, config["weight_floating"].as<double>());
  _dyn_wbc_data->_W_rf =
      DVec::Constant(12, config["weight_reaction_force"].as<double>());

  // for(int i(0); i<3; ++i){
  //   ((BodyPosTask*)_body_pos_task)->_Kp[i] = 10.;
  //   ((BodyPosTask*)_body_pos_task)->_Kd[i] = 3.;

  //   ((BodyOriTask*)_body_ori_task)->_Kp[i] = 10.;
  //   ((BodyOriTask*)_body_ori_task)->_Kd[i] = 3.;

  //   for(int j(0); j<4; ++j){
  //     ((LinkPosTask*)_foot_task[j])->_Kp[i] = 70.0;//70
  //     ((LinkPosTask*)_foot_task[j])->_Kd[i] = 3.0;//3
  //     //((LinkPosTask<T>*)_foot_task[j])->_Kp_kin[i] = 1.5;
  //   }

  //   // _Kp_joint[i] = param->Kp_joint[i];
  //   // _Kd_joint[i] = param->Kd_joint[i];

  //   // WBCtrl::_Kp_joint_swing[i] = param->Kp_joint_swing[i];
  //   // WBCtrl::_Kd_joint_swing[i] = param->Kd_joint_swing[i];
  // }
  // _dyn_wbc_data->_W_floating = DVec::Constant(6, 0.1); // WEIGHT_FLOATING
  // _dyn_wbc_data->_W_rf = DVec::Constant(12, 1); // WEIGHT_FR
}
