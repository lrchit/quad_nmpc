#include "Gait.h"

// Offset - Duration Gait
OffsetDurationGait::OffsetDurationGait(int nSegment,
                                       Eigen::Vector<int, 4> offsets,
                                       Eigen::Vector<int, 4> durations,
                                       const std::string &name)
    : _offsets(offsets.array()), _durations(durations.array()),
      _nIterations(nSegment) {

  _name = name;
  // allocate memory for MPC gait table
  _mpc_table = new int[nSegment * 4];

  _offsetsDouble = offsets.cast<double>() / (double)nSegment;
  _durationsDouble = durations.cast<double>() / (double)nSegment;

  _stance = durations[0];
  _swing = nSegment - durations[0];
}

MixedFrequncyGait::MixedFrequncyGait(int nSegment,
                                     Eigen::Vector<int, 4> periods,
                                     double duty_cycle,
                                     const std::string &name) {
  _name = name;
  _duty_cycle = duty_cycle;
  _mpc_table = new int[nSegment * 4];
  _periods = periods;
  _nIterations = nSegment;
  _iteration = 0;
  _phase.setZero();
}

OffsetDurationGait::~OffsetDurationGait() { delete[] _mpc_table; }

MixedFrequncyGait::~MixedFrequncyGait() { delete[] _mpc_table; }

Eigen::Vector<double, 4> OffsetDurationGait::getContactState() {
  Array4d progress = _phase - _offsetsDouble;

  for (int i = 0; i < 4; i++) {
    if (progress[i] < 0)
      progress[i] += 1.;
    if (progress[i] > _durationsDouble[i]) {
      progress[i] = 0.;
    } else {
      progress[i] = progress[i] / _durationsDouble[i];
    }
  }

  // printf("contact state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1],
  // progress[2], progress[3]);
  return progress.matrix();
}

Eigen::Vector<double, 4> MixedFrequncyGait::getContactState() {
  Array4d progress = _phase;

  for (int i = 0; i < 4; i++) {
    if (progress[i] < 0)
      progress[i] += 1.;
    if (progress[i] > _duty_cycle) {
      progress[i] = 0.;
    } else {
      progress[i] = progress[i] / _duty_cycle;
    }
  }

  // printf("contact state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1],
  // progress[2], progress[3]);
  return progress.matrix();
}

Eigen::Vector<double, 4> OffsetDurationGait::getSwingState() {
  Array4d swing_offset = _offsetsDouble + _durationsDouble;
  for (int i = 0; i < 4; i++)
    if (swing_offset[i] > 1)
      swing_offset[i] -= 1.;
  Array4d swing_duration = 1. - _durationsDouble;

  Array4d progress = _phase - swing_offset;

  for (int i = 0; i < 4; i++) {
    if (progress[i] < 0)
      progress[i] += 1.f;
    if (progress[i] > swing_duration[i]) {
      progress[i] = 0.;
    } else {
      progress[i] = progress[i] / swing_duration[i];
    }
  }

  // printf("swing state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1],
  // progress[2], progress[3]);
  return progress.matrix();
}

Eigen::Vector<double, 4> MixedFrequncyGait::getSwingState() {

  double swing_duration = 1.f - _duty_cycle;
  Array4d progress = _phase - _duty_cycle;
  for (int i = 0; i < 4; i++) {
    if (progress[i] < 0) {
      progress[i] = 0;
    } else {
      progress[i] = progress[i] / swing_duration;
    }
  }

  // printf("swing state: %.3f %.3f %.3f %.3f\n", progress[0], progress[1],
  // progress[2], progress[3]);
  return progress.matrix();
}

int *OffsetDurationGait::getMpcTable() {

  // printf("MPC table:\n");
  for (int i = 0; i < _nIterations; i++) {
    // int iter = (i + _iteration + 1) % _nIterations;
    int iter = (i + _iteration) % _nIterations;
    Array4i progress = iter - _offsets;
    for (int j = 0; j < 4; j++) {
      if (progress[j] < 0)
        progress[j] += _nIterations;
      if (progress[j] < _durations[j])
        _mpc_table[i * 4 + j] = 1;
      else
        _mpc_table[i * 4 + j] = 0;

      // printf("%d ", _mpc_table[i*4 + j]);
    }
    // printf("\n");
  }

  return _mpc_table;
}

int *MixedFrequncyGait::getMpcTable() {
  // printf("MPC table (%d):\n", _iteration);
  for (int i = 0; i < _nIterations; i++) {
    for (int j = 0; j < 4; j++) {
      int progress = (i + _iteration + 1) % _periods[j]; // progress
      if (progress < (_periods[j] * _duty_cycle)) {
        _mpc_table[i * 4 + j] = 1;
      } else {
        _mpc_table[i * 4 + j] = 0;
      }
      // printf("%d %d (%d %d) | ", _mpc_table[i*4 + j], progress, _periods[j],
      // (int)(_periods[j] * _duty_cycle));
    }

    // printf("%d %d %d %d (%.3f %.3f %.3f %.3f)\n", _mpc_table[i*4],
    // _mpc_table[i*4 + 1], _mpc_table[i*4 + ]) printf("\n");
  }
  return _mpc_table;
}

void OffsetDurationGait::setIterations(int iterationsPerMPC,
                                       int currentIteration) {
  _iteration = (currentIteration / iterationsPerMPC) % _nIterations;
  _phase = (double)(currentIteration % (iterationsPerMPC * _nIterations)) /
           (double)(iterationsPerMPC * _nIterations);
}

void MixedFrequncyGait::setIterations(int iterationsBetweenMPC,
                                      int currentIteration) {
  _iteration = (currentIteration / iterationsBetweenMPC); // % _nIterations;
  for (int i = 0; i < 4; i++) {
    int progress_mult = currentIteration % (iterationsBetweenMPC * _periods[i]);
    _phase[i] =
        ((double)progress_mult) / ((double)iterationsBetweenMPC * _periods[i]);
    //_phase[i] = (double)(currentIteration % (iterationsBetweenMPC *
    //_periods[i])) / (double) (iterationsBetweenMPC * _periods[i]);
  }

  // printf("phase: %.3f %.3f %.3f %.3f\n", _phase[0], _phase[1], _phase[2],
  // _phase[3]);
}

int OffsetDurationGait::getCurrentGaitPhase() { return _iteration; }

int MixedFrequncyGait::getCurrentGaitPhase() { return 0; }

double OffsetDurationGait::getCurrentSwingTime(double dtMPC, int leg) {
  (void)leg;
  return dtMPC * _swing;
}

double MixedFrequncyGait::getCurrentSwingTime(double dtMPC, int leg) {
  return dtMPC * (1. - _duty_cycle) * _periods[leg];
}

double OffsetDurationGait::getCurrentStanceTime(double dtMPC, int leg) {
  (void)leg;
  return dtMPC * _stance;
}

double MixedFrequncyGait::getCurrentStanceTime(double dtMPC, int leg) {
  return dtMPC * _duty_cycle * _periods[leg];
}

void OffsetDurationGait::debugPrint() {}

void MixedFrequncyGait::debugPrint() {}