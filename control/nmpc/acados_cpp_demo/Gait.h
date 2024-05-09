
#pragma once

#include <Eigen/Eigen>
#include <queue>
#include <string>

using Eigen::Array4d;
using Eigen::Array4i;

class Gait {
public:
  virtual ~Gait() = default;

  virtual Eigen::Vector<double, 4> getContactState() = 0;
  virtual Eigen::Vector<double, 4> getSwingState() = 0;
  virtual int *getMpcTable() = 0;
  virtual void setIterations(int iterationsBetweenMPC,
                             int currentIteration) = 0;
  virtual double getCurrentStanceTime(double dtMPC, int leg) = 0;
  virtual double getCurrentSwingTime(double dtMPC, int leg) = 0;
  virtual int getCurrentGaitPhase() = 0;
  virtual void debugPrint() {}

protected:
  std::string _name;
};

class OffsetDurationGait : public Gait {
public:
  OffsetDurationGait(int nSegment, Eigen::Vector<int, 4> offset,
                     Eigen::Vector<int, 4> durations, const std::string &name);
  ~OffsetDurationGait();
  Eigen::Vector<double, 4> getContactState();
  Eigen::Vector<double, 4> getSwingState();
  int *getMpcTable();
  void setIterations(int iterationsBetweenMPC, int currentIteration);
  double getCurrentStanceTime(double dtMPC, int leg);
  double getCurrentSwingTime(double dtMPC, int leg);
  int getCurrentGaitPhase();
  void debugPrint();

private:
  int *_mpc_table;
  Array4i _offsets;         // offset in mpc segments
  Array4i _durations;       // duration of step in mpc segments
  Array4d _offsetsDouble;   // offsets in phase (0 to 1)
  Array4d _durationsDouble; // durations in phase (0 to 1)
  int _stance;
  int _swing;
  int _iteration;
  int _nIterations;
  double _phase;
};

class MixedFrequncyGait : public Gait {
public:
  MixedFrequncyGait(int nSegment, Eigen::Vector<int, 4> periods,
                    double duty_cycle, const std::string &name);
  ~MixedFrequncyGait();
  Eigen::Vector<double, 4> getContactState();
  Eigen::Vector<double, 4> getSwingState();
  int *getMpcTable();
  void setIterations(int iterationsBetweenMPC, int currentIteration);
  double getCurrentStanceTime(double dtMPC, int leg);
  double getCurrentSwingTime(double dtMPC, int leg);
  int getCurrentGaitPhase();
  void debugPrint();

private:
  double _duty_cycle;
  int *_mpc_table;
  Array4i _periods;
  Array4d _phase;
  int _iteration;
  int _nIterations;
};
