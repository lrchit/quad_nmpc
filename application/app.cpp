#include "app.h"

#include <mujoco/mujoco.h>

#include <fstream>
#include <future>

#include "array_safety.h"
#include "estimate.h"
#include "glfw_adapter.h"
#include "iostream"
#include "kinematics.h"
#include "simulate.h"
#include "threadpool.h"
#include "yaml-cpp/yaml.h"

namespace mj = ::mujoco;
namespace mju = ::mujoco::sample_util;

std::recursive_mutex g_mutex;

// constants
const double syncMisalign =
    0.1; // maximum mis-alignment before re-sync (simulation seconds)
const double simRefreshFraction =
    0.7;                       // fraction of refresh available for simulation
const int kErrorLength = 1024; // load error string length

// model and data
mjModel *m = nullptr;
mjData *d = nullptr;

// control noise variables
mjtNum *ctrlnoise = nullptr;

using Seconds = std::chrono::duration<double>;

std::shared_ptr<mujoco::Simulate> sim;

// app inst
QuadApp *QuadApp::inst_ptr = nullptr;

/**
 * @brief 加载模型
 *
 * @param file
 * @param sim
 * @return
 */
mjModel *LoadModel(const char *file, mj::Simulate &sim) {
  // this copy is needed so that the mju::strlen call below compiles
  char filename[mj::Simulate::kMaxFilenameLength];
  mju::strcpy_arr(filename, file);

  // make sure filename is not empty
  if (!filename[0]) {
    return nullptr;
  }

  // load and compile
  char loadError[kErrorLength] = "";
  mjModel *mnew = 0;
  if (mju::strlen_arr(filename) > 4 &&
      !std::strncmp(filename + mju::strlen_arr(filename) - 4, ".mjb",
                    mju::sizeof_arr(filename) - mju::strlen_arr(filename) +
                        4)) {
    mnew = mj_loadModel(filename, nullptr);
    if (!mnew) {
      mju::strcpy_arr(loadError, "could not load binary model");
    }
  } else {
    mnew = mj_loadXML(filename, nullptr, loadError, kErrorLength);
    // remove trailing newline character from loadError
    if (loadError[0]) {
      int error_length = mju::strlen_arr(loadError);
      if (loadError[error_length - 1] == '\n') {
        loadError[error_length - 1] = '\0';
      }
    }
  }

  mju::strcpy_arr(sim.load_error, loadError);

  if (!mnew) {
    std::printf("%s\n", loadError);
    return nullptr;
  }

  // compiler warning: print and pause
  if (loadError[0]) {
    // mj_forward() below will print the warning message
    std::printf("Model compiled, but simulation warning (paused):\n  %s\n",
                loadError);
    sim.run = 0;
  }

  return mnew;
}

/**
 * @brief 物理引擎线程
 *
 * @param sim
 */
void PhysicsLoop(mj::Simulate &sim) {
  // cpu-sim syncronization point
  std::chrono::time_point<mj::Simulate::Clock> syncCPU;
  mjtNum syncSim = 0;

  // run until asked to exit
  while (!sim.exitrequest.load()) {
    if (sim.droploadrequest.load()) {
    }

    if (sim.uiloadrequest.load()) {
      sim.LoadMessage(sim.filename);
      mjModel *mnew = LoadModel(sim.filename, sim);
      mjData *dnew = nullptr;
      if (mnew)
        dnew = mj_makeData(mnew);
      if (dnew) {
        sim.Load(mnew, dnew, sim.filename);

        mj_deleteData(d);
        mj_deleteModel(m);

        m = mnew;
        d = dnew;
        mj_forward(m, d);

        YAML::Node config = YAML::LoadFile("../application/config.yaml");
        d->qpos[7] = config["p1"].as<double>();
        d->qpos[8] = config["p2"].as<double>();
        d->qpos[9] = config["p3"].as<double>();
        d->qpos[10] = config["p4"].as<double>();
        d->qpos[11] = config["p5"].as<double>();
        d->qpos[12] = config["p6"].as<double>();
        d->qpos[13] = config["p7"].as<double>();
        d->qpos[14] = config["p8"].as<double>();
        d->qpos[15] = config["p9"].as<double>();
        d->qpos[16] = config["p10"].as<double>();
        d->qpos[17] = config["p11"].as<double>();
        d->qpos[18] = config["p12"].as<double>();

        // allocate ctrlnoise
        free(ctrlnoise);
        ctrlnoise = static_cast<mjtNum *>(malloc(sizeof(mjtNum) * m->nu));
        mju_zero(ctrlnoise, m->nu);
      } else {
        sim.LoadMessageClear();
      }
      sim.uiloadrequest.fetch_sub(1);
    }

    // sleep for 1 ms or yield, to let main thread run
    //  yield results in busy wait - which has better timing but kills
    //  battery life
    if (sim.run && sim.busywait) {
      std::this_thread::yield();
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    {
      // lock the sim mutex
      const std::unique_lock<std::recursive_mutex> lock(sim.mtx);

      // run only if model is present
      if (m) {
        // running
        if (sim.run) {
          // record cpu time at start of iteration
          const auto startCPU = mj::Simulate::Clock::now();

          // elapsed CPU and simulation time since last sync
          const auto elapsedCPU = startCPU - syncCPU;
          double elapsedSim = d->time - syncSim;

          // inject noise
          if (sim.ctrl_noise_std) {
            // convert rate and scale to discrete time
            // (Ornstein–Uhlenbeck)
            mjtNum rate = mju_exp(-m->opt.timestep /
                                  mju_max(sim.ctrl_noise_rate, mjMINVAL));
            mjtNum scale = sim.ctrl_noise_std * mju_sqrt(1 - rate * rate);

            for (int i = 0; i < m->nu; i++) {
              // update noise
              ctrlnoise[i] =
                  rate * ctrlnoise[i] + scale * mju_standardNormal(nullptr);

              // apply noise
              d->ctrl[i] = ctrlnoise[i];
            }
          }

          // requested slow-down factor
          double slowdown = 100 / sim.percentRealTime[sim.real_time_index];

          // misalignment condition: distance from target sim time is
          // bigger than syncmisalign
          bool misaligned = mju_abs(Seconds(elapsedCPU).count() / slowdown -
                                    elapsedSim) > syncMisalign;

          // out-of-sync (for any reason): reset sync times, step
          if (elapsedSim < 0 || elapsedCPU.count() < 0 ||
              syncCPU.time_since_epoch().count() == 0 || misaligned ||
              sim.speed_changed) {
            // re-sync
            syncCPU = startCPU;
            syncSim = d->time;
            sim.speed_changed = false;

            // run single step, let next iteration deal with timing
            mj_step(m, d);
          }

          // in-sync: step until ahead of cpu
          else {
            bool measured = false;
            mjtNum prevSim = d->time;

            double refreshTime = simRefreshFraction / sim.refresh_rate;

            // step while sim lags behind cpu and within refreshTime
            while (Seconds((d->time - syncSim) * slowdown) <
                       mj::Simulate::Clock::now() - syncCPU &&
                   mj::Simulate::Clock::now() - startCPU <
                       Seconds(refreshTime)) {
              // measure slowdown before first step
              if (!measured && elapsedSim) {
                sim.measured_slowdown =
                    std::chrono::duration<double>(elapsedCPU).count() /
                    elapsedSim;
                measured = true;
              }

              // call mj_step
              mj_step(m, d);

              // break if reset
              if (d->time < prevSim) {
                break;
              }
            }
          }
        }

        // paused
        else {
          // run mj_forward, to update rendering and joint sliders
          mj_forward(m, d);
          syncCPU = mj::Simulate::Clock::now();
        }
      }
    } // release std::lock_guard<std::mutex>
  }
}

void NMpcLoop(mj::Simulate &sim) {
  while (!sim.exitrequest.load()) {
    if (!sim.uiloadrequest.load()) {
      {
        // std::chrono::time_point<std::chrono::system_clock> t_start =
        //     std::chrono::system_clock::now();
        if (QuadApp::instance()->FSM.nmpcUpdateNeeded) {
          QuadApp::instance()->FSM.computeNmpc();

          // std::chrono::time_point<std::chrono::system_clock> t_end =
          //     std::chrono::system_clock::now();
          // double time_record =
          //     std::chrono::duration_cast<std::chrono::milliseconds>(t_end -
          //                                                           t_start)
          //         .count();
          // std::cout << "nmpc_time: " << time_record / 1000 << "\n";
        } else
          std::this_thread::sleep_for(std::chrono::milliseconds(2));
      }
    }
  }
}

void MainLoop(mj::Simulate &sim) {
  while (!sim.exitrequest.load()) {
    if (!sim.uiloadrequest.load()) {
      // std::chrono::time_point<std::chrono::system_clock> t_start =
      //     std::chrono::system_clock::now();

      {
        const std::unique_lock<std::recursive_mutex> lock(g_mutex);

        QuadApp::instance()->FSM.stateEstimate(d->sensordata);
        QuadApp::instance()->FSM.mainProgram();

        const std::unique_lock<std::recursive_mutex> unlock(g_mutex);
      }

      std::this_thread::sleep_for(std::chrono::milliseconds(2));

      // std::chrono::time_point<std::chrono::system_clock> t_end =
      //     std::chrono::system_clock::now();
      // double time_record =
      //     std::chrono::duration_cast<std::chrono::milliseconds>(t_end -
      //     t_start)
      //         .count();
      // std::cout << "main loop time: " << time_record / 1000 << "\n";
    }
  }
}

void WbcLoop(mj::Simulate &sim) {
  while (!sim.exitrequest.load()) {
    if (!sim.uiloadrequest.load()) {
      // std::chrono::time_point<std::chrono::system_clock> t_start =
      //     std::chrono::system_clock::now();

      {
        if (QuadApp::instance()->FSM.wbcUpdateNeeded) {
          const std::unique_lock<std::recursive_mutex> lock(g_mutex);

          QuadApp::instance()->FSM.computeWbc();

          const std::unique_lock<std::recursive_mutex> unlock(g_mutex);
        }
      }

      std::this_thread::sleep_for(std::chrono::milliseconds(2));

      // std::chrono::time_point<std::chrono::system_clock> t_end =
      //     std::chrono::system_clock::now();
      // double time_record =
      //     std::chrono::duration_cast<std::chrono::milliseconds>(t_end -
      //     t_start)
      //         .count();
      // std::cout << "wbc_time: " << time_record / 1000 << "\n";
    }
  }
}

void app_controller(const mjModel *m, mjData *d) {
  {
    for (int i = 0; i < 4; ++i) {
      d->ctrl[i * 3 + 0] = QuadApp::instance()->FSM.getJointTorque()(i * 3 + 0);
      d->ctrl[i * 3 + 1] = QuadApp::instance()->FSM.getJointTorque()(i * 3 + 1);
      d->ctrl[i * 3 + 2] = QuadApp::instance()->FSM.getJointTorque()(i * 3 + 2);
    }
  }
}

void QuadApp::Start() {
  // --- init sim ---
  mjvScene scn;
  mjv_defaultScene(&scn);

  mjvCamera cam;
  mjv_defaultCamera(&cam);

  mjvOption opt;
  mjv_defaultOption(&opt);

  mjvPerturb pert;
  mjv_defaultPerturb(&pert);

  // --- simulate object encapsulates the UI ---
  char filename[] = "../unitree_a1/scene.xml";
  sim = std::make_unique<mj::Simulate>(std::make_unique<mj::GlfwAdapter>(),
                                       &scn, &cam, &opt, &pert, true);
  mju::strcpy_arr(sim->filename, filename);

  //  --- setup task threads ---
  ThreadPool main_program_pool(1);
  main_program_pool.Schedule([]() { MainLoop(*sim); });

  ThreadPool nmpc_pool(1);
  nmpc_pool.Schedule([]() { NMpcLoop(*sim); });

  ThreadPool wbc_pool(1);
  wbc_pool.Schedule([]() { WbcLoop(*sim); });

  mjcb_control = app_controller;

  ThreadPool physics_pool(1);
  physics_pool.Schedule([]() { PhysicsLoop(*sim); });

  sim->RenderLoop();
}
