#include "ode_solver.hpp"
#include <vector>
#include <eigen3/Eigen/Core>
#include <iostream>

using namespace std;

/*
 * ODE Definitions
 */ 

// solve using class params
Eigen::VectorXf ODESolver::solve(Eigen::VectorXf args){
  return ODESolver::solve(this->m_y0,
                          {this->m_t0, this->m_tF},
                          args);
}

// solve using specified params 
Eigen::VectorXf ODESolver::solve(Eigen::VectorXf state,
                                 std::vector<float> tspan,
                                 Eigen::VectorXf args){
 
  // unpack data
  float t = tspan[0];
  float dt = this->m_tstep;
  Eigen::VectorXf cur_state = state;

  // iterate through timespan 
  while(t < tspan[1]){
    // step through based on integration scheme 
    cur_state = step(t, cur_state, args);

    // increment t such that tF is never violated 
    dt = (dt < (this->m_tF - t)) ? dt : (this->m_tF - t);
    t = t + dt;
  }

  // return final calculated state
  return cur_state;
}

// solve hist method 
vector<Eigen::VectorXf> ODESolver::solve_hist(Eigen::VectorXf args){
 
  int t_slots = this->m_tF / this->m_tstep + 2;
  int t_idx = 0;
 
  vector<Eigen::VectorXf> hist;
  
   // unpack data
  float t = (double)this->m_t0;
  float dt = this->m_tstep;
  Eigen::VectorXf cur_state = this->m_y0;
  hist.push_back(cur_state);

  // iterate through timespan 
  while(t <= this->m_tF){
    // step through based on integration scheme 
    cur_state = step(t, cur_state, args);

    // store in eigen matrix
    for(int i = 0; i < this->m_num_state; i++){
      hist.push_back(cur_state);
    }
    
    // increment t such that tF is never violated 
    if(t == this->m_tF){
      return hist;
    }
    else{
      dt = (dt < (this->m_tF - t)) ? dt : (this->m_tF - t);
      t += dt;
      t_idx++; 
    }
  }
  return hist;
}

/*
 * Forward Euler Definitions
 */

// Constructor
ForwardEuler::ForwardEuler(Eigen::VectorXf(*func)(float, Eigen::VectorXf, Eigen::VectorXf),
      vector<float> tspan,
      Eigen::VectorXf y0,
      float t_step=1e-3f,
      double tol=1e-8){

  this->m_func = func;
  this->m_y0 = y0;
  this->m_t0 = tspan[0];
  this->m_tF = tspan[1];
  this->m_tstep = t_step;
  this->m_tol = m_tol;
  this->m_num_state = y0.rows();

}

// step method 
Eigen::VectorXf ForwardEuler::step(float time, Eigen::VectorXf state, Eigen::VectorXf args){
  /* Forward Euler (Naive) method of integration:
   *    x_(t+1) = x(t) + dt * xdot(t)
   */
  Eigen::VectorXf dstate = this->m_func(time, state, args);
  return (dstate * this->m_tstep) + state;

}


/*
 * Runge-Kutta-Fehlberg (RKF) Definitions
 */

RKF45::RKF45(Eigen::VectorXf(*func)(float, Eigen::VectorXf, Eigen::VectorXf),
      vector<float> tspan,
      Eigen::VectorXf y0,
      float t_step=1e-3f,
      double tol=1e-8){

  this->m_func = func;
  this->m_y0 = y0;
  this->m_t0 = tspan[0];
  this->m_tF = tspan[1];
  this->m_tstep = t_step;
  this->m_tol = m_tol;
  this->m_num_state = y0.rows();

}

// step method 
Eigen::VectorXf RKF45::step(float time, Eigen::VectorXf state, Eigen::VectorXf args){
  /* RKF45 (ode45) method of integration for 5th order error with variable step size
   * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods    
   */

  // shorten type name
  using vec = Eigen::VectorXf;
  Eigen::VectorXf(*f)(float, Eigen::VectorXf, Eigen::VectorXf) = this->m_func;
  float dt = this->m_tstep;

  // calc constituent parameters
  vec s_1 = f(time,
              state,
              args);
  vec s_2 = f(time + 0.25f*dt,
              state + 0.25f * dt * s_1,
              args);
  vec s_3 = f(time + 0.375f*dt,
              state + 3.0f/32.0f*dt*s_1 + 9.0f/32.0f*dt*s_2,
              args);
  vec s_4 = f(time + 12.0f/13.0f*dt,
              state + 1932.0f/2197.0f*dt*s_1 - 7200.0f/2197.0f*dt*s_2 + 7296.0f/2197.0f*dt*s_3,
              args);
  vec s_5 = f(time + dt,
              state + 439.0f/216.0f*dt*s_1 - 8.0f*dt*s_2 + 3680.0f/513.0f*dt*s_3 - 845.0f/4104.0f*dt*s_4,
              args);
  vec s_6 = f(time + 0.5f*dt,
              state - 8.0f/27.0f*dt*s_1 + 2.0f*dt*s_2 - 3544.0f/2565.0f*dt*s_3 + 1859.0f/4104.0f*dt*s_4 - 11.0f/40.0f*dt*s_5,
              args);

  // calc x_(t+1)
  vec new_state = state + dt * (16.0f/135.0f*s_1 + 
                                6656.0f/12825.0f*s_3 + 
                                28561.0f/56430.0f*s_4 + 
                                -9.0f/50.0f*s_5 +
                                2.0f/55.0f*s_6);

  float err = (1.0f/360.0f*s_1 + 
             -128.0f/4275.0f*s_3 + 
             -2197.0f/75240.0f*s_4 + 
             1.0f/50.0f * s_5 + 
             2.0f / 55.0f * s_6).norm() * dt; 

  // std::cout << "ERROR: " << err << endl;

  return new_state;

}

