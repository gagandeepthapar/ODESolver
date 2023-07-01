#include "odesolvers.hpp"
#include <vector>
#include <eigen3/Eigen/Core>
#include <iostream>

using namespace std;

/*
 * ODE Definitions
 */ 

Eigen::VectorXd ODESolver::solve(Eigen::VectorXd args){
 
  // unpack data
  double t = (double)this->m_t0;
  double dt = this->m_tstep;
  Eigen::VectorXd cur_state = this->m_y0;

  // iterate through timespan 
  while(t < this->m_tF){
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

}

// 

/*
* Runge-Kutta-Fehlberg (RKF) Definitions
*/
