#include "odesolvers.hpp"
#include <vector>
#include <eigen3/Eigen/Core>

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

/*
 * Forward Euler Definitions
 */

// Constructor
ForwardEuler::ForwardEuler(Eigen::VectorXd(*func)(float, Eigen::VectorXd, Eigen::VectorXd),
      vector<float> tspan,
      Eigen::VectorXd y0,
      float t_step=1e-3f,
      double tol){

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
