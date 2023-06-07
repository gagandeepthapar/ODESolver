#include "odesolvers.hpp"

#include <iostream>
#include <vector>

using namespace std;

/*
 * ForwardEuler Definitions 
*/ 
ForwardEuler::ForwardEuler(vector<double>(*func)(float, vector<double>),
        vector<float> tspan,
        float t_step,
        vector<double> y0){

  m_func = func;
  m_t0 = tspan[0];
  m_tF = tspan[1];
  m_tstep = t_step;
  m_y0 = y0;

  solve();

}

void ForwardEuler::solve(){
  /* 
   * Step through in time until final time met
   */

  double t = m_t0;
  vector<double> newstate;
  vector<double> prevstate = m_y0;

  t_hist.push_back(t);
  state_hist.push_back(prevstate);
  
  while(t <= m_tF){
    newstate = step(t, prevstate, m_tstep);

    if(t == m_tF){
      return;
    }

    // check to ensure final time not exceeded
    if(t + m_tstep > m_tF){
      m_tstep = m_tF - t;
    }
    
    t = t + m_tstep;

    t_hist.push_back(t);
    state_hist.push_back(newstate);
    prevstate = newstate;

  }

}

vector<double> ForwardEuler::step(float t, vector<double> state, float dt){

  vector<double> dstate = m_func(t, state);
  vector<double> nstate;

  for(int i = 0; i < dstate.size(); i ++){
    nstate.push_back(dstate[i]*dt + state[i]);
  }

  return nstate;
}

/*
 * RKF45 Definitions 
*/
RKF45::RKF45(vector<double>(*func)(float, vector<double>),
        vector<float> tspan,
        float t_step,
        vector<double> y0){

  m_func = func;
  m_t0 = tspan[0];
  m_tF = tspan[1];
  m_tstep = t_step;
  m_y0 = y0;

}

void RKF45::solve(){
 // TODO: Implement 
}

vector<double> RKF45::step(float t, vector<double> state){
  // TODO: Implement
  return {0.0f, 0.0f};
}
