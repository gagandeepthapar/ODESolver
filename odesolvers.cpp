#include "odesolvers.hpp"

#include <iostream>
#include <vector>

using namespace std;

// RKF-45 Definitions
RKF45::RKF45(vector<double>(*func)(float, vector<double>),
             vector<float> tspan,
             vector<double> y0):
            
             m_func(func),
             m_t0(tspan[0]),
             m_tF(tspan[1]),
             m_y0(y0){
  solve();
}

void RKF45::solve(){
  for(int i = 0; i < (int)m_tF; i++){
    t_hist.push_back(i);
  }
}

vector<double> RKF45::step(){

  vector<double> data = {1.0f, 2.0f};
  return data;
}
