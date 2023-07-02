#include "ode_solver.hpp"
#include "eigen3/Eigen/Core"
#include "matplot/matplot.h"

#include <iostream>

using namespace std;
namespace plt = matplot;

Eigen::VectorXf vanderpol(float time, Eigen::VectorXf state, Eigen::VectorXf args){

  // Vander Pol Equation:
  // x_ddot = mu(1 - x^2)x_dot - x

  float x = state[0];
  float xdot = state[1];

  float d_x = xdot;
  float d_xdot = args[0] * (1.0 - x*x)*xdot - x;

  Eigen::VectorXf dstate {{d_x, d_xdot}};

  return dstate;

}

Eigen::VectorXf lorenz(float time, Eigen::VectorXf state, Eigen::VectorXf args){

  // Lorenz Attractor Equations 
  // xDot = mu * (y - x)
  // yDot = x * (rho - z) - y 
  // zDot = x*y - beta * z 
  
  float mu, rho, beta;
  mu = args[0];
  rho = args[1];
  beta = args[2];

  float x, y, z;
  x = state[0];
  y = state[1];
  z = state[2];

  float d_x, d_y, d_z;
  d_x = mu * (y - x);
  d_y = x * (rho - z) - y;
  d_z = x * y - beta * z;

  Eigen::VectorXf d_state {{d_x, d_y, d_z}};
  return d_state;

}

int main(){
  // sandbox
  
  // initial conditions
  Eigen::VectorXf y0{{1.0f, 2.0f, 3.0f}};

  // instantiate ode solver
  // ForwardEuler ode(lorenz, {0.0f, 50.0f}, y0, 0.01f, 1e-8);
  RKF45 ode(lorenz, {0.0f, 50.0f}, y0, 0.01f, 1e-8);

  // solve ode with arg
  Eigen::VectorXf arg{{10.0f, 28.0f, 8.0f/3.0f, 0.0f}};
  vector<Eigen::VectorXf> hist = ode.solve_hist(arg);

  // plot
  vector<float> x, y, z;
  for(int i = 0; i < hist.size(); i++){
    x.push_back(hist[i][0]);
    y.push_back(hist[i][1]);
    z.push_back(hist[i][2]);
  }

 
  auto p = plt::plot3(x, y, z);
  plt::title(ode.name);
  plt::show();
  
  return 0;
}
