#ifndef ODESOLVERS_HPP 
#define ODESOLVERS_HPP 

#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

using namespace std;

// Abstract base class
class ODESolver{

// functions
protected:
  // step method computes x_(t+1) from x_(t), args
  virtual Eigen::VectorXf step(float time, Eigen::VectorXf state, Eigen::VectorXf args) = 0; 

public:
  // solve method solves the system until t_f
  Eigen::VectorXf solve(Eigen::VectorXf args);

  // get state history if length of integration is known
  // not preferred due to memory requirements
  vector<Eigen::VectorXf> solve_hist(Eigen::VectorXf args);

// member variables
protected:
  Eigen::VectorXf(*m_func)(float, Eigen::VectorXf, Eigen::VectorXf);
  Eigen::VectorXf m_y0; // initial state

  float m_t0; // initial time 
  float m_tF; // final time
  float m_tstep; // time-step (default timestep)
  double m_tol; // tolerance 

public:
  int m_num_state; // number of states in system
  const string name = "Abstract ODE Solver"; // name of solver

};


class ForwardEuler:public ODESolver{
public:
  // Constructor
  ForwardEuler(Eigen::VectorXd(*func)(float, Eigen::VectorXd, Eigen::VectorXd),
      vector<float> tspan,
      Eigen::VectorXd y0,
      float t_step,
      double tol=1e-8);

};


class RKF45:public ODESolver{

public:
  // Constructor
  RKF45(Eigen::VectorXd(*func)(float, Eigen::VectorXd, Eigen::VectorXd),
      vector<float> tspan,
      float t_step,
      Eigen::VectorXd y0,
      double tol=1e-8);

};

#endif // !ODESOLVERS_HPP

