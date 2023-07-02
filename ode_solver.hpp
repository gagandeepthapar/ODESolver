#ifndef ODE_SOLVER_HPP 
#define ODE_SOLVER_HPP 

#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

// Abstract base class
class ODESolver{

// functions
protected:
  // step method computes x_(t+1) from x_(t), args
  virtual Eigen::VectorXf step(float time, Eigen::VectorXf state, Eigen::VectorXf args) = 0; 

public:
  // solve method solves the system until t_f
  Eigen::VectorXf solve(Eigen::VectorXf args);
  Eigen::VectorXf solve(Eigen::VectorXf state,
                        std::vector<float> tspan,
                        Eigen::VectorXf args);

  // get state history if length of integration is known
  // not preferred due to memory requirements
  std::vector<Eigen::VectorXf> solve_hist(Eigen::VectorXf args);

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
  const std::string name = "Abstract ODE Solver"; // name of solver

};


class ForwardEuler:public ODESolver{
public:
  // Constructor
  ForwardEuler(Eigen::VectorXf(*func)(float, Eigen::VectorXf, Eigen::VectorXf),
      std::vector<float> tspan,
      Eigen::VectorXf y0,
      float t_step,
      double tol);

  // override step method
  Eigen::VectorXf step(float time, Eigen::VectorXf state, Eigen::VectorXf args) override;

  // override name
  const std::string name = "Forward Euler Solver";

};


class RKF45:public ODESolver{

public:
  // Constructor
  RKF45(Eigen::VectorXf(*func)(float, Eigen::VectorXf, Eigen::VectorXf),
      std::vector<float> tspan,
      Eigen::VectorXf y0,
      float t_step,
      double tol);

  // override step method
  Eigen::VectorXf step(float time, Eigen::VectorXf state, Eigen::VectorXf args) override; 

  // override name 
  const std::string name = "RKF45 Solver";

};

#endif // !ODESOLVERS_HPP

