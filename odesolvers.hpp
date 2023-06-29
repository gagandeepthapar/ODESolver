#ifndef ODESOLVERS_HPP 
#define ODESOLVERS_HPP 

#include <vector>
#include <eigen3/Eigen/Core>

using namespace std;

// Abstract base class
class ODESolver{
  //makes class abstract  
  virtual void f() = 0;

// functions
protected:
  Eigen::VectorXd step(float time, Eigen::VectorXd state, Eigen::VectorXd args);

public:
  Eigen::VectorXd solve(Eigen::VectorXd args);

// member variables
protected:
  Eigen::VectorXd(*m_func)(float, Eigen::VectorXd, Eigen::VectorXd);
  Eigen::VectorXd m_y0;
  int m_num_state;

  float m_t0;
  float m_tF;
  float m_tstep;
  double m_tol;

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

