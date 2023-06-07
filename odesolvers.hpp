#ifndef ODESOLVERS_HPP 
#define ODESOLVERS_HPP 

#include <vector>

using namespace std;

// Abstract Base Class 
class ODESolver{

  virtual void solve() = 0;  
  virtual vector<double> step() = 0;

};

class RKF45:ODESolver{
 
  // variables
  public:
    vector<float> t_hist;
    vector< vector<double> > state_hist;

  private:
    vector<double>(*m_func)(float, vector<double>);
    vector<double> m_y0;
    float m_t0;
    float m_tF;

  // methods
  public:
    // constructor
    RKF45(vector<double>(*func)(float, vector<double>),
        vector<float> tspan,
        vector<double> y0);

    // destructor

  private:
    // run solver for given parameters
    void solve() override;

    // step
    vector<double> step() override;

};

#endif // !ODESOLVERS_HPP

