#ifndef ODESOLVERS_HPP 
#define ODESOLVERS_HPP 

#include <vector>

using namespace std;

// Abstract Base Class 
class ODESolver{

  void solve();  
  vector<double> step(float t, vector<double> state);

  // variables
  public:
    vector<float> t_hist;
    vector< vector<double> > state_hist;

  protected:
    vector<double>(*m_func)(float, vector<double>);
    vector<double> m_y0;
    float m_t0;
    float m_tF;
    float m_tstep;

};

class ForwardEuler:public ODESolver{
 
  public:
    // Constructor
    ForwardEuler(vector<double>(*func)(float, vector<double>),
        vector<float> tspan,
        float t_step,
        vector<double> y0);

  private:
    // run solver for given parameters
    void solve();

    // step
    vector<double> step(float t, vector<double> state, float dt);

};


class RKF45:public ODESolver{

  public:
    // Constructor
    RKF45(vector<double>(*func)(float, vector<double>),
        vector<float> tspan,
        float t_step,
        vector<double> y0);

  private:
    // run solver for given parameters
    void solve();

    // step
    vector<double> step(float t, vector<double> state);

};


#endif // !ODESOLVERS_HPP

