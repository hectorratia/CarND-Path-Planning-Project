#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <vector>


using namespace std;

class trajectory {
  public:
    double dt;
    // x, xd, xdd
    vector<double> start;
    // x, xd, xdd
    vector<double> goal;

    vector<double> coefficients;

    double cost;

    //trajectory:

    trajectory();
    virtual ~trajectory();
    void set(vector<double> my_start, vector<double> my_goal, double my_dt);
    void solve();
    void calculate_cost();
    double get_point(double t);
};

#endif /* TRAJECTORY_H */
