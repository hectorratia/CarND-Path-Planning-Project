#ifndef PATHPLANNER_H
#define PATHPLANNER_H

#include <vector>
#include "constants.h"
#include "trajectories.h"
#include "Eigen-3.3/Eigen/Core"


using namespace std;

class pathplanner {
public:
    double dt=1;
    // s, sd, sdd, d, dd, ddd
    vector<double> goal;
    // s, sd, sdd, d, dd, ddd
    vector<double> ego_state;
    // id, s, s_d, d, d_d, lane
    int ego_lane;

    vector<vector<double> > sensor_fusion;
    //vector
    string state = "KL";

    //trajectory:
    trajectory s_traj;
    trajectory d_traj;

    pathplanner();
    virtual ~pathplanner();
    void update_telemetry(vector<double> ego_s, vector<vector<double> > sensors);
    void update_state();
    void set_goal();
    void get_trajectory();
    vector<double> get_point(double t);
    void realize_keep_lane();
    double s_distance_to_ego(int i_sensor);
    //update_state(vector<double> ego_state, vector< vector<double>> predictions);
//  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, Eigen::VectorXd weights);
};

#endif /* PATHPLANNER_H */
