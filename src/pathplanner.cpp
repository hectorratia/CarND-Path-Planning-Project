#include "pathplanner.h"
#include <iostream>

using namespace std;

pathplanner::pathplanner() {}
pathplanner::~pathplanner() {}

void pathplanner::update_telemetry(vector<double> ego_s, vector<vector<double> > sensors){
  ego_state = ego_s;
  sensor_fusion = sensors;
  ego_lane = ego_state[3] / 4;
}

void pathplanner::update_state(){
  vector<string> possible_states;
  possible_states.push_back("KL");

  //TODO evaluate possible states and choose the best
  state = possible_states[0];
}

void pathplanner::set_goal(){
  if (state=="KL") {
    realize_keep_lane();
  }
}

void pathplanner::get_trajectory(){
//TODO: calculate perturbations to the goal and choose the best trajectory
//use min_element
//TODO: for now just solving for the one goal
  vector<double> s_start(ego_state.begin(),ego_state.begin()+3);
  vector<double> s_goal(goal.begin(),goal.begin()+3);
  vector<double> d_start(ego_state.begin()+3,ego_state.begin()+6);
  vector<double> d_goal(goal.begin()+3,goal.begin()+6);

  s_traj.set(s_start,s_goal,dt);
  s_traj.solve();
  s_traj.calculate_cost();

  d_traj.set(d_start,d_goal,dt);
  d_traj.solve();
  d_traj.calculate_cost();
}

vector<double> pathplanner::get_point(double t){
  return {s_traj.get_point(t),ego_state[3]};
//      return {ego_state[0], ego_state[3]};
//    return {ego_state[0]+t*v_max, ego_state[3]};
//    return {s_traj.get_point(t),d_traj.get_point(t)};
}

void pathplanner::realize_keep_lane(){

  // goal.push_back(ego_state[0] + dt * v_max);
  // goal.push_back(v_max);
  // goal.push_back(0.0);
  // goal.push_back(ego_lane*4+2);
  // goal.push_back(0.0);
  // goal.push_back(0.0);

//Find nearest car in my lane driving in front of me
  int i_nearest=-1;
  double s_nearest=100000;
  for(int i=0; i<sensor_fusion.size(); i++){
    if((int)sensor_fusion[i][3]/4 == ego_lane){
      double s_difference = s_distance_to_ego(i);
      if(s_difference > 0 and s_difference < s_nearest){
        i_nearest = i;
        s_nearest = s_difference;
      }
    }
  }
  cout << "Ego is at s:" << ego_state[0] <<
          " vs:" << ego_state[1] <<
          " as:" << ego_state[2] <<
          " d:" << ego_state[3] <<
          " vd:" << ego_state[4] <<
          " ad:" << ego_state[5] <<
          " lane:" << ego_lane << endl;
  //no car in front of us...drive freely!
  if(i_nearest == -1){
      cout << "None ahead" << endl;
      double t_min_for_vmax = (v_max - ego_state[0]) / (coeff_a_max * a_max);
      if(t_min_for_vmax>dt) dt=t_min_for_vmax;
      double v_avg = (ego_state[1] + v_max)/2;

      goal.push_back(ego_state[0] + dt * v_avg);
      goal.push_back(v_max);
      goal.push_back(0.0);
      goal.push_back(ego_lane*4+2);
      goal.push_back(0.0);
      goal.push_back(0.0);
  }else{ //car in front...maybe need to adapt driving behaviour
    cout << "Other car id:" << sensor_fusion[i_nearest][0] <<
            " s:" << sensor_fusion[i_nearest][1] <<
            " vs:" << sensor_fusion[i_nearest][2] <<
            " d:" << sensor_fusion[i_nearest][3] <<
            " vd:" << sensor_fusion[i_nearest][4] << endl;
    cout << "Cars ahead " << s_nearest << endl;
      double safety_distance = max(safety_minimum,safety_seconds*sensor_fusion[i_nearest][2]);
      double s_diff = s_nearest-safety_distance;
      double v_diff = ego_state[1] - sensor_fusion[i_nearest][2];
      double t_min = s_diff / (coeff_v_diff_max * v_diff);
      cout << "safety distance:" << safety_distance <<
              " s_diff:" << s_diff <<
              " v_diff:" << v_diff <<
              " t_min:" << t_min << endl;

      if(t_min >0){
        cout << "We faster" << endl;
        //easy case: we drive faster than the other car and we have to catch up to be at
        //the safety distance using the N seconds rule
        // OR the other car is faster and the safety distance is not respected
        //The main factor in these calculations is time to respect the safety distnace
        if(t_min > dt) dt = t_min;
        goal.push_back(sensor_fusion[i_nearest][1] + dt * sensor_fusion[i_nearest][2] - safety_distance);
        goal.push_back(min(sensor_fusion[i_nearest][2],v_max));
        goal.push_back(0.0);
        goal.push_back(ego_lane*4+2);
        goal.push_back(0.0);
        goal.push_back(0.0);
      }else if(v_diff <= 0){
        //The main factor in the last two cases is acceleration. Either acceeleration to
        //reach max speed as in this case...
        cout << "Them faster" << endl;

        double v_goal = min(v_max, sensor_fusion[i_nearest][2]);
        double t_min_for_vmax = (v_goal - ego_state[0]) / (coeff_a_max * a_max);
        if(t_min_for_vmax>dt) dt=t_min_for_vmax;
        double v_avg = (ego_state[1] + v_goal)/2;
        goal.push_back(ego_state[0] + dt * v_avg);
        goal.push_back(v_goal);
        goal.push_back(0.0);
        goal.push_back(ego_lane*4+2);
        goal.push_back(0.0);
        goal.push_back(0.0);
        cout << "dt:" << dt <<
                " goal s:" << goal[0] <<
                " goal vs:" << goal[1] <<
                " goal as:" << goal[2] <<
                " goal d:" << goal[3] <<
                " goal vd:" << goal[4] <<
                " goal ad:" << goal[5] << endl;

      }else{ //s_diff < 0
        //or acceleration to adapt to the other car speed and distance.
        cout << "Safety distance!!!" << endl;

        double v_goal = sensor_fusion[i_nearest][2];
        double t_min_for_vgoal = (ego_state[0]-v_goal) / (coeff_a_max * a_max);
        double t_min_for_safety = s_diff / (coeff_v_fallback * v_goal);
        double t_min = t_min_for_vgoal + t_min_for_safety;
        if(t_min>dt) dt=t_min;
        goal.push_back(sensor_fusion[i_nearest][1] + dt * sensor_fusion[i_nearest][2] - safety_distance);
        goal.push_back(v_goal);
        goal.push_back(0.0);
        goal.push_back(ego_lane*4+2);
        goal.push_back(0.0);
        goal.push_back(0.0);
      }
  }
}

double pathplanner::s_distance_to_ego(int i_sensor){
  double max_s = 6945.554;

  double difference = sensor_fusion[i_sensor][1] - ego_state[0];
  if(sensor_fusion[i_sensor][1] < ego_state[0]){
    difference += max_s;
  }
  return difference;
}
