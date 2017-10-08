#include "trajectories.h"
#include <iostream>

#include "Eigen-3.3/Eigen/Dense"

trajectory::trajectory() {}
trajectory::~trajectory() {}

using namespace std;

void trajectory::set(vector<double> my_start, vector<double> my_goal, double my_dt){
  start = my_start;
  goal = my_goal;
  dt = my_dt;
  coefficients.clear();
  cost=0;
}

void trajectory::solve(){
  double a_0 = start[0];
  double a_1 = start[1];
  double a_2 = start[2]/2.0;

  double c_0 = a_0 + a_1 * dt + a_2 * pow(dt,2);
  double c_1 = a_1 + 2* a_2 * dt;
  double c_2 = 2 * a_2;

  Eigen::MatrixXd A(3,3);
  Eigen::VectorXd b(3);
  Eigen::VectorXd solution(3);
  A <<   pow(dt,3),    pow(dt,4),    pow(dt,5),
       3*pow(dt,2),  4*pow(dt,3),  5*pow(dt,4),
              6*dt, 12*pow(dt,2), 20*pow(dt,3);

  b << goal[0] - c_0,
       goal[1] - c_1,
       goal[2] - c_2;

  solution = A.colPivHouseholderQr().solve(b);

  coefficients.push_back(a_0);
  coefficients.push_back(a_1);
  coefficients.push_back(a_2);
  coefficients.push_back(solution(0));
  coefficients.push_back(solution(1));
  coefficients.push_back(solution(2));

  cout << "c0:" << coefficients[0] <<
          " c1:" << coefficients[1] <<
          " c2:" << coefficients[2] <<
          " c3:" << coefficients[3] <<
          " c4:" << coefficients[4] <<
          " c5:" << coefficients[5] << endl;
}

void trajectory::calculate_cost(){
  //TODO
  cost = 0;
}

double trajectory::get_point(double t){
  double x=0;

  for(int i=0; i<coefficients.size(); i++){
      x += coefficients[i] * pow(t,i);
  }

  return x;
}
