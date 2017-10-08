#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "trajectories.h"
#include "spline.h"
#include "constants.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

  if(closestWaypoint==maps_x.size()) closestWaypoint = 0;

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenetS(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y, vector<double> maps_s)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double frenet_s = maps_s[prev_wp] + proj_norm * (maps_s[next_wp]-maps_s[prev_wp]);
  if(next_wp == 0)
  {
    frenet_s  = (maps_s[prev_wp]) + proj_norm * (maps_s[next_wp]+s_max-maps_s[prev_wp]);
  }
	//double proj_x = proj_norm*n_x;
	//double proj_y = proj_norm*n_y;

  vector<double> pts;
  vector<double> ptsx;
  vector<double> ptsy;
  int wp0 = (prev_wp-2)%maps_x.size();
  int wp1 = (prev_wp-1)%maps_x.size();
  int wp2 = prev_wp;
  int wp3 = (prev_wp+1)%maps_x.size();
  int wp4 = (prev_wp+2)%maps_x.size();
  int wp5 = (prev_wp+3)%maps_x.size();

  double add_s = 0;
  pts.push_back(maps_s[wp0]);
  if(wp1 < 1){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp1]+add_s);
  if(wp2 < 2){add_s = s_max; frenet_s += s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp2]+add_s);
  if(wp3 < 3){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp3]+add_s);
  if(wp4 < 4){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp4]+add_s);
  if(wp5 < 5){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp5]+add_s);

  ptsx.push_back(maps_x[wp0]);
  ptsx.push_back(maps_x[wp1]);
  ptsx.push_back(maps_x[wp2]);
  ptsx.push_back(maps_x[wp3]);
  ptsx.push_back(maps_x[wp4]);
  ptsx.push_back(maps_x[wp5]);

  ptsy.push_back(maps_y[wp0]);
  ptsy.push_back(maps_y[wp1]);
  ptsy.push_back(maps_y[wp2]);
  ptsy.push_back(maps_y[wp3]);
  ptsy.push_back(maps_y[wp4]);
  ptsy.push_back(maps_y[wp5]);

  tk::spline spline_x_s;
  tk::spline spline_y_s;
  spline_x_s.set_points(pts, ptsx);
  spline_y_s.set_points(pts, ptsy);

  double proj_x = spline_x_s(frenet_s)-maps_x[prev_wp];
  double proj_y = spline_y_s(frenet_s)-maps_y[prev_wp];

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

  if(frenet_s>s_max) frenet_s-=s_max;
	return {frenet_s,frenet_d};
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenetSpeed(double x, double y, double theta, double speed, vector<double> maps_x, vector<double> maps_y, vector<double> maps_dx, vector<double> maps_dy)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double prev_yaw_d = atan2(maps_dy[prev_wp],maps_dx[prev_wp]);
  double next_yaw_d = atan2(maps_dy[next_wp],maps_dx[next_wp]);

  double prev_dist = distance(x,y,maps_x[prev_wp],maps_y[prev_wp]);
  double next_dist = distance(x,y,maps_x[next_wp],maps_y[next_wp]);

  double yaw_d = prev_yaw_d + (next_yaw_d-prev_yaw_d)*prev_dist/(prev_dist+next_dist);
  double yaw_s = yaw_d + pi()/2;

  double d_dot = speed * cos(theta - yaw_d);
  double s_dot = speed * cos(theta - yaw_s);

	return {s_dot,d_dot};
}


// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXYspline(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y, vector<double> maps_dx, vector<double> maps_dy)
{

	int prev_wp = -1;
  if(s>s_max) s-=s_max;
  if(s<0) s+= s_max;
	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp0 = (prev_wp-2)%maps_x.size();
  int wp1 = (prev_wp-1)%maps_x.size();
  int wp2 = prev_wp;
  int wp3 = (prev_wp+1)%maps_x.size();
  int wp4 = (prev_wp+2)%maps_x.size();
  int wp5 = (prev_wp+3)%maps_x.size();

  vector<double> pts;

  vector<double> ptsx;
  vector<double> ptsy;

  vector<double> ptsdx;
  vector<double> ptsdy;

  double add_s = 0;
  pts.push_back(maps_s[wp1]);
  if(wp2 < 1){add_s = s_max;s+=s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp2]+add_s);
  if(wp3 < 2){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp3]+add_s);
  if(wp4 < 3){add_s = s_max;}else{add_s=0;}
  pts.push_back(maps_s[wp4]+add_s);

  ptsx.push_back(maps_x[wp1]);
  ptsx.push_back(maps_x[wp2]);
  ptsx.push_back(maps_x[wp3]);
  ptsx.push_back(maps_x[wp4]);

  ptsy.push_back(maps_y[wp1]);
  ptsy.push_back(maps_y[wp2]);
  ptsy.push_back(maps_y[wp3]);
  ptsy.push_back(maps_y[wp4]);

  ptsdx.push_back(maps_dx[wp1]);
  ptsdx.push_back(maps_dx[wp2]);
  ptsdx.push_back(maps_dx[wp3]);
  ptsdx.push_back(maps_dx[wp4]);

  ptsdy.push_back(maps_dy[wp1]);
  ptsdy.push_back(maps_dy[wp2]);
  ptsdy.push_back(maps_dy[wp3]);
  ptsdy.push_back(maps_dy[wp4]);

  tk::spline spline_x_s;
  tk::spline spline_dx_s;
  tk::spline spline_y_s;
  tk::spline spline_dy_s;
  spline_x_s.set_points(pts, ptsx);
  spline_dx_s.set_points(pts, ptsdx);
  spline_y_s.set_points(pts, ptsy);
  spline_dy_s.set_points(pts, ptsdy);

  double x = spline_x_s(s) + d * spline_dx_s(s);
  double y = spline_y_s(s) + d * spline_dy_s(s);

	return {x,y};
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object
          /*double temp = 6850;
          while(temp<7100){
            vector<double> next_point=getXYspline(temp,6.0,map_waypoints_s,map_waypoints_x,map_waypoints_y,map_waypoints_dx,map_waypoints_dy);
            cout << temp << " " << next_point[0] << " " << next_point[1] << endl;
            temp += 20;
          }*/
          /*vector<double> temp_point=getFrenetS(747.018,1130.97,-0.0978821,map_waypoints_x,map_waypoints_y,map_waypoints_s);
          cout << "temp1 " << temp_point[0] << " " << temp_point[1] << endl;
          vector<double> temp_point2=getFrenetS(745.266,1131.14,-0.104734,map_waypoints_x,map_waypoints_y,map_waypoints_s);
          cout << "temp2 " << temp_point2[0] << " " << temp_point2[1] << endl;*/

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            double car_acc = 0;
            int car_lane = car_d / 4;
            if(car_lane<0) car_lane = 0;
            if(car_lane>2) car_lane = 2;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            // The data format for each car is:
            // [ id, x, y, vx, vy, s, d].
            // The id is a unique identifier for that car.
            // The x, y values are in global map coordinates,
            // and the vx, vy values are the velocity components,
            // also in reference to the global map.
            // Finally s and d are the Frenet coordinates for that car.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            car_yaw = deg2rad(car_yaw);
            car_speed *= 0.44704;
          	json msgJson;

            // TODO: prepare data for ego
            vector<double> ego_state;
            int previous_size = previous_path_x.size();
            //take maximum of 5 states...we are giving ourselves a maximum of 100ms
            //to reply with new trajectory

            double elapsed_time = 0.0;

            bool allow_lane_change = true;
            if(previous_size > 2){
              for(int i=1; i<previous_path_x.size(); i++){
                double path_x_0 = previous_path_x[i-1];
                double path_y_0 = previous_path_y[i-1];
                double path_x_1 = previous_path_x[i];
                double path_y_1 = previous_path_y[i];
                double path_yaw_1 = atan2(path_y_1-path_y_0,path_x_1-path_x_0);
                vector<double> car_frenet = getFrenetS(path_x_1,path_y_1,path_yaw_1,map_waypoints_x,map_waypoints_y,map_waypoints_s);
                int path_lane = car_frenet[1]/4.0;
                if(path_lane != car_lane) allow_lane_change = false;
              }
            }

            if(previous_size > 2){//no data, we just have position and speed
              elapsed_time = previous_size * 0.02;

              double path_x_0 = previous_path_x[previous_size-3];
              double path_y_0 = previous_path_y[previous_size-3];
              double path_x_1 = previous_path_x[previous_size-2];
              double path_y_1 = previous_path_y[previous_size-2];
              double path_x_2 = previous_path_x[previous_size-1];
              double path_y_2 = previous_path_y[previous_size-1];

              double path_yaw_1 = atan2(path_y_1-path_y_0,path_x_1-path_x_0);
              double path_distance_1 = sqrt(pow(path_x_1-path_x_0,2)+pow(path_y_1-path_y_0,2));
              double path_speed_1 = path_distance_1 / 0.02;

              double path_yaw_2 = atan2(path_y_2-path_y_1,path_x_2-path_x_1);
              double path_distance_2 = sqrt(pow(path_x_2-path_x_1,2)+pow(path_y_2-path_y_1,2));
              double path_speed_2 = path_distance_2 / 0.02;

              car_x = path_x_2;
              car_y = path_y_2;
              car_yaw = path_yaw_2;
              vector<double> car_frenet = getFrenetS(car_x,car_y,car_yaw,map_waypoints_x,map_waypoints_y,map_waypoints_s);
              car_s = car_frenet[0];
              car_d = car_frenet[1];
              car_speed = path_speed_2;
              car_acc = (path_speed_2 - path_speed_1) / 0.02;
            }

            //Prepare telemetry data
            vector<vector<double> > sensors;
            for(int i=0; i<sensor_fusion.size(); i++){
              double other_car_vx = sensor_fusion[i][3];
              double other_car_vy = sensor_fusion[i][4];

              //double other_car_x = sensor_fusion[i][1];
              //other_car_x += 0.02 * other_car_vx * previous_to_copy;
              //double other_car_y = sensor_fusion[i][2];

              //other_car_y += 0.02 * other_car_vy * previous_to_copy;
              double other_car_speed = sqrt(other_car_vx*other_car_vx+other_car_vy*other_car_vy);
              //double other_car_yaw = atan2(other_car_vy,other_car_vx);

              //vector<double> other_car_frenet = getFrenetS(other_car_x,
              //other_car_y,other_car_yaw,map_waypoints_x,map_waypoints_y,map_waypoints_s);
              double other_car_s = sensor_fusion[i][5];
              other_car_s += other_car_speed * elapsed_time;
              double other_car_d = sensor_fusion[i][6];//OVERWRITE with value from sim
              //vector<double> other_car_frenet_speed = getFrenetSpeed(other_car_x,
              //other_car_y,other_car_yaw,other_car_speed,map_waypoints_x,map_waypoints_y,
              //map_waypoints_dx,map_waypoints_dy);

              vector<double> other_car_sensors = {sensor_fusion[i][0],
                other_car_s,other_car_speed,
                other_car_d};
              sensors.push_back(other_car_sensors);
            }

            //Find nearest car in this lane
            vector<int> i_nearest = {-1,-1,-1};
            vector<double> s_nearest = {6945.554,6945.554,6945.554};
            vector<double> s_behind = {6945.554,6945.554,6945.554};

            for(int i=0; i<sensors.size(); i++)
            {
              //cout << "Other car in lane " << (sensors[i][3]/4 << endl;
              int lane = sensors[i][3]/4;
              if(lane < 0) lane = 0;
              if(lane > 2) lane = 2;
              double s_diff = sensors[i][1] - car_s;
              double behind = car_s - sensors[i][1];
              if(s_diff<0) s_diff += 6945.554;
              if(behind<0) behind += 6945.554;
              if(s_diff<s_nearest[lane]){
                  s_nearest[lane]=s_diff;
                  i_nearest[lane]=i;
              }
              if(behind<s_behind[lane]){
                  s_behind[lane]=behind;
              }
            }
            cout << "In front: ";
            for(int i=0; i<s_nearest.size(); i++) cout << s_nearest[i] << " ";
            cout << endl;
            cout << "Behind: ";
            for(int i=0; i<s_behind.size(); i++) cout << s_behind[i] << " ";
            cout << endl;

            //TODO Behaviour logic
            double advance_s = 100;

            double lane_speed[] = {0,0,0};
            for(int i=0; i<3; i++){
              if(i_nearest[i]==-1){
                lane_speed[i]=v_max;
              }else{
                double safe_nearest = max(sensors[i_nearest[i]][2]*safety_seconds,safety_minimum);
                lane_speed[i] = sensors[i_nearest[i]][2] * (s_nearest[i]/safe_nearest);
                lane_speed[i] = min(v_max,lane_speed[i]);
              }
            }


            int fastest_lane = car_lane;
            double fastest_speed = lane_speed[car_lane];
            int left_lane = car_lane - 1;
            int right_lane = car_lane + 1;
            bool right_lane_ok = false;
            bool left_lane_ok = false;
            if(left_lane>=0 && lane_speed[left_lane]>fastest_speed
              && s_nearest[left_lane]>sensors[i_nearest[left_lane]][2]*safety_seconds && s_behind[left_lane]>safety_minimum){
//              fastest_speed=lane_speed[left_lane];
              fastest_lane=left_lane;
              left_lane_ok = true;
            }
            if(right_lane<=2 && lane_speed[right_lane]>fastest_speed
              && s_nearest[right_lane]>sensors[i_nearest[right_lane]][2]*safety_seconds && s_behind[right_lane]>safety_minimum){
//              fastest_speed=lane_speed[right_lane];
              fastest_lane=right_lane;
              right_lane_ok=true;
            }
            if(left_lane_ok && right_lane_ok){
              if(s_nearest[left_lane]>s_nearest[right_lane]){//} or lane_speed[left_lane]>lane_speed[right_lane]){
                fastest_lane = left_lane;
              }else{
                fastest_lane = right_lane;
              }
            }
            fastest_speed = lane_speed[fastest_lane];

            double v_ref = lane_speed[car_lane];
            double advance_d = 0;
            double offset = 2.0;
            if(allow_lane_change){
              if(car_s>4800 and car_s < 5200 and fastest_lane==2) offset = 1.6;
              advance_d = fastest_lane*4.0+offset-car_d;
            }else{
              if(car_s>4800 and car_s < 5200 and car_lane==2) offset = 1.6;
              advance_d = car_lane*4.0+offset-car_d;
            }
            cout << "s: " << car_s << " d: " << car_d <<" lane: " << car_lane << " speed: " << car_speed << " Fast speed: " << fastest_speed << " ALLOW: " << allow_lane_change << endl;

            // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            vector<double> next_x_vals;
          	vector<double> next_y_vals;

            for(int i=0; i<previous_path_x.size(); i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }
            bool changing_lanes = false;
            if(previous_path_x.size()<25){
              if(fastest_lane!=car_lane && allow_lane_change){
                cout << "CHANGING LANES" << endl;
                cout << "Space in front: " << s_nearest[fastest_lane] << " safe: " << sensors[i_nearest[fastest_lane]][2]*safety_seconds << endl;
                cout << "Space behind: " << s_behind[fastest_lane] << " safe: " << safety_minimum << endl;
                changing_lanes = true;
              }

              vector<double> ptsx;
              vector<double> ptsy;

              ptsx.push_back(-5);
              ptsy.push_back(0);

              ptsx.push_back(0);
              ptsy.push_back(0);
              //cout << "car_x: " << car_x << " car_y: " << car_y << endl;
              int n_points_advance = 8;

              int n_points_lateral = 3;

              double p_lateral[] = {0.25,0.65,0.9,1.0,1.0,1.0,1.0,1.0};
              if(!changing_lanes){
                p_lateral[0] = 1.0;
                p_lateral[1] = 1.0;
                p_lateral[2] = 1.0;
              }
              /*
              if(fastest_lane!=car_lane){
                double t_traj = (advance_s/n_points_advance) / car_speed;
                trajectory traj_d;
                traj_d.set({car_d,0,0},{car_d+advance_d,0,0},advance_s/n_points_advance);
                traj_d.solve();
                for(int i=1; i<n_points_advance; i++){
                  double n_advance_s = (advance_s/n_points_advance) * i / n_points_advance;
                  double next_s = car_s + n_advance_s;
                  double next_d = traj_d.get_point(n_advance_s);

                  vector<double> next_point=getXY(next_s,next_d,map_waypoints_s,map_waypoints_x,map_waypoints_y);
                  //cout << "next s: " << next_s << " next_d: " << next_d << " next_x: " << next_point[0] << " next_y: " << next_point[1] << endl;

                  next_point[0] -= car_x;
                  next_point[1] -= car_y;
                  ptsx.push_back(next_point[0]*cos(-car_yaw) - next_point[1] * sin (-car_yaw));
                  ptsy.push_back(next_point[0]*sin(-car_yaw) + next_point[1] * cos (-car_yaw));
                }
              }*/
              cout << "car x: " << car_x << "car y: " << car_y << " car yaw: " << car_yaw << endl;
              for(int i=1; i<=n_points_advance; i++){
                double next_s = car_s + advance_s * i / n_points_advance;
                double add_d = advance_d * p_lateral[i-1];
                //double add_d = advance_d * min(1.0,((double)i)/n_points_lateral);
                //double add_d = 0;
                //if(i>= n_points_lateral) add_d = advance_d;
                //if(fastest_lane=car_lane) add_d = advance_d;
                double next_d = car_d + add_d;

                vector<double> next_point=getXYspline(next_s,next_d,map_waypoints_s,map_waypoints_x,map_waypoints_y,map_waypoints_dx,map_waypoints_dy);
                cout << "next s: " << next_s << " next_d: " << next_d << " next_x: " << next_point[0] << " next_y: " << next_point[1] << endl;

                next_point[0] -= car_x;
                next_point[1] -= car_y;
                ptsx.push_back(next_point[0]*cos(-car_yaw) - next_point[1] * sin (-car_yaw));
                ptsy.push_back(next_point[0]*sin(-car_yaw) + next_point[1] * cos (-car_yaw));
              }
              /*for(int i=0; i<ptsx.size(); i++){
                cout << ptsx[i] << " " << ptsy[i] << endl;
              }*/

              tk::spline s;
              for(int i=0; i<ptsx.size()-1; i++){
                cout << "Traj " << i << " " << ptsx[i] << " " << ptsy[i] << endl;
              }
              s.set_points(ptsx,ptsy);
              double x_ref = 0;
              double y_ref = 0;
              double a_change = car_acc;
              bool continue_traj = true;
              while(continue_traj){
                double x_ahead = x_ref + 1;
                double y_ahead = s(x_ahead);

                double distance_ahead = sqrt(pow(x_ahead-x_ref,2)+pow(y_ahead-y_ref,2));

                double v_change = v_ref - car_speed;
                int sign = 1;
                if(v_ref < car_speed) sign = -1;

                //control de velocidad

                a_change = a_max;///2 *(1 + min(fabs(fabs(v_change)-1)/5,1.0));
                if(changing_lanes) a_change = a_max/2;
                car_speed += sign * min(fabs(v_change),fabs(a_change * 0.02));

                //control de Jerk
                /*
                double t_for_a_0 = fabs(a_change) / j_max;
                double attainable_v_ref = car_speed + a_change * t_for_a_0 - 0.5 * sign * j_max * pow(t_for_a_0,2);
                double apply_jerk = 0;
                if((sign == 1 && attainable_v_ref < v_ref) or (sign == -1 && attainable_v_ref > v_ref)){
                  if(fabs(a_change) > a_max){
                    apply_jerk = 0;
                  }else{
                    apply_jerk = sign * j_max;
                  }
                }else{
                  apply_jerk = - sign * j_max;
                }
                if(sign == 1 && a_change < 0) apply_jerk = j_max;
                if(sign == -1 && a_change > 0) apply_jerk = -j_max;
                cout << "vfp: " << attainable_v_ref << " car: " << car_speed << " acc: " << a_change << " jerk: " << apply_jerk << endl;
                a_change += apply_jerk*0.02;
                car_speed += a_change * 0.02;*/

                double d_change = car_speed * 0.02;
                x_ahead = x_ref + 1 * d_change / distance_ahead;
                y_ahead = s(x_ahead);

                //cout << "Add point with t-speed: " << car_speed  << " real speed: " << sqrt(pow(x_ahead-x_ref,2)+pow(y_ahead-y_ref,2))/0.02 << endl;

                x_ref = x_ahead;
                y_ref = y_ahead;

                next_x_vals.push_back(x_ahead*cos(car_yaw) - y_ahead * sin (car_yaw)+car_x);
                next_y_vals.push_back(x_ahead*sin(car_yaw) + y_ahead * cos (car_yaw)+car_y);

                //cout << "local x: " << x_ref << " local y: " << y_ref << " map x: " << next_x_vals[next_x_vals.size()-1] << " map y: " << next_y_vals[next_y_vals.size()-1] << endl;
                if(fastest_lane!=car_lane && allow_lane_change){
                  if(next_x_vals.size()>2){
                    double x_1 = next_x_vals[next_x_vals.size()-2];
                    double x_2 = next_x_vals[next_x_vals.size()-1];
                    double y_1 = next_y_vals[next_y_vals.size()-2];
                    double y_2 = next_y_vals[next_y_vals.size()-1];
                    double last_yaw = atan2(y_2-y_1,x_2-x_1);
                    vector<double> frenet_2 = getFrenetS(x_2,y_2,last_yaw,map_waypoints_x,map_waypoints_y,map_waypoints_s);
                    int current_lane = frenet_2[1]/4;
                    v_ref = lane_speed[current_lane];
                    double last_lane = frenet_2[1]/4.0-fastest_lane-0.5;
//                    if(fabs(last_lane)<0.05 or x_ahead>ptsx[n_points_advance]){// && next_x_vals.size()>=50){
//                    if(x_ahead>ptsx[n_points_lateral]){// && next_x_vals.size()>=50){
                    if(fabs(last_lane)<0.1){// && next_x_vals.size()>=50){
                      continue_traj=false;
                    }
                  }
                }else{
                  if(next_x_vals.size()>=25) continue_traj=false;
                }
              }
              /*if(previous_path_x.size() == 0){
                next_x_vals.clear();
                next_y_vals.clear();
                for(int i=0; i<50; i++){
                  vector<double> frenet_end = getXYspline((s_max-500.0)*i/50,6.0,map_waypoints_s, map_waypoints_x,map_waypoints_y,map_waypoints_dx,map_waypoints_dy);
                  next_x_vals.push_back(frenet_end[0]);
                  next_y_vals.push_back(frenet_end[1]);
                }
              }*/

            }


          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
