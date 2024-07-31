#ifndef BIT_STAR_H
#define BIT_STAR_H

#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include "geometry.h"
// #include "kdtree.h"
#include "vector.h"
#include <Eigen/Dense>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

const int MAX_NODES = 400;      // 最多节点
const int MAX_WAYPTS = 100;     // 最多waypoint数量
const float NEAR_DIST = 20.0;   // 当距离小于20cm时就可以直接过去.


// Node class
class Node {
public:
    double x, y;
    Node* parent;
    Node(double _x, double _y);
};

struct state{
  vector2f pos;
  state *parent;
  state *next;
};

// Tree class
class Tree {
public:
    state *start;
    state *goal;
    double r;//搜索半径
    std::set<state*> V;//顶点集合
    std::set<std::pair<state*, state*>> E;//存储树中的搜索边
    std::set<std::pair<state*, state*>> QE;//存储待扩展的边或待处理的边
    std::set<state*> QV;//存储候选顶点和待处理的顶点
    std::set<state*> V_old;//用于存储已处理过的顶点

    Tree(state *start, state *goal);
};

class Env {
public:
    std::pair<int, int> x_range;
    std::pair<int, int> y_range;
    std::vector<std::tuple<float, float, float, float>> obs_boundary;
    std::vector<std::tuple<float, float, float, float>> obs_rectangle;
    std::vector<std::tuple<float, float, float>> obs_circle;

    Env();

    static std::vector<std::tuple<float, float, float, float>> create_obs_boundary();
    static std::vector<std::tuple<float, float, float, float>> create_obs_rectangle();
    static std::vector<std::tuple<float, float, float>> create_obs_circle();
};

class Utils {
public:
    Env env;
    float delta;
    std::vector<std::tuple<float, float, float>> obs_circle;
    std::vector<std::tuple<float, float, float, float>> obs_rectangle;
    std::vector<std::tuple<float, float, float, float>> obs_boundary;

    Utils();

    void update_obs(std::vector<std::tuple<float, float, float>> obs_cir,
                    std::vector<std::tuple<float, float, float, float>> obs_bound,
                    std::vector<std::tuple<float, float, float, float>> obs_rec);

    std::vector<std::vector<std::vector<float>>> get_obs_vertex();

    bool is_intersect_rec(state *start, state *end, std::vector<float> o, std::vector<float> d,
                          std::vector<float> a, std::vector<float> b);

    bool is_intersect_circle(std::vector<float> o, std::vector<float> d, std::vector<float> a, float r);

    bool is_collision(state *start, state *end);

    bool is_inside_obs(state *node);

    static std::tuple<std::vector<float>, std::vector<float>> get_ray(state *start, state *end);

    static float get_dist(state *start, state *end);
};

// BITStar class
class BITStar {
public:
    state node[MAX_NODES];
	state waypoint[MAX_WAYPTS];
	state last_pathpoint[MAX_WAYPTS];  // 记录路径点
	state pathpoint[MAX_WAYPTS];
    state *start;
    state *goal;
    vector2f out_vel; // 最终输出的速度(大小，方向)
    
    // KDTree<state>tree;
    Tree *tree;
    int num_nodes;
	int max_nodes;
	int num_waypoints;
	int last_num_pathpoint;
	int num_pathpoint;
	static const int max_extend = 1;

    float goal_target_prob;
	float waypoint_target_prob;
	float waypoint_target_prob_old;
	float step_size;

    float last_pathlength;    //记录路径长度
	float pathlength;
	float path_target_length;       //记录nearestpoint与障碍物之间的距离
	float last_path_target_length;
	bool path_target_collision;     //记录nearestpoint是否与终点有障碍物
	bool last_path_target_collision;
	int still_count;          //路径没有变化的计数器

    bool init_;
    bool searchInCircle_;
    double circleRadius;
    vector2f circleCenter;

    double eta;
    int iter_max;    
    std::set<state*> X_sample;
    std::map<state*, double> g_T;//start到某个点的距离

public:
    // BITStar(std::pair<double, double> start, std::pair<double, double> goal, double _eta, int _iter_max);
    BITStar(std::pair<double, double> _start, std::pair<double, double> _goal, double _eta, int _iter_max);
    ~BITStar();
    Utils utils;

    // void planning();
    // void init(int _max_nodes, int _num_waypoints, float _goal_target_prob, float _waypoint_target_prob, float _step_size, state _goal, bool searchInCircle = false, vector2f circleCenter = vector2f(0.0, 0.0), double searchCircleRadius = 0.0);
    std::tuple<double, double, Eigen::Vector3d, Eigen::Matrix3d> init();
    float distance(state &state0, state &state1);
    // state random_state();
    // state *add_node(state newExtend, state *parent);
    // state choose_target(int &targ);
    // state *find_nearest(state target);
    // int extend(state *nearest, state target);
    vector<state> ebit(state *initial, state *_goal, bool drawtreeflag);
    vector<state> plan(/*obstacles *_obs, int obs_mask, state initial, vector2f _vel, state _goal*/);
//     void save_newpath();
// 	void cal_vel();
// 	vector2f get_vel();

    std::pair<std::vector<double>, std::vector<double>> ExtractPath(state *node);
    void Prune(double cBest);
    double cost(state *start, state *end);
    bool is_in_range(state *node, double delta);
    std::set<state*> Sample(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C);
    std::set<state*> SampleEllipsoid(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C);
    std::set<state*> SampleFreeSpace(int m);
//     double radius(int q);
    void ExpandVertex(state* v);
    double BestVertexQueueValue();
    double BestEdgeQueueValue();
    state* BestInVertexQueue();
    std::pair<state*, state*> BestInEdgeQueue();
    Eigen::Vector3d SampleUnitNBall();
    Eigen::Matrix3d RotationToWorldFrame(state *x_start, state *x_goal, double L);
    double calc_dist(vector2f start, vector2f end);
//     std::pair<double, double> calc_dist_and_angle(Node* node_start, Node* node_end);
    void animation(const Eigen::Vector3d& xCenter, double cMax, double cMin, double theta);
    void plot_grid(const std::string& name);
    void draw_ellipse(const Eigen::Vector3d& xCenter, double c_best, double dist, double theta);

private:
    double random_double(double min, double max);
    vector<state*> GetNearbyNodes(state node, double radius, bool from_tree = false);
};

#endif // BIT_STAR_H