#ifndef BIT_STAR_H
#define BIT_STAR_H

#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <geometry.h>
#include <Eigen/Dense>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

// Node class
class Node {
public:
    double x, y;
    Node* parent;

    Node(double _x, double _y);
};

// Tree class
class Tree {
public:
    Node* x_start;
    Node* x_goal;
    double r;
    std::set<Node*> V;
    std::set<std::pair<Node*, Node*>> E;
    std::set<std::pair<Node*, Node*>> QE;
    std::set<Node*> QV;
    std::set<Node*> V_old;

    Tree(Node* start, Node* goal);
};

// BITStar class
class BITStar {
public:
    Node* x_start;
    Node* x_goal;
    double eta;
    int iter_max;
    Tree* tree;
    std::set<Node*> X_sample;
    std::map<Node*, double> g_T;

    BITStar(std::pair<double, double> start, std::pair<double, double> goal, double _eta, int _iter_max);
    ~BITStar();

    void planning();
    std::tuple<double, double, Eigen::Vector3d, Eigen::Matrix3d> init();
    std::pair<std::vector<double>, std::vector<double>> ExtractPath();
    void Prune(double cBest);
    double cost(Node* start, Node* end);
    double f_estimated(Node* node);
    double g_estimated(Node* node);
    double h_estimated(Node* node);
    std::set<Node*> Sample(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C);
    std::set<Node*> SampleEllipsoid(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C);
    std::set<Node*> SampleFreeSpace(int m);
    double radius(int q);
    void ExpandVertex(Node* v);
    double BestVertexQueueValue();
    double BestEdgeQueueValue();
    Node* BestInVertexQueue();
    std::pair<Node*, Node*> BestInEdgeQueue();
    Eigen::Vector3d SampleUnitNBall();
    Eigen::Matrix3d RotationToWorldFrame(state x_start, state x_goal, double L);
    double calc_dist(Node* start, Node* end);
    std::pair<double, double> calc_dist_and_angle(Node* node_start, Node* node_end);
    void animation(const Eigen::Vector3d& xCenter, double cMax, double cMin, double theta);
    void plot_grid(const std::string& name);
    void draw_ellipse(const Eigen::Vector3d& xCenter, double c_best, double dist, double theta);

private:
    double random_double(double min, double max);
    std::vector<Node*> GetNearbyNodes(Node* node, double radius, bool from_tree = false);
};

#endif // BIT_STAR_H