#include "BITStarforssl.h"
#include <random>
#include <stdio.h>
#include <time.h>
#include "vector.h"
// #include "obstacle.h"

namespace{
    inline int boost_lrand48() {
        return rand();
    }

    inline float boost_drand48() { // 返回 [0, 1]
        return float(boost_lrand48() % 10000) / 10000;
    }

    inline float sdrand48() { return(2 * boost_drand48() - 1); }  //返回 [-1, 1]
}

// Node::Node(double _x, double _y) : x(_x), y(_y), parent(nullptr) {}

Tree::Tree(state *start, state *goal) : start(start), goal(goal), r(4.0) {}

BITStar::BITStar(std::pair<double, double> _start, std::pair<double, double> _goal, double _eta, int _iter_max)
    : eta(_eta), iter_max(_iter_max) {
    state *s_temp=new state();
    state *g_temp=new state();
    s_temp->pos = vector2f(_start.first, _start.second);
    s_temp->parent = nullptr;
    s_temp->next = nullptr;
    start =  s_temp;
    g_temp->pos = vector2f(_goal.first, _goal.second);
    g_temp->parent = nullptr;
    g_temp->next = nullptr;
    goal =  g_temp;
    tree = new Tree(start, goal);
    srand(time(0));
    init_ = false;
}

BITStar::~BITStar() {
    delete tree;
}
// // state BITStar::random_state(){//生成随即状态
// //     state s;
// //     if(!searchInCircle_)
// //         s.pos = vector2f((PARAM::Field::PITCH_LENGTH / 2 + PARAM::Field::GOAL_DEPTH) * sdrand48(), PARAM::Field::PITCH_LENGTH / 2 * sdrand48());
// //     else {
// //         double r = boost_drand48() * circleRadius;
// //         double angle = sdrand48() * PARAM::Math::PI;
// //         s.pos = vector2f(circleCenter.x + r * cos(angle), circleCenter.y + r * sin(angle));
// //     }
// //     s.parent = NULL;
// //     return(s);
// // }

// // 初始化bit信息，包括最大节点数，waypoints中储存的节点数，以goal为目标的概率，以waypoint中的点为目标的概率，步长.
// // 在调用时采用如下参数path_planner[i].init(150, 80, 0.15, 0.65, PARAM::Field::MAX_PLAYER_SIZE, initial);
// // void BITStar::init()(int _max_nodes, int _num_waypoints, float _goal_target_prob, float _waypoint_target_prob, float _step_size, state _goal, bool searchInCircle, vector2f circleCenter, double searchCircleRadius) {
// //     if ( init_ && !searchInCircle) {
// //         return;
// //     }
// //     vector2f minv;
// //     vector2f maxv;

// //     max_nodes = _max_nodes;
// //     num_waypoints = _num_waypoints;
// //     goal_target_prob = _goal_target_prob;
// //     waypoint_target_prob = _waypoint_target_prob;
// //     waypoint_target_prob_old = _waypoint_target_prob;
// //     step_size = _step_size;
// //     num_pathpoint = 0;
// //     last_num_pathpoint = 0;
// //     pathlength = 0;
// //     last_pathlength = 99999999.0f;
// //     path_target_length = 0;
// //     last_path_target_length = 999999999.0f;
// //     last_path_target_collision = false;
// //     path_target_collision = false;
// //     still_count = 0;
// //     out_vel = vector2f(0.0, 0.0);
// //     searchInCircle_ = searchInCircle;
// //     this->circleCenter = circleCenter;
// //     this->circleRadius = searchCircleRadius;

// //     // 将waypoint数组初始化为随机点
// //     for(int i = 0; i < num_waypoints; i++) {
// //         waypoint[i] = random_state();
// //     }

// //     //长方形(minv, maxv), 场地左下角和右上角
// //     if(!searchInCircle) {
// //         minv = vector2f(-PARAM::Field::PITCH_LENGTH / 2 - PARAM::Field::GOAL_DEPTH, -PARAM::Field::PITCH_WIDTH / 2);
// //         maxv = vector2f(PARAM::Field::PITCH_LENGTH / 2 + PARAM::Field::GOAL_DEPTH / 2, PARAM::Field::PITCH_WIDTH / 2);
// //     } else {
// //         minv = vector2f(circleCenter.x - searchCircleRadius, circleCenter.y - searchCircleRadius);
// //         maxv = vector2f(circleCenter.x + searchCircleRadius, circleCenter.y + searchCircleRadius);
// //         if(searchCircleRadius < step_size)
// //             cout << "Oh shit!!! Search radius is smaller than search step!!! ---path_planner.cpp" << endl;
// //         if((circleCenter - _goal->pos).length() > searchCircleRadius)
// //             cout << "Oh shit!!! Target is out of circle!!! ---path_planner.cpp" << endl;
// //     }
// //     // tree.setdim(minv, maxv, 16, 8); // 重置一棵树
// //     tree.V.add(start);
// //     X_sample.add(goal);

// //     g_T[start] = 0.0;
// //     g_T[goal] = std::numeric_limits<double>::infinity();//无限大

// //     double cMin = distance(start,goal);
// //     init_ = true;
// // }

Eigen::Matrix3d BITStar::RotationToWorldFrame(state *start, state *goal, double L) {
    Eigen::Vector3d a1((goal->pos.x - start->pos.x) / L, (goal->pos.y - start->pos.y) / L, 0.0);//起始点到目标点的单位向量
    Eigen::Vector3d e1(1.0, 0.0, 0.0);//x轴单位向量
    Eigen::Matrix3d M = a1 * e1.transpose();//矩阵乘积
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);//矩阵奇异值分解
    Eigen::Matrix3d C = svd.matrixU() * Eigen::Matrix3d::Identity() * svd.matrixV().transpose();
    return C;
}

std::tuple<double, double, Eigen::Vector3d, Eigen::Matrix3d> BITStar::init(){
    utils = Utils();
    tree->V.insert(start);
    X_sample.insert(goal);
    max_nodes = MAX_NODES;

    g_T[start] = 0.0;
    g_T[goal] = std::numeric_limits<double>::infinity();//无限大
    // std::cout << "g_T[start] = " << g_T[start] << std::endl;
    // std::cout<<start->pos.x<<std::endl;

    double cMin = calc_dist(start->pos,goal->pos);;//起始点和终点之间的距离,也是路径规划的最小距离
    double angle = atan2(goal->pos.y-start->pos.y,goal->pos.x-start->pos.x);//起始点和终点之间的角度
    Eigen::Matrix3d C = RotationToWorldFrame(start, goal, cMin); //计算旋转矩阵
    Eigen::Vector3d xCenter((start->pos.x + goal->pos.x) / 2.0, (start->pos.y + goal->pos.y) / 2.0,0.0);//计算起点和终点的中心点

    return {angle, cMin, xCenter, C};
}

//返回两点间的距离
inline float BITStar::distance(state &state0, state &state1) {
    return(Vector::distance(state0.pos, state1.pos));
}

double BITStar::calc_dist(vector2f start, vector2f end) {
    return std::hypot(start.x - end.x, start.y - end.y);
}

// 将节点n加入kdtree,并返回新加入的节点
// state *path_planner::add_node(state newExtend, state *parent) {
//     // 传入的n是扩展后的节点也就是nearest+step之后的节点
//     // parent是nearest节点
//     if(num_nodes >= max_nodes) {
//         return(NULL);
//     }
//     newExtend.parent = parent;
//     node[num_nodes] = newExtend;
//     tree.add(&node[num_nodes]); //state和cmu的node的转换
//     num_nodes++;
//     return(&node[num_nodes - 1]);
// }

void BITStar::Prune(double cBest) {
    std::set<state*> new_X_sample;
    for (auto& x : X_sample) {
        if ((Vector::distance(start->pos,x->pos)+Vector::distance(x->pos,goal->pos)) < cBest) {
            new_X_sample.insert(x);
        }
    }
    X_sample = new_X_sample;

    std::set<state*> new_V;
    for (auto& v : tree->V) {
        if ((Vector::distance(start->pos,v->pos)+Vector::distance(v->pos,goal->pos)) <= cBest) {
            new_V.insert(v);
        }
    }
    tree->V = new_V;

    std::set<std::pair<state*, state*>> new_E;
    for (auto& [v, w] : tree->E) {
        if ((Vector::distance(start->pos,v->pos)+Vector::distance(v->pos,goal->pos)) <= cBest && (Vector::distance(start->pos,w->pos)+Vector::distance(w->pos,goal->pos)) <= cBest) {
            new_E.insert(std::make_pair(v, w));
        }
    }
    tree->E = new_E;

    for (auto& v : tree->V) {
        if (g_T[v] == std::numeric_limits<double>::infinity()) {
            X_sample.insert(v);
        }
    }
    std::set<state*> temp_V;
    for (auto& v : tree->V){
        if (g_T[v] < std::numeric_limits<double>::infinity()) {
            temp_V.insert(v);
        }
        // std::cout<<g_T[start]<<std::endl;
    }
    tree->V = temp_V;
}

std::set<state*> BITStar::Sample(int m, double cMax, double cMin, const Eigen::Vector3d&  xCenter, const Eigen::Matrix3d& C) {
    if (cMax < std::numeric_limits<double>::infinity()) {
        return SampleEllipsoid(m, cMax, cMin, xCenter, C);
    }
    return SampleFreeSpace(m);
}


Eigen::Vector3d BITStar::SampleUnitNBall() {
     // 创建一个随机数生成器
    std::random_device rd;// 用于生成种子
    std::mt19937 gen(rd());// 以 rd() 作为种子初始化 Mersenne Twister 生成器
      // 定义一个分布范围
    std::uniform_real_distribution<> dis(-1, 1);//定义了一个均匀分布，用于生成范围在 -1.0 到 1.0 之间的实数
    while (true) {
        double x = dis(gen);
        double y = dis(gen);
        if (x * x + y * y < 1) {
            return Eigen::Vector3d(x, y, 0.0);
        }
    }
}

bool BITStar::is_in_range(state *node, double delta) {
    double x_range[2] = {0.0, 50.0};
    double y_range[2] = {0.0, 30.0};
    return (x_range[0] + delta <= node->pos.x && node->pos.x <= x_range[1] - delta &&
            y_range[0] + delta <= node->pos.y && node->pos.y <= y_range[1] - delta);
}

//在椭圆内规划
std::set<state*> BITStar::SampleEllipsoid(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C) {
    double r[3] = {cMax / 2.0, std::sqrt(cMax * cMax - cMin * cMin) / 2.0, std::sqrt(cMax * cMax - cMin * cMin) / 2.0};
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L.diagonal() << r[0], r[1], r[2];

    int ind = 0;
    double delta = 0.5;
    std::set<state*> sample_set;

    while (ind < m) {
        Eigen::Vector3d xBall = SampleUnitNBall();
        Eigen::Vector3d x_rand = C * L * xBall + xCenter;
        state *node=new state();
        node->pos = vector2f(x_rand(0, 0),x_rand(1, 0));
        node->parent = nullptr;
        node->next = nullptr;

        if (!utils.is_inside_obs(node) && is_in_range(node, delta)) {
            sample_set.insert(node);
            ind++;
        } else {
            delete node;
        }
    }
    return sample_set;
}

//一直在场地上任意找点
std::set<state*> BITStar::SampleFreeSpace(int m) {
    double delta = 0.5;
    std::set<state*> sample_set;

    int ind = 0;
    while (ind < m) {
        state *node=new state();
        node->parent = nullptr;
        node->next = nullptr;
        node->pos=vector2f(random_double(0 + delta, 50 - delta), random_double(0 + delta, 30 - delta));
        if (!utils.is_inside_obs(node)) {
            sample_set.insert(node);
            ind++;
        } else {
            delete node;// 如果节点在障碍物内，删除节点以防止内存泄漏
        }
    }
    return sample_set; 
}

double BITStar::BestEdgeQueueValue() {
    if (tree->QE.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& [v, x] : tree->QE) {
        min_value = std::min(min_value, g_T[v] + Vector::distance(x->pos,v->pos) + Vector::distance(x->pos,goal->pos));
    }
    return min_value;
}

double BITStar::BestVertexQueueValue() {
    if (tree->QV.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& v : tree->QV) {
        min_value = std::min(min_value, g_T[v] + Vector::distance(v->pos,goal->pos));
    }
    return min_value;
}

std::vector<state*> BITStar::GetNearbyNodes(state node, double radius, bool from_tree) {
    std::vector<state*> nearby_nodes;
    const std::set<state*>& nodes = from_tree ? tree->V : X_sample;
    for (auto& n : nodes) {
        if (Vector::distance(n->pos, node.pos)<= radius) {
            nearby_nodes.push_back(n);
        }
    }
    return nearby_nodes;
}

void BITStar::ExpandVertex(state* v) {
    tree->QV.erase(v);
    int count=0;
    auto X_near = GetNearbyNodes(*v, tree->r);//找顶点
    for (auto& x : X_near) {
        if (Vector::distance(start->pos,v->pos) + Vector::distance(x->pos,v->pos) +Vector::distance(x->pos,goal->pos) < g_T[goal]) {
            g_T[x] = std::numeric_limits<double>::infinity();
            tree->QE.insert(std::make_pair(v, x));
        }
        count++;
    }

    for(auto& [y,x] : tree->QE){
        // std::cout<<x->pos.x<<" " << x->pos.y<<std::endl;
    }


    if (tree->V_old.find(v) == tree->V_old.end()) {
        auto V_near = GetNearbyNodes(*v, tree->r, true);//在树内找

        for (auto& w : V_near) {
            if (tree->E.find(std::make_pair(v, w)) == tree->E.end() &&
                Vector::distance(start->pos,v->pos) + Vector::distance(w->pos,v->pos) +Vector::distance(w->pos,goal->pos) < g_T[goal] &&
                g_T[v] + Vector::distance(w->pos,v->pos) < g_T[w]) {
                tree->QE.insert(std::make_pair(v, w));
                if (g_T.find(w) == g_T.end()) {
                    g_T[w] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }
}

state* BITStar::BestInVertexQueue() {
    if (tree->QV.empty()) {
        std::cerr << "QV is Empty!" << std::endl;
        return nullptr;
    }
    state* best_node;
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& v : tree->QV) {
        double value = g_T[v] + Vector::distance(v->pos, goal->pos);
        if (value < min_value) {
            min_value = value;
            best_node = v;
        }
    }
    return best_node;
}

std::pair<state*, state*> BITStar::BestInEdgeQueue() {
    if (tree->QE.empty()) {
        std::cerr << "QE is Empty!" << std::endl;
        return {nullptr, nullptr};
    }
    std::pair<state*, state*> best_edge = {nullptr, nullptr};
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& [v, x] : tree->QE) {
        double value = g_T[v] + Vector::distance(v->pos,x->pos)  + Vector::distance(goal->pos,x->pos);
        if (value < min_value) {
            min_value = value;
            best_edge = std::make_pair(v, x);
        }
    }
    return best_edge;
}

double BITStar::cost(state *start, state *end) {
    if (utils.is_collision(start, end)) {
        return std::numeric_limits<double>::infinity();
    }
    return Vector::distance(start->pos, end->pos);
}

void BITStar::plot_grid(const std::string& name) {
    for (const auto& obs : utils.obs_boundary) {
        double ox, oy, w, h;
        std::tie(ox, oy, w, h) = obs;
        std::vector<double> x = {ox, ox + w, ox + w, ox, ox};
        std::vector<double> y = {oy, oy, oy + h, oy + h, oy};
        plt::plot(x, y, "black");
    }

    for (const auto& obs : utils.obs_rectangle) {
        double ox, oy, w, h;
        std::tie(ox, oy, w, h) = obs;
        std::vector<double> x = {ox, ox + w, ox + w, ox, ox};
        std::vector<double> y = {oy, oy, oy + h, oy + h, oy};
        plt::plot(x, y, "gray");
    }

    for (const auto& obs : utils.obs_circle) {
        double ox, oy, r;
        std::tie(ox, oy, r) = obs;
        const int num_points = 100; // 圆的点数
        std::vector<double> x(num_points), y(num_points);
        for (int i = 0; i < num_points; ++i) {
            double angle = 2 * M_PI * i / num_points;
            x[i] = ox + r * cos(angle);
            y[i] = oy + r * sin(angle);
        }
        plt::plot(x, y, "gray");
    }
    std::vector<double> x_start = {start->pos.x};
    std::vector<double> y_start = {start->pos.y};
    std::vector<double> x_goal = {goal->pos.x};
    std::vector<double> y_goal = {goal->pos.y};
    
    // plt::plot(start->pos.x, start->pos.y, "bs", {{"linewidth", 3}});
    plt::scatter(x_start, y_start, 32.0, {{"color", "red"}});
    plt::scatter(x_goal, y_goal, 32.0, {{"color", "blue"}});
    // plt::plot(goal->pos.x, goal->pos.y, "rs", {{"linewidth", 3}});

    plt::title(name);
    plt::axis("equal");

    // plt::show();
}

void BITStar::draw_ellipse(const Eigen::Vector3d& xCenter, double c_best, double dist, double theta) {
    double a = std::sqrt(c_best * c_best - dist * dist) / 2.0;
    double b = c_best / 2.0;
    double angle = M_PI / 2.0 - theta;
    double cx = xCenter(0);
    double cy = xCenter(1);
    std::vector<double> t;
    for (double i = 0; i <= 2 * M_PI + 0.1; i += 0.2) {
        t.push_back(i);
    }
    std::vector<double> x, y;
    for (auto& it : t) {
        x.push_back(a * std::cos(it));
        y.push_back(b * std::sin(it));
    }
    Eigen::Matrix2d rot;
    rot << std::cos(-angle), -std::sin(-angle), std::sin(-angle), std::cos(-angle);
    Eigen::MatrixXd fx = rot * Eigen::Map<const Eigen::MatrixXd>(x.data(), 2, x.size());
    std::vector<double> px(fx.row(0).data(), fx.row(0).data() + fx.row(0).size());
    std::vector<double> py(fx.row(1).data(), fx.row(1).data() + fx.row(1).size());
    for (size_t i = 0; i < px.size(); ++i) {
        px[i] += cx;
        py[i] += cy;
    }
    plt::plot(px, py, "b-");
}

void BITStar::animation(const Eigen::Vector3d& xCenter, double cMax, double cMin, double theta) {
    plt::clf();
    plot_grid("Batch Informed Trees (BIT*)");

    for (auto& v : X_sample) {
        std::vector<double> x_vals = {v->pos.x};
        std::vector<double> y_vals = {v->pos.y};
        plt::scatter(x_vals, y_vals, 2.0, {{"color", "lightgrey"}});  // 使用scatter函数设置颜色和标记大小
    }

    if (cMax < std::numeric_limits<double>::infinity()) {
        // draw_ellipse(xCenter, cMax, cMin, theta);
    }

    // for (auto& v : tree->QV) {
    //     std::vector<double> x_vals = {v->pos.x};
    //     std::vector<double> y_vals = {v->pos.y};
    //     plt::scatter(x_vals, y_vals, 32.0, {{"color", "green"}});  // 使用scatter函数设置颜色和标记大小
    // }

    for (auto& [v, w] : tree->E) {
        plt::plot({v->pos.x, w->pos.x}, {v->pos.y, w->pos.y}, "g-");
    }

    // for (auto& [v, w] : tree->QE) {
    //     plt::plot({v->pos.x, w->pos.x}, {v->pos.y, w->pos.y}, "r-");
    // }

    plt::pause(0.001);
}

vector<state> BITStar::ebit(state *initial, state *_goal, bool drawtreeflag){
    auto [angle, cMin, xCenter, C] = init();
/* 0、如果初始点已经距离goal比较近，采用简化方法。否则执行bit的plan
   1、
*/
    vector<state> targetList;
    // state target,*nearest,*nearest_goal;
    float tempDist;
    float newDist;
    goal=_goal;
    // 步骤0: 如果initial与goal足够近时的简便方法: 直接考虑target在initial与goal连线上,尽量靠近goal,不避障碍
    tempDist = Vector::distance(initial->pos, goal->pos);

    // if(obs->check(initial,goal)){
    //     // 起点与终点连线之间没有碰撞
    //     pathlength =  0;
    //     path_target_length = 0;
    //     path_target_collision = 1;
    //     num_pathpoint = 2;
    //     pathpoint[0] = initial;
    //     pathpoint[1] = goal;
    //     if (isnan(goal->pos.x) || isnan(goal->pos.y)) {
    //         printf("ERROR INIT3\n");
    //     }
    //     targetList.push_back(goal);
    //     return(targetList);
    // }else{//如果有碰撞，则进行路径规划
        // tree.clear();
        int i,iter_limit;
        int targ_type;
        i = num_nodes = 0;
        // nearest = nearest_goal = add_node(initial, NULL);
        // tempDist = distance(*nearest, goal);
        // plan bit规划,直到距离目标点足够近,或者循环次数达到上限,足够近时就可能碰到障碍物
        iter_limit = max_nodes; // 尝试的次数等于的结点数
        // i到达数目或者num_nodes到达数目且nearest_goal离目标点大于设定值时会跳出while
        // i到了,num_nodes到了,或者找到离目标点很近的点了,就会跳出while
        while(i < iter_limit){
            if(tree->QE.empty() && tree->QV.empty()){//如果树节点为空
                if(i==0){
                    num_nodes=350;
                }else{
                    num_nodes=200;
                }
                
                if(goal->parent!=nullptr){
                    auto [path_x, path_y] = ExtractPath(goal);
                    std::cout<<"plan success"<<std::endl;
                    plt::plot(path_x, path_y, "r-");
                    plt::pause(0.5);
                }
                
                Prune(g_T[goal]);
                X_sample.merge(Sample(num_nodes, g_T[goal], cMin, xCenter, C));//采样并合并新节点
                tree->V_old = tree->V;
                tree->QV = tree->V;
            
            }
            // plt::clf();
            // for (auto& v : X_sample) {
            //     std::vector<double> x_vals = {v->pos.x};
            //     std::vector<double> y_vals = {v->pos.x};
            //     plt::scatter(x_vals, y_vals, 2.0, {{"color", "lightgrey"}});  // 使用scatter函数设置颜色和标记大小
            // }
            // plt::show();
            while (BestVertexQueueValue() <= BestEdgeQueueValue()) {
                ExpandVertex(BestInVertexQueue());
                // std::cout<<BestVertexQueueValue() <<std::endl;
                // std::cout<<BestEdgeQueueValue() <<std::endl;
            }   

            // for(auto& [y,x] : tree->QE){
                // std::cout<<x->pos.x<<" " << x->pos.y<<std::endl;
                // std::cout<<y->pos.x<<" " << y->pos.y<<std::endl;
            // }
            auto [vm, xm] = BestInEdgeQueue();
            tree->QE.erase(std::make_pair(vm, xm));

            // double actual_cost = cost(vm, xm);
            // g_T[xm] = g_T[vm] /*+ actual_cost*/;
            // tree->E.insert(std::make_pair(vm, xm));
            // xm->parent = vm;

            if (g_T[vm] + Vector::distance(vm->pos,xm->pos) + Vector::distance(vm->pos,xm->pos) < g_T[goal]) {
                double actual_cost = cost(vm, xm);
                if (Vector::distance(start->pos,vm->pos) + actual_cost + Vector::distance(goal->pos,xm->pos) < g_T[goal]) {
                    if (g_T[vm] + actual_cost < g_T[xm]) {
                        if (tree->V.find(xm) != tree->V.end()) {
                            std::set<std::pair<state*, state*>> edge_delete;
                            for (auto& [v, x] : tree->E) {
                                if (x == xm) edge_delete.insert(std::make_pair(v, x));
                            }
                            for (auto& edge : edge_delete) tree->E.erase(edge);
                        } else {
                            X_sample.erase(xm);
                            tree->V.insert(xm);
                            tree->QV.insert(xm);
                        }
                        g_T[xm] = g_T[vm] + actual_cost;
                        tree->E.insert(std::make_pair(vm, xm));
                        xm->parent = vm;
                        // if()
                        // std::cout<<xm->parent->pos.x<<std::endl;
                        

                        std::set<std::pair<state*, state*>> set_delete;
                        for (auto& [v, x] : tree->QE) {
                            if (x == xm && g_T[v] +  Vector::distance(v->pos,xm->pos) >= g_T[xm]) {
                                set_delete.insert(std::make_pair(v, x));
                            }
                        }
                        for (auto& edge : set_delete) tree->QE.erase(edge);
                    }
                }
            } else {
                // std::cout<<(int)(goal->parent!=nullptr)<<std::endl;
                std::cout << "null" <<std::endl;
                tree->QE.clear();
                tree->QV.clear();
                // std::cout<<tree->QE.empty()<<std::endl;
                // std::cout<<tree->QV.empty()<<std::endl;
                // std::cout<<i<<std::endl;
            }

            if (i % 5 == 0) {
                animation(xCenter, g_T[goal], cMin, angle);
            }
        i++;
    }
        auto [path_x, path_y] = ExtractPath(goal);
        plt::show();
    // }
    return targetList;
}

vector<state> BITStar::plan(){
    vector<state> target_now; // 要返回的目标点，每次只返回一个点，下一次再进来返回另外一个点，planner是全局的
    target_now = ebit(start,goal,1);
    return target_now;
}

std::pair<std::vector<double>, std::vector<double>> BITStar::ExtractPath(state *node) {
    state *_node=node;
    std::vector<double> path_x = {_node->pos.x};
    std::vector<double> path_y = {_node->pos.y};

    while (_node->parent) {
        _node = _node->parent;
        path_x.push_back(_node->pos.x);
        path_y.push_back(_node->pos.y);
        // std::cout<<"extract posx:"<<_node->pos.x<<std::endl;
    }

    return {path_x, path_y};
}

double BITStar::random_double(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

Utils::Utils() : delta(0.5) {
    env = Env();
    obs_circle = env.obs_circle;
    obs_rectangle = env.obs_rectangle;
    obs_boundary = env.obs_boundary;
}

void Utils::update_obs(std::vector<std::tuple<float, float, float>> obs_cir,
                       std::vector<std::tuple<float, float, float, float>> obs_bound,
                       std::vector<std::tuple<float, float, float, float>> obs_rec) {
    obs_circle = obs_cir;
    obs_boundary = obs_bound;
    obs_rectangle = obs_rec;
}

std::vector<std::vector<std::vector<float>>> Utils::get_obs_vertex() {
    float delta = this->delta;
    std::vector<std::vector<std::vector<float>>> obs_list;

    for (auto& obs : obs_rectangle) {
        float ox = std::get<0>(obs);
        float oy = std::get<1>(obs);
        float w = std::get<2>(obs);
        float h = std::get<3>(obs);

        std::vector<std::vector<float>> vertex_list = {
            {ox - delta, oy - delta},
            {ox + w + delta, oy - delta},
            {ox + w + delta, oy + h + delta},
            {ox - delta, oy + h + delta}
        };
        obs_list.push_back(vertex_list);
    }

    return obs_list;
}

bool Utils::is_intersect_rec(state *start, state *end, std::vector<float> o, std::vector<float> d,
                             std::vector<float> a, std::vector<float> b) {
    std::vector<float> v1 = {o[0] - a[0], o[1] - a[1]};
    std::vector<float> v2 = {b[0] - a[0], b[1] - a[1]};
    std::vector<float> v3 = {-d[1], d[0]};

    float div = v2[0] * v3[0] + v2[1] * v3[1];

    if (div == 0) {
        return false;
    }

    float t1 = std::abs(v2[0] * v1[1] - v2[1] * v1[0]) / div;
    float t2 = (v1[0] * v3[0] + v1[1] * v3[1]) / div;

    if (t1 >= 0 && 0 <= t2 && t2 <= 1) {
        state *shot;
        shot->pos = vector2f(o[0] + t1 * d[0], o[1] + t1 * d[1]);
        float dist_obs = get_dist(start, shot);
        float dist_seg = get_dist(start, end);
        if (dist_obs <= dist_seg) {
            return true;
        }
    }

    return false;
}

bool Utils::is_intersect_circle(std::vector<float> o, std::vector<float> d, std::vector<float> a, float r) {
    float d2 = d[0] * d[0] + d[1] * d[1];
    float delta = this->delta;

    if (d2 == 0) {
        return false;
    }

    float t = ((a[0] - o[0]) * d[0] + (a[1] - o[1]) * d[1]) / d2;

    if (0 <= t && t <= 1) {
        state *shot=new state(),*node=new state();
        shot->pos = vector2f(o[0] + t * d[0], o[1] + t * d[1]);
        node->pos = vector2f(a[0], a[1]);
        if (get_dist(shot, node) <= r + delta) {
            return true;
        }
    }

    return false;
}

bool Utils::is_collision(state* start, state* end) {
    if (is_inside_obs(start) || is_inside_obs(end)) {
        return true;
    }

    std::vector<float> o, d;
    std::tie(o, d) = get_ray(start, end);
    auto obs_vertex = get_obs_vertex();

    for (auto& vertex : obs_vertex) {
        if (is_intersect_rec(start, end, o, d, vertex[0], vertex[1]) ||
            is_intersect_rec(start, end, o, d, vertex[1], vertex[2]) ||
            is_intersect_rec(start, end, o, d, vertex[2], vertex[3]) ||
            is_intersect_rec(start, end, o, d, vertex[3], vertex[0])) {
            return true;
        }
    }

    for (auto& circle : obs_circle) {
        if (is_intersect_circle(o, d, {std::get<0>(circle), std::get<1>(circle)}, std::get<2>(circle))) {
            return true;
        }
    }

    return false;
}

bool Utils::is_inside_obs(state *node) {
    float delta = this->delta;

    for (auto& circle : obs_circle) {
        if (std::hypot(node->pos.x - std::get<0>(circle), node->pos.y - std::get<1>(circle)) <= std::get<2>(circle) + delta) {
            return true;
        }
    }

    for (auto& rect : obs_rectangle) {
        float x = std::get<0>(rect);
        float y = std::get<1>(rect);
        float w = std::get<2>(rect);
        float h = std::get<3>(rect);
        if (0 <= node->pos.x - (x - delta) && node->pos.x - (x - delta) <= w + 2 * delta &&
            0 <= node->pos.y - (y - delta) && node->pos.y - (y - delta) <= h + 2 * delta) {
            return true;
        }
    }

    for (auto& bound : obs_boundary) {
        float x = std::get<0>(bound);
        float y = std::get<1>(bound);
        float w = std::get<2>(bound);
        float h = std::get<3>(bound);
        if (0 <= node->pos.x- (x - delta) && node->pos.x - (x - delta) <= w + 2 * delta &&
            0 <= node->pos.y - (y - delta) && node->pos.y - (y - delta) <= h + 2 * delta) {
            return true;
        }
    }

    return false;
}

std::tuple<std::vector<float>, std::vector<float>> Utils::get_ray(state *start, state *end) {
    std::vector<float> orig = {start->pos.x, start->pos.y};
    std::vector<float> direc = {end->pos.x - start->pos.x, end->pos.y - start->pos.y};
    return std::make_tuple(orig, direc);
}

float Utils::get_dist(state *start, state *end) {
    return std::hypot(end->pos.x - start->pos.x, end->pos.y - start->pos.y);
}

Env::Env() {
    x_range = std::make_pair(0, 50);
    y_range = std::make_pair(0, 30);
    obs_boundary = create_obs_boundary();
    obs_rectangle = create_obs_rectangle();
    obs_circle = create_obs_circle();
}

std::vector<std::tuple<float, float, float, float>> Env::create_obs_boundary() {
    std::vector<std::tuple<float, float, float, float>> obs_boundary = {
        std::make_tuple(0, 0, 1, 30),
        std::make_tuple(0, 30, 50, 1),
        std::make_tuple(1, 0, 50, 1),
        std::make_tuple(50, 1, 1, 30)
    };
    return obs_boundary;
}

std::vector<std::tuple<float, float, float, float>> Env::create_obs_rectangle() {
    std::vector<std::tuple<float, float, float, float>> obs_rectangle = {
        std::make_tuple(14, 12, 8, 2),
        std::make_tuple(18, 22, 8, 3),
        std::make_tuple(26, 7, 2, 12),
        std::make_tuple(32, 14, 10, 2)
    };
    return obs_rectangle;
}

std::vector<std::tuple<float, float, float>> Env::create_obs_circle() {
    std::vector<std::tuple<float, float, float>> obs_circle = {
        std::make_tuple(7, 12, 3),
        std::make_tuple(46, 20, 2),
        std::make_tuple(15, 5, 2),
        std::make_tuple(37, 7, 3),
        std::make_tuple(37, 23, 3)
    };
    return obs_circle;
}
