#include "BITStar.h"
#include <random>

Node::Node(double _x, double _y) : x(_x), y(_y), parent(nullptr) {}

Tree::Tree(Node* start, Node* goal) : x_start(start), x_goal(goal), r(4.0) {}

BITStar::BITStar(std::pair<double, double> start, std::pair<double, double> goal, double _eta, int _iter_max)
    : eta(_eta), iter_max(_iter_max) {
    x_start = new Node(start.first, start.second);
    x_goal = new Node(goal.first, goal.second);
    tree = new Tree(x_start, x_goal);
}

BITStar::~BITStar() {
    delete x_start;
    delete x_goal;
    delete tree;
}

void BITStar::planning() {
    auto [theta, cMin, xCenter, C] = init();//初始化

    for (int k = 0; k < 500; ++k) {
        if (tree->QE.empty() && tree->QV.empty()) {//树上两节点为空
            int m = (k == 0) ? 350 : 200;
            if (x_goal->parent != nullptr) { //如果有最初节点
                auto [path_x, path_y] = ExtractPath();//提取最初节点路径
                plt::plot(path_x, path_y, "r-");
                plt::pause(0.5);
            }

            Prune(g_T[x_goal]);//裁减路径
            X_sample.merge(Sample(m, g_T[x_goal], cMin, xCenter, C));//采样并合并新节点
            tree->V_old = {tree->V.begin(), tree->V.end()};
            tree->QV = {tree->V.begin(), tree->V.end()};
        }

        while (BestVertexQueueValue() <= BestEdgeQueueValue()) {
            ExpandVertex(BestInVertexQueue());
        }

        auto [vm, xm] = BestInEdgeQueue();
        tree->QE.erase(std::make_pair(vm, xm));

        if (g_T[vm] + calc_dist(vm, xm) + h_estimated(xm) < g_T[x_goal]) {
            double actual_cost = cost(vm, xm);
            if (g_estimated(vm) + actual_cost + h_estimated(xm) < g_T[x_goal]) {
                if (g_T[vm] + actual_cost < g_T[xm]) {
                    if (tree->V.find(xm) != tree->V.end()) {
                        std::set<std::pair<Node*, Node*>> edge_delete;
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

                    std::set<std::pair<Node*, Node*>> set_delete;
                    for (auto& [v, x] : tree->QE) {
                        if (x == xm && g_T[v] + calc_dist(v, xm) >= g_T[xm]) {
                            set_delete.insert(std::make_pair(v, x));
                        }
                    }
                    for (auto& edge : set_delete) tree->QE.erase(edge);
                }
            }
        } else {
            tree->QE.clear();
            tree->QV.clear();
        }

        if (k % 5 == 0) {
            animation(xCenter, g_T[x_goal], cMin, theta);
        }
    }

    auto [path_x, path_y] = ExtractPath();
    plt::plot(path_x, path_y, "r-");
    plt::pause(0.01);
    plt::show();
}

std::tuple<double, double, Eigen::Vector3d, Eigen::Matrix3d> BITStar::init() {
    tree->V.insert(x_start);
    X_sample.insert(x_goal);

    g_T[x_start] = 0.0;
    g_T[x_goal] = std::numeric_limits<double>::infinity();//无限大

    double cMin;//起始点和终点之间的距离
    double theta;//起始点和终点之间的角度
    std::tie(cMin, theta) = calc_dist_and_angle(x_start, x_goal);
    Eigen::Matrix3d C = RotationToWorldFrame(x_start, x_goal, cMin); //计算旋转矩阵
    Eigen::Vector3d xCenter((x_start->x + x_goal->x) / 2.0, (x_start->y + x_goal->y) / 2.0, 0.0);//计算中心点

    return {theta, cMin, xCenter, C};
}

std::pair<std::vector<double>, std::vector<double>> BITStar::ExtractPath() {
    Node* node = x_goal;
    std::vector<double> path_x = {node->x};
    std::vector<double> path_y = {node->y};

    while (node->parent) {
        node = node->parent;
        path_x.push_back(node->x);
        path_y.push_back(node->y);
    }

    return {path_x, path_y};
}

void BITStar::Prune(double cBest) {
    std::set<Node*> new_X_sample;
    for (auto& x : X_sample) {
        if (f_estimated(x) < cBest) {
            new_X_sample.insert(x);
        }
    }
    X_sample = new_X_sample;

    std::set<Node*> new_V;
    for (auto& v : tree->V) {
        if (f_estimated(v) <= cBest) {
            new_V.insert(v);
        }
    }
    tree->V = new_V;

    std::set<std::pair<Node*, Node*>> new_E;
    for (auto& [v, w] : tree->E) {
        if (f_estimated(v) <= cBest && f_estimated(w) <= cBest) {
            new_E.insert(std::make_pair(v, w));
        }
    }
    tree->E = new_E;

    for (auto& v : tree->V) {
        if (g_T[v] == std::numeric_limits<double>::infinity()) {
            X_sample.insert(v);
        }
    }
    std::set<Node*> temp_V;
    for (auto& v : tree->V) {
        if (g_T[v] < std::numeric_limits<double>::infinity()) {
            temp_V.insert(v);
        }
    }
    tree->V = temp_V;
}

double BITStar::cost(Node* start, Node* end) {
    if (is_collision(start, end)) {
        return std::numeric_limits<double>::infinity();
    }
    return calc_dist(start, end);
}

double BITStar::f_estimated(Node* node) {
    return g_estimated(node) + h_estimated(node);
}

//计算出发点和node之间的距离
double BITStar::g_estimated(Node* node) { 
    return calc_dist(x_start, node);
}

//计算node到目标点之间的距离
double BITStar::h_estimated(Node* node) {
    return calc_dist(node, x_goal);
}

std::set<Node*> BITStar::Sample(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C) {
    if (cMax < std::numeric_limits<double>::infinity()) {
        return SampleEllipsoid(m, cMax, cMin, xCenter, C);
    }
    return SampleFreeSpace(m);
}

std::set<Node*> BITStar::SampleEllipsoid(int m, double cMax, double cMin, const Eigen::Vector3d& xCenter, const Eigen::Matrix3d& C) {
    double r[3] = {cMax / 2.0, std::sqrt(cMax * cMax - cMin * cMin) / 2.0, std::sqrt(cMax * cMax - cMin * cMin) / 2.0};
    Eigen::Matrix3d L = Eigen::Matrix3d::Zero();
    L.diagonal() << r[0], r[1], r[2];

    int ind = 0;
    double delta = 0.1;
    std::set<Node*> sample_set;

    while (ind < m) {
        Eigen::Vector3d xBall = SampleUnitNBall();
        Eigen::Vector3d x_rand = C * L * xBall + xCenter;
        Node* node = new Node(x_rand(0, 0), x_rand(1, 0));

        if (!is_inside_obs(node) && is_in_range(node, delta)) {
            sample_set.insert(node);
            ind++;
        } else {
            delete node;
        }
    }
    return sample_set;
}

std::set<Node*> BITStar::SampleFreeSpace(int m) {
    double delta = 0.1;
    std::set<Node*> sample_set;

    int ind = 0;
    while (ind < m) {
        Node* node = new Node(random_double(0 + delta, 50 - delta), random_double(0 + delta, 50 - delta));
        if (!is_inside_obs(node)){
            sample_set.insert(node);
            ind++;
        }else{
            delete node;
        }
    }
    return sample_set;
}

double BITStar::radius(int q) {
    double cBest = g_T[x_goal];
    int lambda_X = std::count_if(tree->V.begin(), tree->V.end(), [&](Node* v) { return f_estimated(v) <= cBest; });
    return 2 * eta * std::sqrt((1.5 * lambda_X / M_PI * std::log(q)) / q);
}

void BITStar::ExpandVertex(Node* v) {
    tree->QV.erase(v);
    auto X_near = GetNearbyNodes(v, tree->r);

    for (auto& x : X_near) {
        if (g_estimated(v) + calc_dist(v, x) + h_estimated(x) < g_T[x_goal]) {
            g_T[x] = std::numeric_limits<double>::infinity();
            tree->QE.insert(std::make_pair(v, x));
        }
    }

    if (tree->V_old.find(v) == tree->V_old.end()) {
        auto V_near = GetNearbyNodes(v, tree->r, true);

        for (auto& w : V_near) {
            if (tree->E.find(std::make_pair(v, w)) == tree->E.end() &&
                g_estimated(v) + calc_dist(v, w) + h_estimated(w) < g_T[x_goal] &&
                g_T[v] + calc_dist(v, w) < g_T[w]) {
                tree->QE.insert(std::make_pair(v, w));
                if (g_T.find(w) == g_T.end()) {
                    g_T[w] = std::numeric_limits<double>::infinity();
                }
            }
        }
    }
}

double BITStar::BestVertexQueueValue() {
    if (tree->QV.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& v : tree->QV) {
        min_value = std::min(min_value, g_T[v] + h_estimated(v));
    }
    return min_value;
}

double BITStar::BestEdgeQueueValue() {
    if (tree->QE.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& [v, x] : tree->QE) {
        min_value = std::min(min_value, g_T[v] + calc_dist(v, x) + h_estimated(x));
    }
    return min_value;
}

Node* BITStar::BestInVertexQueue() {
    if (tree->QV.empty()) {
        std::cerr << "QV is Empty!" << std::endl;
        return nullptr;
    }
    Node* best_node = nullptr;
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& v : tree->QV) {
        double value = g_T[v] + h_estimated(v);
        if (value < min_value) {
            min_value = value;
            best_node = v;
        }
    }
    return best_node;
}

std::pair<Node*, Node*> BITStar::BestInEdgeQueue() {
    if (tree->QE.empty()) {
        std::cerr << "QE is Empty!" << std::endl;
        return {nullptr, nullptr};
    }
    std::pair<Node*, Node*> best_edge = {nullptr, nullptr};
    double min_value = std::numeric_limits<double>::infinity();
    for (auto& [v, x] : tree->QE) {
        double value = g_T[v] + calc_dist(v, x) + h_estimated(x);
        if (value < min_value) {
            min_value = value;
            best_edge = std::make_pair(v, x);
        }
    }
    return best_edge;
}

Eigen::Vector3d BITStar::SampleUnitNBall() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1, 1);
    while (true) {
        double x = dis(gen);
        double y = dis(gen);
        if (x * x + y * y < 1) {
            return Eigen::Vector3d(x, y, 0.0);
        }
    }
}

Eigen::Matrix3d BITStar::RotationToWorldFrame(state start, state goal, double L) {
    Eigen::Vector3d a1((goal.pos.x - start.pos.x) / L, (goal.pos.y - start.pos.y) / L, 0.0);//起始点到目标点的单位向量
    Eigen::Vector3d e1(1.0, 0.0, 0.0);//x轴单位向量
    Eigen::Matrix3d M = a1 * e1.transpose();//矩阵乘积
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);//矩阵奇异值分解
    Eigen::Matrix3d C = svd.matrixU() * Eigen::Matrix3d::Identity() * svd.matrixV().transpose();
    return C;
}

bool BITStar::is_collision(Node* start, Node* end) {
    // Implement collision detection logic here
    return false;
}

bool BITStar::is_inside_obs(Node* node) {
    // Implement obstacle detection logic here
    return false;
}

bool BITStar::is_in_range(Node* node, double delta) {
    double x_range[2] = {0.0, 50.0};
    double y_range[2] = {0.0, 30.0};
    return (x_range[0] + delta <= node->x && node->x <= x_range[1] - delta &&
            y_range[0] + delta <= node->y && node->y <= y_range[1] - delta);
}

double BITStar::calc_dist(Node* start, Node* end) {
    return std::hypot(start->x - end->x, start->y - end->y);
}

std::pair<double, double> BITStar::calc_dist_and_angle(Node* node_start, Node* node_end) {
    double dx = node_end->x - node_start->x;
    double dy = node_end->y - node_start->y;
    return {std::hypot(dx, dy), std::atan2(dy, dx)};
}

void BITStar::animation(const Eigen::Vector3d& xCenter, double cMax, double cMin, double theta) {
    plt::clf();
    plot_grid("Batch Informed Trees (BIT*)");

    for (auto& v : X_sample) {
        std::vector<double> x_vals = {v->x};
        std::vector<double> y_vals = {v->y};

        plt::scatter(x_vals, y_vals, 2.0, {{"color", "lightgrey"}});  // 使用scatter函数设置颜色和标记大小
    }

    if (cMax < std::numeric_limits<double>::infinity()) {
        draw_ellipse(xCenter, cMax, cMin, theta);
    }

    for (auto& [v, w] : tree->E) {
        plt::plot({v->x, w->x}, {v->y, w->y}, "g-");
    }

    plt::pause(0.001);
}

void BITStar::plot_grid(const std::string& name) {
    plt::title(name);
    plt::axis("equal");
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

double BITStar::random_double(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

std::vector<Node*> BITStar::GetNearbyNodes(Node* node, double radius, bool from_tree) {
    std::vector<Node*> nearby_nodes;
    const std::set<Node*>& nodes = from_tree ? tree->V : X_sample;
    for (auto& n : nodes) {
        if (calc_dist(node, n) <= radius) {
            nearby_nodes.push_back(n);
        }
    }
    return nearby_nodes;
}