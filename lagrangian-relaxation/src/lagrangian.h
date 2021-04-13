#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include <vector>
#include "kruskal.h"

struct node_info {
    std::vector<std::pair<int, int>> edges_illegal;
    vii edges_illegal_new;;
    vii tree;
    std::vector<double> vec_penalty;
    std::vector<int> subgrad;
    double cost;
    bool fertility;
};

extern double primal_bound;
extern double s_cost_optimal;

void Matrix_penalty_apply(Matrix & cost, std::vector<double> vec);
int node_degree(const vii & edges, int node);
std::vector<int> subgrad_calc(const vii & edges, int dimension);
double step_size_calc(const std::vector<int> & subgrad, double lower_bound, double upper_bound, double epsilon);
void vec_penalty_update(std::vector<double> & vec_penalty, const std::vector<int> & subgrad, double step_size);
bool subgrad_validate(const std::vector<int> & subgrad);
void lagrangian_dual_solve(const Matrix & cost, struct node_info & node_parent);
void cost_restriction(Matrix & cost, const vii & edge_illegal);

#endif
