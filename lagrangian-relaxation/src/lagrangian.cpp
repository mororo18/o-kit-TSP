#include "lagrangian.h"

#define INFINITE 999999

void cost_restriction(Matrix & cost, const vii & edge_illegal){
    for(int i = 0; i < edge_illegal.size(); i++){
        int a = edge_illegal[i].first;
        int b = edge_illegal[i].second;
        cost[a][b] = INFINITE;
        cost[b][a] = INFINITE;
    }
}

void Matrix_penalty_apply(Matrix & cost, std::vector<double> vec){
    /*
    for(int j = 1; j < cost.size(); j++){
        
        cost[0][j] -= vec[j];
        cost[j][0] -= vec[j];
    }
    */

    for(int i = 0; i < cost.size(); i++){
        double vec_i = vec[i];
        for(int j = i+1; j < cost.size(); j++){
            cost[i][j] -= vec_i + vec[j];
            cost[j][i] -= vec_i + vec[j];
        }
    }
}

int node_degree(const vii & edges, int node){  
    int degree = 0;
    for(int i = 0; i < edges.size(); i++){
        if(edges[i].first == node || edges[i].second == node)
            degree++;
    }

    return degree;
}

std::vector<int> subgrad_calc(const vii & edges, int dimension){
    std::vector<int> vec (dimension);

    for(int i = 0; i < dimension; i++){
        vec[i] = 2 - node_degree(edges, i); 
    }

    return vec;
}

double step_size_calc(const std::vector<int> & subgrad, double lower_bound, double upper_bound, double epsilon){
    double step_size;
    double dif = upper_bound - lower_bound;
    int d = 0;

    for(auto &i: subgrad){
        d += i * i;
        //std::cout << "square  " << (i << 1) << std::endl;
    }
    
    step_size = epsilon * dif / d; 
    
    return step_size;
}

void vec_penalty_update(std::vector<double> & vec_penalty, const std::vector<int> & subgrad, double step_size){
    double a;
    for(int i = 0; i < vec_penalty.size(); i++){
        a = step_size * subgrad[i];
        //double diff = vec_penalty[i] - a;
        //a = (diff > 0 ? diff : 0.0f);
        vec_penalty[i] += a;
    }
}

bool subgrad_validate(const std::vector<int> & subgrad){
    for(int i = 0; i < subgrad.size(); i++)
        if(subgrad[i] != 0)
            return false;

    return true;
}

void lagrangian_dual_solve(const Matrix & cost, struct node_info & node){

    const int dimension = cost.size();

    Matrix cost_cpy = cost;
    vii restriction = node.edges_illegal;
    cost_restriction(cost_cpy, restriction);
    Matrix cost_lagrange = cost_cpy;

    std::vector<double> vec_penalty_new = node.vec_penalty;
    std::vector<double> vec_penalty_best;
    std::vector<int> subgrad (dimension);
    std::vector<int> subgrad_max (dimension);
    vii tree_edges;
    vii best_edges;
    double step_size;
    double epsilon = 1.0f;
    if(dimension > 77)
        epsilon = 2.0f;
    double obj_value;
    double obj_value_max = 0;
    int count = 0 ;
    int iteration_sum = 0;

    while(true){
        Matrix_penalty_apply(cost_lagrange, vec_penalty_new);
        Kruskal tree(cost_lagrange);

        // solve/generate the optimal 1-tree with the current penalties
        obj_value = tree.MST(dimension);
        tree_edges = tree.getEdges();

        // subgradient 
        subgrad = subgrad_calc(tree_edges, dimension);

        if(obj_value_max < obj_value){
            obj_value_max = obj_value;
            best_edges = tree_edges;
            vec_penalty_best = vec_penalty_new;
            subgrad_max  = subgrad;
            //std::cout << obj_value << std::endl;
        }

        //vector_print_int(subgrad, "subgrad");

        // break condition
        if(subgrad_validate(subgrad)){
            obj_value_max = obj_value;
            best_edges = tree_edges;
            vec_penalty_best = vec_penalty_new;
            subgrad_max  = subgrad;
            break;
        }

        // break condition
        //if(primal_bound <= obj_value)
            //break;
        //if(obj_value > s_cost_optimal){
            //break;
        //}

        if(obj_value <= obj_value_max ){
            iteration_sum++;

            // if the value not increases after 30 iterations, epsilon is decreased
            if(iteration_sum >= 30){
                iteration_sum = 0;
                epsilon /= 2;

                if(epsilon < 0.0005f)
                    break;
            }

        }else
            iteration_sum = 0;

        if(primal_bound - obj_value_max <= TOL){
            //std::cout << "break yx" << std::endl;
            //std::cout << obj_value << std::endl;
            //exit(0);
            break;
        }

        step_size = step_size_calc(subgrad, obj_value, primal_bound, epsilon);
        //std::cout << "step size   " << step_size << "\n";

        vec_penalty_update(vec_penalty_new, subgrad, step_size);
        //vector_print_dbl(vec_penalty_new, "penalty");

        cost_lagrange = cost_cpy;
        count++;
    }

    //struct node_info node = node_parent;

    node.tree = best_edges;
    node.vec_penalty = vec_penalty_best;
    node.cost = obj_value_max;
    //node.edges_illegal = node_parent.edges_illegal;
    node.subgrad = subgrad_max;

}
