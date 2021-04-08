#include <iostream>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <list>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "data.h"
#include "hungarian.h"
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


int cost_optimal;
double s_cost_optimal = INFINITE;
double upper_bd;
std::vector<int> s;

int cost_optimal_get(char * instance_name){
    string path = "benchmark/target_data";
    int value;
    string file;
    ifstream inTSP(path, ios::in);

    inTSP >> file; while ( file.find((string)instance_name+":") != 0) {
        inTSP >> file;
    }

    int pos = file.find(":") + 1;
    string value_str = file.substr(pos);

    value = stoi(value_str);

    return value;
}

void matrix_print(const Matrix &cost){
    for(int i = 0; i < cost.size(); i++){
        std::cout << i << " | ";
        for(int j = 0; j < cost.size(); j++){
            std::cout  << cost[i][j] << " ";

        }
        std::cout << std::endl;
    }
}

void cost_matrix_alloc(int *** matrix, int dimension){
    int ** ptr =  (int**)calloc(dimension, sizeof(int*));
    hungarian_test_alloc(ptr);
    for(int i = 0; i < dimension; i++){
        ptr[i] =  (int*)calloc(dimension, sizeof(int));
        hungarian_test_alloc(ptr[i]);
    }

    *matrix = ptr;
}

inline void cost_matrix_free(int *** ptr, int dimension){
    int ** ar = *ptr;
    for(int i = 0; i < dimension; i++)
        free(ar[i]);
    free(ar);

    *ptr = NULL;
}

void cost_restriction(Matrix & cost, const vii & edge_illegal){
    for(int i = 0; i < edge_illegal.size(); i++){
        int a = edge_illegal[i].first;
        int b = edge_illegal[i].second;
        cost[a][b] = INFINITE;
        cost[b][a] = INFINITE;
    }
}

void vector_print_pair(const vii &vec, std::string name){

    std::cout << "Vector " << name << " :  ";
    for(int i = 0; i < vec.size(); i++){
        std::cout << vec[i].first << " "  << vec[i].second << " |";
    }
    std::cout << std::endl;
}

void vector_print_int(const vector<int> &vec, std::string name){

    std::cout << "Vector " << name << " :  ";
    for(int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void vector_print_dbl(const vector<double> &vec, std::string name){

    std::cout << "Vector " << name << " :  ";
    for(int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

inline ii array_item_max(int * arr, int size){
    int max = arr[0];
    int index = 0;
    ii node_and_degree;
    for(int i = 1; i < size; i++){
        if(arr[i] > max){
            max = arr[i];
            index = i;
        }
    }
    node_and_degree.first = index;
    node_and_degree.second = max;
    return node_and_degree;
}

ii tree_node_degree_max_find(vii edge, int dimension){
    int node_degree[dimension] = {0};
    ii edge_current;

    for(int i = 0; i < edge.size(); i++){
        edge_current = edge[i];
        node_degree[edge_current.first]++;
        node_degree[edge_current.second]++;
    }

    return array_item_max(node_degree, dimension);
}

vii tree_node_children(vii edges, int parent){
    vii children;
        for(int i = 0; i < edges.size(); i++){
            if(edges[i].first == parent || edges[i].second == parent)
                children.push_back(edges[i]);
        }

    return children;
}

void vector_pair_sort(vii & vec, const Matrix & cost){
    for(int i = 0; i < vec.size() -1; i++)
        for(int j = i+1; j < vec.size(); j++){
            double i_cost = cost[vec[i].first][vec[i].second];
            double j_cost = cost[vec[j].first][vec[j].second];

            if(j_cost > i_cost){
                ii aux = vec[j];
                vec[j] = vec[i];
                vec[i] = aux;
            }
            
        }
}

vii tree_nodes_children(const Matrix & cost, vii edges, std::vector<int> parent){
//vii tree_node_children(vii edges, int parent){
    std::vector<vii> kids;
    int sz = parent.size();
    double parent_costs[sz] = {0};

    for(int j = 0; j < parent.size(); j++){
        vii children;
        for(int i = 0; i < edges.size(); i++){
            if(edges[i].first == parent[j] || edges[i].second == parent[j]){
                children.push_back(edges[i]);
                parent_costs[j] += cost[edges[i].first][edges[i].second];
            }
        }
        //vector_print_pair(children, "childs");
        //std::cout << "cost "<< parent_costs[j] << std::endl;
        //vector_pair_sort(children, cost);
        kids.push_back(children);
    }

    int i_min;
    double cost_min = INFINITE;

    /*
    for(int i = 0; i < parent.size(); i++){
        if(cost[kids[i][0].first][kids[i][0].second] < cost_min){
            cost_min =cost[kids[i][0].first][kids[i][0].second]; 
            i_min = i;
        }
    }
    */
    for(int i = 0; i < parent.size(); i++){
        if(parent_costs[i] < cost_min){
            cost_min = parent_costs[i];
            i_min = i;
        }
    }

    //vector_print_pair(kids[i_min], "selected");
    //std::cout << "cost selected "<< parent_costs[i_min] << std::endl;

    return kids[i_min];
}

void Matrix_penalty_apply(Matrix & cost, std::vector<double> vec){
    for(int i = 0; i < cost.size(); i++){
        double vec_i = vec[i];
        for(int j = i; j < cost.size(); j++){
            cost[i][j] -= vec_i + vec[j];
            cost[j][i] -= vec_i + vec[j];
        }
    }
}

inline int node_degree(const vii & edges, int node){  
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

std::vector<int> tree_nodes_degree_max(std::vector<int> subgrad){
    std::vector<int> parents;
    int min = 0;
    //find the min subgrad (max degree)
    for(int i = 0; i < subgrad.size(); i++){
        if(subgrad[i] < min){
            min = subgrad[i];
        }
    }

    parents.push_back(2 - min);

    for(int i = 0; i < subgrad.size(); i++){
        if(subgrad[i] == min){
            parents.push_back(i);
        }
    }

    return parents;

}

double vec_sum(const std::vector<double> vec){
    double sum = 0;
    for(auto & i : vec)
        sum += i;

    return sum;
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

bool cycle_search_path(vii edges, int v, bool visited[], int parent)
{

    // Mark the current node as visited
    visited[v] = true;

    // Recur for all the vertices
    // adjacent to this vertex
    vii node_adjacent = tree_node_children(edges, v) ;
    for (int i = 0; i < node_adjacent.size(); ++i)
    {
        int adj;
        if(node_adjacent[i].first == v)
            adj = node_adjacent[i].second;
        if(node_adjacent[i].second == v)
            adj = node_adjacent[i].first;

        // If an adjacent vertex is not visited,
        //then recur for that adjacent
        if (!visited[adj])
        {
            if (cycle_search_path(edges, adj, visited, v))
                return true;
        }

        // If an adjacent vertex is visited and
        // is not parent of current vertex,
        // then there exists a cycle in the graph.
        else if (adj != parent)
            return true;
    }
    return false;
}

bool cycle_search(vii edges, int dimension)
{

    // Mark all the vertices as not
    // visited and not part of recursion
    // stack
    bool visited[dimension];
    for (int i = 0; i < dimension; i++)
        visited[i] = false;

    // Call the recursive helper
    // function to detect cycle in different
    // DFS trees
    for (int u = 0; u < dimension; u++)
    {

        // Don't recur for u if
        // it is already visited
        if (!visited[u])
            if (cycle_search_path(edges, u, visited, -1))
                return true;
    }
    return false;
}

double primal_bound(Matrix cost){
    const int dimension = cost.size();
    std::priority_queue<std::pair<double, ii>> edges;

    for(int i = 0; i < dimension-1; i++){
        for(int j = 0; j < dimension; j++){
           edges.push( make_pair(-cost[i][j], make_pair(i, j)) ); 
        }
    }

    vii graph;
    double cost_total = 0;
    //graph.size() < dimension
    while(graph.size() < dimension){
        vii graph_cpy = graph; 
        std::pair<double, ii> edge_current = edges.top();
        edges.pop();

        graph_cpy.push_back(edge_current.second);
        int a = edge_current.second.first;
        int b = edge_current.second.second;
        double edge_cost = -edge_current.first;
        //std::cout << a << " " << b << " | " << "cost  " << edge_cost << std::endl;

        if(node_degree(graph, a) < 2 
        && node_degree(graph, b) < 2 
        && (!cycle_search(graph_cpy, dimension) || graph.size() == dimension - 1)){
            graph = graph_cpy;
            cost_total += edge_cost;
        }

    }
    
    /*
    std::vector<int> subgrad (dimension);
    subgrad = subgrad_calc(graph, dimension);
    vector_print_int(subgrad, "subgrad");
    */

    return cost_total;
}

bool subgrad_validate(const std::vector<int> & subgrad){
    for(int i = 0; i < subgrad.size(); i++)
        if(subgrad[i] != 0)
            return false;

    return true;
}

vii result;

struct node_info lagrangian_dual_solve(const Matrix & cost, struct node_info node_parent){

    const int dimension = cost.size();

    Matrix cost_cpy = cost;
    vii restriction = node_parent.edges_illegal;
    cost_restriction(cost_cpy, restriction);
    Matrix cost_lagrange = cost_cpy;

    std::vector<double> vec_penalty_new = node_parent.vec_penalty;
    std::vector<double> vec_penalty_best;
    std::vector<int> subgrad (dimension);
    std::vector<int> subgrad_max (dimension);
    vii tree_edges;
    vii best_edges;
    double step_size;
    double epsilon = 2.0f;
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
        if(upper_bd  <= obj_value)
            break;
        if(obj_value > s_cost_optimal){
            break;
        }

        if(obj_value <= obj_value_max ){
            iteration_sum++;

            // if the value not increases after 20 iterations, epsilon is decreased
            if(iteration_sum >= 30){
                iteration_sum = 0;
                epsilon /= 2;

                if(epsilon < 0.0005f)
                    break;
            }

        }else
            iteration_sum = 0;

        if((upper_bd - obj_value_max)/upper_bd < 0.0001)
            break;

        step_size = step_size_calc(subgrad, obj_value, upper_bd, epsilon);
        //std::cout << "step size   " << step_size << "\n";

        vec_penalty_update(vec_penalty_new, subgrad, step_size);
        //vector_print_dbl(vec_penalty_new, "penalty");

        cost_lagrange = cost_cpy;
        count++;
    }

    struct node_info node;

    node.tree = best_edges;
    node.vec_penalty = vec_penalty_best;
    node.cost = obj_value_max;
    node.edges_illegal = node_parent.edges_illegal;
    node.subgrad = subgrad_max;

    return node;
}
// =========  Depth Search Function BEGIN ========= 

void branch_and_bound_depth(const Matrix cost, struct node_info node_parent, int gen){

    const int dimension = cost.size();
    
    struct node_info node; 
    node = lagrangian_dual_solve(cost, node_parent);

    vii tree_edges = node.tree;
    double obj_value = node.cost;

    /*
    std::cout << "\n" << "GEN   " << gen << "\n";
    std::cout << "Limitante primal " << upper_bd << std::endl;
    std::cout << "current cost   " << obj_value << "\n";
    */
    //vector_print_int(subgrad, "Subgrad");
    //vector_print_pair(restriction, "restriction");


    //exit(0);

    if(obj_value + DBL_EPSILON < s_cost_optimal - 1){
        //vector_print(restriction, "restriction");

        ii tree_node_degree_max = tree_node_degree_max_find(tree_edges, dimension);
        std::vector<int> parents = tree_nodes_degree_max(node.subgrad);
        //int parent = tree_node_degree_max.first;
        int parent_degree = parents[0];
        parents.erase(parents.begin());
        //int parent_degree = tree_node_degree_max.second;

        // valid solution
        if(parent_degree == 2){ 
            
            s_cost_optimal = obj_value;
            upper_bd = obj_value;
            std::cout << "current cost  solution  " << obj_value << "\n";
            //vector_print_pair(tree_edges, "solution ");
            result = tree_edges;
        }else if(parent_degree > 2){
            //vii node_children = tree_node_children(tree_edges, parent);
            /*
            vector_print_pair(tree_edges, " | ");
            for(int i = 0; i < parents.size(); i++){
                
                std::cout << parents[i] << " ";
                vector_print_pair(tree_node_children(tree_edges, parents[i]), " | ");
            }
            */
            vii node_children = tree_nodes_children(cost, tree_edges, parents);
            //exit(0);
            vector_pair_sort(node_children, cost);
            //vector_print_pair(node_children, " | ");
            //vector_print_int(node.subgrad, "subgrad");
            //exit(0);

            int gen_next = gen + 1;
            // invalid solution -> new branchs for each new restriction
            for(int i = 0; i < node_children.size(); i++){
                struct node_info node_child = node;

                node_child.edges_illegal.push_back(node_children[i]);

                //vector_print_pair(node_child.edges_illegal, " | ");
                //if(gen == 3) exit(0);

                branch_and_bound_depth(cost, node_child, gen_next);
                //node_children.erase(node_children.begin() + i);

            }
        }
    }

}

// =========  Depth Search Function END ========= 


// =========  Breadth Search Functions BEGIN ========= 

inline void node_solve(struct node_info &node_parent, const Matrix & cost){

    const int dimension = cost.size();
    
    struct node_info node; 
    node = lagrangian_dual_solve(cost, node_parent);

    vii tree_edges = node.tree;
    double obj_value = node.cost;

    //std::cout << obj_value << std::endl;

    if(ceil(obj_value + DBL_EPSILON)  < s_cost_optimal - 1){

        ii tree_node_degree_max = tree_node_degree_max_find(tree_edges, dimension);
        std::vector<int> parents = tree_nodes_degree_max(node.subgrad);
        int parent_degree = parents[0];
        parents.erase(parents.begin());

        // valid solution
        if(parent_degree == 2){ 
            
            s_cost_optimal = obj_value;
            upper_bd = obj_value;
            std::cout << "current cost  solution  " << obj_value << "\n";
            //vector_print_pair(tree_edges, "solution ");
            result = tree_edges;
            node_parent.fertility = false;
        }else if(parent_degree > 2){   
            vii node_children = tree_nodes_children(cost, tree_edges, parents);
            vector_pair_sort(node_children, cost);
            //vector_pair_sort(node_children, cost);
            node_parent.edges_illegal_new = node_children;
            node_parent.fertility = true;
        }

    }else
        node_parent.fertility = false;
    
}

inline void node_procreate(std::list<struct node_info> & recipient, struct node_info & node_parent){
    for(int i = 0; i < node_parent.edges_illegal_new.size(); i++){
        struct node_info node_son;

        node_son.vec_penalty = node_parent.vec_penalty;
        node_son.edges_illegal = node_parent.edges_illegal;

        std::pair<int, int> edge = node_parent.edges_illegal_new[i];

        node_son.edges_illegal.push_back(edge);
        recipient.push_back(node_son);
    }
}

void branch_and_bound_breadth(const Matrix cost, int gen){

    const int dimension = cost.size();
    std::vector<double> vec_penalty_default (dimension, 0);
    struct node_info node_root;
    node_root.vec_penalty = vec_penalty_default;
    std::list<struct node_info> layer_A;
    std::list<struct node_info> layer_B;

    layer_A.push_back(node_root);
    int count = 0;

    while(!layer_A.empty() || !layer_B.empty()){

        std::cout << "Gen " << count++ << " size " << layer_A.size() << std::endl;
        while(!layer_A.empty()){
            node_solve(layer_A.front(), cost);

            if(layer_A.front().fertility == true)
                // generate the sons in the layer B
                node_procreate(layer_B, layer_A.front());
            
            // kill/erase the node
            layer_A.pop_front();
            
        }

        std::cout << "Gen " << count++ << " size " << layer_B.size() << std::endl;

        while(!layer_B.empty()){
            node_solve(layer_B.front(), cost);

            if(layer_B.front().fertility == true)
                // generate the sons in the layer A
                node_procreate(layer_A, layer_B.front());
            
            // kill/erase the node
            layer_B.pop_front();
            
        }

    }

}

// =========  Breadth Search Functions END ========= 

double cost_calc(vii & edges, Matrix & cost){
    double total = 0;

    for(auto i : edges){
        int a  = i.first;
        int b  = i.second;

        total += cost[a][b];
    }

    return total;
}

int main(int argc, char** argv) {

    const char * breadth_flag = "--breadth";
    const char * depth_flag = "--depth";
    const char * bound_opt = "-o";

    //srand(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    Data * data = new Data(argc, argv[1]);
    data->readData();
    int dimension = data->getDimension();
    
    Matrix cost (dimension, std::vector<double>(dimension));

    for (int i = 0; i < data->getDimension(); i++){
        for (int j = 0; j < data->getDimension(); j++){
            cost[i][j] = data->getDistance(i, j);
        }
    }

    std::vector<double> penalty (dimension, 0); 
    struct node_info node_root;
    node_root.vec_penalty = penalty;
    upper_bd = primal_bound(cost);

    auto t1 = std::chrono::high_resolution_clock::now();

    if(!strcmp(breadth_flag, argv[2]))
        branch_and_bound_breadth(cost, 0);
    else if(!strcmp(depth_flag, argv[2]))
        branch_and_bound_depth(cost, node_root, 0);
    else{
        std::cout << "[!!] Error: Search Mode not identified\n" << std::endl;
        exit(-1);
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "COST: " << s_cost_optimal << std::endl;

    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "TIME: " << exec_time / 10e2  << std::endl;

    return 0;
}
