#include <iostream>
#include <chrono>
#include <vector>
#include <list>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include "data.h"
#include "hungarian.h"
#include "kruskal.h"

struct node_info {
    std::vector<std::pair<int, int>> edges_illegal;
    //std::vector<std::vector<int>> subtour;
    //double bound_lower;
    //int subtour_selected;
    //bool node_fertility;
};

struct node_info2 {
    std::vector<std::pair<int, int>> edges_illegal;
    std::vector<int> subtour;
    bool fertility;
};

int cost_optimal;
int s_cost_optimal = INFINITE;
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

void cost_restriction(Matrix &cost, vii &edge_illegal){
    for(int i = 0; i < edge_illegal.size(); i++){
        int a = edge_illegal[i].first;
        int b = edge_illegal[i].second;
        cost[a][b] = INFINITE;
        //cost[b][a] = INFINITE;
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
    int index;
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

void Matrix_penalty_apply(Matrix & cost, std::vector<double> subgrad){
    for(int i = 0; i < cost.size(); i++)
        for(int j = 0; j < cost.size(); j++){
            cost[i][j] -= subgrad[i] + subgrad[j];
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

std::vector<int> subgrad_calc(vii edges, int dimension){
    std::vector<int> vec (dimension);

    for(int i = 0; i < dimension; i++){
        vec[i] = 2 - node_degree(edges, i); 
    }

    return vec;
}

double vec_sum(const std::vector<double> vec){
    double sum = 0;
    for(auto & i : vec)
        sum += i;

    return sum;
}

double step_size_calc(std::vector<int> subgrad, double lower_bound, double upper_bound, double epsilon){
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

void vec_penalty_update(std::vector<double> & vec_penalty, std::vector<int> subgrad, double epsilon){
    double a;
    for(int i = 0; i < vec_penalty.size(); i++){
        a = epsilon * subgrad[i];
        //a = ( a > 0 ? a : 0.0f);
        vec_penalty[i] += a;
    }
}

// =========  Depth Search Function BEGIN ========= 

void branch_and_bound_depth(Matrix cost, const std::vector<double> &vec_penalty, vii restriction, int gen){
    
    cost_restriction(cost, restriction);
    const int dimension = cost.size();
    Matrix cost_lagrange = cost;
    //std::cout << "after restrsiction   "  << "\n";

    std::vector<double> vec_penalty_new = vec_penalty;
    std::vector<int> subgrad (dimension);
    vii tree_edges;
    double step_size;
    double epsilon = 1.0f;
    double upper_bound = 148.0f;
    double obj_value;
    double obj_value_before;
    int count = 0 ;
    bool flag = false;

    int iteration_sum = 0;
    while(true){
        //if (count == 5)
            //exit(0);
        Matrix_penalty_apply(cost_lagrange, vec_penalty_new);
        matrix_print(cost_lagrange);
        Kruskal tree(cost_lagrange);

        // solve/generate the optimal 1-tree for the current penalties
        obj_value = tree.MST(dimension);
        obj_value += 2*vec_sum(vec_penalty_new);

        std::cout << "current cost   " << obj_value << "\n";

        if(upper_bound - obj_value <= 1)
            break;

        // if the value not increases after 30 iterations, epsilon is decreased
        if(flag && obj_value - obj_value_before - DBL_EPSILON <= 0){
            iteration_sum++;

            if(iteration_sum >= 30){
                iteration_sum = 0;
                epsilon /= 2;
            }
        }
        
        tree_edges = tree.getEdges();

        // subgradient 
        subgrad = subgrad_calc(tree_edges, dimension);
        vector_print_int(subgrad, "subgrad");

        step_size = step_size_calc(subgrad, obj_value, upper_bound, epsilon);
        std::cout << "step size   " << step_size << "\n";

        vec_penalty_update(vec_penalty_new, subgrad, step_size);
        vector_print_dbl(vec_penalty_new, "penalty");

        obj_value_before = obj_value;
        cost_lagrange = cost;
        flag = true;
        count++;
    }

    vector_print_pair(tree_edges, "result");

    std::cout << "current cost   " << obj_value << "\n";

    exit(0);


    //std::cout << "after treee   "  << "\n";

    //std::cout << "after solving   "  << "\n";


    /*
    std::cout << "\n\nGen   " << gen << "\n"; //"  Restriction  " << res.first << "   " <<  res.second << "\n";
    if(!restriction.empty()){
        vii copy = restriction;
        ii par = copy.back();
        int a = par.first;
        int b = par.second; 
        std::cout << "Restriction  " << a << "   " <<  b << "\n";
        matrix_print(cost_new);
    }
    vector_print(edges, " arvore");
    if(obj_value <= s_cost_optimal){
        //vector_print(restriction, "restriction");

        ii tree_node_degree_max = tree_node_degree_max_find(edges, dimension);
        int parent = tree_node_degree_max.first;
        //bool free = node_edges_possible(cost_new);

        //if(!free)
            //return;

        //std::cout << "degree max   " <<  tree_node_degree_max.second << "\n";
        if(tree_node_degree_max.second == 2){ 
            //std::cout << "degree max  yes  " << "\n";
            // valid solution
            s_cost_optimal = obj_value;
            exit(0);
        }else {
            //std::cout << "degree max noo   " << "\n";
            vii node_children = tree_node_children(edges, parent);
            //node_children_disp(node_children, parent, cost_new);
            
            std::sort(node_children.begin(), node_children.end(), 
                [](const ii &left, const ii &right){
                    return (left.first + left.second) < (right.first + right.second);
                }
            );
            
            // Debug ;;
                std::cout << "Pai :" << parent << "\nGrau: " << tree_node_degree_max.second << std::endl;
                vector_print(node_children, "filhos");


            //if(gen == 10)
                //exit(0);

            int gen_next = gen + 1;
            // invalid solution - new searches with distinct restritions 
            for(int i = 0; i < node_children.size(); i++){
                vii restriction_new = restriction;

                restriction_new.push_back(node_children[i]);

                branch_and_bound_depth(cost, restriction_new, gen_next);
            }
        }
    }//else
            */

}

// =========  Depth Search Function END ========= 


// =========  Breadth Search Functions BEGIN ========= 
/*
inline void node_solve(struct node_info2 &node, int **cost, int dimension){
    
    const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_problem_t p;

    int ** cost_new;
    cost_matrix_alloc(&cost_new, dimension);
    cost_restriction(node.edges_illegal, cost, cost_new, dimension);

    hungarian_init(&p, cost_new, dimension, dimension, mode); 
    cost_matrix_free(&cost_new, dimension);

    int node_cost = hungarian_solve(&p);

    if(node_cost <= s_cost_optimal){

        subtour_lower_get(node.subtour, &p, dimension);
        hungarian_free(&p);

        if(node.subtour.size() >= dimension){ 
            // valid solution
            s_cost_optimal = node_cost;
            s = node.subtour;
            node.fertility = false;
        }else   
            node.fertility = true;

    }else{
        hungarian_free(&p);
        node.fertility = false;
    }
}

inline void node_procreate(std::list<struct node_info2> & recipient, struct node_info2 & node_parent){
    for(int i = 0; i < node_parent.subtour.size() - 1; i++){
        struct node_info2 node_son;

        node_son.edges_illegal = node_parent.edges_illegal;

        std::pair<int, int> edge;
        edge.first = node_parent.subtour[i];
        edge.second = node_parent.subtour[i+1];

        node_son.edges_illegal.push_back(edge);
        recipient.push_back(node_son);
    }
}

void branch_and_bound_breadth(Data * data, int ** cost){

    const int dimension = data->getDimension();
    struct node_info2 node_root;
    std::list<struct node_info2> layer_A;
    std::list<struct node_info2> layer_B;

    layer_A.push_back(node_root);

    while(!layer_A.empty() || !layer_B.empty()){

        while(!layer_A.empty()){

            node_solve(layer_A.front(), cost, dimension);

            if(layer_A.front().fertility == true)
                // generate the sons in the layer B
                node_procreate(layer_B, layer_A.front());
            
            // kill/erase the node
            layer_A.pop_front();
            
        }

        while(!layer_B.empty()){

            node_solve(layer_B.front(), cost, dimension);

            if(layer_B.front().fertility == true)
                // generate the sons in the layer A
                node_procreate(layer_A, layer_B.front());
            
            // kill/erase the node
            layer_B.pop_front();
            
        }

    }

}

// =========  Breadth Search Functions END ========= 


*/



int main(int argc, char** argv) {

    const char * breadth_flag = "--breadth";
    const char * depth_flag = "--depth";
    const char * bound_opt = "-o";

    /*
    Data * data = new Data(argc, argv[1]);
    data->readData();
    int dimension = data->getDimension();
    */
    int dimension = 5;

    //Matrix cost (dimension, std::vector<double>(dimension));
    Matrix cost =   {{ 99999 , 30 , 26 , 50 , 40},
             { 30 , 99999 , 24 , 40 , 50},
             { 26 , 24 , 99999 , 24 , 26},
             { 50 , 40 , 24 , 99999 , 30},
             { 40 , 50 , 26 , 30 , 99999}};

    /*
    for (int i = 0; i < data->getDimension(); i++){
        //cost[i] = new int[data->getDimension()];
        for (int j = 0; j < data->getDimension(); j++){
            cost[i][j] = data->getDistance(i, j);
        }
    }
    */


    /*
    Kruskal tree(cost);


    std::cout << "Cost: " << tree.MST(dimension) << std::endl;
    vii arvore = tree.getEdges();

    for(int i =0; i < arvore.size(); i++){
        std::cout << arvore[i].first << " " << arvore[i].second << std::endl;
    }
    std::cout << std::endl;
    */
    std::vector<double> penalty = {0,0,0,0,0};
    vii restriction;
    branch_and_bound_depth(cost, penalty, restriction, 0);

    //std::cout << "COST: " << s_cost_optimal << std::endl;

/*
    cost_optimal = cost_optimal_get(argv[1]);

    if(argc == 4 && !strcmp(bound_opt, argv[3]))
        s_cost_optimal = cost_optimal;

    struct node_info node_initial;

    auto t1 = std::chrono::high_resolution_clock::now();

    if(!strcmp(breadth_flag, argv[2]))
        branch_and_bound_breadth(data, cost);
    else if(!strcmp(depth_flag, argv[2]))
        branch_and_bound_depth(data, cost, node_initial, 0);
    else{
        std::cout << "[!!] Error: Search Mode not identified\n" << std::endl;
        exit(-1);
    }

    auto t2 = std::chrono::high_resolution_clock::now();


    for(int i = 0; i < s.size(); i++){
        std::cout << s[i] << " ";
    }
    std::cout << std::endl;

    if(!s.empty())
        std::cout << "COST: " << s_cost_optimal << std::endl;
    else
        std::cout << "Couldnt find a solution." << std::endl;

    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "TIME: " << exec_time / 10e2  << std::endl;
    */

    return 0;
}
