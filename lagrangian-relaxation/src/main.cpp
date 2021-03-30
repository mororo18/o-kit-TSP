#include <iostream>
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

void Matrix_penalty_apply(Matrix & cost, std::vector<double> vec){
    for(int i = 0; i < cost.size(); i++)
        for(int j = 0; j < cost.size(); j++){
            cost[i][j] -= vec[i] + vec[j];
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

void vec_penalty_update(std::vector<double> & vec_penalty, std::vector<int> subgrad, double step_size){
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

    //vector_print_pair(graph, "result");

    return cost_total;
}

// =========  Depth Search Function BEGIN ========= 

void branch_and_bound_depth(Matrix cost, const std::vector<double> &vec_penalty, vii restriction, int gen){
    
    const int dimension = cost.size();
    cost_restriction(cost, restriction);
    Matrix cost_lagrange = cost;
    //std::cout << "after restrsiction   "  << "\n";

    std::vector<double> vec_penalty_new = vec_penalty;
    std::vector<double> vec_penalty_best;
    std::vector<int> subgrad (dimension);
    vii tree_edges;
    double step_size;
    double epsilon = 2.0f;
    //upper_bd = primal_bound(cost);
    double obj_value;
    double obj_value_before;
    double obj_value_max;
    int count = 0 ;
    bool flag = false;
    vii best_edges;

    int iteration_sum = 0;
    while(true){
        Matrix_penalty_apply(cost_lagrange, vec_penalty_new);
        Kruskal tree(cost_lagrange);
        //matrix_print(cost_lagrange);

        // solve/generate the optimal 1-tree for the current penalties
        obj_value = tree.MST(dimension);

        /*
        std::cout << "iteration   " << iteration_sum << "\n";
        std::cout << "epsilon   " << epsilon << "\n";
        std::cout << "value   " << obj_value << "\n";
        */

        // break condition
        if(upper_bd - obj_value <= 1)
            break;
        if(obj_value > s_cost_optimal)
            break;

        tree_edges = tree.getEdges();
        ii tree_node_degree_max = tree_node_degree_max_find(tree_edges, dimension);

        if(!flag){
            obj_value_max = obj_value;
            best_edges = tree_edges;
            vec_penalty_best = vec_penalty_new;
        }else if(obj_value_max < obj_value){
            obj_value_max = obj_value;
            best_edges = tree_edges;
            vec_penalty_best = vec_penalty_new;
        }

        if(flag && obj_value <= obj_value_max )
            iteration_sum++;
        else
            iteration_sum = 0;

        // if the value not increases after 30 iterations, epsilon is decreased
        if(iteration_sum >= 20){
            iteration_sum = 0;
            epsilon /= 2;
        }

        // break condition
        if(epsilon < 0.0005f)
            break;
        

        // break condition
        if(tree_node_degree_max.second == 2)
            break;

        // subgradient 
        subgrad = subgrad_calc(tree_edges, dimension);
        //vector_print_int(subgrad, "subgrad");

        step_size = step_size_calc(subgrad, obj_value, upper_bd, epsilon);
        //std::cout << "step size   " << step_size << "\n";

        vec_penalty_update(vec_penalty_new, subgrad, step_size);
        //vector_print_dbl(vec_penalty_new, "penalty");

        cost_lagrange = cost;
        count++;
        obj_value_before = obj_value;
        flag = true;
    }

    //vector_print_pair(tree_edges, "result");


    std::cout << "\n" << "GEN   " << gen << "\n";
    std::cout << "Limitante primal " << upper_bd << std::endl;
    std::cout << "Iter max   " << count << "\n";
    std::cout << "current cost   " << obj_value << "\n";
    vector_print_int(subgrad, "Subgrad");


    //exit(0);

    tree_edges = best_edges;
    vec_penalty_new = vec_penalty_best;

    if(obj_value < s_cost_optimal){
        //vector_print(restriction, "restriction");

        //upper_bd = ceil(obj_value);
        ii tree_node_degree_max = tree_node_degree_max_find(tree_edges, dimension);
        int parent = tree_node_degree_max.first;
        int parent_degree = tree_node_degree_max.second;
        //bool free = node_edges_possible(cost_new);

        //if(!free)
            //return;

        //std::cout << "degree max   " <<  tree_node_degree_max.second << "\n";
        if(parent_degree == 2){ 
            //std::cout << "degree max  yes  " << "\n";
            // valid solution
            s_cost_optimal = obj_value;
            std::cout << "current cost   " << obj_value << "\n";
            vector_print_pair(tree_edges, "solution ");
            //exit(0);
        }else if(parent_degree > 2){
            //std::cout << "degree max noo   " << "\n";
            vii node_children = tree_node_children(tree_edges, parent);
            //node_children_disp(node_children, parent, cost_new);

            
               /*
               std::partial_sort(node_children.begin(), node_children.begin() + 2, node_children.end(), 
               [](const ii &left, const ii &right){
               return left.first + left.second < right.first + right.second;
               }
               );
            
            */

            // Debug ;;
            /*
                std::cout << "Pai :" << parent << "\nGrau: " << parent_degree << std::endl;
                //vector_print(node_children, "filhos");
                */


            //if(gen == 10)
                //exit(0);

            int gen_next = gen + 1;
            // invalid solution - new searches with distinct restritions 
            for(int i = 0; i < node_children.size(); i++){
                vii restriction_new = restriction;

                restriction_new.push_back(node_children[i]);

                branch_and_bound_depth(cost, vec_penalty_new, restriction_new, gen_next);
            }
        }
    }//else

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

    


    Data * data = new Data(argc, argv[1]);
    data->readData();
    int dimension = data->getDimension();
    
    Matrix cost (dimension, std::vector<double>(dimension));

    
    for (int i = 0; i < data->getDimension(); i++){
        //cost[i] = new int[data->getDimension()];
        for (int j = 0; j < data->getDimension(); j++){
            cost[i][j] = data->getDistance(i, j);
        }
    }
    
    /*

    int dimension = 5;

    Matrix cost =   {{ 99999 , 30 , 26 , 50 , 40},
             { 30 , 99999 , 24 , 40 , 50},
             { 26 , 24 , 99999 , 24 , 26},
             { 50 , 40 , 24 , 99999 , 30},
             { 40 , 50 , 26 , 30 , 99999}};
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
    //std::vector<double> penalty = {0,0,0,0,0};
    std::vector<double> penalty (dimension, 0); //= {0,0,0,0,0};
    vii restriction;
    upper_bd = primal_bound(cost);
    branch_and_bound_depth(cost, penalty, restriction, 0);

    std::cout << "COST: " << s_cost_optimal << std::endl;
    //vii test = {make_pair(0,1), make_pair(1,2), make_pair(2,3), make_pair(3,2)};
    //std::cout << "primal bound: " << primal_bound(cost) << std::endl;

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
