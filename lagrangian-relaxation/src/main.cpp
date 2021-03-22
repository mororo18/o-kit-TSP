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
        for(int j = 0; j < cost.size(); j++){
            std::cout << cost[i][j] << " ";

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

void cost_restriction(Matrix &cost, ii edge_illegal){
    // cpy the old to the new matrix of costs
    int a = edge_illegal.first;
    int b = edge_illegal.second;
    cost[a][b] = INFINITE;
}


bool subtour_cmp(std::vector<int> &a, std::vector<int> &b){
    return (a.size() < b.size());// && a[0] < b[0];
}

void subtour_lower_get(std::vector<int> &subtour, hungarian_problem_t * solution, int dimension){
    std::vector<std::vector<int>> tour(1, std::vector<int> (1, 1));
    bool nodes_known[dimension] = {0};
    bool flag;
    int count = 0;
    int index = 0;
    int tour_node_first = 1;
    nodes_known[index] = true;

    while(true){

        //std::cout << index+1 << " ";
        // identify the edges
        for(int i = 0; i < dimension; i++){
            if(solution->assignment[index][i]){
                nodes_known[index] = true;
                index = i;
                break;
            }
        }

        tour[count].push_back(index+1);

        // identify the end of the subtour 
        if(tour_node_first == index+1){
            //std::cout << index+1 << " ";
            //std::cout << std::endl;
            count++;
            flag = false;
            // take the next free node and creates a new subtour
            for(int i = 1; i < dimension; i++)
                if(!nodes_known[i]){
                    tour_node_first = i+1;
                    tour.push_back(std::vector<int> (1, tour_node_first));
                    //tour[count].push_back(tour_node_first);
                    index = i;
                    flag = true;
                    break;
                }

            // when there is no avaliable node for a new tour
            if(!flag)
                break;
        }
    }

    //std::sort(tour.begin(),tour.end(), subtour_cmp); 
    std::partial_sort(tour.begin(), tour.begin()+1,tour.end(), subtour_cmp); 

    /*
    std::cout << std::endl;
    for(int i = 0; i < tour.size(); i++){
        for(int j = 0; j < tour[i].size(); j++)
            std::cout << tour[i][j] << " ";
        std::cout << std::endl;
    }
    */

    subtour = tour[0];
}

void vector_print(const std::vector<int> &vec, std::string name){

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

std::vector<int> tree_node_children(vii edges, int parent){
    std::vector<int> children;
    for(int i = 0; i < edges.size(); i++){
        if(edges[i].first == parent)
            children.push_back(edges[i].second);
        if(edges[i].second == parent)
            children.push_back(edges[i].first);
    }

    return children;
}

// =========  Depth Search Function BEGIN ========= 

void branch_and_bound_depth(Matrix cost, ii restriction, int gen){
    std::cout << "Gen   " << gen << "\n";

    cost_restriction(cost, restriction);
    std::cout << "after restrsiction   "  << "\n";
    Kruskal tree(cost);
    std::cout << "after treee   "  << "\n";
    const int dimension = cost.size();

    int obj_value = tree.MST(dimension);
    std::cout << "after solving   "  << "\n";

    std::cout << "current cost   " << obj_value << "\n";
    if(obj_value <= s_cost_optimal){

        
        ii tree_node_degree_max = tree_node_degree_max_find(tree.getEdges(), dimension);
        std::cout << "degree max   " <<  tree_node_degree_max.second << "\n";
        if(tree_node_degree_max.second == 2){ 
        std::cout << "degree max  yes  " << "\n";
            // valid solution
            s_cost_optimal = obj_value;
        }else {
            std::cout << "degree max noo   " << "\n";
            int parent = tree_node_degree_max.first;
            std::vector<int> node_children = tree_node_children(tree.getEdges(), parent);

            // Debug ;;
            std::cout << "Pai :" << parent << "\nGrau: " << tree_node_degree_max.second << std::endl;
            vector_print(node_children, "filhos");

            matrix_print(cost);

            //exit(0);
            if(gen == 10)
                exit(0);

            int gen_next = gen + 1;
            // invalid solution - new searches with distinct restritions 
            for(int i = 0; i < node_children.size(); i++){
                ii edge_illegal;
                edge_illegal.first = parent;
                edge_illegal.second = node_children[i];
                branch_and_bound_depth(cost, edge_illegal, gen_next);
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

    Kruskal tree(cost);


    std::cout << "Cost: " << tree.MST(dimension) << std::endl;
    vii arvore = tree.getEdges();

    for(int i =0; i < arvore.size(); i++){
        std::cout << arvore[i].first << " " << arvore[i].second << std::endl;
    }
    std::cout << std::endl;

    ii restriction;
    branch_and_bound_depth(cost, restriction, 0);

    std::cout << "COST: " << s_cost_optimal << std::endl;

/*
    for (int i = 0; i < data->getDimension(); i++) {
        for(int j = 0; j < dimension; j++)
            std::cout << cost[i][j] << " ";
        std::cout << std::endl;
    }
    */

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
