#include <iostream>
#include <chrono>
#include <vector>
#include <list>
#include <queue>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include "data.h"
#include "hungarian.h"

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
    double s_cost;
    bool fertility;
};

double cost_optimal;
double s_cost_optimal = INFINITE;
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

void matrix_print(double ** ar, int dimension){
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << ar[i][j] << " ";

        }
        std::cout << std::endl;
    }
}

void cost_matrix_alloc(double *** matrix, int dimension){
    double ** ptr =  (double**)calloc(dimension, sizeof(double*));
    hungarian_test_alloc(ptr);
    for(int i = 0; i < dimension; i++){
        ptr[i] =  (double*)calloc(dimension, sizeof(double));
        hungarian_test_alloc(ptr[i]);
    }

    *matrix = ptr;
}

inline void cost_matrix_free(double *** ptr, int dimension){
    double ** ar = *ptr;
    for(int i = 0; i < dimension; i++)
        free(ar[i]);
    free(ar);

    *ptr = NULL;
}

void cost_restriction(std::vector<std::pair<int, int>> &edges, double ** cost_old, double ** cost_new, int dimension){
    // cpy the old to the new matrix of costs
    for(int i = 0; i < dimension; i++)
        memcpy(cost_new[i], cost_old[i], dimension*sizeof(double));  

    // apply the restrictions
    for(int i = 0; i < edges.size(); i++){
        int a = edges[i].first - 1;
        int b = edges[i].second - 1;

        cost_new[a][b] = INFINITE;
    }
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

// =========  Depth Search Function BEGIN ========= 

void branch_and_bound_depth(Data *data, double ** cost, struct node_info node_root, int gen){

    hungarian_problem_t p;
    const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    const int dimension = data->getDimension();

    double ** cost_new;
    cost_matrix_alloc(&cost_new, dimension);
    cost_restriction(node_root.edges_illegal, cost, cost_new, dimension);

    hungarian_init(&p, cost_new, dimension, dimension, mode); 
    cost_matrix_free(&cost_new, dimension);

    double obj_value = hungarian_solve(&p);

    if(obj_value <= s_cost_optimal){

        std::vector<int> subtour;
        subtour_lower_get(subtour, &p, dimension);
        hungarian_free(&p);

        if(subtour.size() >= dimension){ 
            // valid solution
            s_cost_optimal = obj_value;
            std::cout << "BestLB  " << s_cost_optimal << std::endl;
            s = subtour;
        }else {
            int gen_next = gen + 1;
            // invalid solution - new searches with distinct restritions 
            for(int i = 0; i < subtour.size() - 1; i++){
                struct node_info node_son;

                node_son.edges_illegal = node_root.edges_illegal;

                std::pair<int, int> edge;
                edge.first = subtour[i];
                edge.second = subtour[i+1];

                node_son.edges_illegal.push_back(edge);
                branch_and_bound_depth(data, cost, node_son, gen_next);

            }
        }
    }else
        hungarian_free(&p);

}

// =========  Depth Search Function END ========= 


// =========  Breadth Search Functions BEGIN ========= 

inline void node_solve(struct node_info2 & node, double **cost, int dimension){
    
    const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_problem_t p;

    double ** cost_new;
    cost_matrix_alloc(&cost_new, dimension);
    cost_restriction(node.edges_illegal, cost, cost_new, dimension);

    hungarian_init(&p, cost_new, dimension, dimension, mode); 
    cost_matrix_free(&cost_new, dimension);

    double node_cost = hungarian_solve(&p);


    node.s_cost = node_cost;


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

void branch_and_bound_breadth(Data * data, double ** cost){

    const int dimension = data->getDimension();
    struct node_info2 node_root;
    std::list<struct node_info2> layer_A;
    std::list<struct node_info2> layer_B;

    layer_A.push_back(node_root);

    while(!layer_A.empty() || !layer_B.empty()){
        std::cout << "Left  " << layer_A.size() << std::endl;

        while(!layer_A.empty()){

            node_solve(layer_A.front(), cost, dimension);

            if(layer_A.front().fertility == true)
                // generate the sons in the layer B
                node_procreate(layer_B, layer_A.front());
            
            // kill/erase the node
            layer_A.pop_front();
            
        }
        std::cout << "Left  " << layer_B.size() << std::endl;

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

/*
struct node_info2 {
    std::vector<std::pair<int, int>> edges_illegal;
    std::vector<int> subtour;
    double s_cost;
    bool fertility;
};

*/

// =========  Best-First Search Functions BEGIN ========= 

bool operator<(const node_info2 & a, const node_info2 & b){
    return a.s_cost > b.s_cost - DBL_EPSILON;
}

bool still_fertil(struct node_info2 & node){
    return node.s_cost < (s_cost_optimal - DBL_EPSILON);
}

void branch_and_bound_best(Data * data, double ** cost){

    const int dimension = data->getDimension();
    struct node_info2 node_root;
    std::priority_queue<struct node_info2> tree;

    node_solve(node_root, cost, dimension);
    tree.push(node_root);

    while(!tree.empty()){
        struct node_info2 senior = tree.top();
        tree.pop();

        //std::cout << "fertilidade  " << senior.fertility <<  std::endl << "Custo  " << senior.s_cost << std::endl;
        if(!still_fertil(senior)){
            continue;
        }

        std::cout << "Left  " << tree.size() << std::endl;

        std::list<struct node_info2> progeny;
        node_procreate(progeny, senior);

        while(!progeny.empty()){
            struct node_info2 junior = progeny.front();
            progeny.pop_front();

            node_solve(junior, cost, dimension);

            if(junior.fertility == true)
                tree.push(junior);
        }

    }


}

// =========  Best-First Search Functions END ========= 

int main(int argc, char** argv) {

    const char * breadth_flag = "--breadth";
    const char * depth_flag = "--depth";
    const char * best_flag = "--best";
    const char * bound_opt = "-o";

    Data * data = new Data(argc, argv[1]);
    data->readData();

    double **cost = new double*[data->getDimension()];
    for (int i = 0; i < data->getDimension(); i++){
        cost[i] = new double[data->getDimension()];
        for (int j = 0; j < data->getDimension(); j++){
            cost[i][j] = data->getDistance(i, j);
        }
    }

    cost_optimal = cost_optimal_get(argv[1]);

    if(argc == 4 && !strcmp(bound_opt, argv[3]))
        s_cost_optimal = cost_optimal;

    struct node_info node_initial;

    auto t1 = std::chrono::high_resolution_clock::now();

    if(!strcmp(breadth_flag, argv[2]))
        branch_and_bound_breadth(data, cost);
    else if(!strcmp(depth_flag, argv[2]))
        branch_and_bound_depth(data, cost, node_initial, 0);
    else if(!strcmp(best_flag, argv[2]))
        branch_and_bound_best(data, cost);
    else{
        std::cout << "[!!] Error: Search Mode not identified\n" << std::endl;
        exit(-1);
    }


    auto t2 = std::chrono::high_resolution_clock::now();


    /*
    for(int i = 0; i < s.size(); i++){
        std::cout << s[i] << " ";
    }
    std::cout << std::endl;
    */

    if(!s.empty())
        std::cout << "COST: " << s_cost_optimal << std::endl;
    else
        std::cout << "Couldnt find a solution." << std::endl;

    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "TIME: " << exec_time / 10e2  << std::endl;

    for (int i = 0; i < data->getDimension(); i++) 
        delete [] cost[i];
    delete [] cost;
    delete data;

    return 0;
}
