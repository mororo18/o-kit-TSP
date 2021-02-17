#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include "data.h"
#include "hungarian.h"

struct node_info {
    std::vector<std::pair<int, int>> edges_illegal;
    std::vector<std::vector<int>> subtour;
    double bound_lower;
    int subtour_selected;
    bool node_fertility;
};

int cost_optimal;
std::vector<int> s;

void cost_matrix_print(int ** ar, int dimension){
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << ar[i][j] << " ";

        }
        std::cout << std::endl;
    }
}

bool subtour_cmp(std::vector<int> &a, std::vector<int> &b){
    return (a.size() < b.size());// && a[0] < b[0];
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

void cost_restriction(std::vector<std::pair<int, int>> &edges, int ** cost_old, int ** cost_new, int dimension){
    // cpy the old to the new matrix of costs
    for(int i = 0; i < dimension; i++)
        memcpy(cost_new[i], cost_old[i], dimension*sizeof(int));  

    // apply the restrictions
    for(int i = 0; i < edges.size(); i++){
        int a = edges[i].first - 1;
        int b = edges[i].second - 1;

        cost_new[a][b] = INFINITE;
    }
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

    std::partial_sort(tour.begin(), tour.begin()+1,tour.end(), subtour_cmp); 

    subtour = tour[0];
}

void branch_AND_bound(Data *data, int ** cost, struct node_info node_root, int gen){

    hungarian_problem_t p;
    const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    const int dimension = data->getDimension();

    int ** cost_new;
    cost_matrix_alloc(&cost_new, dimension);
    cost_restriction(node_root.edges_illegal, cost, cost_new, dimension);

    hungarian_init(&p, cost_new, dimension, dimension, mode); 
    cost_matrix_free(&cost_new, dimension);

    int obj_value = hungarian_solve(&p);

    if(obj_value <= cost_optimal){

        std::vector<int> subtour;
        subtour_lower_get(subtour, &p, dimension);
        hungarian_free(&p);

        if(subtour.size() >= dimension){ 
            // valid solution
            cost_optimal = obj_value;
            s = subtour;
        }else {
            // invalid solution - new searches with distinct restritions 
            for(int i = 0; i < subtour.size() - 1; i++){
                struct node_info node_son;

                node_son.edges_illegal = node_root.edges_illegal;

                std::pair<int, int> edge;
                edge.first = subtour[i];
                edge.second = subtour[i+1];

                node_son.edges_illegal.push_back(edge);
                branch_AND_bound(data, cost, node_son, gen+1);

            }
        }
    }else
        hungarian_free(&p);

}

int cost_optimal_get(char * instance_name){
    string path = "benchmark/target_data";
    int value;
    string file;
    ifstream inTSP(path, ios::in);

    inTSP >> file;
    while ( file.find((string)instance_name+":") != 0) {
        inTSP >> file;
    }

    int pos = file.find(":") + 1;
    string value_str = file.substr(pos);

    value = stoi(value_str);

    return value;
}

int main(int argc, char** argv) {

    Data * data = new Data(argc, argv[1]);
    data->readData();

    int **cost = new int*[data->getDimension()];
    for (int i = 0; i < data->getDimension(); i++){
        cost[i] = new int[data->getDimension()];
        for (int j = 0; j < data->getDimension(); j++){
            cost[i][j] = data->getDistance(i, j);
        }
    }

    cost_optimal = cost_optimal_get(argv[1]);

    struct node_info node_initial;

    auto t1 = std::chrono::high_resolution_clock::now();
    branch_AND_bound(data, cost, node_initial, 0);
    auto t2 = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < s.size(); i++){
        std::cout << s[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Cost: " << cost_optimal << std::endl;
    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Execution time: " << exec_time / 10e2  << std::endl;

    for (int i = 0; i < data->getDimension(); i++) 
        delete [] cost[i];
    delete [] cost;
    delete data;

    return 0;
}
