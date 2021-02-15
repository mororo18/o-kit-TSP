#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cfloat>
#include "data.h"
#include "hungarian.h"

//using namespace std;

struct node_info {
    std::vector<std::pair<int, int>> edges_illegal;
    std::vector<std::vector<int>> subtour;
    double bound_lower;
    int subtour_selected;
    bool node_fertility;
};

bool subtour_cmp(std::vector<int> &a, std::vector<int> &b){
    return a.size() < b.size();
}

void matrix_alloc(double *** matrix, int dimension){
	//double ** ptr = new double*[dimension];
	double ** ptr =  (double**)calloc(dimension, sizeof(double*));
    hungarian_test_alloc(ptr);
    for(int i = 0; i < dimension; i++){
        //ptr[i] = new double[dimension];
        ptr[i] =  (double*)calloc(dimension, sizeof(double));
        hungarian_test_alloc(ptr[i]);
    }

    *matrix = ptr;
}

inline void matrix_free(double *** ptr, int dimension){
    double ** ar = *ptr;
    for(int i = 0; i < dimension; i++)
        //printf("%p ", ptr[i]);
        free(ar[i]);
        //delete [] ar[i];
	//delete [] ar;
    free(ar);
    //printf("\n%p \n", *ptr);

    *ptr = NULL;

}

void cost_restriction(std::vector<std::pair<int, int>> &edges, double ** cost_old, double ** cost_new, int dimension){
    for(int i = 0; i < dimension; i++)
        memcpy(cost_new[i], cost_old[i], dimension*sizeof(double));  

    //can be reduced for just load the last illegal edge
    for(int i = 0; i < edges.size(); i++){
        int a = edges[i].first - 1;
        int b = edges[i].second - 1;

        //std::cout << "EDGE SZ :" << edges.size() << std::endl;
        //std::cout << b << std::endl;
        cost_new[a][b] = INFINITE;// DBL_MAX;
        cost_new[b][a] = INFINITE;//DBL_MAX;
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
    //tour[count].push_back(index+1);
    std::cout << "yes"  << std::endl;

    while(true){

        // identify the edges
        for(int i = 0; i < dimension; i++){
            //std::cout << solution->assignment[index][i];// << std::endl;
            if(solution->assignment[index][i]){
                nodes_known[index] = true;
                index = i;
                //std::cout << i << std::endl;
                break;
            }
        }
        //std::cout << std::endl;

        tour[count].push_back(index+1);

        // identify the end of the subtour 
        if(tour_node_first == index+1){
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

            // when there is no node avaliable for a new tour
            if(!flag)
                break;
        }
    }

    std::sort(tour.begin(), tour.end(), subtour_cmp);

    /*
    for(int i = 0; i< tour.size(); i++){
        for(int j = 0; j < tour[i].size(); j++)
            std::cout << tour[i][j] << " ";

        std::cout << std::endl;
    }
    */

    subtour = tour[0];
}

void matrix_print(double ** ar, int dimension){
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << ar[i][j] << " ";

        }
        std::cout << std::endl;
    }
}

void branch_AND_bound(Data *data, double ** cost, struct node_info node_root){

	hungarian_problem_t p;
	const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    const int dimension = data->getDimension();

    double ** cost_new;
    matrix_alloc(&cost_new, dimension);
    cost_restriction(node_root.edges_illegal, cost, cost_new, dimension);
    //std::cout << "old :" << std::endl;
    //matrix_print(cost, dimension);
    //std::cout << "new :" << std::endl;
    matrix_print(cost_new, dimension);

    std::cout << "aqui1" << std::endl;
	hungarian_init(&p, cost_new, dimension, dimension, mode); // Carregando o problema
    matrix_free(&cost_new, dimension);
    std::cout << "aqui2" << std::endl;

	double obj_value = hungarian_solve(&p);
    std::cout << "aqui3" << std::endl;
	std::cout << "Obj. value: " << obj_value << std::endl;
	hungarian_print_assignment(&p);
    
    std::vector<int> subtour;
    subtour_lower_get(subtour, &p, dimension);
	hungarian_free(&p);

    //if(obj_value > 3300) return;
    if(subtour.size() == dimension)
        return;
    /*
    for(int j = 0; j < subtour.size(); j++)
        std::cout << subtour[j] << " ";
    
    std::cout << std::endl;
    */

    for(int i = 0; i < subtour.size() - 2; i++){
        struct node_info node_son;

        node_son.edges_illegal = node_root.edges_illegal;

        std::pair<int, int> edge;
        edge.first = subtour[i];
        edge.second = subtour[i+1];

        node_son.edges_illegal.push_back(edge);
        branch_AND_bound(data, cost, node_son);
        
    }

}

int main(int argc, char** argv) {

	Data * data = new Data(argc, argv[1]);
	data->readData();

	double **cost = new double*[data->getDimension()];
	for (int i = 0; i < data->getDimension(); i++){
		cost[i] = new double[data->getDimension()];
		for (int j = 0; j < data->getDimension(); j++){
			cost[i][j] = data->getDistance(i, j);
		}
	}

    /*
	hungarian_problem_t p;
	const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
	hungarian_init(&p, cost, data->getDimension(), data->getDimension(), mode); // Carregando o problema

	double obj_value = hungarian_solve(&p);
	std::cout << "Obj. value: " << obj_value << std::endl;

	std::cout << "Assignment" << std::endl;
	hungarian_print_assignment(&p);

	hungarian_free(&p);
    */

    struct node_info node_initial;
    branch_AND_bound(data, cost, node_initial);
	for (int i = 0; i < data->getDimension(); i++) 
        delete [] cost[i];
	delete [] cost;
	delete data;

	return 0;
}
