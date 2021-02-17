#include <iostream>
#include <chrono>
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

int cost_optimal;

bool subtour_cmp(std::vector<int> &a, std::vector<int> &b){
    return (a.size() < b.size());// && a[0] < b[0];
}

bool subtour_cmp2(std::vector<int> &a, std::vector<int> &b){
    return a[0] < b[0];
}

void matrix_alloc(int *** matrix, int dimension){
	//double ** ptr = new double*[dimension];
	int ** ptr =  (int**)calloc(dimension, sizeof(int*));
    hungarian_test_alloc(ptr);
    for(int i = 0; i < dimension; i++){
        //ptr[i] = new double[dimension];
        ptr[i] =  (int*)calloc(dimension, sizeof(int));
        hungarian_test_alloc(ptr[i]);
    }

    *matrix = ptr;
}

inline void matrix_free(int *** ptr, int dimension){
    int ** ar = *ptr;
    for(int i = 0; i < dimension; i++)
        //printf("%p ", ptr[i]);
        free(ar[i]);
        //delete [] ar[i];
	//delete [] ar;
    free(ar);
    //printf("\n%p \n", *ptr);

    *ptr = NULL;

}

void cost_restriction(std::vector<std::pair<int, int>> &edges, int ** cost_old, int ** cost_new, int dimension){
    for(int i = 0; i < dimension; i++)
        memcpy(cost_new[i], cost_old[i], dimension*sizeof(int));  

    //can be reduced for just load the last illegal edge
    for(int i = 0; i < edges.size(); i++){
        int a = edges[i].first - 1;
        int b = edges[i].second - 1;

        //std::cout << "EDGE SZ :" << edges.size() << std::endl;
        //std::cout << b << std::endl;
        cost_new[a][b] = INFINITE;// DBL_MAX;
        //cost_new[b][a] = INFINITE;//DBL_MAX;
    }
}

void subtour_lower_get(std::vector<std::vector<int>> &subtour, hungarian_problem_t * solution, int dimension){
    std::vector<std::vector<int>> tour(1, std::vector<int> (1, 1));
    bool nodes_known[dimension] = {0};
    bool flag;
    int count = 0;
    int index = 0;
    int tour_node_first = 1;
    nodes_known[index] = true;
    //tour[count].push_back(index+1);
    //std::cout << "yes"  << std::endl;

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

            // when there is no node avaliable for a new tour
            if(!flag)
                break;
        }
    }
    /*
    for(int i = 0; i < tour.size(); i++){
        for(int j = 0; j < tour[i].size(); j++)
            std::cout << tour[i][j] << " ";

        std::cout << std::endl;
    }

    int sz = tour.size();
    for(int i = sz-1; i >= 0 ; i--)
        if(tour[i].size() <= 3){
            tour.erase(tour.begin() + i);
        }

    if(tour.empty()){
        //std::vector<int> vec_empty();
        //subtour = vec_empty;
        subtour.clear();
        return;
    }

    */

    std::sort(tour.begin()/*, tour.begin()+1*/,tour.end(), subtour_cmp); 
    //std::partial_sort(tour.begin(), tour.begin(),tour.begin()+1, subtour_cmp2);

    /*
    for(int i = 0; i < tour[0].size(); i++){
        //for(int j = 0; j < tour[i].size(); j++)
            std::cout << tour[0][i] << " ";

        //std::cout << std::endl;
    }
    std::cout << std::endl;

    for(int i = 0; i < tour.size(); i++){
        for(int j = 0; j < tour[i].size(); j++)
            std::cout << tour[i][j] << " ";

        std::cout << std::endl;
    }
    */
    subtour = tour;
}

void matrix_print(int ** ar, int dimension){
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << ar[i][j] << " ";

        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<int>> s;
std::vector<int> bigest;

inline void branch_AND_bound(Data *data, int ** cost, struct node_info node_root, int gen){

    //int count;
    //std::cout << std::endl << gen << "th generation" << std::endl;
    //count = gen+1;

	hungarian_problem_t p;
	const int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    const int dimension = data->getDimension();

    int ** cost_new;
    matrix_alloc(&cost_new, dimension);
    cost_restriction(node_root.edges_illegal, cost, cost_new, dimension);
    //std::cout << "old :" << std::endl;
    //matrix_print(cost, dimension);
    //std::cout << "new :" << std::endl;
    //matrix_print(cost_new, dimension);

    //std::cout << "aqui1" << std::endl;
	hungarian_init(&p, cost_new, dimension, dimension, mode); // Carregando o problema
    matrix_free(&cost_new, dimension);
    //std::cout << "aqui2" << std::endl;

	int obj_value = hungarian_solve(&p);
    //std::cout << "aqui3" << std::endl;
	//std::cout << "Obj. value: " << obj_value << std::endl;
	//hungarian_print_assignment(&p);
    
    std::vector<std::vector<int>> subtour;
    subtour_lower_get(subtour, &p, dimension);
	hungarian_free(&p);
    //for(int j = 0; j < subtour.size(); j++){
    //for(int j = 0; j < 1; j++){
        //exit(0);
        //std::cout << " aqui4 " << subtour[0] << subtour[1]<< std::endl;


        //if(obj_value > cost_optimal) return;

        //int sol_sz = s.size();
        int j = 0;
        
        //std::cout << subtour.empty() << std::endl;
        if(subtour[j].size() >= dimension && obj_value <= cost_optimal){ 
            /*
               for(int i = 0; i< dimension+1; i++)
               std::cout << subtour[i] << " ";
               std::cout << std::endl;
             */
            //std::cout << obj_value << std::endl;
            subtour[j].push_back(obj_value);
            s.push_back(subtour[j]);
            //exit(0);
            //return;
        }else if(obj_value <= cost_optimal && !subtour[j].empty()){

            for(int i = 0; i < subtour[j].size() - 1; i++){
                struct node_info node_son;

                node_son.edges_illegal = node_root.edges_illegal;

                std::pair<int, int> edge;
                edge.first = subtour[j][i];
                edge.second = subtour[j][i+1];

                node_son.edges_illegal.push_back(edge);
                branch_AND_bound(data, cost, node_son, gen+1);


            }
        }

        //if(sol_sz != s.size())
            //break;
    //}

}

int getOptimalCost(char * instanceName){
     string path = "benchmark/target_data";
     int value;
     string file;
     ifstream inTSP(path, ios::in);
 
     //string file;   //They are used into the reader
 
     inTSP >> file;
     while ( file.find((string)instanceName+":") != 0) {
         inTSP >> file;
         //std::cout << file<< std::endl;
     }
     int pos = file.find(":") + 1;
     string value_str = file.substr(pos);
     //if ( file.compare((string)instanceName ) == 0 )  inTSP >> file;
     //std::cout << pos<< std::endl;
 
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

    cost_optimal = getOptimalCost(argv[1]);
    //exit(0);
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

    auto t1 = std::chrono::high_resolution_clock::now();
    branch_AND_bound(data, cost, node_initial, 0);
    auto t2 = std::chrono::high_resolution_clock::now();

    for(int i = 0; i< s.size(); i++){
        for(int j = 0; j< s[i].size()-1; j++)
            std::cout << s[i][j] << " ";

        std::cout << "Cost :" << s[i][s[i].size()-1] << std::endl;
    }

    auto exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

	std::cout << "Execution time: " << exec_time / 10e2  << std::endl;
	for (int i = 0; i < data->getDimension(); i++) 
        delete [] cost[i];
	delete [] cost;
	delete data;

	return 0;
}
