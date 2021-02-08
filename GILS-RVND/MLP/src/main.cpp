#include <iostream>
#include <cstdint>
#include <cfloat>
#include <new>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "readData.h"
#include <cassert>

#define REINSERTION 1
#define OR_OPT2 	2
#define OR_OPT3 	3
#define SWAP 		4
#define TWO_OPT		5

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::system_clock;
using std::chrono::high_resolution_clock;

struct alignas(4) insertion_info {
    double cost;
    int edge_removed;
    int node_new;
};

struct alignas(alignof(int)) subseq {
    double T;
    double C;
    int W;
};

struct sub_info{
    std::vector<int> seq;
    short prev;
    short next;
    short pos_next;
};

double **c;
int dimension;

std::chrono::_V2::system_clock::time_point t3;
std::chrono::_V2::system_clock::time_point t4;


bool state;
bool flag;
int sum_t = 0;
long swap_t  = 0;
long two_opt_t  = 0;
long reinsertion_t  = 0;
long opt2_t = 0;
long opt3_t = 0;
long construct_t = 0;
long search_t; 
long search_t_average = 0;

#define TABLE_SZ 26

double R_table(int i){
    static const double table[] = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    //assert((unsigned) i >= (unsigned) TABLE_SZ);

    return table[i];
}

void after(){
    if(!flag)
        return;
    t3 = high_resolution_clock::now();
}

void before(int a){
    if(!flag)
        return;
    t4 = high_resolution_clock::now();

    switch(a){
        case SWAP:
            swap_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;
        case TWO_OPT:
            two_opt_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;				
        case REINSERTION:
            reinsertion_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;				
        case OR_OPT2:
            opt2_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;				
        case OR_OPT3:
            opt3_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;				
        case 6:
            construct_t += std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            break;				
        case 7:
            search_t = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
            search_t_average += search_t;
            break;				

    }
}

bool cost_compare(const insertion_info &a, const insertion_info &b){
    return a.cost < b.cost - DBL_EPSILON;
}

void candidates_load(std::vector<int> &cand, int dim){
    for(int i = 2; i <= dim; ++i){
        cand.push_back(i);
    }
}	

inline void binary_search(std::vector<int> &vec, int node, short from, short to){
    //short middle = (from + to) / 2; 
    #define MIDDLE (from + to) / 2

    if(vec[MIDDLE] == node){
        vec.erase(vec.begin() + MIDDLE);
        return;
    }else if(vec[MIDDLE] < node){
        binary_search(vec, node, MIDDLE + 1, to);
    }else {
        binary_search(vec, node, from, MIDDLE);
    }
}

void construct(std::vector<int> &s, const double alpha){

    //s.clear();
    s = {1, 1};
    std::vector<int> candidates;
    candidates.reserve(dimension);
    candidates_load(candidates, dimension);

    int subtour_inicial = 3;

    /*
    for(int i = subtour_inicial; --i;){
        
        int node_rand_index = (unsigned)rand() % candidates.size(); 
        int node_rand = candidates[node_rand_index];

        s.insert(s.begin() + 1, node_rand);
        candidates.erase(candidates.begin() + node_rand_index);
    }
    */

    int r = 1;

    while(!candidates.empty()){

        //int dim = s.size() - 1;

        //int insertion_possiblts = (dim * candidates.size());
        std::vector<struct insertion_info> insertion_cost (candidates.size());

        for(int i = 0; i < candidates.size(); ++i){   
            insertion_cost[i].cost = c[r][candidates[i]];  
            insertion_cost[i].node_new = candidates[i];
        }
        //printf("valor de l  =  %d\n", l);


        int node_rand_range = alpha * insertion_cost.size() + 1;
        int node_rand_index = (unsigned)rand() % node_rand_range;
        std::partial_sort(insertion_cost.begin(), insertion_cost.begin() + node_rand_range, insertion_cost.end(), cost_compare);

        int node_rand = insertion_cost[node_rand_index].node_new;
        s.insert(s.end() - 1, node_rand); 
        r = node_rand;
        //printf("inserted node: %d\n", node_rand);

        /*
        for(int i = 0; i < candidates.size(); ++i){
            if(candidates[i] ==  node_rand){
                candidates.erase(candidates.begin() + i);
                break;
            }
        }
        */

        binary_search(candidates, node_rand, 0, candidates.size() );

    }

    /*  
    for(int i = 0; i< s.size(); i++)
        std::cout << s[i] << " ";

    std::cout << std::endl;
    */

    //exit(-1);

}	

inline void swap_2(std::vector<int> &vec, int i, int j){
    std::iter_swap(vec.begin() + i, vec.begin() + j);
}

inline void two_opt(std::vector<int> &vec, int i, int j){
    std::reverse(vec.begin() + i, vec.begin() + j+1);
    //std::cout << "yas\n";
}

inline void reinsert(std::vector<int> &vec, int i, int j, int pos){
    if(pos < i){
        std::vector<int> copy;
        copy.reserve(j - i);
        std::vector<int>::iterator t1 = vec.begin() + i;
        std::vector<int>::iterator t2 = vec.begin() + j; 
        copy.insert(copy.begin(),t1 ,t2 );
        vec.erase(t1, t2);
        vec.insert(vec.begin() + pos, copy.begin(), copy.end());
    }else{
        std::vector<int> copy;
        copy.reserve(j - i );
        copy.insert(copy.begin(), vec.begin() + i, vec.begin() + j);
        vec.insert(vec.begin() + pos, copy.begin(), copy.end());
        vec.erase(vec.begin() + i, vec.begin() + j);
    }

}

inline void load_subseq_info(std::vector<std::vector<struct subseq>> &seq, std::vector<int> &s){
    
    //std::cout << "cost_1: " << std::endl;
    /*
    for(int i = 0; i < dimension+1; i++){
        for(int j = i; j < dimension+1; j++){
            seq[i][j].T = time_calc(s, i, j);
            seq[i][j].C = cost_acumulated_calc(s, i, j);
            seq[i][j].W = j - i + 1 -(!i) - (j == dimension-1);
            //std::cout << seq[s[i]][s[j]].W << " " ;

            seq[j][i].C = seq[i][j].C;
            seq[j][i].T = seq[i][j].T;
            seq[j][i].W = seq[i][j].W;

        }
        //std::cout << std::endl;
    }
    //std::cout << seq[0][dimension].C << "\n" ;
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << seq[s[i]][s[j]].C << " " ;

        }
        std::cout << std::endl;
    }

    //std::cout << "cost_2: " << std::endl;
    for(int i = 0; i < dimension+1; i++){
        for(int j = i; j < dimension+1; j++){
            seq[i][j].T = time_calc(s, i, j);
            seq[i][j].C = time_calc(s, i, j) + seq[i][j-1].C*(j!=i);
            //std::cout << seq[s[i]][s[j]].C << " " ;
            seq[i][j].W = j - i + 1 -(!i) - (j == dimension-1);

        }
        //std::cout << std::endl;
    }

*/
    //std::cout << "cost_3: " << std::endl;
    for(int i = 0; i < dimension+1; i++){
        for(int j = i; j < dimension+1; j++){
            seq[i][j].T = (c[s[j-1]][s[j]] + seq[i][j-1].T)*(j!=i);
            seq[i][j].C = seq[i][j].T + seq[i][j-1].C*(j!=i);
            //std::cout << seq[s[i]][s[j]].C << " " ;
            seq[i][j].W = j - i + 1 -(!i) ;//- (j == dimension);

            seq[j][i].C = seq[i][j].C;
            seq[j][i].T = seq[i][j].T;
            seq[j][i].W = seq[i][j].W;
        }
        //std::cout << std::endl;
    }

    /*
    for(int i = 0; i < dimension+1; i++){
        for(int j = 0; j < dimension+1; j++){
            std::cout << seq[i][j].C << " " ;

        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */
}
inline void load_subseq_info2(std::vector<std::vector<struct subseq>> &seq, std::vector<int> &s, int i = 0){
    
    //std::cout << "cost_1: " << std::endl;
    /*
    for(int i = 0; i < dimension+1; i++){
        for(int j = i; j < dimension+1; j++){
            seq[i][j].T = time_calc(s, i, j);
            seq[i][j].C = cost_acumulated_calc(s, i, j);
            seq[i][j].W = j - i + 1 -(!i) - (j == dimension-1);
            //std::cout << seq[s[i]][s[j]].W << " " ;

            seq[j][i].C = seq[i][j].C;
            seq[j][i].T = seq[i][j].T;
            seq[j][i].W = seq[i][j].W;

        }
        //std::cout << std::endl;
    }
    //std::cout << seq[0][dimension].C << "\n" ;
    for(int i = 0; i < dimension; i++){
        for(int j = 0; j < dimension; j++){
            std::cout << seq[s[i]][s[j]].C << " " ;

        }
        std::cout << std::endl;
    }

    //std::cout << "cost_2: " << std::endl;
    for(int i = 0; i < dimension+1; i++){
        for(int j = i; j < dimension+1; j++){
            seq[i][j].T = time_calc(s, i, j);
            seq[i][j].C = time_calc(s, i, j) + seq[i][j-1].C*(j!=i);
            //std::cout << seq[s[i]][s[j]].C << " " ;
            seq[i][j].W = j - i + 1 -(!i) - (j == dimension-1);

        }
        //std::cout << std::endl;
    }

*/
    //std::cout << "cost_3: " << std::endl;
    int from = i;
    for(int i = 0; i < dimension+1; i++){
        for(int j = (i>=from?from++:from); j < dimension+1; j++){
            seq[i][j].T = (c[s[j-1]][s[j]] + seq[i][j-1].T)*(j!=i);
            seq[i][j].C = seq[i][j].T + seq[i][j-1].C*(j!=i);
            //std::cout << seq[s[i]][s[j]].C << " " ;
            seq[i][j].W = j - i + 1 -(!i) ;//- (j == dimension);

            seq[j][i].C = seq[i][j].C;
            seq[j][i].T = seq[i][j].T;
            seq[j][i].W = seq[i][j].W;
        }
        //std::cout << std::endl;
    }

    /*
    for(int i = 0; i < dimension+1; i++){
        for(int j = 0; j < dimension+1; j++){
            std::cout << seq[i][j].C << " " ;

        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */
}


#define DBL_SZ 8
#define INT_SZ 4

inline void cost_calc(std::vector<int> &s, double *cost_best){
    *cost_best = 0;
    double sum = 0;
    for(unsigned i = 0; i < (unsigned)dimension; ++i)
        sum += c[s[i]][s[i+1]];

    *cost_best = sum;
}

inline double time_calc(std::vector<int> &s, int from, int to){
    double sum = 0;
    for(int i = from; i < to; ++i)
        sum += c[s[i]][s[i+1]];

    return sum;
}

inline double cost_acumulated_calc(std::vector<int> &s, int from, int to){
    double sum = 0;
    for(int i = to; i >= from; --i)
        for(int j = from; j < i; j++)
            sum += c[s[j]][s[j+1]];
        

    return sum;
}

inline double cost_reverse_calc(std::vector<int> &s, int from, int to){
    double sum = 0;
    for(int i = from; i < to; ++i)
        for(int j = i; j < to; j++)
            sum += c[s[j]][s[j+1]]; 
    return sum;
}

inline void neighbor_swap_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){
    alignas(DBL_SZ) double dif;
    alignas(DBL_SZ) double dif1;
    alignas(DBL_SZ) double dif3;
    alignas(DBL_SZ) double dif_lower;
    alignas(INT_SZ) int i;
    alignas(INT_SZ) int j;
    alignas(INT_SZ) int e;
    alignas(INT_SZ) int b;
    alignas(INT_SZ) int d;
    alignas(INT_SZ) int a;
    alignas(INT_SZ) int o;
    alignas(INT_SZ) int p;
    alignas(INT_SZ) int q;
    alignas(INT_SZ) int dim = dimension - 2;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool t = 1;
    alignas(1) bool l;

    for(i = 1; i < dimension - 1; ++i){
        //o = s[i];
        d = i - 1;
        a = i + 1;
        q = a + 1;
        l = false;

        //dif1 = - c[s[i]][s[d]];
        //dif3 = dif1 - c[s[i]][s[a]];

        //consecutive nodes
        //dif = (dif1  + c[s[i]][s[q]]) + (c[s[a]][s[d]] - c[s[a]][s[q]]);
        dif = seq[0][i-1].C + seq[a][i].W * (seq[0][i-1].T + c[s[i-1]][s[a]]) +  c[s[a]][s[i]];
        dif = dif + seq[a+1][dimension].W * (seq[0][i-1].T + c[s[i-1]][s[a]] + seq[a][i].T  + c[s[i]][s[a+1]] ) +  seq[a+1][dimension].C;

        if(dif < dif_lower || t){
            dif_lower = dif - DBL_EPSILON;
            i_best = i;
            j_best = a;
            t = 0;
        }

        if(i == dim) continue;

        for(j = q; j < dimension; ++j){
            p = s[j];
            b = j + 1;
            e = j - 1;

            dif = seq[0][i-1].C + seq[j][j].W * (seq[0][i-1].T + c[s[i-1]][s[j]]) +  c[s[j]][s[j]];
            dif = dif + seq[i+1][j-1].W * (seq[0][i-1].T + c[s[i-1]][s[j]] + seq[j][j].T  + c[s[j]][s[i+1]] ) +  seq[i+1][j-1].C;
            dif = dif + seq[i][i].W * (seq[0][i-1].T + c[s[i-1]][s[j]] + seq[j][j].T  + c[s[j]][s[i+1]] + seq[i+1][j-1].T+ c[s[j-1]][s[i]]) + seq[i][i].C;
            dif = dif + seq[j+1][dimension].W * (seq[0][i-1].T + c[s[i-1]][s[j]] + seq[j][j].T  + c[s[j]][s[i+1]] + seq[i+1][j-1].T+ c[s[j-1]][s[i]] + seq[i][i].T + c[s[i]][s[j+1]]) + seq[j+1][dimension].C;

            if(dif < dif_lower){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
            }
        }
    }

    if(dif_lower < seq[0][dimension].C - DBL_EPSILON){
        swap_2(s, i_best, j_best);
        load_subseq_info2(seq, s, i_best);
        /*
        std::cout << dif_lower << std::endl;
        std::cout << seq[0][dimension].C << std::endl;

        std::cout  << std::endl;
        */
        state = true;
    }
}

inline void neighbor_two_opt_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){
    alignas(DBL_SZ) double dif;
    alignas(DBL_SZ) double dif1;
    alignas(DBL_SZ) double dif_lower;
    alignas(INT_SZ) int b;
    alignas(INT_SZ) int a;
    alignas(INT_SZ) int j;
    alignas(INT_SZ) int i;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool t = 1;

    for(i = 1; i < dimension - 1; ++i){
        //b = s[i - 1];
        //dif1 = - c[s[i]][b];
        for(j = i + 2; j < dimension; ++j){
            //a = j + 1;
            a = cost_reverse_calc(s, i, j);

            //dif = (dif1 - c[s[j]][s[a]])  + (c[s[j]][b] + c[s[i]][s[a]]);
            dif = seq[0][i-1].C + seq[j][i].W * (seq[0][i-1].T + c[s[j]][s[i-1]]) +  a;
            dif = dif + seq[j+1][dimension].W * (seq[0][i-1].T + c[s[j]][s[i-1]] + seq[i][j].T  + c[s[j+1]][s[i]] ) +  seq[j+1][dimension].C;
            //std::cout << dif << std::endl;
            

            if(dif < dif_lower || t){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                t = 0;

            }
        }
    }

    if(dif_lower < seq[0][dimension].C - DBL_EPSILON){
        two_opt(s, i_best, j_best);
        load_subseq_info2(seq, s, i_best);
        /*
        std::cout << dif_lower << std::endl;
        std::cout << seq[0][dimension].C << std::endl;
        std::cout << i_best << std::endl;
        std::cout << j_best << std::endl;
        for(int i = 0; i< dimension+1; i++)
            std::cout << s[i] << " ";
        std::cout  << std::endl;
        */
        state = true;
    }
}


inline void neighbor_reinsertion_better(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq, int sz){
    alignas(DBL_SZ) double dif;
    alignas(DBL_SZ) double dif1;
    alignas(DBL_SZ) double dif_lower;
    alignas(INT_SZ) int i;
    alignas(INT_SZ) int j;
    alignas(INT_SZ) int k;
    alignas(INT_SZ) int g;
    alignas(INT_SZ) int b;
    alignas(INT_SZ) int a;
    alignas(INT_SZ) int o;
    alignas(INT_SZ) int e;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(INT_SZ) int pos_new;
    alignas(INT_SZ) int k_lim = dimension - sz - 1;
    alignas(1) bool l;
    alignas(1) bool t = 1;

    for(i = 1, j = sz + i - 1; i < dimension - sz + 1; ++i, ++j){
        //a = s[i];
        //o = s[j];
        //e = j + 1;
        b = i - 1;
        l = 1;

        //dif1 = c[s[b]][s[e]] - c[s[i]][s[b]] - c[s[j]][s[e]];

        //k -> edges 
        for(k = 0; k < b; ++k){


            /*
            if(l && k == b){
                k = k + sz + 1;
                l = 0;
            }
            */

            /*
            moving the 2nd and the 3rd elements to the 6th position,
            for example, is the same as moving the 4th and the 5th 
            elements to	the 2nd position.
             */

            //if(k == a)
                //continue;
            
            //g = k + 1;

            dif = seq[0][k].C + seq[i][j].W*(seq[0][k].T + c[s[k]][s[i]]) + seq[i][j].C; 
            dif = dif + seq[k+1][i-1].W*(seq[0][k].T + c[s[k]][s[i]] + seq[i][j].T + c[s[j]][s[k+1]]) + seq[k+1][i-1].C;
            dif = dif + seq[j+1][dimension].W*(seq[0][k].T + c[s[k]][s[i]] + seq[i][j].T + c[s[j]][s[k+1]] + seq[k+1][i-1].T + c[s[i-1]][s[j+1]]) + seq[j+1][dimension].C;

            if( dif < dif_lower || t){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
                t = 0;

            }
        }

        for(k = i + sz; k < k_lim; ++k){


            /*
            if(l && k == b){
                k = k + sz + 1;
                l = 0;
            }
            */

            /*
            moving the 2nd and the 3rd elements to the 6th position,
            for example, is the same as moving the 4th and the 5th 
            elements to	the 2nd position.
             */

            //if(k == a)
                //continue;
            
            //g = k + 1;

            //dif = (dif1 + c[a][s[k]]) + (c[o][s[g]] - c[s[g]][s[k]]); 

            dif = seq[0][i-1].C + seq[j+1][k].W*(seq[0][i-1].T + c[s[i-1]][s[j+1]]) + seq[j+1][k].C; 
            dif = dif + seq[i][j].W*(seq[0][i-1].T + c[s[i-1]][s[j+1]] + seq[j+1][k].T + c[s[k]][s[i]]) + seq[i][j].C;
            dif = dif + seq[k+1][dimension].W*(seq[0][i-1].T + c[s[i-1]][s[j+1]] + seq[j+1][k].T + c[s[k]][s[i]] + seq[i][j].T + c[s[j]][s[k+1]]) + seq[k+1][dimension].C;

            if( dif < dif_lower){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
            }
        }
    }

    if(dif_lower < seq[0][dimension].C - DBL_EPSILON){
        reinsert(s, i_best, j_best + 1, pos_new + 1);
        load_subseq_info(seq, s);
        state = true;
    }

}

inline void neighbd_list_repopulate(std::vector<int> &list){
    list.clear();
    list = {1,2,3,4,5};
    //printf("REPOPULATE\n");
}

void RVND(std::vector<int> &s, std::vector<std::vector<struct subseq>> &seq){

    //printf("for -> while -> RVND\n");
    alignas(8) std::vector<int> neighbd_list = {1,2,3,4,5};
    alignas(4) int neighbd_rand_index;
    alignas(4) int neighbd_rand;

    while(!neighbd_list.empty()){

        neighbd_rand_index = (unsigned)rand() % neighbd_list.size();
        neighbd_rand = neighbd_list[neighbd_rand_index];

        state = false;

        switch(neighbd_rand){
            case REINSERTION:
                //after();
                neighbor_reinsertion_better(s, seq, REINSERTION);
                //before(REINSERTION);
                break;				
            case OR_OPT2:
                //after();
                neighbor_reinsertion_better(s, seq, OR_OPT2);
                //before(OR_OPT2);
                break;				
            case OR_OPT3:
                //after();
                neighbor_reinsertion_better(s, seq, OR_OPT3);
                //before(OR_OPT3);
                break;				
            case SWAP:
                //after();
                neighbor_swap_better(s, seq);
                //before(SWAP);
                break;
            case TWO_OPT:
                //after();
                neighbor_two_opt_better(s, seq);
                //before(TWO_OPT);
                break;				
        }

        if(state){
            neighbd_list_repopulate(neighbd_list);
            //load_subseq_info(seq, s);
        }else
            neighbd_list.erase(neighbd_list.begin() + neighbd_rand_index);
        

    }
}

void perturb(std::vector<int> &sl, std::vector<int> &s){

    s.clear();
    s = sl;
    int dimension = s.size() - 1;
    int first_size = ((unsigned)rand() % (int)(dimension / 10)) + 2;
    int first_pos = (unsigned)rand() % (dimension - first_size - 1) + 1;

    struct sub_info *sub_seq1 = new(std::nothrow) struct sub_info;
    sub_seq1->prev = sl[first_pos - 1];
    sub_seq1->next = sl[first_pos + first_size];
    sub_seq1->seq.reserve(first_size);
    sub_seq1->seq.insert(sub_seq1->seq.begin(), sl.begin() + first_pos, sl.begin() + first_pos + first_size);

    int second_size = ((unsigned)rand() % ((int)(dimension / 10))) + 2;
    int second_pos;
    // second_pos restrictions
    while(true){
        int pos = (unsigned)rand() % (dimension - second_size - 1) + 1;
        if(pos <( first_pos - second_size) || pos > first_pos + first_size){
            second_pos = pos;
            break;
        }
    }
    struct sub_info *sub_seq2 = new(std::nothrow) struct sub_info;
    sub_seq2->prev = sl[second_pos - 1];
    sub_seq2->next = sl[second_pos + second_size];
    sub_seq2->seq.reserve(second_size);
    sub_seq2->seq.insert(sub_seq2->seq.begin(), sl.begin() + second_pos, sl.begin() + second_pos + second_size);

    if(first_pos > second_pos){
        // elimina sub_seq1
        s.erase(s.begin() + first_pos, s.begin() + first_pos + first_size);
        //elimina sub_seq2	
        s.erase(s.begin() + second_pos, s.begin() + second_pos + second_size);
    }else{
        //elimina sub_seq2	
        s.erase(s.begin() + second_pos, s.begin() + second_pos + second_size);
        // elimina sub_seq1
        s.erase(s.begin() + first_pos, s.begin() + first_pos + first_size);
    }

    dimension = s.size() - 1;
    while(true){
        int pos = (unsigned)rand() % dimension + 1;
        if(s[pos-1] == sub_seq1->prev || s[pos] == sub_seq1->next){
            continue;
        }else{
            sub_seq1->pos_next = pos;
            break;
        }
    }

    //dimension = s.size() - 1;
    while(true){
        int pos = (unsigned)rand() % dimension + 1;
        if(s[pos-1] == sub_seq2->prev || s[pos] == sub_seq2->next){
            continue;
        }else{
            sub_seq2->pos_next = pos;
            break;
        }
    }

    if(sub_seq1->pos_next > sub_seq2->pos_next){
        s.insert(s.begin() + sub_seq1->pos_next, sub_seq1->seq.begin(), sub_seq1->seq.begin() + first_size);
        s.insert(s.begin() + sub_seq2->pos_next, sub_seq2->seq.begin(), sub_seq2->seq.end());
    }else{
        s.insert(s.begin() + sub_seq2->pos_next, sub_seq2->seq.begin(), sub_seq2->seq.end());
        s.insert(s.begin() + sub_seq1->pos_next, sub_seq1->seq.begin(), sub_seq1->seq.begin() + first_size);
    }

    delete sub_seq1;
    delete sub_seq2;
}


void GILS_RVND(int Imax, int Iils){

    alignas(alignof(std::vector<int>)) std::vector<int> s;
    alignas(alignof(std::vector<int>)) std::vector<int> sl;
    alignas(alignof(std::vector<int>)) std::vector<int> s_final;
    alignas(alignof(std::vector<std::vector<struct subseq>>)) std::vector<std::vector<struct subseq>> subseq_info (dimension+1, std::vector<struct subseq>(dimension+1)) ;
    sl.reserve(dimension+1);
    s.reserve(dimension+1);
    s_final.reserve(dimension+1);
    //subseq_info.reserve(dimension+1);
    double cost_rvnd_current;
    double cost_rvnd_best;
    double cost_sl;
    double cost_final;

    for(int i = 0; i < Imax; ++i){

        after();
        int aux = (unsigned)rand() % TABLE_SZ;
        double alpha = R_table(aux);

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	
        //after();
        construct(s, alpha);
        /*
        s= {1};
        candidates_load(s, dimension);
        s.push_back(1);
         */

        load_subseq_info(subseq_info, s);
        
        // inserir de 7 a 12(12+1) no 5 (5+1)(burma14)

        /*
        double dif;
        int j = 12;
        int it = 7;

        dif = subseq_info[0][it-1].C + subseq_info[j][j].W * (subseq_info[0][it-1].T + c[s[it-1]][s[j]]) +  c[s[j]][s[j]];
        dif = dif + subseq_info[it+1][j-1].W * (subseq_info[0][it-1].T + c[s[it-1]][s[j]] + subseq_info[j][j].T  + c[s[j]][s[it+1]] ) +  subseq_info[it+1][j-1].C;
        dif = dif + subseq_info[it][it].W * (subseq_info[0][it-1].T + c[s[it-1]][s[j]] + subseq_info[j][j].T  + c[s[j]][s[it+1]] + subseq_info[it+1][j-1].T+ c[s[j-1]][s[it]]) + subseq_info[it][it].C;
        dif = dif + subseq_info[j+1][dimension].W * (subseq_info[0][it-1].T + c[s[it-1]][s[j]] + subseq_info[j][j].T  + c[s[j]][s[it+1]] + subseq_info[it+1][j-1].T + c[s[j-1]][s[it]] + subseq_info[it][it].T + c[s[it]][s[j+1]]) + subseq_info[j+1][dimension].C;

        std::cout << "manual : " << dif << std::endl;


        swap_2(s, 7, 12);
        load_subseq_info(subseq_info, s, 0);

        //std::cout << "pos : " << subseq_info[0][dimension].C << std::endl;
        //std::cout << "pos : " << cost_acumulated_calc(s, 0, dimension) << std::endl;
        exit(1);
        */

        sl = s;
        int Iterils = 0;
        //before(6);

        printf("\t[+] Looking for the best Neighbor..\n");
        //cost_calc(sl, &cost_rvnd_best);
        cost_rvnd_best = subseq_info[0][dimension].C - DBL_EPSILON;
        //std::cout << cost_rvnd_best << std::endl;
        //std::cout << cost_acumulated_calc(s, 0, dimension) << std::endl;

        while(Iterils < Iils){
            RVND(s, subseq_info);
            load_subseq_info(subseq_info, s);
            cost_rvnd_current = subseq_info[0][dimension].C - DBL_EPSILON;
            //cost_calc(s, &cost_rvnd_current);
            if(cost_rvnd_current < cost_rvnd_best){
                sl.clear();
                sl = s;
                cost_rvnd_best = cost_rvnd_current - DBL_EPSILON;
                //std::cout << cost_rvnd_best << std::endl;
                //cost_calc(sl, &cost_rvnd_best);
                Iterils = 0;
            }
            perturb(sl, s);
            load_subseq_info(subseq_info, s);
            Iterils++;
            //printf("aqui\n");
        }

        //cost_calc(sl, &cost_sl);
        load_subseq_info(subseq_info, sl);
        cost_sl = subseq_info[0][dimension].C - DBL_EPSILON;
        if(cost_sl < cost_final - DBL_EPSILON || i == 0){
            s_final.clear();
            s_final =  s;
            cost_final = cost_sl;
        }

        before(7);

        std::cout << "\tCurrent best cost: "<< cost_final << std::endl;
        std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;

    }
    std::cout << "COST: " << cost_final << std::endl;
}


int main(int argc, char **argv){
    int Imax = 10;
    int Iils;

    flag = true;

    srand(duration_cast<microseconds>(system_clock::now().time_since_epoch()).count());
    readData(argc, argv, &dimension, &c);

    int ar[] = {100, dimension};
    
    Iils = ar[dimension < 100];

    auto t1 = high_resolution_clock::now();

    GILS_RVND(Imax, Iils);

    auto t2 = high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    double res = (double)duration / 10e2;
    std::cout << "TIME: " << res << std::endl;

    if(flag){
        std::cout << "Construction time: " << construct_t/10e5 << std::endl;
        std::cout << "Swap time: " << swap_t/10e5 << std::endl;
        std::cout << "two_opt time: " << two_opt_t/10e5 << std::endl;
        std::cout << "reinsertion time: " << reinsertion_t/10e5 << std::endl;
        std::cout << "or_opt2 time: " << opt2_t/10e5 << std::endl;
        std::cout << "or_opt3 time: " << opt3_t /10e5<< std::endl;
    }
    return 0;
}

