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
#include <sched.h>

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
    float cost;
    int edge_removed;
    int node_new;
};

struct sub_info{
    std::vector<int> seq;
    int prev;
    int next;
    int pos_next;
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


//-----  to measure the execution time BEGIN  -----

void before(){
    if(!flag)
        return;
    t3 = high_resolution_clock::now();
}

void after(int a){
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

//-----  to measure the execution time END  -----

bool cost_compare(const insertion_info &a, const insertion_info &b){
    return a.cost < b.cost - DBL_EPSILON;
}

void candidates_load(std::vector<int> &cand, int dim){
    for(int i = 2; i <= dim; ++i){
        cand.push_back(i);
    }
}	

inline void binary_search(std::vector<int> &vec, int node, short from, short to){
    #define MIDDLE ((from + to) / 2)

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

    for(int i = subtour_inicial; --i;){
        
        int node_rand_index = (unsigned)rand() % candidates.size(); 
        int node_rand = candidates[node_rand_index];

        s.insert(s.begin() + 1, node_rand);
        candidates.erase(candidates.begin() + node_rand_index);
    }


    while(!candidates.empty()){

        int dim = s.size() - 1;

        int insertion_possiblts = (dim * candidates.size());
        std::vector<struct insertion_info> insertion_cost (insertion_possiblts);

        for(int i = 0, j = 1, l = 0; i < dim; ++i, ++j){   
            double a = c[s[i]][s[j]];
            for(auto k : candidates){

                insertion_cost[l].cost = c[s[i]][k] - a + c[s[j]][k] ;  
                insertion_cost[l].node_new = k;
                insertion_cost[l].edge_removed = i; 
                ++l;

            }
        }

        int node_rand_range = ceil(alpha * insertion_cost.size());
        int node_rand_index = (unsigned)rand() % node_rand_range;
        std::partial_sort(insertion_cost.begin(), insertion_cost.begin() + node_rand_range, insertion_cost.end(), cost_compare);

        int node_rand = insertion_cost[node_rand_index].node_new;
        int edge = insertion_cost[node_rand_index].edge_removed;
        s.insert(s.begin() + edge + 1, node_rand); 

        binary_search(candidates, node_rand, 0, candidates.size() );

    }

}	

inline void swap_2(std::vector<int> &vec, int i, int j){
    std::iter_swap(vec.begin() + i, vec.begin() + j);
}

inline void two_opt(std::vector<int> &vec, int i, int j){
    std::reverse(vec.begin() + i, vec.begin() + j+1);
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

#define DBL_SZ 8
#define INT_SZ 4

inline void neighbor_swap_better(std::vector<int> &s){
    alignas(DBL_SZ) double dif, dif1, dif3, dif_lower;
    alignas(INT_SZ) int i, j, e, b, d, a, o, p, q;
    alignas(INT_SZ) int dim = dimension - 2;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool l, t = 1;

    for(i = 1; i < dimension - 1; ++i){
        o = s[i];
        d = i - 1;
        a = i + 1;
        q = a + 1;
        l = false;

        dif1 = - c[s[i]][s[d]];
        dif3 = dif1 - c[s[i]][s[a]];

        //consecutive nodes
        dif = (dif1  + c[s[i]][s[q]]) + (c[s[a]][s[d]] - c[s[a]][s[q]]);

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

            dif = ((dif3 + c[o][s[e]]) + (c[o][s[b]] + c[p][s[d]])) + ((c[p][s[a]] - c[p][s[e]]) - c[p][s[b]]);

            if(dif < dif_lower){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
            }
        }
    }

    if(dif_lower < - DBL_EPSILON){
        swap_2(s, i_best, j_best);
        state = true;
    }
}

inline void neighbor_two_opt_better(std::vector<int> &s){
    alignas(DBL_SZ) double dif, dif1, dif_lower;
    alignas(INT_SZ) int b, a, j, i;
    alignas(INT_SZ) int i_best;
    alignas(INT_SZ) int j_best;
    alignas(1) bool t = 1;

    for(i = 1; i < dimension - 1; ++i){
        b = s[i - 1];
        dif1 = - c[s[i]][b];
        for(j = i + 2; j < dimension; ++j){
            a = j + 1;

            dif = (dif1 - c[s[j]][s[a]])  + (c[s[j]][b] + c[s[i]][s[a]]);

            if(dif < dif_lower || t){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                t = 0;

            }
        }
    }

    if(dif_lower < - DBL_EPSILON){
        two_opt(s, i_best, j_best);
        state = true;
    }
}


inline void neighbor_reinsertion_better(std::vector<int> &s, int sz){
    alignas(DBL_SZ) double dif, dif1, dif_lower;
    alignas(INT_SZ) int i, j, k, g, b, a, o, e;
    alignas(INT_SZ) int i_best, j_best, pos_new;
    alignas(INT_SZ) int k_lim = dimension - sz - 1;
    alignas(1) bool l;
    alignas(1) bool t = 1;

    for(i = 1, j = sz + i - 1; i < dimension - sz + 1; ++i, ++j){
        a = s[i];
        o = s[j];
        e = j + 1;
        b = i - 1;
        l = 1;

        dif1 = c[s[b]][s[e]] - c[s[i]][s[b]] - c[s[j]][s[e]];

        //k -> edges 
        for(k = 0; k < b; ++k){
            
            g = k + 1;

            dif = (dif1 + c[a][s[k]]) + (c[o][s[g]] - c[s[g]][s[k]]); 

            if( dif < dif_lower || t){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
                t = 0;

            }
        }

        for(k = i + sz; k < k_lim; ++k){
            
            g = k + 1;

            dif = (dif1 + c[a][s[k]]) + (c[o][s[g]] - c[s[g]][s[k]]); 

            if( dif < dif_lower){
                dif_lower = dif - DBL_EPSILON;
                i_best = i;
                j_best = j;
                pos_new = k;
            }
        }
    }

    if(dif_lower < - DBL_EPSILON){
        reinsert(s, i_best, j_best + 1, pos_new + 1);
        state = true;
    }

}

inline void cost_calc(std::vector<int> &s, double *cost_best){
    *cost_best = 0;
    double sum = 0;
    for(unsigned i = 0; i < (unsigned)dimension; ++i)
        sum += c[s[i]][s[i+1]];

    *cost_best = sum;
}

inline void neighbd_list_repopulate(std::vector<int> &list){
    list.clear();
    list = {1, 2, 3, 4, 5};
    //printf("REPOPULATE\n");
}

void RVND(std::vector<int> &s){

    //printf("for -> while -> RVND\n");
    alignas(8) std::vector<int> neighbd_list = {1, 2, 3, 4, 5};
    alignas(4) int neighbd_rand_index;
    alignas(4) int neighbd_rand;

    while(!neighbd_list.empty()){

        neighbd_rand_index = (unsigned)rand() % neighbd_list.size();
        neighbd_rand = neighbd_list[neighbd_rand_index];

        state = false;

        switch(neighbd_rand){
            case REINSERTION:
                //before();
                neighbor_reinsertion_better(s, REINSERTION);
                //after(REINSERTION);
                break;				
            case OR_OPT2:
                //before();
                neighbor_reinsertion_better(s, OR_OPT2);
                //after(OR_OPT2);
                break;				
            case OR_OPT3:
                //before();
                neighbor_reinsertion_better(s, OR_OPT3);
                //after(OR_OPT3);
                break;				
            case SWAP:
                //before();
                neighbor_swap_better(s);
                //after(SWAP);
                break;
            case TWO_OPT:
                //before();
                neighbor_two_opt_better(s);
                //after(TWO_OPT);
                break;				
        }

        if(state)
            neighbd_list_repopulate(neighbd_list);
        else
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
        s.erase(s.begin() + first_pos, s.begin() + first_pos + first_size);
        s.erase(s.begin() + second_pos, s.begin() + second_pos + second_size);
    }else{
        s.erase(s.begin() + second_pos, s.begin() + second_pos + second_size);
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
    sl.reserve(dimension+1);
    s.reserve(dimension+1);
    s_final.reserve(dimension+1);
    double cost_rvnd_current;
    double cost_rvnd_best;
    double cost_sl;
    double cost_final;

    for(int i = 0; i < Imax; ++i){

        before();
        int a = 99;
        int aux = (unsigned)rand() % a + 1;
        double alpha = 1.0 / aux;

        printf("[+] Search %d\n", i+1);
        printf("\t[+] Constructing..\n");	
        //before();
        construct(s, alpha);
        sl = s;
        int Iterils = 0;
        //after(6);

        printf("\t[+] Looking for the best Neighbor..\n");
        cost_calc(sl, &cost_rvnd_best);
        cost_rvnd_best -= DBL_EPSILON;

        while(Iterils < Iils){
            RVND(s);
            cost_calc(s, &cost_rvnd_current);
            if(cost_rvnd_current < cost_rvnd_best){
                sl.clear();
                sl = s;
                cost_rvnd_best = cost_rvnd_current - DBL_EPSILON;
                Iterils = 0;
            }
            perturb(sl, s);
            Iterils++;
        }

        cost_calc(sl, &cost_sl);
        if(cost_sl < cost_final - DBL_EPSILON || i == 0){
            s_final.clear();
            s_final =  s;
            cost_final = cost_sl;
        }

        after(7);

        std::cout << "\tCurrent best cost: "<< cost_final << std::endl;
        std::cout << "\tCurrent search time: "<< search_t / 10e5<< std::endl;
        std::cout << "\tCurrent search time average: "<< (search_t_average / (i+1)) / 10e5 << std::endl;

    }
    std::cout << "COST: " << cost_final << std::endl;
}


int main(int argc, char **argv){
    int Imax = 50;
    int Iils;

    flag = true;

    srand(duration_cast<microseconds>(system_clock::now().time_since_epoch()).count());
    readData(argc, argv, &dimension, &c);

    if(dimension >= 150)
        Iils = dimension / 2;
    else
        Iils = dimension;

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

