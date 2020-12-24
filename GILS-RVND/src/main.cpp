#include <iostream>
#include <new>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "readData.h"

#define REINSERTION 1
#define OR_OPT2 	2
#define OR_OPT3 	3
#define SWAP 		4
#define TWO_OPT		5

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::system_clock;
using std::chrono::high_resolution_clock;


struct insertion_info {
  int node_new;
  int edge_removed;
  int cost;
};

struct neighbor_info{
  //std::vector<int> neighbor;
  int i;
  int j;
  int pos_new;
  int cost_dif;
};

struct sub_info{
  std::vector<int> seq;
  int pos_next;
  int prev;
  int next;
};

double **c;
int dimension;
int cost_rvnd_current;
std::vector<int> candidates;

std::chrono::_V2::system_clock::time_point t3;
std::chrono::_V2::system_clock::time_point t4;

long swap_t  = 0;
long two_opt_t  = 0;
long reinsertion_t  = 0;
long opt2_t = 0;
long opt3_t = 0;

void after(){
	t3 = high_resolution_clock::now();
}

void before(int a){
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
	}
}

bool cost_compare(const insertion_info &a, const insertion_info &b){
	return a.cost < b.cost;
}

void candidates_load(std::vector<int> &cand, int dim){
	for(int i = 2; i <= dim; i++){
		cand.insert(cand.begin() + i - 2, i);
	}
	/*
	for(unsigned i = 0; i < cand.size(); i++)
		std::cout << " " << cand[i];

	std::cout << std::endl;
	*/
}	

int sum_t = 0;

void construct(std::vector<int> &s, double alpha){

	//s.clear();
	s = { 1, 1};
	candidates.reserve(dimension);
	candidates_load(candidates, dimension);
	int subtour_inicial = 3;

	/*
	sum_t = 0;
	for(unsigned i = 0; i < candidates.size(); i++)
		sum_t += candidates[i];
		*/

	for(int i = 0; i < subtour_inicial; i++){
		int node_rand_index = rand() % candidates.size();
		int node_rand = candidates[node_rand_index];

		s.insert(s.begin() + 1, node_rand);
		candidates.erase(candidates.begin() + node_rand_index);
	}

	
	while(!candidates.empty()){


		int insertion_possiblts = ((s.size() - 1) * candidates.size());
		//printf("possibilidades  %d\n", insertion_possiblts);
		std::vector<struct insertion_info> insertion_cost (insertion_possiblts);

		int l = 0; // 'insertion_cost' index

		int dim = s.size() - 1;
		for(auto k : candidates){
			for(int i = 0, j = 1; i < dim; i++, j++){   
				
				insertion_cost[l].cost = c[s[i]][k] + c[s[j]][k] - c[s[i]][s[j]];  
				insertion_cost[l].node_new = k;
				insertion_cost[l].edge_removed = i; 
				l++;

			}
		}
		//printf("valor de l  =  %d\n", l);

		
		int node_rand_range = (int) (alpha * insertion_cost.size());
		int node_rand_index = rand() % node_rand_range;
		std::partial_sort(insertion_cost.begin(), insertion_cost.begin() + node_rand_range, insertion_cost.end(), cost_compare);

		int node_rand = insertion_cost[node_rand_index].node_new;
		int edge = insertion_cost[node_rand_index].edge_removed;
		s.insert(s.begin() + edge + 1, node_rand); 
		//printf("inserted node: %d\n", node_rand);

		for(int i = 0; i < candidates.size(); i++){
			if(candidates[i] ==  node_rand){
				candidates.erase(candidates.begin() + i);
				break;
			}
		};

		/*
		for(unsigned i = 0; i < s.size(); i++)
		  std::cout << " " << s[i];

		std::cout << std::endl;
		*/
	}

}	

void swap_2(std::vector<int> &vec, /*struct neighbor_info &cheapest,*/ int i, int j){
	/*
	int aux;

	aux = cheapest.neighbor[i];
	cheapest.neighbor[i] = cheapest.neighbor[j];
	cheapest.neighbor[j] = aux;
	*/
	std::iter_swap(vec.begin() + i, vec.begin() + j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
	  std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

void two_opt(std::vector<int> &vec, /*struct neighbor_info &cheapest,*/ int i, int j){
	std::reverse(vec.begin() + i, vec.begin() + j+1);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
	  std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

void reinsert(std::vector<int> &vec, /*struct neighbor_info &cheapest, */ int i, int j, int pos){
	if(pos < i){
		if(i == j){
			int aux = vec[i];
			vec.erase(vec.begin() + i);
			vec.insert(vec.begin() + pos, aux);
		}else{
			std::vector<int> copy;
			copy.reserve(j - i + 1);
			copy.insert(copy.begin(), vec.begin() + i, vec.begin() + j+1);
			vec.erase(vec.begin() + i, vec.begin() + j+1);
			vec.insert(vec.begin() + pos, copy.begin(), copy.end());
		}
	}else{
		if(i == j){
			vec.insert(vec.begin() + pos, vec[i]);
			vec.erase(vec.begin() + i);
		}else{
			std::vector<int> copy;
			copy.reserve(j - i + 1);
			copy.insert(copy.begin(), vec.begin() + i, vec.begin() + j+1);
			vec.insert(vec.begin() + pos, copy.begin(), copy.end());
			vec.erase(vec.begin() + i, vec.begin() + j+1);
		}
	}

}

void neighbor_swap_better(struct neighbor_info &cheapest, std::vector<int> &vec, std::vector<int> &s){
	int dif_lower;

	for(int i = 1; i < dimension - 1; i++){
		int dif1 = - c[s[i]][s[i-1]] - c[s[i]][s[i+1]];
		for(int j = i + 1; j < dimension; j++){
			int dif = 0;
			// subtract the old distances
			dif +=  - c[s[j]][s[j-1]] - c[s[j]][s[j+1]];
			// add the new distances
			dif += c[s[i]][s[j-1]] + c[s[i]][s[j+1]] + c[s[j]][s[i-1]] + c[s[j]][s[i+1]];

			if(i + 1 == j) //consecutive
				dif += 2 * c[s[i]][s[j]];

			//std::cout << " " << dif;
			dif += dif1;
			if((i == 1 && j == 2) || dif < dif_lower){
				dif_lower = dif;
				cheapest.i = i;
				cheapest.j = j;
				cheapest.cost_dif = dif_lower;

			}
		}
	}
	//printf("\n");

	if(cheapest.cost_dif < 0)
	  swap_2(vec, cheapest.i, cheapest.j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_two_opt_better(struct neighbor_info &cheapest, std::vector<int> &vec, std::vector<int> &s){
	int dif_lower;
	int sub_sz = rand() % (dimension - 1);
	while(sub_sz == 1 || !sub_sz)
		sub_sz = rand() % (dimension - 1);

	
	/*
	printf("subsequence size =  %d\n", sub_sz);
	*/

	for(int i = 1; i < (dimension+1) - sub_sz; i++){
		int dif = 0;
		int j = i + sub_sz-1; // modificacao aqui
		// subtract the old distances
		dif += - c[s[i]][s[i-1]]  - c[s[j]][s[j+1]];
		// add the new distances
		dif +=  (c[s[i]][s[j+1]] + c[s[j]][s[i-1]]);

		//std::cout << " " << dif;

		if((i == 1) || dif < dif_lower){
		  dif_lower = dif;
		  cheapest.i = i;
		  cheapest.j = j;
		  cheapest.cost_dif = dif_lower;

		}
	}
	//printf("\n");

	if(cheapest.cost_dif < 0)
	  two_opt(vec, cheapest.i, cheapest.j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_reinsertion_better(struct neighbor_info &cheapest, std::vector<int> &vec, std::vector<int> &s, int sz){
	int dif_lower;
	int t = 1;

	for(int i = 1; i < (dimension+1) - sz; i++){
		//k : arestas 
		int dif1;
		int j = i + sz - 1;
		// subtract the old distances
		dif1 = (- c[s[i]][s[i-1]] - c[s[j]][s[j+1]] + c[s[i-1]][s[j+1]]);
		for(int k = 0; k < dimension; k++){
			int pos_new;
			int dif;
			int dif2;
			dif = 0;

			if(k < (i - 1) || k > j){
				pos_new = k + 1;
			}else if(k >= (i - 1) && k <= j){
				continue;
			}		
				
			// add the new distances
			dif2 = (c[s[i]][s[pos_new-1]] + c[s[pos_new]][s[j]] - c[s[pos_new]][s[pos_new-1]]); 
			dif = dif1 + dif2;
			//std::cout << " " << dif;

			if( t || dif < dif_lower){
				dif_lower = dif;
				cheapest.i = i;
				cheapest.j = j;
				cheapest.pos_new = pos_new;
				cheapest.cost_dif = dif_lower;
				t = 0;
	
			}
		}
		//printf("\n");
	}
	//printf("\n");

	//printf("\nsize = %d\ni = %d\nj = %d\npos = %d\n", sze, cheapest->i, cheapest->j, cheapest->pos_new);
	if(cheapest.cost_dif < 0)
		reinsert(vec, cheapest.i, cheapest.j, cheapest.pos_new);

	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

void cost_calc(std::vector<int> &s, int *cost_best){
	*cost_best = 0;
	for(unsigned i = 0; i < dimension; i++)
		*cost_best += c[s[i]][s[i+1]];
}

void neighbd_list_repopulate(std::vector<int> &list){
	list.clear();
	list = {1, 2, 3, 4, 5};
	//printf("REPOPULATE\n");
}
		
void RVND(std::vector<int> &s){
	
	//printf("for -> while -> RVND\n");
	std::vector<int> neighbd_list = {1, 2, 3, 4, 5};

	while(!neighbd_list.empty()){
		//printf("for -> while -> RVND %ld vezes ", l);
		
		int neighbd_rand_index = rand() % neighbd_list.size();
		int neighbd_rand = neighbd_list[neighbd_rand_index];


		//struct neighbor_info *cheapest = (struct neighbor_info*)malloc(sizeof(struct neighbor_info)+sizeof(int)*(dimension+1));
		struct neighbor_info cheapest ;//=  new(std::nothrow) struct neighbor_info;
		std::vector<int> cheapest_vec;
		cheapest_vec.reserve(dimension+1);
		cheapest_vec = s;

		
		switch(neighbd_rand){
			case SWAP:
				//printf("SWAP\n");
				after();
				neighbor_swap_better(cheapest, cheapest_vec, s);
				before(SWAP);
				break;
			case TWO_OPT:
				//printf("TWO_OPT \n");
				after();
				neighbor_two_opt_better(cheapest, cheapest_vec, s);
				before(TWO_OPT);
				break;				
			case REINSERTION:
				//printf("REINSERTION \n");
				after();
				neighbor_reinsertion_better(cheapest, cheapest_vec, s, REINSERTION);
				before(REINSERTION);
				break;				
			case OR_OPT2:
				//printf("OR_OPT2 \n");
				after();
				neighbor_reinsertion_better(cheapest, cheapest_vec, s, OR_OPT2);
				before(OR_OPT2);
				break;				
			case OR_OPT3:
				//printf("OR_OPT3 \n");
				after();
				neighbor_reinsertion_better(cheapest, cheapest_vec, s, OR_OPT3);
				before(OR_OPT3);
				break;				
		}
		
		/*
		int sum = 0;
		for(int i = 1; i < cheapest->neighbor.size() -1; i++)
			sum += cheapest->neighbor[i];

		if(sum != sum_t){
			//printf("sum %d\n"
					//"sum_t %d\n", sum, sum_t);
			exit(-1);
		}
		*/

		if(cheapest.cost_dif < 0){
			s.clear();
			s = cheapest_vec;
			cost_calc(s, &cost_rvnd_current);
			neighbd_list_repopulate(neighbd_list);
		}else{
			neighbd_list.erase(neighbd_list.begin() + neighbd_rand_index);
		}

		//delete cheapest;
		//free(cheapest);
		//delete &(cheapest->neighbor);
		/*
		if (cheapest != NULL)     
			operator delete cheapest;
			*/
	}
}

void perturb(std::vector<int> &sl, std::vector<int> &s){
	//printf("for -> perturb\n");

	s.clear();
	s = sl;
	int dimension = s.size() - 1;
	int first_size = (rand() % (int)(dimension / 10)) + 2;
	int first_pos = rand() % (dimension - first_size - 1) + 1;

	struct sub_info *sub_seq1 = new(std::nothrow) struct sub_info;
	sub_seq1->prev = sl[first_pos - 1];
	sub_seq1->next = sl[first_pos + first_size];
	sub_seq1->seq.reserve(first_size);
	sub_seq1->seq.insert(sub_seq1->seq.begin(), sl.begin() + first_pos, sl.begin() + first_pos + first_size);

	int second_size = (rand() % ((int)(dimension / 10))) + 2;
	int second_pos;
	// second_pos restrictions
	while(true){
		int pos = rand() % (dimension - second_size - 1) + 1;
		if((pos <( first_pos - second_size) || pos > first_pos + first_size)){
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
		int pos = rand() % dimension + 1;
		if(s[pos-1] == sub_seq1->prev || s[pos] == sub_seq1->next){
			continue;
		}else{
			sub_seq1->pos_next = pos;
			break;
		}
	}

	dimension = s.size() - 1;
	while(true){
		int pos = rand() % dimension + 1;
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

	/*
	printf("1a pos = %d\n"
			"1a next = %d\n"
			"1o tamanho = %d\n"
			"2a pos = %d\n"
			"2a next = %d\n"
			"2o tamanho = %d\n", first_pos, sub_seq1->pos_next, first_size, second_pos, sub_seq2->pos_next, second_size);
			*/

	/*
	for(unsigned i = 0; i <  copy.size(); i++)
		std::cout << " " << copy[i];

	std::cout << std::endl;
	*/

	delete sub_seq1;
	delete sub_seq2;
}

void GILS_RVND(int Imax, int Iils){
	std::vector<int> s_final;
	s_final.reserve(dimension+1);
	int cost_final;
	for(int i = 0; i < Imax; i++){
		int aux = rand() % 10 + 1;
		double alpha = 1.0 / aux;
		//double alpha = 0.5;
		//printf("alpha  = %lf\n", alpha);

		printf("[+] Search %d\n", i);
		printf("\t[+] Constructing..\n");	
		std::vector<int> s;
		std::vector<int> sl;
		s.reserve(dimension+1);
		sl.reserve(dimension+1);
		construct(s, alpha);
		sl = s;
		int Iterils = 0;

		printf("\t[+] Looking for the best Neighbor..\n");
		int cost_rvnd_best;
		cost_calc(sl, &cost_rvnd_best);
		while(Iterils < Iils){
			RVND(s);
			if(cost_rvnd_current < cost_rvnd_best){
				sl.clear();
				sl = s;
				cost_calc(sl, &cost_rvnd_best);
				Iterils = 0;
			}
			perturb(sl, s);
			Iterils++;
			//printf("iter ils %d\n", Iterils);
		
		}

		int cost_sl;
		cost_calc(sl, &cost_sl);
		if(i == 0 || cost_sl < cost_final){
			s_final.clear();
			s_final =  s;
			cost_final = cost_sl;
		}

		printf("\tCurrent best result: \n");
		for(unsigned i = 0; i <  dimension + 1; i++)
			std::cout << " " << s_final[i];

		std::cout << std::endl;
		printf("\tCurrent best cost: %d\n", cost_final);
	
	}
	printf("better cost = %d\n", cost_final);
}


int main(int argc, char **argv){
	int Imax = 50;
	int Iils;

	srand(duration_cast<microseconds>(system_clock::now().time_since_epoch()).count());
	readData(argc, argv, &dimension, &c);
	/*
	std::vector<int> sa = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
	int sum = 0;
	for(int i = 1; i< sa.size() - 1; i++)
		sum+=sa[i];
	while(1){
	  int sun = 0;
	  sa = perturb(sa);
	  for(int i = 1; i< sa.size() - 1; i++)
		  sun+=sa[i];

	  if(sun != sum)
	  	break;
	}

	*/
	if(dimension >= 150)
		Iils = dimension / 2;
	else
		Iils = dimension;


	auto t1 = high_resolution_clock::now();

	GILS_RVND(Imax, Iils);

	auto t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	double res = (double)duration / 10e2;
	std::cout << "Total time: " << res << std::endl;
	std::cout << "Swap time: " << swap_t/10e5 << std::endl;
	std::cout << "two_opt time: " << two_opt_t/10e5 << std::endl;
	std::cout << "reinsertion time: " << reinsertion_t/10e5 << std::endl;
	std::cout << "or_opt2 time: " << opt2_t/10e5 << std::endl;
	std::cout << "or_opt3 time: " << opt3_t /10e5<< std::endl;
	return 0;
}

