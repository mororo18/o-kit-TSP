#include <iostream>
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

auto swap_time = 0;

struct insertion_info {
  int node_new;
  int edge_removed;
  int cost;
};

struct neighbor_info{
  std::vector<int> neighbor;
  int i;
  int j;
  int pos_new;
  int cost_dif;
};

struct sub_info{
  std::vector<int> seq;
  int pos_next;
  int size;
  int prev;
  int next;
};

double **c;
int dimension;
int cost_rvnd_best;
std::vector<int> candidates;

std::vector<int> s_final;
int cost_final;

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

std::vector<int> construct(double alpha){

	std::vector<int> s = { 1, 1};
	candidates_load(candidates, dimension);
	int subtour_inicial = 3;

	sum_t = 0;
	for(unsigned i = 0; i < candidates.size(); i++)
		sum_t += candidates[i];

	for(int i = 0; i < subtour_inicial; i++){
		int node_rand_index = rand() % candidates.size();
		int node_rand = candidates[node_rand_index];

		s.insert(s.begin() + 1, node_rand);
		candidates.erase(candidates.begin() + node_rand_index);
	}

	
	while(candidates.size() != 0){


		int insertion_possiblts = ((s.size() - 1) * candidates.size());
		//printf("possibilidades  %d\n", insertion_possiblts);
		std::vector<struct insertion_info> insertion_cost (insertion_possiblts);

		int l = 0; // 'insertion_cost' index

		for(auto k : candidates){
			for(int i = 0, j = 1; i < s.size() - 1; i++, j++){   
				
				insertion_cost[l].cost = c[s[i]][k] + c[s[j]][k] - c[s[i]][s[j]];  
				insertion_cost[l].node_new = k;
				insertion_cost[l].edge_removed = i; 
				l++;

			}
		}
		//printf("valor de l  =  %d\n", l);

		
		int node_rand_range = (int) (alpha * insertion_cost.size());
		int node_rand_index = rand() % node_rand_range;
		std::partial_sort (insertion_cost.begin(), insertion_cost.begin() + node_rand_range, insertion_cost.end(), cost_compare);

		/*
		std::cout << "alpha :" << alpha << std::endl;
		std::cout << "size : " << insertion_cost.size()  << std::endl;
		std::cout << node_rand_range << std::endl;
		std::cout << node_rand_index << std::endl;
		*/
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


	return s;
}	

void swap_2(struct neighbor_info &cheapest, int i, int j){
	int aux;

	aux = cheapest.neighbor[i];
	cheapest.neighbor[i] = cheapest.neighbor[j];
	cheapest.neighbor[j] = aux;
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
	  std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

void two_opt(struct neighbor_info &cheapest, int i, int j){
	std::reverse(cheapest.neighbor.begin() + i, cheapest.neighbor.begin() + j+1);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
	  std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

void reinsert(struct neighbor_info &cheapest, int i, int j, int pos){
	int size = j - i + 1;
	std::vector<int> copy;
	if(pos < i){
		if(size == 1){
			int aux = cheapest.neighbor[i];
			cheapest.neighbor.erase(cheapest.neighbor.begin() + i);
			cheapest.neighbor.insert(cheapest.neighbor.begin() + pos, aux);
		}else{
			copy.insert(copy.begin(), cheapest.neighbor.begin() + i, cheapest.neighbor.begin() + j+1);
			cheapest.neighbor.erase(cheapest.neighbor.begin() + i, cheapest.neighbor.begin() + j+1);
			cheapest.neighbor.insert(cheapest.neighbor.begin() + pos, copy.begin(), copy.end());
					/*
					for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
						std::cout << " " << cheapest->neighbor[i];

					std::cout << std::endl;
					*/
			
		}
	}else{
		if(i == j){
			cheapest.neighbor.insert(cheapest.neighbor.begin() + pos, cheapest.neighbor[i]);
			cheapest.neighbor.erase(cheapest.neighbor.begin() + i);
		}else{
			copy.insert(copy.begin(), cheapest.neighbor.begin() + i, cheapest.neighbor.begin() + j+1);
			cheapest.neighbor.insert(cheapest.neighbor.begin() + pos, copy.begin(), copy.end());
			cheapest.neighbor.erase(cheapest.neighbor.begin() + i, cheapest.neighbor.begin() + j+1);
		}
	}

}

void neighbor_swap_better(struct neighbor_info &cheapest, std::vector<int> &s){
	int dif;
	int dif_lower;

	for(int i = 1; i < s.size() - 2; i++){
		for(int j = i + 1; j < s.size() - 1; j++){
			dif = 0;
			// subtract the old distances
			dif += - c[s[i]][s[i-1]] - c[s[i]][s[i+1]] - c[s[j]][s[j-1]] - c[s[j]][s[j+1]];
			// add the new distances
			dif += c[s[i]][s[j-1]] + c[s[i]][s[j+1]] + c[s[j]][s[i-1]] + c[s[j]][s[i+1]];

			if(i + 1 == j) //consecutive
				dif += 2 * c[s[i]][s[j]];

			//std::cout << " " << dif;

			if((i == 1 && j == 2) || dif < dif_lower){
				dif_lower = dif;
				cheapest.i = i;
				cheapest.j = j;
				cheapest.cost_dif = dif_lower;

			}
		}
	}
	//printf("\n");

	swap_2(cheapest, cheapest.i, cheapest.j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_two_opt_better(struct neighbor_info &cheapest, std::vector<int> &s){
	double dif;
	double dif_lower;
	int j;
	int sub_sz = rand() % (s.size() - 2);
	while(sub_sz == 1 || !sub_sz)
		sub_sz = rand() % (s.size() - 2);

	
	/*
	printf("subsequence size =  %d\n", sub_sz);
	*/

	for(int i = 1; i < s.size() - sub_sz; i++){
		dif = 0;
		j = i + sub_sz-1; // modificacao aqui
		// subtract the old distances
		dif += (- c[s[i]][s[i-1]]  - c[s[j]][s[j+1]]);
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

	two_opt(cheapest, cheapest.i, cheapest.j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_reinsertion_better(struct neighbor_info &cheapest, std::vector<int> &s, int sz){
	int dif;
	int sze = sz;
	int dif_lower;
	int j;
	int pos_new;
	int t = 1;

	for(int i = 1; i < s.size() - sze; i++){
		//k : arestas 
		j = i + sze - 1;
		for(int k = 0; k < s.size() - 1; k++){
			dif = 0;

			if(k < (i - 1) || k > j){
				pos_new = k + 1;
			}else if(k >= (i - 1) && k <= j){
				continue;
			}		
				
			// subtract the old distances
			dif += (- c[s[i]][s[i-1]] - c[s[j]][s[j+1]] - c[s[pos_new]][s[pos_new-1]]);
			// add the new distances
			dif += (c[s[i]][s[pos_new-1]] + c[s[pos_new]][s[j]] + c[s[i-1]][s[j+1]]);

			//std::cout << " " << dif;

			if( t || dif < dif_lower){
				dif_lower = dif;
				cheapest.i = i;
				cheapest.j = j;
				cheapest.pos_new = pos_new;
				cheapest.cost_dif = dif_lower;

			}
			t = 0;
		}
		//printf("\n");
	}
	//printf("\n");

	//printf("\nsize = %d\ni = %d\nj = %d\npos = %d\n", sze, cheapest->i, cheapest->j, cheapest->pos_new);
	reinsert(cheapest, cheapest.i, cheapest.j, cheapest.pos_new);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	*/
}

int cost_calc(std::vector<int> &s){
	int cost_best = 0;
	for(unsigned i = 0; i < s.size() - 1; i++)
		cost_best += c[s[i]][s[i+1]];

	return cost_best;
}

void neighbd_list_repopulate(std::vector<int> &list){
	list.clear();
	list = {1, 2, 3, 4, 5};
	//printf("REPOPULATE\n");
}
		

//long l=0;

void RVND(std::vector<int> &s){
	
	//printf("for -> while -> RVND\n");
	std::vector<int> neighbd_list = {1, 2, 3, 4, 5};
	//l++;
	while(!neighbd_list.empty()){
		//printf("for -> while -> RVND %ld vezes ", l);
		//printf("preco atual = %d\n", cost_final);
		
		int neighbd_rand_index = rand() % neighbd_list.size();
		int neighbd_rand = neighbd_list[neighbd_rand_index];


		struct neighbor_info *cheapest  (new struct neighbor_info) ;
		//std::copy(s.begin(), s.end(), cheapest->neighbor.begin());
		cheapest->neighbor = s;

		
		switch(neighbd_rand){
			case SWAP:
				{
				//printf("SWAP\n");
				//auto t1 = high_resolution_clock::now();
				neighbor_swap_better(*(cheapest), s);
				//auto t2 = high_resolution_clock::now();
				//auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
				//swap_time += duration;
				}
				break;
			case TWO_OPT:
				//printf("TWO_OPT ");
				neighbor_two_opt_better(*(cheapest), s);
				break;				
			case REINSERTION:
				//printf("REINSERTION ");
				neighbor_reinsertion_better(*(cheapest), s, REINSERTION);
				break;				
			case OR_OPT2:
				//printf("OR_OPT2 ");
				neighbor_reinsertion_better(*(cheapest), s, OR_OPT2);
				break;				
			case OR_OPT3:
				//printf("OR_OPT3 ");
				neighbor_reinsertion_better(*(cheapest), s, OR_OPT3);
				break;				
		}
		
		int sum = 0;
		for(int i = 1; i < cheapest->neighbor.size() -1; i++)
			sum += cheapest->neighbor[i];

		if(sum != sum_t){
			//printf("sum %d\n"
					//"sum_t %d\n", sum, sum_t);
			exit(-1);
		}

		if(cheapest->cost_dif < 0){
			s.clear();
			s = cheapest->neighbor;
			cost_rvnd_best = cost_calc(s);
			neighbd_list_repopulate(neighbd_list);
		}else{
			neighbd_list.erase(neighbd_list.begin() + neighbd_rand_index);
		}

		delete(cheapest);
	}
}

void perturb(std::vector<int> &sl, std::vector<int> &s){
	//printf("for -> perturb\n");

	s.clear();
	s = sl;
	int dimension = s.size() - 1;
	int first_size = (rand() % (int)(dimension / 10)) + 2;
	int first_pos = rand() % (dimension - first_size - 1) + 1;

	struct sub_info *sub_seq1 (new struct sub_info);
	sub_seq1->size = first_size;
	sub_seq1->prev = sl[first_pos - 1];
	sub_seq1->next = sl[first_pos + first_size];
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
	struct sub_info *sub_seq2 (new struct sub_info);
	sub_seq2->size = second_size;
	sub_seq2->prev = sl[second_pos - 1];
	sub_seq2->next = sl[second_pos + second_size];
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
		s.insert(s.begin() + sub_seq1->pos_next, sub_seq1->seq.begin(), sub_seq1->seq.end());
		s.insert(s.begin() + sub_seq2->pos_next, sub_seq2->seq.begin(), sub_seq2->seq.end());
	}else{
		s.insert(s.begin() + sub_seq2->pos_next, sub_seq2->seq.begin(), sub_seq2->seq.end());
		s.insert(s.begin() + sub_seq1->pos_next, sub_seq1->seq.begin(), sub_seq1->seq.end());
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

	delete(sub_seq1);
	delete(sub_seq2);
}

void GILS_RVND(int Imax, int Iils){
	for(int i = 0; i < Imax; i++){
		int aux = rand() % 10 + 1;
		double alpha = 1.0 /(double) aux;
		//double alpha = 0.5;
		//printf("alpha  = %lf\n", alpha);

		printf("[+] Search %d\n", i);
		printf("\t[+] Constructing..\n");	
		std::vector<int> s = construct(alpha);
		std::vector<int> sl = s;
		int Iterils = 0;

		printf("\t[+] Looking for the best Neighbor..\n");
		int cost_rvnd_current = cost_calc(sl);
		while(Iterils < Iils){
			RVND(s);
			if(cost_rvnd_best < cost_rvnd_current){
				sl.clear();
				sl = s;
				cost_rvnd_current = cost_calc(sl);
				Iterils = 0;
			}
			perturb(sl, s);
			Iterils++;
			//printf("iter ils %d\n", Iterils);
		
		}

		int cost_sl = cost_calc(sl);
		if(i == 0 || cost_sl < cost_final){
			s_final =  s;
			cost_final = cost_sl;
		}

		printf("\tCurrent best result: \n");
		for(unsigned i = 0; i <  s_final.size(); i++)
			std::cout << " " << s_final[i];

		std::cout << std::endl;
		printf("\tCurrent best cost: %d\n", cost_final);
	
	}
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
	printf("better cost = %d\n", cost_final);

	auto t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	double res = (double)duration / 10e2;
	std::cout << res << std::endl;
	return 0;
}

