#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "readData.h"

struct insertion_info {
  int node_new;
  int edge_removed;
  double cost;
};

struct neighbor_info{
  std::vector<int> neighbor;
  int i;
  int j;
  int pos_new;
  double cost_dif;
};

double **c;
int dimension;
double cost_rvnd_best;
std::vector<int> candidatos;

std::vector<int> s_final;
double cost_final;

bool cost_compare(const insertion_info &a, const insertion_info &b){
	return a.cost < b.cost;
}

void candidates(std::vector<int> &cand, int dim){
	for(int i = 2; i <= dim; i++){
		cand.insert(cand.begin() + i - 2, i);
	}
	for(unsigned i = 0; i < cand.size(); i++)
		std::cout << " " << cand[i];

	std::cout << std::endl;
}	

std::vector<int> construcao(double alpha){
	
  	srand(time(NULL));

	std::vector<int> s = { 1, 1};
	candidates(candidatos, dimension);
	int subtour_inicial = 3;

	for(int i = 0; i < subtour_inicial; i++){
		int node_rand_index = rand() % candidatos.size();
		int node_rand = candidatos[node_rand_index];

		s.insert(s.begin() + 1, node_rand);
		candidatos.erase(candidatos.begin() + node_rand_index);
	}

	/*
	for(unsigned i = 0; i < s.size(); i++)
		std::cout << " " << s[i];

	std::cout << std::endl;
	*/
	
	while(candidatos.size() != 0){


		int insertion_possiblts = ((s.size() - 1) * candidatos.size());
		//printf("possibilidades  %d\n", insertion_possiblts);
		std::vector<struct insertion_info> insertion_cost (insertion_possiblts);

		int l = 0; // 'insertion_cost' index

		for(auto k : candidatos){
			for(int i = 0, j = 1; i < s.size() - 1; i++, j++){  //modificacao ('i < s.size()' para 'i <= s.size()') 
				
				insertion_cost[l].cost = c[s[i]][k] + c[s[j]][k] - c[s[i]][s[j]];  //modificacao temporaria
				insertion_cost[l].node_new = k;
				insertion_cost[l].edge_removed = i; 
				l++;

			}
		}
		//printf("valor de l  =  %d\n", l);

		std::sort (insertion_cost.begin(), insertion_cost.end(), cost_compare);
		
		int node_rand_range = (int) (alpha * insertion_cost.size());
		int node_rand_index = rand() % node_rand_range;

		int node_rand = insertion_cost[node_rand_index].node_new;
		int edge = insertion_cost[node_rand_index].edge_removed;
		s.insert(s.begin() + edge + 1, node_rand); 
		//printf("inserted node: %d\n", node_rand);

		for(int i = 0; i < candidatos.size(); i++){
			if(candidatos[i] ==  node_rand){
				candidatos.erase(candidatos.begin() + i);
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

void swap(std::vector<int> &s, int i, int j){
	int aux;

	aux = s[i];
	s[i] = s[j];
	s[j] = aux;
}

static std::vector<int> DEFAULT_VECTOR (0);

void two_opt(struct neighbor_info *cheapest = NULL, int i = 0, int j = 0, std::vector<int> &s = DEFAULT_VECTOR){
	if(cheapest != NULL)
		std::reverse(cheapest->neighbor.begin() + i, cheapest->neighbor.begin() + j);
	else
		std::reverse(s.begin() + i, s.begin() + j);
}

void reinsert(struct neighbor_info *cheapest = NULL, int i = 0, int j = 0, int pos = 0, std::vector<int> &s = DEFAULT_VECTOR){
	int size = j - i + 1;
	if(pos < i){
		if(i == j){
			cheapest->neighbor.insert(cheapest->neighbor.begin() + pos, cheapest->neighbor[i]);
			cheapest->neighbor.erase(cheapest->neighbor.begin() + i + size);
		}else{
			if(cheapest != NULL){
				//structure
				for(int k = i, m = 0; k <= j; k++, m++){
					cheapest->neighbor.insert(cheapest->neighbor.begin() + pos + m, cheapest->neighbor[k]);
					cheapest->neighbor.erase(cheapest->neighbor.begin() + k+1);
					/*
					for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
						std::cout << " " << cheapest->neighbor[i];

					std::cout << std::endl;
					*/
				}
			}else if(!s.empty()){
				//vector
				for(int k = i, m = 0; k <= j; k++, m++){
					s.insert(s.begin() + pos + m, s[k]);
					s.erase(s.begin() + k+1);
					/*
					for(unsigned i = 0; i < s.size(); i++)
						std::cout << " " << s[i];

					std::cout << std::endl;
					*/
				}
			}
		}
	}else{
		if(i == j){
			cheapest->neighbor.insert(cheapest->neighbor.begin() + pos, cheapest->neighbor[i]);
			cheapest->neighbor.erase(cheapest->neighbor.begin() + i);
		}else{
			if(cheapest != NULL){
				for(int k = i, m = 0; k <= j; k++, m++){
					cheapest->neighbor.insert(cheapest->neighbor.begin() + pos, cheapest->neighbor[i]);
					cheapest->neighbor.erase(cheapest->neighbor.begin() + i);

					/*
					std::cout << "Posterior" << std::endl;
					for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
						std::cout << " " << cheapest->neighbor[i];

					std::cout << std::endl;
					*/
				}
			}else if(!s.empty()){
				for(int k = i, m = 0; k <= j; k++, m++){
					s.insert(s.begin() + pos, s[i]);
					s.erase(s.begin() + i);

					//std::cout << "Posterior" << std::endl;
					/*
					for(unsigned i = 0; i < s.size(); i++)
						std::cout << " " << s[i];

					std::cout << std::endl;
					*/
				}

			}
		}
	}

}

void neighbor_swap_better(struct neighbor_info *cheapest, std::vector<int> &s){
	double dif;
	double dif_lower;
	printf("SWAP\n");
	for(int i = 1; i < s.size() - 2; i++){
		for(int j = i + 1; j < s.size() - 1; j++){
			dif = 0;
			// subtract the old distances
			dif += - c[s[i]][s[i-1]] - c[s[i]][s[i+1]] - c[s[j]][s[j-1]] - c[s[j]][s[j+1]];
			// add the new distances
			dif += c[s[i]][s[j-1]] + c[s[i]][s[j+1]] + c[s[j]][s[i-1]] + c[s[j]][s[i+1]];

			if(i + 1 == j) //consecutive
				dif += 2 * c[s[i]-1][s[j]-1];

			//std::cout << " " << dif;

			if(dif < dif_lower || (i == 1 && j == 2)){
				dif_lower = dif;
				cheapest->i = i;
				cheapest->j = j;
				cheapest->cost_dif = dif_lower;

			}
		}
	}
	//printf("\n");

	swap(cheapest->neighbor, cheapest->i, cheapest->j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_two_opt_better(struct neighbor_info *cheapest, std::vector<int> &s){
	srand(time(NULL));
	int dif;
	int dif_lower;
	int j;
	int sub_sz = rand() % (s.size() - 2);
	while(sub_sz == 1 || !sub_sz)
		sub_sz = rand() % (s.size() - 2);

	
	printf("TWO OPT\n");
	/*
	printf("subsequence size =  %d\n", sub_sz);
	*/

	for(int i = 1; i < s.size() - sub_sz; i++){
		dif = 0;
		j = i + sub_sz;
		// subtract the old distances
		dif += - c[s[i]][s[i-1]]  - c[s[j]][s[j+1]];
		// add the new distances
		dif +=  c[s[i]][s[j+1]] + c[s[j]][s[i-1]];

		//std::cout << " " << dif;

		if(dif < dif_lower || (i == 1) ){
		  dif_lower = dif;
		  cheapest->i = i;
		  cheapest->j = j;
		  cheapest->cost_dif = dif_lower;

		}
	}
	//printf("\n");

	two_opt(cheapest, cheapest->i, cheapest->j);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

void neighbor_reinsertion_better(struct neighbor_info *cheapest, std::vector<int> &s, int size){
	double dif;
	double dif_lower;
	int j;
	int pos_new;
	printf("REINSERTION\n");
	for(int i = 1; i < s.size() - size; i++){
		//k : arestas 
		for(int k = 0; k < s.size() - 2; k++){
			dif = 0;
			j = i + size - 1;

			if(k < i - 1){
				pos_new = k + 1;
			}else if(k >= i - 1 && k <= j){
				continue;
			}else if(k > j){
				pos_new = k + 1;
			}		
				
			// subtract the old distances
			dif += - c[s[i]][s[i-1]] - c[s[j]][s[j+1]] - c[s[pos_new]][s[pos_new-1]];
			// add the new distances
			dif += c[s[i]][s[pos_new-1]] + c[s[pos_new]][s[j]] + c[s[i-1]][s[j+1]];

			//std::cout << " " << dif;

			if(dif < dif_lower || (i == 1)){
				dif_lower = dif;
				cheapest->i = i;
				cheapest->j = j;
				cheapest->pos_new = pos_new;
				cheapest->cost_dif = dif_lower;

			}
		}
		//printf("\n");
	}
	//printf("\n");

	reinsert(cheapest, cheapest->i, cheapest->j, cheapest->pos_new);
	/*
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
	*/
}

double cost_calc(std::vector<int> &s){
	double cost_best = 0;
	for(unsigned i = 0; i < s.size() - 1; i++)
		cost_best += c[s[i]][s[i+1]];

	return cost_best;
}

void neighbd_list_repopulate(std::vector<int> &list){
	list.clear();
	list = {1, 2, 3, 4, 5};
}
		
#define REINSERTION 1
#define OR_OPT2 	2
#define OR_OPT3 	3
#define SWAP 		4
#define TWO_OPT		5

long l=0;

void RVND(std::vector<int> &s){

	srand(time(NULL));
	std::vector<int> neighbd_list = {1, 2, 3, 4, 5};
	l++;
	while(!neighbd_list.empty()){
		printf(" %ld vezes\n", l);
		
		int neighbd_rand_index = rand() % neighbd_list.size();
		int neighbd_rand = neighbd_list[neighbd_rand_index];


		struct neighbor_info *cheapest  (new struct neighbor_info) ;
		//std::copy(s.begin(), s.end(), cheapest->neighbor.begin());
		cheapest->neighbor = s;

		
		switch(neighbd_rand){
			case SWAP:
				neighbor_swap_better(cheapest, s);
				break;
			case TWO_OPT:
				neighbor_two_opt_better(cheapest, s);
				break;				
			case REINSERTION:
				neighbor_reinsertion_better(cheapest, s, REINSERTION);
				break;				
			case OR_OPT2:
				neighbor_reinsertion_better(cheapest, s, OR_OPT2);
				break;				
			case OR_OPT3:
				neighbor_reinsertion_better(cheapest, s, OR_OPT3);
				break;				
		}

		if((int)cheapest->cost_dif < 0){
			s.clear();
			s = cheapest->neighbor;
			cost_rvnd_best = cost_calc(s);
			neighbd_list_repopulate(neighbd_list);
		}else
			neighbd_list.erase(neighbd_list.begin() + neighbd_rand_index);

		delete(cheapest);
	}
}

void perturb(std::vector<int> &s){
	srand(time(NULL));
	int dimension = s.size();
	int first_size = (rand() % ((dimension / 10) - 2)) + 2;
	int second_size = (rand() % ((dimension / 10) -2)) + 2;

	int first_pos = rand() % (dimension - first_size - 1) + 1;
	int second_pos;
	int flag = 1;
	// second_pos restrictions
	while(flag){
		int pos = rand() % (dimension - second_size);
		if(pos && (pos <( first_pos - second_size) || pos > (first_pos + first_size))){
			second_pos = pos;
			flag = 0;
		}
	}

	flag = 1;
	std::vector<int> list (dimension);
	list[0] = 0;

	for(int i = first_pos; i < first_pos + first_size; i++)
		list[i] = 1;
	
	for(int i = second_pos; i < second_pos + second_size; i++)
		list[i] = 2;

	for(unsigned i = 0; i <  list.size(); i++)
		std::cout << " " << list[i];

	std::cout << std::endl;

	int first_pos_next;
	while(flag){
		int pos = rand() % (dimension - first_size + 1) + 1;
		if(pos == first_pos)
			continue;
		if(!(list[pos] == 2) != !(list[pos-1] == 2)){
			first_pos_next = pos;
			flag = 0;
		}
	}

	reinsert(NULL, first_pos, first_pos + first_size - 1, first_pos_next, list);


	for(unsigned i = 0; i <  list.size(); i++)
		std::cout << " " << list[i];

	std::cout << std::endl;

	flag = 1;	
	int second_pos_next;
	while(flag){
		int pos = rand() % (dimension - second_size + 1) + 1;
		if(pos == second_pos)
			continue;
		if(!(list[pos] == 1) != !(list[pos-1] == 1)){
			second_pos_next = pos;
			flag = 0;
		}
	}

	if(first_pos > second_pos && first_pos_next <= second_pos){
		second_pos += first_size;
	}else if(first_pos < second_pos && first_pos_next > second_pos){
		second_pos -= first_size;
	}
	
	reinsert(NULL, second_pos, second_pos + second_size - 1, second_pos_next, list);

	for(unsigned i = 0; i <  list.size(); i++)
		std::cout << " " << list[i];

	std::cout << std::endl;
	printf("1a pos = %d\n"
			"1o tamanho = %d\n"
			"2a pos = %d\n"
			"2o tamanho = %d\n"
			"1a next = %d\n"
			"2a next = %d\n", first_pos, first_size, second_pos, second_size, first_pos_next, second_pos_next);
	
	/*
	int nodes = s.size() - 1; // subtract the last
	int sub_seq_sz = (int) nodes / 4;
	int flag;
	int count = 1 + 3*sub_seq_sz;

	if(nodes % 2 != 0)
		flag = 0;
	else
		flag = 1;
	
	reinsert(NULL, 1, sub_seq_sz - flag, s.size()-1, s);

	for(int i = 1; i < s.size() - 1; i += sub_seq_sz){
	  if(i == count)
		  two_opt(NULL, i, i + sub_seq_sz - flag, s);
	  else
		  two_opt(NULL, i, i + sub_seq_sz, s);
	}
	printf("perturb\n");
*/	
	
}

void GILS_RVND(int Imax, int Iils){
	for(int i = 0; i < Imax; i++){
		srand(time(NULL));
		int aux = rand() % 10;
		double alpha =(double)( 1.0 / aux);
		//double alpha = 0.5;
		printf("alpha  = %lf\n", alpha);

		std::vector<int> s = construcao(alpha);
		std::vector<int> sl = s;
		int Iterils = 0;

		while(Iterils < Iils){
			RVND(s);
			if((int)cost_rvnd_best <(int) cost_calc(sl)){
				sl = s;
				Iterils = 0;
			}
			perturb(sl);
			s.clear();
			s = sl;
			Iterils++;
			printf("iter ils %d\n", Iterils);
		
		}

		double cost_sl = cost_calc(sl);
		if(i == 1 || cost_sl < cost_final){
			s_final =  s;
			cost_final = cost_sl;
		}
	}
}


int main(int argc, char **argv){
	int Imax = 50;
	int Iils;

	srand(time(NULL));
	std::vector<int> sa = {1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20};
	
	  perturb(sa);
	/*
	readData(argc, argv, &dimension, &c);

	if(dimension >= 150)
		Iils = dimension / 2;
	else
		Iils = dimension;

	GILS_RVND(Imax, Iils);
	printf("better cost = %lf\n", cost_final);
	return 0;
	*/
}

