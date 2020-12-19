#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "readData.h"


int c[6][6] = {
	{0, 60, 50, 40, 30, 20},
	{60, 0, 55, 45, 35, 25},
	{50, 55, 0, 50, 40, 30},
	{40, 45, 50, 0, 45, 35},
	{30, 35, 40, 45, 0, 40},
	{20, 25, 30, 35, 40, 0}
};

struct insertion_info {
  int node_new;
  int edge_removed;
  double cost;
};

struct neighbor_info{
  std::vector<int> neighbor;
  int i;
  int j;
  int cost_dif;
};

double **matrix;
int dimension;

bool cost_compare(const insertion_info &a, const insertion_info &b){
	return a.cost < b.cost;
}

std::vector<int> construcao(float alpha){
	
  	srand(time(NULL));

	std::vector<int> s = { 1, 1};
	std::vector<int> candidatos = { 2, 3, 4, 5, 6};
	int subtour_inicial = 3;

	for(int i = 0; i < subtour_inicial; i++){
		int node_rand_index = rand() % candidatos.size();
		int node_rand = candidatos[node_rand_index];

		s.insert(s.begin() + 1, node_rand);
		candidatos.erase(candidatos.begin() + node_rand_index);
	}

	for(unsigned i = 0; i < s.size(); i++)
		std::cout << " " << s[i];

	std::cout << std::endl;
	
	while(candidatos.size() != 0){


		int insertion_possiblts = ((s.size() - 1) * candidatos.size());
		printf("possibilidades  %d\n", insertion_possiblts);
		std::vector<struct insertion_info> insertion_cost (insertion_possiblts);

		int l = 0; // 'insertion_cost' index

		for(auto k : candidatos){
			for(int i = 0, j = 1; i < s.size() - 1; i++, j++){  //modificacao ('i < s.size()' para 'i <= s.size()') 
				
				insertion_cost[l].cost = c[s[i]-1][k-1] + c[s[j]-1][k-1] - c[s[i]-1][s[j]-1];  //modificacao temporaria
				insertion_cost[l].node_new = k;
				insertion_cost[l].edge_removed = i; 
				l++;

			}
		}
		printf("valor de l  =  %d\n", l);

		std::sort (insertion_cost.begin(), insertion_cost.end(), cost_compare);
		
		int node_rand_range = (int) (alpha * insertion_cost.size());
		int node_rand_index = rand() % node_rand_range;

		int node_rand = insertion_cost[node_rand_index].node_new;
		int edge = insertion_cost[node_rand_index].edge_removed;
		s.insert(s.begin() + edge + 1, node_rand); 
		printf("inserted node: %d\n", node_rand);

		for(int i = 0; i < candidatos.size(); i++){
			if(candidatos[i] ==  node_rand){
				candidatos.erase(candidatos.begin() + i);
				break;
			}
		};

		for(unsigned i = 0; i < s.size(); i++)
		  std::cout << " " << s[i];

		std::cout << std::endl;
	}


	return s;
}	

void swap(std::vector<int> &s, int i, int j){
	int aux;

	aux = s[i];
	s[i] = s[j];
	s[j] = aux;
}

void two_opt(struct neighbor_info *cheapest, int i, int j){
	std::reverse(cheapest->neighbor.begin() + i, cheapest->neighbor.begin() + j);
}

void neighbor_swap_better(struct neighbor_info *cheapest, std::vector<int> &s){
	int dif = 0;
	int dif_lower;
	for(int i = 1; i < s.size() - 2; i++){
		for(int j = i + 1; j < s.size() - 1; j++){
			// subtract the old distances
			dif += - c[s[i]-1][s[i-1]-1] - c[s[i]-1][s[i+1]-1] - c[s[j]-1][s[j-1]-1] - c[s[j]-1][s[j+1]-1];
			// add the new distances
			dif += c[s[i]-1][s[j-1]-1] + c[s[i]-1][s[j+1]-1] + c[s[j]-1][s[i-1]-1] + c[s[j]-1][s[i+1]-1];

			if(i + 1 == j) //consecutive
				dif += 2 * c[s[i]-1][s[j]-1];

			std::cout << " " << dif;

			if(dif < dif_lower || (i == 1 && j == 2)){
				dif_lower = dif;
				cheapest->i = i;
				cheapest->j = j;
				cheapest->cost_dif = dif_lower;

			}
		}
		dif = 0;
	}
	printf("\n");

	swap(cheapest->neighbor, cheapest->i, cheapest->j);
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
}

void neighbor_two_opt_better(struct neighbor_info *cheapest, std::vector<int> &s){
	int dif = 0;
	int dif_lower;
	int j;
	int sub_sz = rand() % (s.size() - 2);
	while(sub_sz == 1 || !sub_sz)
		sub_sz = rand() % (s.size() - 2);

	printf("subsequence size =  %d\n", sub_sz);

	for(int i = 1; i < s.size() - sub_sz; i++){
		j = i + sub_sz - 1;
		// subtract the old distances
		dif += - c[s[i]-1][s[i-1]-1]  - c[s[j]-1][s[j+1]-1];
		// add the new distances
		dif +=  c[s[i]-1][s[j+1]-1] + c[s[j]-1][s[i-1]-1];

		std::cout << " " << dif;

		if(dif < dif_lower || (i == 1) ){
		  dif_lower = dif;
		  cheapest->i = i;
		  cheapest->j = j;
		  cheapest->cost_dif = dif_lower;

		}
		dif = 0;
	}
	printf("\n");

	two_opt(cheapest, cheapest->i, cheapest->j);
	for(unsigned i = 0; i < cheapest->neighbor.size(); i++)
		std::cout << " " << cheapest->neighbor[i];

	std::cout << std::endl;
	std::cout << "cost-dif =  " << cheapest->cost_dif << std::endl;
}

#define SWAP 		1
#define TWO_OPT		2
#define REINSERTION 3
#define OR_OPT2 	4
#define OR_OPT3 	5

void RVND(std::vector<int> &s){

	std::vector<int> neighbd_list = {2/*, 2, 3, 4, 5*/};
	/* 
	list of 5 neighborhood structures

	1 : swap
	2 : 2-opt
	3 : reinsertion
	4 : or-opt2
	5 : or-opt3
	*/
	int n =1;
	while(n--/*!neighbd_list.empty()*/){
		
		int neighbd_rand_index = rand() % neighbd_list.size();
		int neighbd_rand = neighbd_list[neighbd_rand_index];


		struct neighbor_info *cheapest = (struct neighbor_info *) malloc(sizeof(struct neighbor_info));
		//std::copy(s.begin(), s.end(), cheapest->neighbor.begin());
		cheapest->neighbor = s;

		
		switch(neighbd_rand){
			case SWAP:
				neighbor_swap_better(cheapest, s);
				break;
			case TWO_OPT:
				neighbor_two_opt_better(cheapest, s);
				break;				
		}
		free(cheapest);
		

	}
}


void GILS_RVND(int Imax, int Iils){
	for(int i = 0; i < Imax; i++){
		float alpha = 1 / rand();

		std::vector<int> s = construcao(alpha);
		RVND(s);

	}
}


int main(){

	std::vector<int> s = construcao(0.5);
	RVND(s);
	//readData(argc, argv, &dimension, &matrix);
	return 0;
}

