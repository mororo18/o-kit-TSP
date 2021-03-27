#include "kruskal.h"

Kruskal::Kruskal(Matrix dist){
    for(int i = 0; i < dist.size(); i++){
        costs.push( make_pair(-dist[0][i], make_pair(0, i)) );
        costs.push( make_pair(-dist[i][0], make_pair(i, 0)) );
    }

    for(int i = 1; i < dist.size(); ++i){
        for(int j = 1; j < dist[i].size(); ++j){
            graph.push( make_pair(-dist[i][j], make_pair(i, j)) );
            //graph.push_back( make_pair(-dist[i][j], make_pair(i, j)) );
        }	
    }
}

void Kruskal::initDisjoint(int n){
    pset.resize(n);
    for (int i = 0; i < n; ++i){
        pset[i] = i;
    }
}

int Kruskal::findSet(int i){
    return (pset[i] == i) ? i : (pset[i] = findSet(pset[i]));
}

void Kruskal::unionSet(int i, int j){
    pset[findSet(i)] = findSet(j);
}

bool Kruskal::isSameSet(int i, int j){
    return (findSet(i) == findSet(j))? true:false;
}

vii Kruskal::getEdges(){
    return edges;
}

double Kruskal::MST(int nodes){
    initDisjoint(nodes);
    double cost = 0;

    /*
    std::sort(graph.begin(), graph.end(), 
        [](const std::pair<double,ii> &left, const std::pair<double,ii> &right) { 
            return left.first < right.first;
        }
    );
    */

    while(!graph.empty()){
        pair<double, ii> p = graph.top();
        //pair<double, ii> p = graph.back();

        //graph.pop_back();
        graph.pop();

        if(!isSameSet(p.second.first, p.second.second)){
            edges.push_back(make_pair(p.second.first, p.second.second));
            cost += (-p.first);
            unionSet(p.second.first, p.second.second);
        }
    }


    /*
    std::partial_sort(costs.begin() , costs.begin() + 2 ,costs.end(), 
        [](const std::pair<double,ii> &left, const std::pair<double,ii> &right) { 
            return -left.first < -right.first;
        }
    );
    for(int i = 0; i < costs.size(); i++){

        std::cout << costs[i].first << " ";
    }
    std::cout << std::endl;
    std::cout << costs[0].second.first << " ";
    std::cout << costs[0].second.second << " ";
    std::cout << costs[1].second.first << " ";
    std::cout << costs[1].second.second << " ";
    std::cout << std::endl;

    std::cout << edge_first.first << " " <<  edge_first.second;
    std::cout << std::endl;

    std::cout << edge_second.first << " " <<  edge_second.second;
    std::cout << std::endl;
    */



    pair<double, ii> edge_significant = costs.top();
    ii edge_first =  edge_significant.second;
    ii edge_second;//=  edge_significant.second;
    edges.push_back(edge_first);
    cost += -edge_significant.first;
    costs.pop();

    edge_second =  edge_significant.second;
    int sum_first = edge_first.first + edge_first.second;
    int sum_second = edge_second.first + edge_second.second;

    if(sum_first == sum_second){
        costs.pop();
        
        edge_significant = costs.top();
        edge_second =  edge_significant.second;
    }

    cost += -edge_significant.first;
    edges.push_back(edge_second);
    costs.pop();

    //cost += -costs[0].first -costs[1].first;

    return cost;
}
