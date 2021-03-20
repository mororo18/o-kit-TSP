#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>

using namespace std;

typedef pair<int, int> ii;
typedef vector <vector<double> > Matrix;
typedef vector<ii> vii;

class Kruskal{
    public:
        Kruskal(Matrix dist);

        double MST(int nodes);
        vii getEdges();


    private:
        priority_queue <pair<double,ii> > graph;
        vector <pair<double, ii> >costs;
        vector <int> pset;
        vii edges;

        void initDisjoint(int n);
        int findSet(int i);
        void unionSet(int i, int j);
        bool isSameSet(int i, int j);
        bool compare(const std::pair<double,ii> &left, const std::pair<double,ii> &right);
};

#endif
