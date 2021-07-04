#ifndef MASTER_H
#define MASTER_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include "auxFunctions.h"

using namespace std;

class Master {
    public:
        Master(int);
        ~Master();
        void addColumn(vector<int>);
        //vector<int> getDuals();
        IloNumArray* getDuals();
        void solve();
        void printResult();
        
    private:
        IloEnv          env;
        IloModel        model;
        IloObjective    objF;
        IloRangeArray   cstr;
        vector<IloNumVar>  solution;
        IloNumArray *   duals;
        int             initial_size;
        int             col_qnt;


        void build();
        void getColumn();
        void getSolution(IloCplex);
        vector<int> d_vec;
        vector<int> getColumn(IloNumVar, IloRangeArray, int);
        
};
#endif 
