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
        double solve();
        void printResult();
        vector<int> getColumn(IloNumVar &);
        vector<IloNumVar> getSolution();
        vector<double> getSolutionValues();
        
    private:
        IloEnv          env;
        IloModel        model;
        IloObjective    objF;
        IloRangeArray   cstr;
        vector<IloNumVar>  solution;
        vector<double>  solution_values;
        IloNumArray *   duals;
        int             initial_size;
        int             col_qnt;


        void build();
        void getSolution(IloCplex);
        vector<int> d_vec;
        
};
#endif 
