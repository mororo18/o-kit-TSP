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
        virtual ~Master();
        void addColumn(vector<bool>);
        //vector<int> getDuals();
        IloNumArray getDuals();
        double solve();
        vector<vector<int>> getBins();
        vector<bool> getColumn(IloNumVar &);
        vector<IloNumVar> getSolution();
        vector<double> getSolutionValues();
        void enforce(pair<int, int>);
        void exclude(pair<int, int>);
        void reinsert(pair<int, int>);
        bool getStatus();
        
    private:
        IloEnv              env;
        IloModel            model;
        IloObjective        objF;
        IloRangeArray       cstr;
        IloCplex            BPP;
        vector<IloNumVar>   var_vec;
        vector<IloNumVar>   solution;
        vector<vector<bool>> A;
        vector<bool> var_standby_flag;
        vector<pair<pair<int, int>, int>> var_standby;   
        vector<double>      solution_values;
        IloNumArray         duals;
        int                 initial_size;
        int                 col_qnt;
        bool                feasibility;

        void build();
        void getSolution(IloCplex &);
        //void removeVar(int);
        vector<int> d_vec;
        
};
#endif 
