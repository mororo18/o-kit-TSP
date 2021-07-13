#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <cfloat>
#include "auxFunctions.h"

using namespace std;

class SubProblem {
    public:
        SubProblem(int*, int, int, vector<pair<int, int>>, vector<pair<int, int>>);
        ~SubProblem();
        //vector<int> solve(vector<int>);
        double solve(IloNumArray);
        vector<bool> getColumn();
    private:
        int *   weights;
        int     dimension;
        int     C;

        IloEnv          env;
        IloModel model;
        IloExpr LHS;
        IloBoolVarArray x ;

        vector<bool> column;

        vector<pair<int, int>> cstr_enforce;
        vector<pair<int, int>> cstr_exclude;
};

#endif
