#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include "auxFunctions.h"

using namespace std;

class SubProblem {
    public:
        SubProblem(double*, double, int);
        ~SubProblem();
        //vector<int> solve(vector<int>);
        double solve(IloNumArray*);
        vector<int> getColumn();
    private:
        double *        weights;
        int             dimension;
        double          C;

        IloEnv          env;
        IloModel model;
        IloExpr LHS;
        IloBoolVarArray x ;

        vector<int> column;
};

#endif
