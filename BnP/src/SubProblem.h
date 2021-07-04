#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>

using namespace std;

class SubProblem {
    public:
        SubProblem(double*, double, int);
        ~SubProblem();
        //vector<int> solve(vector<int>);
        vector<int> solve(IloNumArray*);
    private:
        double *        weights;
        int             dimension;
        double          C;

        IloEnv          env;
        IloModel model;
        IloExpr LHS;
        IloBoolVarArray x ;
};

#endif
