#include "SubProblem.h"

SubProblem::SubProblem(double * w, double binC, int n): env(), model(env), LHS(env), x(env, n), column(n){

    //this->weights = w;
    weights = (double*)calloc(n, sizeof(double));
    memcpy(weights, w, n*sizeof(double));
    this->dimension = n;
    this->C = binC;

    // add variables
    for(int i = 0; i < this->dimension; i++){
        char name[100];
        sprintf(name, "x_%d", i);
        x[i].setName(name);
        model.add(x[i]);
    }

    // add cstrs
    IloRange cstr;
    for(int i = 0; i < this->dimension; i++){
        LHS += this->weights[i] * x[i];
    }
    cstr = (LHS <= this->C);
    model.add(cstr);
    
}

SubProblem::~SubProblem(){
    free(weights);
    env.end();
}

//vector<int> SubProblem::solve(vector<int> duals) {
double SubProblem::solve(IloNumArray * duals){


    // add OF
    IloExpr obj(env);

    for(int i = 0; i < this->dimension; ++i){
        //obj += duals[i] * x[i];
        obj += (*duals)[i] * x[i];
    }

    IloObjective f = IloMaximize(this->env, obj);
    this->model.add(f);

    IloCplex KP(model);
    KP.setParam(IloCplex::Param::TimeLimit, 60);
    KP.setParam(IloCplex::Param::Threads, 1);


    double after, before;

    try{ 
        before = KP.getTime();
        KP.solve();
        after = KP.getTime();
    }
    catch(IloException& e){
        std::cout << e;
    }

    for(int i = 0; i < this->dimension; i++){
        column[i] = KP.getValue(x[i]);
        cout << column[i] << " ";
    }
    cout << endl;

    printResults(KP, "sacola", after-before);
    double obj_value = -(KP.getObjValue() - 1);

    model.remove(f);
    f.end();

    return obj_value;
}

vector<int> SubProblem::getColumn(){
    return column;
}
