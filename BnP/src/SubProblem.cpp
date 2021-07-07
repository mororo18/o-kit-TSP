#include "SubProblem.h"

SubProblem::SubProblem(int * w, int binC, int n): env(), model(env), LHS(env), x(env, n), column(n){

    //this->weights = w;
    weights = (int*)calloc(n, sizeof(int));
    memcpy(weights, w, n*sizeof(int));
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
        //cout << (*duals)[i] << " ";
    }
    //cout << endl;

    IloObjective f = IloMaximize(this->env, obj);
    this->model.add(f);

    IloCplex KP(model);
    KP.setParam(IloCplex::Param::TimeLimit, 60);
    KP.setParam(IloCplex::Param::Threads, 1);
    KP.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-08);
    KP.setOut(env.getNullStream());
    //KP.exportModel("kp.lp");

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
        double value = KP.getValue(x[i]);
        if(value >= 0.9)
            column[i] = 1;
        else
            column[i] = 0;

        //cout << column[i] << " ";
    }
    //cout << endl;

    //printResults(KP, "sacola", after-before);
    double obj_value = 1.0 - KP.getObjValue();

    //K.exportModel("kp.lp");
    obj.end();
    f.end();

    return obj_value;
}

vector<int> SubProblem::getColumn(){
    return column;
}
