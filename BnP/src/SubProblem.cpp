#include "SubProblem.h"
#include "minknap.c"

SubProblem::SubProblem(int * w, int binC, int n, vector<pair<int, int>> exclude, vector<pair<int, int>> enforce): env(), model(env), LHS(env), x(env, n), column(n){

    //this->weights = w;
    weights = (int*)calloc(n, sizeof(int));
    memcpy(weights, w, n*sizeof(int));
    this->dimension = n;
    this->C = binC;

    this->cstr_exclude = exclude;
    this->cstr_enforce = enforce;
    if(exclude.empty() && enforce.empty()){
        this->mode = PISINGER;
        //this->mode = BY_MIP;
        //MIP_build();
    }else{
        this->mode = BY_MIP;
        MIP_build();
    }

}

SubProblem::~SubProblem(){
    free(weights);
    cstr_enforce.clear();
    cstr_exclude.clear();

    if(this->mode == BY_MIP){
        this->model.end();
        this->LHS.end();
        this->x.end();
        this->env.end();
    }
}

double SubProblem::solve(IloNumArray duals){
    if(this->mode == BY_MIP){
        return 1.0 - MIP_solve(duals);
    }else if(this->mode == PISINGER){
        return 1.0 - PISINGER_solve(duals);
    }
}

double SubProblem::PISINGER_solve(IloNumArray duals){
    int p[this->dimension] = {};
    int x[this->dimension];

    for(int i = 0; i < this->dimension; i++){
        p[i] = duals[i] * BIG_FAC;
    }

    double obj_value = minknap(this->dimension, p, this->weights, x, this->C);
    for(int i = 0; i < this->dimension; i++){
        column[i] = (bool)x[i];
    }

    obj_value *= BIG_FAC_inv;

    return obj_value;
}

void SubProblem::MIP_build(){
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

    cstr = (LHS <= (this->C ));
    model.add(cstr);
    
    // branching cstrs
    IloExpr exp_a (env);
    IloExpr exp_b (env);
    for(int i = 0; i < cstr_exclude.size(); i++){
        IloRange cstr_branch_a;
        int item_a = cstr_exclude[i].first;
        int item_b = cstr_exclude[i].second;
        exp_a += x[item_a] + x[item_b];
        cstr_branch_a = (exp_a <= 1);
        model.add(cstr_branch_a);
        exp_a.clear();
    }

    for(int i = 0; i < cstr_enforce.size(); i++){
        IloRange cstr_branch_a;
        int item_a = cstr_enforce[i].first;
        int item_b = cstr_enforce[i].second;
        exp_a += x[item_a] - x[item_b];
        cstr_branch_a = (exp_a == 0);
        model.add(cstr_branch_a);
        exp_a.clear();
    }

    exp_a.end();
    exp_b.end();
}

double SubProblem::MIP_solve(IloNumArray duals){
    // add OF
    IloExpr obj(env);

    for(int i = 0; i < this->dimension; ++i){
        obj += duals[i] * x[i];
    }

    IloObjective f = IloMaximize(this->env, obj);
    this->model.add(f);

    IloCplex KP(model);
    KP.setParam(IloCplex::Param::TimeLimit, 60);
    KP.setParam(IloCplex::Param::Threads, 1);
    KP.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-08);
    KP.setOut(env.getNullStream());
    //KP.exportModel("kp.lp");

    try{ 
        KP.solve();
    }
    catch(const IloException& e){
        cerr << e;
    }

    double value;
    bool arr_branch[] = {0, 1};
    for(int i = 0; i < this->dimension; i++){
        value = KP.getValue(x[i]);
        column[i] = arr_branch[value >= 0.9];
    }

    double obj_value = KP.getObjValue();
    KP.end();
    obj.end();
    f.end();

    return obj_value;
}

vector<bool> SubProblem::getColumn(){
    return column;
}
