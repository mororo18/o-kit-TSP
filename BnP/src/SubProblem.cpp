#include "SubProblem.h"

SubProblem::SubProblem(int * w, int binC, int n, vector<pair<int, int>> exclude, vector<pair<int, int>> enforce): env(), model(env), LHS(env), x(env, n), column(n){

    this->weights = w;
    weights = (int*)calloc(n, sizeof(int));
    memcpy(weights, w, n*sizeof(int));
    this->dimension = n;
    this->C = binC;

    this->cstr_exclude = exclude;
    this->cstr_enforce = enforce;

    // add variables
    for(int i = 0; i < this->dimension; i++){
        char name[100];
        sprintf(name, "x_%d", i);
        x[i].setName(name);
        model.add(x[i]);
    }

    /*
    for(int i = 0; i < cstr_exclude.size(); i++){
        int item_a = cstr_exclude[i].first;
        int item_b = cstr_exclude[i].second;

        this->weights[item_a] = this->C + 1;
        this->weights[item_b] = this->C + 1;
    }

    for(int i = 0; i < cstr_enforce.size(); i++){
        int item_a = cstr_enforce[i].first;
        int item_b = cstr_enforce[i].second;

        this->weights[item_a] = 0;
        this->weights[item_b] = 0;
    }
    */

    // add cstrs
    IloRange cstr;
    for(int i = 0; i < this->dimension; i++){
        LHS += this->weights[i] * x[i];
    }

    /*
    int sum = 0;
    for(int i = 0; i < cstr_enforce.size(); i++){
        int a = cstr_enforce[i].first;
        int b = cstr_enforce[i].second;
        sum += this->weights[a] + this->weights[b];
        cout << "A " << a << " e B " << b << endl;
    }
    */

    //cstr = (LHS <= (this->C - sum));
    cstr = (LHS <= (this->C ));
    model.add(cstr);
    
    // branching cstrs
        IloExpr exp_a (env);
        IloExpr exp_b (env);
    for(int i = 0; i < cstr_exclude.size(); i++){
        IloRange cstr_branch_a;
        IloRange cstr_branch_b;
        int item_a = cstr_exclude[i].first;
        int item_b = cstr_exclude[i].second;
        exp_a += x[item_a] + x[item_b];
        //exp_b += x[item_b];
        cstr_branch_a = (exp_a <= 1);
        //cstr_branch_b = (exp_b == 0);
        model.add(cstr_branch_a);
        //model.add(cstr_branch_b);
        exp_a.clear();
        exp_b.clear();
    }

    for(int i = 0; i < cstr_enforce.size(); i++){
        IloRange cstr_branch_a;
        IloRange cstr_branch_b;
        int item_a = cstr_enforce[i].first;
        int item_b = cstr_enforce[i].second;
        exp_a += x[item_a] + x[item_b];
        //exp_b += x[item_b];
        cstr_branch_a = (exp_a == 2);
        //cstr_branch_b = (exp_b == 1);
        model.add(cstr_branch_a);
        //model.add(cstr_branch_b);
        exp_a.clear();
        exp_b.clear();
    }
}

SubProblem::~SubProblem(){
    free(weights);
    this->env.end();
    cstr_enforce.clear();
    cstr_exclude.clear();
}

//vector<int> SubProblem::solve(vector<int> duals) {
double SubProblem::solve(IloNumArray duals){


    // add OF
    IloExpr obj(env);

    for(int i = 0; i < this->dimension; ++i){
        //obj += duals[i] * x[i];
        obj += duals[i] * x[i];
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
    catch(const IloException& e){
        cerr << e;
    }

    //cout << "opa" << endl;
    for(int i = 0; i < this->dimension; i++){
        double value = 0;
        try{ 
            value = KP.getValue(x[i]);
        }
        catch(const IloException& e){
            cerr << e;
        }

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

vector<bool> SubProblem::getColumn(){
    return column;
}
