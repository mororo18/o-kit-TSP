#include "Master.h"

Master::Master(int size) : env(), model(env), duals(new IloNumArray(env)), objF(), cstr(){
    this->initial_size = size;
    this->col_qnt = size;

    this->env.setName("Column Generation");
    this->model.setName("Bin Packing Problem");
    
    build();
}

Master::~Master(){
    this->env.end();
    delete duals;
}

void Master::build(){
    this->objF = IloAdd(this->model, IloMinimize(this->env));
    this->cstr = IloAdd(this->model, IloRangeArray(this->env, initial_size, 1, 1));

    for(int i = 0; i < this->initial_size; i++){
        IloNumVar var(this->objF(1) + this->cstr[i](1), 0, IloInfinity);
        char var_name[100];
        sprintf(var_name, "b_%u", i);

        var.setName(var_name);
        this->model.add(var);
    }

}

void Master::solve(){

    IloCplex BPP(this->model);
    BPP.setParam(IloCplex::Param::TimeLimit, 1*60*60);
    BPP.setParam(IloCplex::Param::Threads, 1);

    BPP.exportModel("model.lp");

    double after;
    double before;

    try{
        before = BPP.getTime();
        BPP.solve();
        after = BPP.getTime();
    }catch(IloException& e){
        cout << e;
    }

    printResults(BPP, "instancia", before-after);
    //getSolution(BPP);
    BPP.getDuals((*this->duals), this->cstr);

    for(int i = 0; i < this->initial_size; i++){
        //cout << (*this->duals)[i] << " ";
    }
    cout << endl;
}

void Master::addColumn(vector<int> col_vec){
    IloNumColumn col = this->objF(1);
    
    for(int i = 0; i < initial_size; i++){
        if(col_vec[i] > 0.00000005){
            col += this->cstr[i](1);
        }
    }

    char var_name[100];
    sprintf(var_name, "b_%u", this->col_qnt++);

    IloNumVar var_new (col, 0, IloInfinity);
    var_new.setName(var_name);

    this->model.add(var_new);
}

void Master::getSolution(IloCplex cplex){
    IloExpr vars = objF.getExpr();

    this->solution.clear();
    for(IloExpr::LinearIterator it = vars.getLinearIterator(); it.ok(); ++it){
        double var_value = cplex.getValue(it.getVar());   
        if(var_value > 0.0005){
            solution.push_back(it.getVar());
            cout << /*it.getVar() << " " <<*/ var_value << " ";
        }
    }
    cout << endl;
    //exit(1);

}

IloNumArray * Master::getDuals(){
    return duals;
}


vector<int> Master::getColumn(IloNumVar var, IloRangeArray cstr, int n){
    vector<int> column;

    for(int i = 0; i < n; i++){
        int coef = 0;
        for (IloExpr::LinearIterator it = IloExpr(cstr[i].getExpr()).getLinearIterator(); it.ok();++it)
            if(it.getVar().getName() == var.getName()){
                //cout << it.getCoef() << endl;
                coef = 1;
                break;
            }

        column.push_back(coef);
        //cout << coef << endl;    
    }

    return column;
}

void Master::printResult(){
    vector<vector<int>> bins;
    for(auto & var: solution){
        vector<int> column;// = getColumn(var, cstr, initial_size);
        //if(
        column = getColumn(var, cstr, initial_size);

        vector<int> bin;
        for(int i = 0; i < column.size(); i++){
            if(column[i] >= 1 - 0.3){
                bin.push_back(i+1); 
                //cout << i+1 << " ";
            }
        }
        //cout << endl;
        //exit(1);

        bins.push_back(bin);
    }

    //cout << "Total bins: " << bins.size() << endl;
}
