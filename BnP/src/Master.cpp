#include "Master.h"

Master::Master(int size) : env(), model(env), duals(env), objF(), cstr(){
    this->initial_size = size;
    this->col_qnt = size;

    this->env.setName("Column Generation");
    this->model.setName("Bin Packing Problem");
    
    build();
}

Master::~Master(){
    this->env.end();
    //delete duals;
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

        this->var_vec.push_back(var);
        this->var_standby_flag.push_back(0);
        vector<bool> col(initial_size, 0);
        col[i] = 1;
        this->A.push_back(col);
    }

    /*
    for(int i = 0; i < initial_size; i++){
        for(int j = 0; j < initial_size; j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    exit(1);
    */

}

double Master::solve(){

    IloCplex BPP(this->model);
    BPP.setParam(IloCplex::Param::TimeLimit, 1*60*60);
    BPP.setParam(IloCplex::Param::Threads, 1);
    BPP.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-08);
    BPP.setOut(env.getNullStream());

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

    //printResults(BPP, "instancia", before-after);
    getSolution(BPP);
    this->duals.clear();
    BPP.getDuals(this->duals, this->cstr);

    double obj = BPP.getObjValue();
    BPP.end();
    return obj;
}

void Master::addColumn(vector<bool> col_vec){
    IloNumColumn col = this->objF(1);
    
    for(int i = 0; i < initial_size; i++){
        /*
        if(col_vec[i] > 0.0005){
            col += this->cstr[i](1);
        }
        */

        col += this->cstr[i](col_vec[i]);
    }

    char var_name[100];
    sprintf(var_name, "b_%u", this->col_qnt++);

    //IloNumVar var_new (col, 0, 1);
    IloNumVar var_new (col, 0, IloInfinity);
    var_new.setName(var_name);

    this->model.add(var_new);

    this->var_vec.push_back(var_new);
    this->A.push_back(col_vec);
    this->var_standby_flag.push_back(0);

    col.end();
}

void Master::getSolution(IloCplex & cplex){
    IloExpr vars = objF.getExpr();

    this->solution.clear();
    this->solution_values.clear();
    for(IloExpr::LinearIterator it = vars.getLinearIterator(); it.ok(); ++it){
        double var_value = cplex.getValue(it.getVar());   
        if(var_value > 0.00000005){
            solution.push_back(it.getVar());
            solution_values.push_back(var_value);
            //cout << /*it.getVar() << " " <<*/ var_value << " ";
        }
    }
    //cout << endl;
    //exit(1);
    //
    vars.end();

}

IloNumArray Master::getDuals(){
    return duals;
}

vector<bool> Master::getColumn(IloNumVar & var){
    vector<bool> column;

    for(int i = 0; i < initial_size; i++){
        bool coef = 0;
        for (IloExpr::LinearIterator it = IloExpr(this->cstr[i].getExpr()).getLinearIterator(); it.ok();++it)
            if(it.getVar().getName() == var.getName()){
                coef = 1;
                break;
            }

        column.push_back(coef);
    }

    return column;
}

vector<double> Master::getSolutionValues(){
    return solution_values;
}

vector<IloNumVar> Master::getSolution(){
    return solution;
}

void Master::enforce(pair<int, int> m_pair){
    int a = m_pair.first; 
    int b = m_pair.second; 

    cout << "A size " << A.size() << endl;

    for(int j = initial_size; j < A.size(); j++){
        if(this->var_standby_flag[j] == 0 && (A[j][a] == 0 && A[j][b] == 0)){
            //removeVar(var_vec[j]);
            cout << "aqui 1" << endl;
            //this->model.remove(var_vec[j]);
            this->var_standby_flag[j] = 1;
            cout << var_vec[j] << endl;    
            var_vec[j].end();
            var_standby.push_back(make_pair(m_pair, j));
        }
    }
}

/*
void Master::removeVar(IloNumVar var){
    this->model.remove(var);
    var_standbydk
}
*/

void Master::exclude(pair<int, int> m_pair){
    int a = m_pair.first; 
    int b = m_pair.second; 

    cout << "A size " << A.size() << endl;
    for(int j = initial_size; j < A.size(); j++){
        if(this->var_standby_flag[j] == 0 && A[j][a] == 1 && A[j][b] == 1){
            //removeVar(var_vec[j]);
            cout << "aqui 1" << endl;
            //this->model.remove(var_vec[j]);
            this->var_standby_flag[j] = 1;
            cout << var_vec[j] << endl;    
            var_vec[j].end();
            cout << "aqui 2" << endl;
            var_standby.push_back(make_pair(m_pair, j));
        }
    }
    
}

inline bool operator==(const pair<int, int> & A, const pair<int, int> & B){ 
    return (A.first == B.first) && (A.second == B.second);
}

void Master::reinsert(pair<int, int> m_pair){

    for(int i = 0; i < var_standby.size(); i++){
        if(m_pair == var_standby[i].first){
            int var_index = var_standby[i].second;
            IloNumColumn col = this->objF(1);

            for(int j = 0; j < initial_size; j++){
                col += this->cstr[i](this->A[var_index][i]);
            }

            char var_name[100];
            sprintf(var_name, "b_%u", var_index);

            IloNumVar var_new (col, 0, IloInfinity);
            var_new.setName(var_name);
            cout << var_name << endl;

            this->model.add(var_new);

            var_vec[var_index] = var_new;
            var_standby_flag[var_index] = 0;
        }
    }
}

void Master::printResult(){
    vector<vector<int>> bins;
    for(auto & var: solution){
        vector<bool> column;// = getColumn(var, cstr, initial_size);
        //if(
        column = getColumn(var);

        vector<int> bin;
        for(int i = 0; i < column.size(); i++){
            if(column[i] == 1){
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
