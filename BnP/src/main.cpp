#include "auxFunctions.h"
#include "Data.h"
#include "Master.h"
#include <ilcplex/ilocplex.h>

int main(int argc, char * argv[]){
    Data data(argc, argv[1]);

    data.loadData();

    Master M(data.getItemQnt());
    M.solve();



    /*
    cout << data.getItemQnt() << endl;

    cout << data.getBinCapacity() << endl;

    for(int i=0; i < data.getItemQnt(); i++){
        cout << data.getWeight(i) << endl;
    }

    IloEnv env;
    IloModel model(env);

    env.setName("Column Generation");
    model.setName("Bin Packing Problem");

    IloObjective objF;
    IloRangeArray cstr;
    vector<IloNumVar> vars;
    vector<vector<int>> A;

    objF = IloAdd(model, IloMinimize(env));
    cstr = IloAdd(model, IloRangeArray(env, data.getItemQnt(), 1, 1));

    for(int i = 0; i < data.getItemQnt(); i++){
        IloNumVar var(objF(1) + cstr[i](1), 0, IloInfinity);
        char var_name[100];
        sprintf(var_name, "q_%u", i);

        var.setName(var_name);

        vars.push_back(model.add(var).asVariable());

        vector<int> col (data.getItemQnt(), 0);
        col[i] = 1;

        A.push_back(col);
    }

    IloNumColumn coluna = objF(1) + cstr[0](1) + cstr[47](1);

    IloNumVar var (coluna, 0, IloInfinity);
    char var_name[100];
    sprintf(var_name, "q_%u", 50);

    var.setName(var_name);

    model.add(var);

    IloCplex BPP(model);
    BPP.setParam(IloCplex::TiLim, 1*60*60);
    BPP.setParam(IloCplex::Threads, 1);

    BPP.exportModel("model.lp");

    double after;
    double before;

    try{
        after = BPP.getTime();
        BPP.solve();
        before = BPP.getTime();
    }catch(IloException& e){
        cout << e;
    }

    printResults(BPP, argv[1], after-before);
    IloNumArray duals(env);
    BPP.getDuals(duals, cstr);
    for(int i = 0; i < data.getItemQnt(); i++){
        cout << duals[i] << " ";
    }

    getColumn(var, cstr, data.getItemQnt());
    cout << endl;
    env.end();
    
    */

    return 0;

}
