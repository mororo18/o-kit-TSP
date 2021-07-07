#include "auxFunctions.h"
#include "Data.h"
#include "Master.h"
#include "SubProblem.h"
#include <ilcplex/ilocplex.h>

int main(int argc, char * argv[]){
    Data data(argc, argv[1]);

    data.loadData();

    /*
    cout << data.getItemQnt() << endl;
    cout << data.getBinCapacity() << endl;
    for(int i = 0; i < data.getItemQnt(); i++){
        cout << i << " " << data.getWeight(i) << endl;
    }
    cout << endl;
    */
    //exit(1);
    

    Master      M  (data.getItemQnt());
    SubProblem  kp (data.getWeights(), data.getBinCapacity(), data.getItemQnt());

    clock_t begin = clock();
    double obj_value;
    double pricing;
    while(true){
        obj_value = M.solve();
        pricing = kp.solve(M.getDuals());

        if(pricing >= -0.000001)
            break;

        vector<int> column = kp.getColumn();
        M.addColumn(column);
    }
    clock_t end = clock();

    //M.printResult();
    cout << "ObjValue: " << obj_value << endl;
    cout << "Pricing Value " << pricing << endl;
    cout << "Time: " << (double)(end-begin)/CLOCKS_PER_SEC << endl;;



    /*IloNumArray * duals = M.getDuals();

    int dimension = data.getItemQnt();

    IloEnv          env;
    IloModel model(env);

    IloBoolVarArray x (env, dimension);
    // add variables
    for(int i = 0; i < dimension; i++){
        char name[100];
        sprintf(name, "x_%d", i);
        x[i].setName(name);
        model.add(x[i]);
    }

    // add OF
    IloExpr obj(env);

    for(int i = 0; i < dimension; ++i){
        obj += (*duals)[i] * x[i];
    }
    model.add(IloMaximize(env, obj));

    // add cstrs
    {
    IloRange cstr;
    IloExpr LHS(env);
    for(int i = 0; i < dimension; i++){
        LHS += data.getWeight(i) * x[i];
    }
    cstr = (LHS <= data.getBinCapacity());
    model.add(cstr);
    }
    
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
    for(int i = 0; i < dimension; i++){
        cout << KP.getValue(x[i]) << " ";
    }
    cout << "opa" << endl;
    env.end();
    cout << "opa" << endl;

    */
    /*
    cout << data.getItemQnt() << endl;

    cout << data.getBinCapacity() << endl;


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
