#include "auxFunctions.h"
#include "Data.h"
#include "Master.h"
#include "SubProblem.h"
#include "measure.h"
#include <cfloat>
#include <ilcplex/ilocplex.h>

#define DEPTH   0
#define BREADTH 1
#define BEST_BD 2

MEASURE_BLOCK enforce;
MEASURE_BLOCK exclude;
MEASURE_BLOCK GC;
MEASURE_BLOCK masterP;
MEASURE_BLOCK subP;
MEASURE_BLOCK add_col;
MEASURE_BLOCK s_load;

double upper_bd = DBL_MAX; 
int colunas = 0;

struct node_info {
    vector<double> lambda_solution;
    vector<vector<bool>> columns;
    vector<pair<int, int>> exclude;
    vector<pair<int, int>> enforce;
    double value;
    bool feasibility;
};

struct node_info node_best;

void solution_load(struct node_info & node, Master & M){
    vector<IloNumVar> solution = M.getSolution();
    node.lambda_solution = M.getSolutionValues();

    for(int i = 0; i < solution.size(); i++){
        vector<bool> column = M.getColumn(solution[i]);
        node.columns.push_back(column);
    }

    cout << node.value << " " << solution.size() << endl;

    if(node.value >= solution.size() - 2*DBL_EPSILON)
        node.feasibility = true;
    else
        node.feasibility = false;

}

void column_generation(struct node_info & node, Master & BPP, Data & data){
    //vector<pair<int, int>> test = {make_pair(18, 39)};

    //SubProblem  KP (data.getWeights(), data.getBinCapacity(), data.getItemQnt(), node.exclude, test);
    SubProblem  KP (data.getWeights(), data.getBinCapacity(), data.getItemQnt(), node.exclude, node.enforce);

    double obj_value;
    double pricing;

    cout << "ANTES" << endl;
    while(true){
        measure_before(masterP);
        obj_value = BPP.solve();
        measure_after(masterP);
        //cout << "crnt obj " << obj_value << endl;
        measure_before(subP);
        pricing = KP.solve(BPP.getDuals());
        measure_after(subP);
        //cout << "crnt prcing " << pricing << endl;
        //cout << endl;

        if(pricing >= -0.000001)
            break;

        colunas++;

        measure_before(add_col);
        vector<bool> column = KP.getColumn();
        BPP.addColumn(column);
        measure_after(add_col);

    }
    cout << "DEPOIS" << endl;

    cout << obj_value << endl;
    cout << colunas << " Colunas" << endl;

    //if(node.exclude.size() == 5)
        //exit(1);

    //exit(1);
    //if(node.exclude.empty() && !node.enforce.empty())
    //exit(1);

    cout << "Exclude ";
    for(int i = 0; i < node.exclude.size(); i++){
        cout << node.exclude[i].first << "-" << node.exclude[i].second << " ";
    }
    cout << endl;


    cout << "Enforce ";
    for(int i = 0; i < node.enforce.size(); i++){
        cout << node.enforce[i].first << "-" << node.enforce[i].second << " ";
    }
    cout << endl;

    node.value = obj_value;
    measure_before(s_load);
    solution_load(node, BPP);
    measure_after(s_load);
    cout << "auq" << endl;
}

pair<int, int> get_most_fractional_pair(vector<double> var_value, vector<vector<bool>> columns){
    int n = columns[0].size();
    double pair_value[n][n];
    memset(pair_value, 0.0, sizeof(double)*n*n);

    for(int k = 0; k < var_value.size(); k++){
        if(var_value[k] < 1.0 - DBL_EPSILON && var_value[k] > DBL_EPSILON){
            vector<int> index;
            index.reserve(n);
            for(int i = 0; i < n; i++){
                if(columns[k][i] == 1)
                    index.push_back(i);
            }

            for(int i = 0; i < index.size()-1; i++){
                for(int j = i+1; j < index.size(); j++){
                    pair_value[index[i]][index[j]] += var_value[k];
                }
            }
        }
    }

    int item_a;
    int item_b;
    double diff_min = DBL_MAX;

    for(int i = 0; i < n-1; i++){
        for(int j = i+1; j < n; j++){
            double diff = abs(pair_value[i][j] - 0.5);
            if(diff < diff_min - DBL_EPSILON){
                diff_min = diff;
                item_a = i;
                item_b = j;
            }
            //cout << pair_value[i][j] << " ";
        }
        //cout << endl;
    }
    cout << diff_min << " com  " << item_a << " e " << item_b << endl;
    return make_pair(item_a, item_b);
}

void search_depth(struct node_info & node, Master & M, Data & data, int gen){
    measure_before(GC);
    column_generation(node, M, data);
    measure_after(GC);

    //if(gen == 1)
        //exit(1);

    if(node.value - DBL_EPSILON <= upper_bd - 1 ){
        if(node.feasibility == true){
            //cout << "aue"<< endl;
            node_best = node;
            upper_bd = node.value;
            //exit(1);
        }else{
            // branching rule
            //cout << "pee" << endl;
            pair<int, int> m_frac = get_most_fractional_pair(node.lambda_solution, node.columns);
            for(int i = 0; i < node.columns.size(); i++){
                node.columns[i].clear(); 
            }
            node.columns.clear(); 
            node.lambda_solution.clear(); 
            //cout << "paee" << endl;


            cout << endl << " BRANCH Exclude "<<endl;
            // enforce

            struct node_info node_son_B;
            node_son_B.exclude = node.exclude;
            node_son_B.enforce = node.enforce;

            node_son_B.enforce.push_back(m_frac);


            node.exclude.clear();

            measure_before(enforce);
            M.enforce(m_frac);
            measure_after(enforce);

            search_depth(node_son_B, M, data, gen+1);
            M.reinsert(m_frac);



            // se a diferenca entre o valor do pai e o upper bound ainda continua maior q 1
            if(node.value - DBL_EPSILON <= upper_bd - 1.0f){
                
            // exclude 
            struct node_info node_son_A;
            node_son_A.exclude = node_son_B.exclude;
            node_son_A.enforce = node.enforce;

            node.enforce.clear();

            //cout << "paee2" << endl;
            node_son_A.exclude.push_back(m_frac);

            measure_before(exclude);
            M.exclude(m_frac);
            measure_after(exclude);

            search_depth(node_son_A, M, data, gen+1);
            M.reinsert(m_frac);

                cout << endl << " BRANCH Enforce "<<endl;
            }
        }
    }
}

void branch_and_price(Master & M, Data & data, int opt, int gen){
    struct node_info node_root;
    switch (opt){
        case DEPTH:
            search_depth(node_root, M, data, gen);
            break;
        case BREADTH:
            break;
        case BEST_BD:
            break;
    }
}

int main(int argc, char * argv[]){
    Data data(argc, argv[1]);

    data.loadData();

    measure_block_init(&enforce);
    measure_block_init(&exclude);
    measure_block_init(&GC);
    measure_block_init(&masterP);
    measure_block_init(&subP);
    measure_block_init(&add_col);
    measure_block_init(&s_load);


    double before = now();
    Master  BPP (data.getItemQnt());
    branch_and_price(BPP, data, DEPTH, 0);
    double after = now();

    cout << "Enforce\t" << measure_total(enforce) << endl;
    cout << "Exclude\t" << measure_total(exclude) << endl;
    cout << "Geracao de Col\t" << measure_total(GC) << endl;
    cout << "Mestre\t" << measure_total(masterP) << endl;
    cout << "SubProb\t" << measure_total(subP) << endl;
    cout << "Solution Load\t" << measure_total(s_load) << endl;
    cout << "Add Col\t" << measure_total(add_col) << endl;


    cout << "Time: " << (after- before) << endl;;
    cout << "Solution:  " << node_best.value << endl;
    

    /*
    vector<vector<int>> bins = BPP.getBins();
    for(int i = 0; i < node_best.columns.size(); i++){
        vector<int> bin;
        for(int j = 0; j < node_best.columns[i].size(); j++){
            if(node_best.columns[i][j]){
                bin.push_back(j);
                cout << j << " ";
            }
        }
        cout << endl;
        bins.push_back(bin);
    }
    */

    return 0;

}
