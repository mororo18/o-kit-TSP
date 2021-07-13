#include "auxFunctions.h"
#include "Data.h"
#include "Master.h"
#include "SubProblem.h"
#include <cfloat>
#include <ilcplex/ilocplex.h>

#define DEPTH   0
#define BREADTH 1
#define BEST_BD 2

double upper_bd = DBL_MAX; 
int colunas = 0;

struct node_info {
    //vector<vector<int>> solution;
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

    double count = 0;
    for(int i = 0; i < solution.size(); i++){
        count += node.lambda_solution[i];
        vector<bool> column = M.getColumn(solution[i]);
        node.columns.push_back(column);
    }

    if(count >= solution.size() - 2*DBL_EPSILON)
        node.feasibility = true;
    else
        node.feasibility = false;

}

void column_generation(struct node_info & node, Master & BPP, Data & data){
    vector<pair<int, int>> test = {make_pair(18, 39)};

    //SubProblem  KP (data.getWeights(), data.getBinCapacity(), data.getItemQnt(), node.exclude, test);
    SubProblem  KP (data.getWeights(), data.getBinCapacity(), data.getItemQnt(), node.exclude, node.enforce);

    double obj_value;
    double pricing;

    while(true){
        obj_value = BPP.solve();
        //cout << "crnt obj " << obj_value << endl;
        pricing = KP.solve(BPP.getDuals());
        //cout << "crnt prcing " << pricing << endl;
        IloNumArray dual = BPP.getDuals();
        for(int i = 0; i < data.getItemQnt(); i++){
            //cout << dual[i] << " ";
        }
        //cout << endl;

        if(pricing >= -0.000001)
            break;

        colunas++;

        vector<bool> column = KP.getColumn();
        BPP.addColumn(column);
    }

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
    solution_load(node, BPP);
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

            /*
            for(int i = 0; i < n-1; i++){
                if(column[k][i] == 1)
                    for(int j = i+1; j < n; j++){
                        if(columns[k][j] == 1)
                            pair_value[i][j] += var_value[k];
                    }
            }
            */
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
    //cout << diff_min << " com  " << item_a << " e " << item_b << endl;
    return make_pair(item_a, item_b);
}

void search_depth(struct node_info & node, Master & M, Data & data){
    column_generation(node, M, data);

    if(node.value  < upper_bd - DBL_EPSILON){
        if(node.feasibility == true){
            //cout << "aue"<< endl;
            node_best = node;
            upper_bd = node.value;
            //exit(1);
        }else{
            // branching rule
            //cout << "pee" << endl;
            pair<int, int> m_frac = get_most_fractional_pair(node.lambda_solution, node.columns);
            node.columns.clear(); 
            node.lambda_solution.clear(); 
            //cout << "paee" << endl;

            // exclude 
            struct node_info node_son_A;
            node_son_A.exclude = node.exclude;
            node_son_A.enforce = node.enforce;


            //cout << "paee2" << endl;
            node_son_A.exclude.push_back(m_frac);

            cout << endl << " BRANCH Exclude "<<endl;
            M.exclude(m_frac);
            search_depth(node_son_A, M, data);
            M.reinsert(m_frac);

            // enforce

            struct node_info node_son_B;
            node_son_B.exclude = node.exclude;
            node_son_B.enforce = node.enforce;

            node_son_B.enforce.push_back(m_frac);

            cout << endl << " BRANCH Enforce "<<endl;
            //M.enforce(m_frac);
            search_depth(node_son_B, M, data);
            //exit(1);
            //M.reinsert(m_frac);
        }
    }
}

void branch_and_price(Master & M, Data & data, int opt){
    struct node_info node_root;
    switch (opt){
        case DEPTH:
            search_depth(node_root, M, data);
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

    clock_t begin = clock();
    Master  BPP (data.getItemQnt());
    branch_and_price(BPP, data, DEPTH);
    clock_t end = clock();

    cout << "Solution  " << node_best.value << endl;

    cout << "Time: " << (double)(end-begin)/CLOCKS_PER_SEC << endl;;

    return 0;

}
