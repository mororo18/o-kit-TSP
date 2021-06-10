#include "separation.h"

double weight_sum(int vertex, double ** x, int n){
    double sum = 0;
    
    // column
    for(int j = vertex + 1; j < n; j++){
        if(x[vertex][j] > DBL_EPSILON)
            sum += x[vertex][j];
    }

    // row
    for(int i = 0; i < vertex; i++){
        if(x[i][vertex] > DBL_EPSILON)
            sum += x[i][vertex];
    }

    return sum;
}

double max_back_calc(int vertex, const vector<int> S, double ** x){
    double sum = 0;
    for(int u : S){
        if(vertex < u)
            sum += x[vertex][u];
        else
            sum += x[u][vertex];
    }

    return sum;
}

bool belong_to(int vertex, vector<int> S){
    for(int & u : S){
        if(u == vertex)
            return true;
    }

    return false;
}

void max_back_vec_init(vector<double> & vec, vector<int> S, double ** x){
    for(int i = 0; i < vec.size(); i++){
        if(!belong_to(i, S)){
            vec[i] = max_back_calc(i, S, x);
        }
    }
}

int max_index_get(vector<double> vec, vector<int> S){
    double max = -999999;
    double max_i;
    for(int i = 0; i < vec.size(); i++){
        if(vec[i] > max && !belong_to(i, S)){
            max = vec[i];
            max_i = i;
        }

    }

    return max_i;
}

void max_back_update(vector<double> & vec, int vertex, vector<int> S, double ** x){
    for(int i = 0; i < vec.size(); i++){
        if(!belong_to(i, S)){
            if(i < vertex)
                vec[i] += x[i][vertex];
            else
                vec[i] += x[vertex][i];
        }
    }
}

void vec_print_dbl(vector<int> vec, string title){
    cout << title;

    for(int i = 0; i < vec.size(); i++){
        cout << vec[i] << " ";
    }

    cout << endl;
}

double min_cut_of(vector<int> & S_min, int set_size, double ** x, int n){

    vector<int> S = S_min;
    double cut_min = weight_sum(S[0], x, n);
    double cut_value = cut_min;
    vector<double> max_back_vec (n);

    max_back_vec_init(max_back_vec, S, x);

    while(S.size() < set_size){
        int v = max_index_get(max_back_vec, S); 
        S.push_back(v);

        cut_value += 2 - 2*max_back_vec[v];

        max_back_update(max_back_vec, v, S, x);

        if(cut_value < cut_min - DBL_EPSILON){
            cut_min = cut_value;
            S_min = S;
        }
    }

    return cut_min;

    /*
       for(int i = 0 ; i < S_min.size(); i++){
       cout << S_min[i] << " ";
       }
       cout << "\n";
       */

}

bool all_known(bool found[], int n){
    for(int i = 0; i < n; i++){
        if(found[i] == false)
            return false;
    }

    return true;
}

vector<vector<int>> MaxBack(double ** x, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(x[i][j] < DBL_EPSILON &&  x[i][j] > - DBL_EPSILON)
                x[i][j] = 0.0f;

        }
    }

    vector<vector<int>> set_pool;
    bool found[n] = {false};
    int vertex;

    while(!all_known(found, n)){
        for(int i = 0; i < n; i++){
            if(found[i] != true){
                vertex = i;
                break;
            }
        }
        
        vector<int> S_min = {vertex};
        min_cut_of(S_min, n, x, n);

        if(S_min.size() == n)
            break;

        for(int i = 0; i < S_min.size(); i++)
            found[S_min[i]] = true;

        set_pool.push_back(S_min);
    }

    return set_pool;
}

void V_set_init(vector<vector<int>> & vec, int n){
    vec.clear();

    for(int i = 0; i < n; i++){
        vector<int> vertex = {i};
        vec.push_back(vertex);
    }
}

double min_cut_phase(vector<int> & A, int set_size, double ** x, int n){
    int v;
    vector<double> weight_vec(n);

    max_back_update(weight_vec, A[0], A, x);

    while(A.size() < set_size){
        v = max_index_get(weight_vec, A);
        max_back_update(weight_vec, v, A, x);
        A.push_back(v);
    }

    return weight_sum(v, x, n);
}

bool edge_from_to(int a, int b, double ** x){
    if(a < b && x[a][b] > DBL_EPSILON){
        return true; 
    }else if(x[b][a] > DBL_EPSILON){
        return true;
    }else
        return false;
}

void V_set_merge(vector<vector<int>> & set, int s, int t, double ** x, int n){
    int a, b;
    int a_i;

    if(s < t){
        a = s;
        b = t;
    }else{
        a = t;
        b = s;
    }

    // merges the vertices
    // the one with smaller index represents the group
    for(int i = 0; i < set.size(); i++){
        if(a == set[i][0]){
            a_i = i;
        }else if(b == set[i][0]){
            set[a_i].insert(set[a_i].end(), set[i].begin(), set[i].end());
            set.erase(set.begin() + i);
        }
    }

    // find common vertex connected to s and t

    x[a][b] = 0;
    for(int i = 0; i < set.size(); i++){
        int v = set[i][0];
        if(v == s)
            continue;

        // merges the edges on the smaller index vertex
        if(edge_from_to(v, s, x) || edge_from_to(v, t, x)){

            if(v < a){
                x[v][a] += x[v][b];
                x[v][b] = 0;
            }else if(v > a && v < b){
                x[a][v] += x[v][b];
                x[v][b] = 0;
            }else if(v > b){
                x[a][v] += x[b][v];
                x[b][v] = 0;
            }
        }
    }

}

vector<int> set_complement(vector<int> set, int n){
    vector<int> comp;
    bool mark[n] = {false};

    for(int i = 0; i < set.size(); i++){
        mark[set[i]] = true;
    }

    for(int i = 0; i < n; i++){
        if(mark[i] != true){
            comp.push_back(i);
        }
    }

    return comp;
}

vector<vector<int>> MinCut(double ** x, int n){
    vector<vector<int>> V_set;
    vector<vector<int>> set_pool;

    srand(time(NULL));

    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            if(x[i][j] < DBL_EPSILON &&  x[i][j] > - DBL_EPSILON)
                x[i][j] = 0.0f;

        }
    }

    int count = 0;
    const double violation_bound = 2.0f - EPSILON;

    vector<int> min_cut_vec;
    vector<int> cut_vec;
    double cut_min = 999999;
    int cut_i;

    int a = rand() % n;

    V_set_init(V_set, n);

    while(V_set.size() > 1){
        vector<int> A = {a};

        double cut_value = min_cut_phase(A, V_set.size(), x, n); 

        int t = A[A.size() - 1];
        int s = A[A.size() - 2];

            /*
        //cout << "minimum  " << cut_value << endl;
        if(cut_value < cut_min ){
            cut_min = cut_value;

            // stores the minCut
            cut_i = t;

            for(int i = 0; i < V_set.size(); i++){
                if(cut_i == V_set[i][0]){
                    min_cut_vec = V_set[i];
                    break;
                }
            }
        }
            */

        for(int i = 0; i < V_set.size(); i++){
            if(t == V_set[i][0]){
                cut_vec = V_set[i];
                break;
            }
        }

        // verifies if the current cut represents a violation of subtour constraint 
        if(cut_vec.size() <= n / 2 && cut_value < violation_bound){ 
            set_pool.push_back(cut_vec);

        }else if(cut_vec.size() > n / 2 && cut_value < violation_bound){
            vector<int> comp = set_complement(cut_vec, n);
            set_pool.push_back(comp);
        }

        V_set_merge(V_set, s, t, x, n);
    }

    /*
    if(cut_min >= 2.0f - DBL_EPSILON){
        cout << "(vazio)\nsem violacoes\n";
    }
    */

    /*
    set_pool.clear();

    if(!min_cut_vec.empty()){
        set_pool.push_back(min_cut_vec);
        //if(!min_cut_vec.empty());
        set_pool.push_back(set_complement(min_cut_vec, n));

        for(int i = 0; i < set_pool.size(); i++){
            for(int j = 0; j < set_pool[i].size(); j++){
                cout << set_pool[i][j] << " ";
            }
            cout << endl;
        }

    }
    */

    return set_pool;
}
