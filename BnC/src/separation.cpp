#include "separation.h"

double weight_sum(int vertex, double ** x, int n){
    double sum = 0;
    for(int j = vertex + 1; j < n; j++){
        if(x[vertex][j] > DBL_EPSILON)
            sum += x[vertex][j];
    }

    return sum;
}

double max_back_calc(int vertex, const vector<int> S, double ** x, int n){
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

void max_back_vec_init(vector<double> & vec, vector<int> S, double ** x, int n){
    for(int i = 0; i < vec.size(); i++){
        if(!belong_to(i, S)){
            vec[i] = max_back_calc(i, S, x, n);
        }
    }
}

int max_index_get(vector<double> vec, vector<int> S){
    double max = -999999;
    double max_i;
    for(int i = 0; i < vec.size(); i++){
        if(vec[i] - DBL_EPSILON >= max && !belong_to(i, S)){
            max = vec[i];
            max_i = i;
        }

    }

    return max_i;
}

void max_back_update(vector<double> & vec, int vertex, vector<int> S, double ** x, int n){
    for(int i = 0; i < vec.size(); i++){
        if(!belong_to(i, S)){
            if(i < vertex)
                vec[i] += x[i][vertex];
            else
                vec[i] += x[vertex][i];
        }
    }
}

void vec_print_dbl(vector<double> vec, string title){
    cout << title;

    for(int i = 0; i < vec.size(); i++){
        cout << vec[i] << " ";
    }

    cout << endl;
}

double min_cut_of(vector<int> & S_min, double ** x, int n){

    vector<int> S = S_min;
    double cut_min = weight_sum(S[0], x, n);
    double cut_value = cut_min;
    vector<double> max_back_vec (n);

    max_back_vec_init(max_back_vec, S, x, n);

    while(S.size() < n){
        int v = max_index_get(max_back_vec, S); 
        S.push_back(v);

        cut_value += 2 - 2*max_back_vec[v];

        max_back_update(max_back_vec, v, S, x, n);

        if(cut_value < cut_min){
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
    bool found[n] = {};
    int vertex;
    
    while(!all_known(found, n)){
        for(int i = 0; i < n; i++){
            if(found[i] != true){
                vertex = i;
                break;
            }
        }

        vector<int> S_min = {vertex};
        min_cut_of(S_min, x, n);

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

void V_set_merge(vector<vector<int>> & set, int s, int t, double ** x, int n){
    int a, b;

    if(s < t){
        a = s;
        b = t;
    }

    for(int i = 0; i < set.size(); i++){
        if(a == set[i][0]){
            a_i = i;
        }else if(b == set[i][0]){
            set.insert(set[a_i].end(), set[i].begin(), set[i].end());
            set.erase(set.begin() + i);
        }
    }

}

vector<vector<int>> MinCut(double ** x, int n){
    vector<vector<int>> V_set;


    while(V_set.size() > 1){
        vector<int> S_min  = {0};
        double cut_value = min_cut_of(S_min, x, n);

        int t = S_min[S_min.size() - 1];
        int s = S_min[S_min.size() - 2];
        
    }

    return empty;
}
