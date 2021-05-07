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
    for(int u : S){
        if(u == vertex)
            return true;
    }

    return false;
}

void max_back_vec_init(vector<double> & vec, vector<int> S, double ** x, int n){
    for(int i = 0; i < vec.size(); i++){
        if(!belong_to(i, S))
            vec[i] += max_back_calc(i, S, x, n);
    }
}

int max_index_get(vector<double> vec, vector<int> S){
    double max = -999999;
    double max_i;
    for(int i = 0; i < vec.size(); i++){
        if(vec[i] - DBL_EPSILON > max && !belong_to(i, S)){
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

vector<vector<int>> MaxBack(double ** x, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(x[i][j] < DBL_EPSILON &&  x[i][j] > - DBL_EPSILON)
                x[i][j] = 0.0f;

            if(j < j)
                cout << "0 ";
            else
                cout << x[i][j] << " ";
        }
        cout << endl;
    }

    
    
    vector<int> S = {3};
    vector<int> S_min;
    double cut_min = weight_sum(S[0], x, n);
    double cut_value = cut_min;
    vector<double> max_back_vec (n);

    max_back_vec_init(max_back_vec, S, x, n);

    while(S.size() < n){
        int v = max_index_get(max_back_vec, S); 
        //cout << "V " << v << endl;
        S.push_back(v);

        //cout << "cut before " << cut_value << endl;
        cut_value += 2 - 2*max_back_vec[v];
        //cout << "cut after " << cut_value << endl;

        //vec_print_dbl(max_back_vec, "vec before  ");
        max_back_update(max_back_vec, v, S, x, n);
        //vec_print_dbl(max_back_vec, "vec after  ");

        if(cut_value < cut_min){
            cut_min = cut_value;
            S_min = S;
        }

        //cout << "min cut " << cut_min << endl;

    }

    for(int i = 0 ; i < S_min.size(); i++){
        cout << S_min[i] << " ";
    }
    cout << "\n";

    exit(0);
}

vector<vector<int>> MinCut(double ** x, int n){
    vector<vector<int>> empty;

    return empty;
}
