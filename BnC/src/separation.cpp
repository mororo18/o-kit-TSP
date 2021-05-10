#include "separation.h"

double weight_sum(int vertex, double ** x, int n){
    double sum = 0;
    for(int j = vertex + 1; j < n; j++){
        if(x[vertex][j] > DBL_EPSILON)
            sum += x[vertex][j];
    }

    // test
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

    for(int i = 0; i < set.size(); i++){
        if(a == set[i][0]){
            a_i = i;
        }else if(b == set[i][0]){
            set[a_i].insert(set[a_i].end(), set[i].begin(), set[i].end());
            set.erase(set.begin() + i);
        }
    }

    // find common connected vertex to s and t

    x[a][b] = 0;
    for(int i = 0; i < set.size(); i++){
        int v = set[i][0];
        if(v == s)
            continue;

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

vector<int> minimum_cut(int a, double ** x, int n){
    vector<vector<int>> V_set;
    V_set_init(V_set, n);

    int count = 0;

    vector<int> min_cut_vec ;
    double cut_min = 999999;
    int cut;

    while(V_set.size() > 1){
        //if(count++ == 7)
            //break;
        vector<int> A = {a};

        double cut_value = min_cut_phase(A, V_set.size(), x, n); 
        //vec_print_dbl(A, "A group ");


        int t = A[A.size() - 1];
        int s = A[A.size() - 2];

        //cout << "minimum  " << cut_value << endl;
        if(cut_value < cut_min){
            cut_min = cut_value;
            cut = t;

            for(int i =0; i < V_set.size(); i++){
                if(cut == V_set[i][0]){
                    min_cut_vec = V_set[i];
                    break;
                }
            }

            if(cut_min < DBL_EPSILON)
                break;
        }


        V_set_merge(V_set, s, t, x, n);

    }


    return min_cut_vec;
}

void matrix_alloc(double *** matrix, int dimension){
    double ** ptr =  (double**)calloc(dimension, sizeof(double*));
    test_alloc(ptr);
    for(int i = 0; i < dimension; i++){
        ptr[i] =  (double*)calloc(dimension, sizeof(double));
        test_alloc(ptr[i]);
    }

    *matrix = ptr;
}

void matrix_free(double *** ptr, int dimension){
    double ** ar = *ptr;
    for(int i = 0; i < dimension; i++)
        free(ar[i]);
    free(ar);

    *ptr = NULL;
}

void matrix_cpy(double ** ptr, double ** src, int n){
    for(int i = 0; i < n; i++)
        memcpy(ptr[i], src[i], n * sizeof(double));
}

vector<vector<int>> MinCut(double ** x, int n){
    vector<vector<int>> set_pool;

    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            if(x[i][j] < DBL_EPSILON &&  x[i][j] > - DBL_EPSILON)
                x[i][j] = 0.0f;

        }
    }

    vector<int> cut;

    /*
    double ** x_cpy;

    matrix_alloc(&x_cpy, n);
    matrix_cpy(x_cpy, x, n);
    int v = 0;
    bool find[n] = {false};

    while(v < n){
        cut = minimum_cut(v++, x_cpy, n);
        //vec_print_dbl(cut, "subtour ");

        if(find[cut[0]] == false && cut.size() != 1){
            find[cut[0]] = true;
            set_pool.push_back(cut);
        }
        matrix_cpy(x_cpy, x, n);
    }

    //cout << "aqui\n";
    matrix_free(&x_cpy, n);


    bool mark[n] = {false};
    for(int i = 0; i < set_pool.size(); i++){
        for(int j = 0; j < set_pool[i].size(); j++){
            mark[set_pool[i][j]] = true;
        }
    }

    vector<int> vec;

    for(int i = 0; i < n; i++){
        if(mark[i] != true){
            vec.push_back(i);
        }
    }

    set_pool.push_back(vec);
    */

    
    cut = minimum_cut((int)(n/2), x, n);

    if(cut.size() == 1){
        cout << "vazio" <<endl;
        return set_pool;
    }

    bool mark[n] = {false};


    for(int i = 0; i < cut.size(); i++){
        mark[cut[i]] = true;
    }

    vector<int> vec;

    for(int i = 0; i < n; i++){
        if(mark[i] != true){
            vec.push_back(i);
        }
    }

    set_pool.push_back(cut);
    set_pool.push_back(vec);

    /*
    for(int i = 0; i < set_pool.size(); i++){
        for(int j = 0; j < set_pool[i].size(); j++){
            cout << set_pool[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout << endl;
    */

    //exit(0);

    return set_pool;
}
