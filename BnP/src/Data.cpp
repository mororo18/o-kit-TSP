#include "Data.h"

Data::Data(int param_qnt, char * file_name){

    if (param_qnt < 2) {
        cout << "Missing parameters\n";
        cout << " ./exe [Instance]"  << endl;
        exit(1);        
    }

    if (param_qnt > 2) {
        cout << "Too many parameters\n";
        cout << " ./exe [Instance]"  << endl;
        exit(1);
    }

    instance_name = file_name;
}

Data::~Data(){
    delete [] weight;
}

void Data::loadData(){
    ifstream file;

    file.open(instance_name);

    string line;

    // get the qnty n of itens
    std::getline(file, line);
    n = stoi(line);

    // get the capacity of the bins
    getline(file, line);
    bin_C = stod(line);

    // get the weights of the itens
    weight = new int[n];

    int count = 0;
    while(std::getline(file, line)){
        weight[count] = stod(line); 
        count++;
    }

    if(count != n){
        cout << "Invalid instancen\n";
    }

    file.close();

}

int Data::getItemQnt(){
    return n;
}

int Data::getBinCapacity(){
    return bin_C;
}

int * Data::getWeights(){
    int * p = (int*) calloc(n, sizeof(int));
    memcpy(p, weight, sizeof(int)*n);
    return p;
}

int Data::getWeight(int i){
    return weight[i];
}
