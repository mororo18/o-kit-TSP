#ifndef DATA_H
#define DATA_H


#include <fstream> 
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

class Data {
    public:
        Data(int, char *);
        ~Data();
        int getItemQnt();
        int getBinCapacity();
        int getWeight(int);
        int * getWeights();
        void loadData();


    private:
        int n;
        int * weight;
        int bin_C;

        std::string instance_name;
};


#endif
