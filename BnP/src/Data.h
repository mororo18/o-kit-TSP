#ifndef DATA_H
#define DATA_H


#include <fstream> 
#include <iostream>
#include <string>

using namespace std;

class Data {
    public:
        Data(int, char *);
        ~Data();
        int getItemQnt();
        double getBinCapacity();
        double getWeight(int);
        void loadData();


    private:
        int n;
        double * weight;
        double bin_C;

        std::string instance_name;
};


#endif
