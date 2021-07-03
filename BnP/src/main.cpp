#include "Data.h"
#include <ilcplex/ilocplex.h>

int main(int argc, char * argv[]){
    Data data(argc, argv[1]);

    data.loadData();

    cout << data.getItemQnt() << endl;

    cout << data.getBinCapacity() << endl;

    for(int i=0; i < data.getItemQnt(); i++){
        cout << data.getWeight(i) << endl;
    }

    return 0;

}
