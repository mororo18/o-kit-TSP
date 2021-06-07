//---------------------------------------------------------------------------

/***************************************************/
/* Functions prototypes by Prof. Anand Subramanian */
/***************************************************/

#ifndef Separation_H
#define Separation_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <cfloat>
#include <time.h>

#define EPSILON 0.00000001

using namespace std;

extern vector <vector<int> > MaxBack(double** x, int n);
extern vector <vector<int> > MinCut(double** x, int n);

#endif

//---------------------------------------------------------------------------
