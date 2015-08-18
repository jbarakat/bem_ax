#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "/home/jmb/Documents/install/lapack-3.5.0/lapacke/include/lapacke.h"
#include "/home/jmb/Documents/install/gsl-1.16/gsl_types.h"

using namespace std;

int main(){
	vector<double> x(5);

	cout << x.size() << endl;

	int i, N;
	N = 5;

	for (i = 0; i < N; i++){
		x.push_back(i);
		printf("%.4f\n", x[i]);
	}
	
}
