#ifndef DFT_INCLUDED
#define DFT_INCLUDED

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
using namespace std;

vector<complex<double>> dft(const vector<complex<double>> &X) {
	int N = X.size();
	vector<complex<double>> Y(N, polar(0.0,0.0));
	vector<complex<double>> omega(N);

	for (int i=0; i<N; i++) {
		omega[i] = polar(1.0, -2*i*M_PI/N);
	}

	for (int k=0; k<N; k++) {
		for (int j=0; j<N; j++) {
			Y[k] += X[j] * omega[(k*j) % N];
		}
	}

	return Y;
}

#endif //DFT_INCLUDED
