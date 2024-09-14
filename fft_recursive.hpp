#ifndef FFT_REC_INCLUDED
#define FFT_REC_INCLUDED

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
using namespace std;

vector<complex<double>> fft_rec(const vector<complex<double>> &X) {
	if (X.size()==1) {
		return X;
	}

	size_t m = X.size()/2;
	size_t n = X.size();

	vector<complex<double>> Xe, Xo;
	for (size_t i=0; i<n; i+=2) {
		Xe.push_back(X[i]);
	}

	for (size_t i=1; i<n; i+=2) {
		Xo.push_back(X[i]);
	}

	vector<complex<double>> Ye = fft_rec(Xe);
	vector<complex<double>> Yo = fft_rec(Xo);

	vector<complex<double>> Y(n);
	for (size_t k=0; k<m; k++) {
		complex<double> omega = polar(1.0, k * (-2*M_PI/n));
		Y[k] = Ye[k] + omega*Yo[k];
		Y[m + k] = Ye[k] - omega*Yo[k];
	}

	return Y;
}

/*
vector<complex<double>> fft(const vector<complex<double>> &X) {
	size_t size = 1;
	while (size<X.size()) {
		size <<= 1;
	}

	vector<complex<double>> Xc(size, polar(0.0, 0.0));
	for (size_t i=0; i<X.size(); i++) {
		Xc[i] = X[i];
	}

	vector<complex<double>> Y = fft_rec(Xc);

	while (Y.size()!=X.size()) {
		Y.pop_back();
	}
	
	return Y;
}
*/
#endif //FFT_REC_INCLUDED
