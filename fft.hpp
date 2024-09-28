#ifndef FFT_INCLUDED
#define FFT_INCLUDED

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
using namespace std;

vector<complex<double>> fft(const vector<complex<double>> &X) {
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

	vector<complex<double>> Ye = fft(Xe);
	vector<complex<double>> Yo = fft(Xo);

	vector<complex<double>> Y(n);
	for (size_t k=0; k<m; k++) {
		complex<double> omega = polar(1.0, k * (-2*M_PI/n));
		Y[k] = Ye[k] + omega*Yo[k];
		Y[m + k] = Ye[k] - omega*Yo[k];
	}

	return Y;
}

vector<complex<double>> ifft_rec(const vector<complex<double>> &X) {
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

	vector<complex<double>> Ye = ifft_rec(Xe);
	vector<complex<double>> Yo = ifft_rec(Xo);

	vector<complex<double>> Y(n);
	for (size_t k=0; k<m; k++) {
		complex<double> omega = polar(1.0, k * (2*M_PI/n));
		Y[k] = Ye[k] + omega*Yo[k];
		Y[m + k] = Ye[k] - omega*Yo[k];
	}

	return Y;
}

vector<complex<double>> ifft(const vector<complex<double>> &X) {
	vector<complex<double>> Y = ifft_rec(X);
	for (size_t i=0; i<Y.size(); i++) {
		Y[i] /= X.size();
	}

	return Y;
}

/*
int main() {
	vector<complex<double>> X = {{1, 0}, {1, 0}, {1, 0}, {0, 0}};
	vector<complex<double>> IFFT = ifft(X);
	for (auto x : IFFT) {
		cout << x << endl;
	}
}
*/
#endif //FFT_INCLUDED
