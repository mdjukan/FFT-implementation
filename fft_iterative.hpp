#ifndef FFT_ITER_INCLUDED
#define FFT_ITER_INCLUDED

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
using namespace std;

size_t flip_binary(size_t x, size_t len) {
	size_t nd = 0;
	while (len!=1) {
		nd ++;
		len >>= 1;
	}

	for (size_t i=0; i<nd/2; i++) {
		size_t bit_left = (x&(1<<i))>>i;
		size_t bit_right = (x&(1<<(nd-1-i)))>>(nd-1-i);
		x &= ~(1<<i);
		x &= ~(1<<(nd-1-i));
		x |= bit_left<<(nd-1-i);
		x |= bit_right<<i;
	}

	return x;
}

vector<complex<double>> fft_iter(const vector<complex<double>> &X) {
	size_t size = 1;
	while (size<X.size()) {
		size <<= 1;
	}

	vector<complex<double>> Xc(size, polar(0.0, 0.0));
	for (size_t i=0; i<X.size(); i++) {
		Xc[i] = X[i];
	}

	vector<complex<double>> Y(Xc.size());
	vector<complex<double>> Ytmp(Xc.size());

	for (size_t i=0; i<Xc.size(); i++) {
		Y[flip_binary(i, Xc.size())] = Xc[i];
	}

	size_t chunk = 1;
	while (chunk!=Xc.size()) {
		for (size_t i=0; i<Xc.size(); i += 2*chunk) {
			size_t bot = i;
			size_t top = chunk + i;
			for (size_t j=0; j<chunk; j++) {
				complex<double> omega = polar(1.0, - M_PI / chunk * j);
				Ytmp[bot+j] = Y[bot+j] + omega * Y[top+j];
				Ytmp[top+j] = Y[bot+j] - omega * Y[top+j];
			}
		}

		chunk *= 2;
		Y = Ytmp;
	}

	return Y;
}
#endif //FFT_ITER_INCLUDED
