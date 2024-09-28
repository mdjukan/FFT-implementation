#include "fft.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>
#include <complex>

typedef vector<double> POLY;
typedef vector<complex<double>> SIG;
#define EPS 1e-5

SIG poly_to_sig(const POLY &p, size_t length) {
	SIG X(length, 0);
	for (size_t i=0; i<p.size(); i++) {
		X[i] = p[i];
	}

	return X;
}

bool all_imaginary_zero(const SIG &X) {
	for (auto x : X) {
		if (abs(imag(x))>EPS) {
			return false;
		}
	}

	return true;
}

POLY sig_to_poly(const SIG &X) {
	assert(all_imaginary_zero(X));
	POLY p;
	for (auto x : X) {
		p.push_back(real(x));
	}

	auto it = p.end()-1;
	while (abs(*it)<EPS) {
		it --;
		p.pop_back();
	}
	return p;
}

POLY poly_mul(const POLY &P, const POLY &Q) {
	size_t length = P.size() + Q.size() + 1;
	size_t length_pow2 = 1;
	while (length>length_pow2) {
		length_pow2 <<= 1;
	}
	
	SIG p = poly_to_sig(P, length_pow2);
	SIG q = poly_to_sig(Q, length_pow2);

	SIG pt = fft(p);
	SIG qt = fft(q);

	SIG rt;
	for (size_t i=0; i<pt.size(); i++) {
		rt.push_back(pt[i] * qt[i]);
	}

 SIG r = ifft(rt);
	POLY R = sig_to_poly(r);
	return R;
}

int main() {
	POLY P = {1, 0, 1};
	POLY Q = {0, 1, 1, 1};
	POLY R = poly_mul(P, Q);
	for (double c : R) {
		cout << c << " ";
	}
}
