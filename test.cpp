#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;

#include "dft.hpp"
#include "fft_recursive.hpp"
#include "fft_iterative.hpp"


const double EPSILON = 1e-9;

bool compare_complex_vectors(const vector<complex<double>> &a, const vector<complex<double>> &b, double epsilon = EPSILON) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (abs(a[i] - b[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

void test_single_element() {
    vector<complex<double>> input = { {1.0, 0.0} };
    auto result_dft = dft(input);
    auto result_fft_rec = fft_rec(input);
    auto result_fft_iter = fft_iter(input);

    assert(compare_complex_vectors(result_dft, input));
    assert(compare_complex_vectors(result_fft_rec, input));
    assert(compare_complex_vectors(result_fft_iter, input));

    std::cout << "Test Case 1 Passed: Single element vector.\n";
}

void test_two_elements() {
    vector<complex<double>> input = { {1.0, 0.0}, {0.0, 0.0} };
    vector<complex<double>> expected = { {1.0, 0.0}, {1.0, 0.0} };

    auto result_dft = dft(input);
    auto result_fft_rec = fft_rec(input);
    auto result_fft_iter = fft_iter(input);

    assert(compare_complex_vectors(result_dft, expected));
    assert(compare_complex_vectors(result_fft_rec, expected));
    assert(compare_complex_vectors(result_fft_iter, expected));

    std::cout << "Test Case 2 Passed: Two-element vector.\n";
}

void test_four_elements() {
    vector<complex<double>> input = { {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
    vector<complex<double>> expected = {
        {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}
    };

    auto result_dft = dft(input);
    auto result_fft_rec = fft_rec(input);
    auto result_fft_iter = fft_iter(input);

    assert(compare_complex_vectors(result_dft, expected));
    assert(compare_complex_vectors(result_fft_rec, expected));
    assert(compare_complex_vectors(result_fft_iter, expected));

    std::cout << "Test Case 3 Passed: Four-element vector.\n";
}

void test_random_complex_input() {
    vector<complex<double>> input = {
        {1.0, -1.0}, {0.0, 1.0}, {0.0, 0.0}, {1.0, 0.0}
    };

    auto result_dft = dft(input);
    auto result_fft_rec = fft_rec(input);
    auto result_fft_iter = fft_iter(input);

    assert(compare_complex_vectors(result_fft_rec, result_dft));
    assert(compare_complex_vectors(result_fft_iter, result_dft));

    std::cout << "Test Case 4 Passed: Random complex input.\n";
}

int main() {
    test_single_element();
    test_two_elements();
    test_four_elements();
    test_random_complex_input();

    std::cout << "All test cases passed successfully!\n";
    return 0;
}
