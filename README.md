# Fast Fourier Transform (FFT) Implementation

This repository contains C++ implementations of the **Fast Fourier Transform (FFT)** algorithm in both **iterative** and **recursive** versions. FFT is a widely used algorithm to compute the discrete Fourier transform (DFT) and its inverse, which has applications in signal processing, image processing, and more.

## Overview
This project provides two implementations of FFT:
- **Recursive FFT**: A divide-and-conquer algorithm that recursively breaks down a DFT of any composite size.
- **Iterative FFT**: A non-recursive, in-place computation using the Cooley-Tukey algorithm.

Both implementations are written in C and tested with example datasets.

## Features
- **Recursive and Iterative Implementations**.
- Supports **complex number** inputs for real-world applications.
- Includes **test cases** for correctness and performance comparison.
- Outputs frequency domain data for any 1D input signal.

## How FFT Works
FFT is an efficient algorithm to compute the Discrete Fourier Transform (DFT), reducing the time complexity from \( O(n^2) \) to \( O(n \log n) \). It decomposes a signal into its constituent frequencies.

- **Recursive Version**: Uses a divide-and-conquer approach.
- **Iterative Version**: An in-place computation, leveraging the Cooley-Tukey algorithm for efficiency.
