#ifndef FOURIER_SERIES_PARALLEL_H
#define FOURIER_SERIES_PARALLEL_H

// Discrete Fourier Transform
// First step for a function approximation, used to get DFT coefficients 
// returns an array of DFT coefficients
double* DFTSampleParallel(double** pointsArray, const int arraySize, const int k); // one should not use that func outside of DFTSamples()
double** DFTSamplesParallel(double** pointsArray, const int arraySize, const int numberOfElementsPerCPU, const int procRank);

double* Linearize2DTo1D(double** pointsArray, const int dimensions, const int arraySize);

// Inverse Discrete Fourier Transform
// Second step for a function approximation, used to get an array of points {x,y} of an approximated function
double IDFTSampleParallel(double* fourierCoefficients, const int approximationAccuracy, const int k, const int numberOfPoints);
double* IDFTSamplesParallel(double* fourierCoefficients, const int samplesSize, const int approximationAccuracy,
							const double leftBound, const double rightBound, const int numberOfPoints, const int pointsPerCPU, const int procRank);

#endif // !FOURIER_SERIES_PARALLEL_H
