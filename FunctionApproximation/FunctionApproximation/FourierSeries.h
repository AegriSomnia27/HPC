#ifndef FOURIER_SERIES_H
#define FOURIER_SERIES_H

// Discrete Fourier Transform
// First step for a function approximation, used to get DFT coefficients 
// returns an array of DFT coefficients
double* DFTSample(double** pointsArray, const int arraySize, const int k); // one should not use that func outside of DFTSamples()
double** DFTSamples(double** pointsArray, const int arraySize);

// Inverse Discrete Fourier Transform
// Second step for a function approximation, used to get an array of points {x,y} of an approximated function
double IDFTSample(double** samples, const int samplesSize, const int approximationAccuracy,
				  const int k, const int numberOfPoints);
double** IDFTSamples(double** samples, const int samplesSize, const int approximationAccuracy,
					 const double leftBound, const double rightBound, const int numberOfPoints);

#endif // !FOURIER_SERIES_H