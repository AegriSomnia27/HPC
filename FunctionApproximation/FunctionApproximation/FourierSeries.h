#ifndef FOURIER_SERIES_H
#define FOURIER_SERIES_H

double* DFTSample(double** pointsArray, const int arraySize, const int numberOfIterations = 20);
double** DFTSamples(double** pointsArray, const int arraySize);
double* IDFTSample(double** samples, const int arraySize);
double** IDFTSamples(double** samples, const int arraySize);

double** ApproximateFourierSeries(double** pointsArray, const double** fourierCoefficients);

#endif // !FOURIER_SERIES_H