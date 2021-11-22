#include "FourierSeries.h"
#include <cmath>
#include <iostream>

const double PI = 3.14159265358979;
const int dimensions = 2;
const int x = 0;
const int y = 1;
const int real = 0;
const int im = 1;

double* DFTSample(double** pointsArray, const int arraySize, const int k){
	double* coefficientsDFT = new double[dimensions];
	coefficientsDFT[real] = 0; coefficientsDFT[im] = 0;

	double standartSample = 2 * PI * (pointsArray[k][x] - pointsArray[0][x]) / static_cast<double>(arraySize);

	for (int n = 0; n < arraySize-1; ++n) {
		double realPart = pointsArray[n][y] * std::cos(2 * PI * k * n / arraySize);
		double imPart = pointsArray[n][y] * std::sin(2 * PI * k * n / arraySize);

		coefficientsDFT[real] += realPart;
		coefficientsDFT[im] += imPart;
		
		std::cout << "	n = " << n << ": (" << coefficientsDFT[real] << ", " << coefficientsDFT[im] << ")\n";
	}

	coefficientsDFT[real] /= (static_cast<double>(arraySize));
	coefficientsDFT[im] /= (static_cast<double>(arraySize));
	
	std::cout << "  Sum of the above and divided by " << arraySize << "\n";

	return coefficientsDFT;
}

double** DFTSamples(double** pointsArray, const int arraySize) {
	double** samples = new double* [arraySize];

	for (int k = 0; k < arraySize; ++k) {
		std::cout << "  k = " << k << "\n";
		samples[k] = DFTSample(pointsArray, arraySize, k);
		std::cout << "  X_" << k << " = (" << samples[k][real] << ", " << samples[k][im] << ")\n\n";
	}

	return samples;
}

double* IDFTSample(double** samples, const int arraySize){

}

double** IDFTSamples(double** samples, const int arraySize){

}

double** ApproximateFourierSeries(double** pointsArray, double** fourierCoefficients, const int arraySize){
	return nullptr;
}