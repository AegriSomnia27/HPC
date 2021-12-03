#include "fourier_series_serial.h"
#include <cmath>
#include <iostream>

const double PI = 3.14159265358979;
const int dimensions = 2;
const int x = 0; const int y = 1;
const int real = 0; const int im = 1;

double* DFTSample(double** pointsArray, const int arraySize, const int k){
	double* coefficientsDFT = new double[dimensions];
	coefficientsDFT[real] = 0; coefficientsDFT[im] = 0;

	// Find A(i) and B(i) coefficients 
	// A(i) is a real part and B(i) is an imaginary part
	for (int n = 0; n < arraySize-1; ++n) {
		double realPart = pointsArray[n][y] * std::cos(2 * PI * k * n / static_cast<double>(arraySize));
		double imPart = pointsArray[n][y] * std::sin(2 * PI * k * n / static_cast<double>(arraySize));

		coefficientsDFT[real] += realPart;
		coefficientsDFT[im] += imPart;
		
		//std::cout << "	n = " << n << ": (" << coefficientsDFT[real] << ", " << coefficientsDFT[im] << ")\n";
	}

	coefficientsDFT[real] *= (2.0/(static_cast<double>(arraySize)));
	coefficientsDFT[im] *= (2.0/(static_cast<double>(arraySize)));
	
	//std::cout << "  Sum of the above and divided by " << arraySize << "\n";

	return coefficientsDFT;
}

double** DFTSamples(double** pointsArray, const int arraySize, const int approximationAccuracy) {
	// Calculate coefficients for k iteration
	if (approximationAccuracy > arraySize) {
		throw std::length_error("Approximation accuracy must not be bigger than sample size!\n");
		return nullptr;
	}

	double** samples = new double* [approximationAccuracy];

	for (int k = 0; k < approximationAccuracy; ++k) {
		//std::cout << "  k = " << k << "\n";
		samples[k] = DFTSample(pointsArray, approximationAccuracy, k);
		//std::cout << "  X_" << k << " = (" << samples[k][real] << ", " << samples[k][im] << ")\n\n";
	}

	return samples;
}

double IDFTSample(double** fourierCoefficients, const int samplesSize, const int approximationAccuracy,
				  const int k, const int numberOfPoints){
	// Calculate y coordinate for k iteration
	if (approximationAccuracy > samplesSize) {
		throw std::length_error("Approximation accuracy must not be bigger than the size of samples!\n");
		return 0.0;
	}

	double yCoordinate = 0.0;
	double summ = 0.0;

	for (int i = 1; i <= approximationAccuracy; ++i) {
		summ = fourierCoefficients[i][real] * std::cos(i * 2.0 * PI * k / static_cast<double>(numberOfPoints))
			+ fourierCoefficients[i][im] * std::sin(i * 2.0 * PI * k / static_cast<double>(numberOfPoints));
		yCoordinate += summ;
	}
	yCoordinate += fourierCoefficients[0][real] / 2;

	return yCoordinate;
}

double** IDFTSamples(double** fourierCoefficients, const int samplesSize, const int approximationAccuracy, 
					const double leftBound, const double rightBound, const int numberOfPoints){
	// Calculate new coordinates for approximated function using DFT coefficients
	double** approximatePointsArray = new double* [numberOfPoints];

	double offset = (std::abs(rightBound - leftBound)) / static_cast<double>(numberOfPoints);
	double xCoordinate = leftBound;

	for (int k = 0; k < numberOfPoints; ++k) {
		approximatePointsArray[k] = new double[dimensions];

		approximatePointsArray[k][x] = xCoordinate;
		approximatePointsArray[k][y] = IDFTSample(fourierCoefficients, samplesSize, approximationAccuracy, k, numberOfPoints);
		//std::cout << "  Coordinates_" << k << " = (" << approximatePointsArray[k][x] << ", " << approximatePointsArray[k][y] << ")\n\n";

		xCoordinate += offset;
	}

	return approximatePointsArray;
}