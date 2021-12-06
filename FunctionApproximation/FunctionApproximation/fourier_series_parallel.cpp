#include "fourier_series_parallel.h"

#include <iostream>
#include <cmath>
#include <mpi.h>

const double PI = 3.14159265358979;
const int dimensions = 2;
const int x = 0; const int y = 1;
const int real = 0; const int im = 1;

double* DFTSampleParallel(double** pointsArray, const int arraySize, const int k){
	double* coefficientsDFT = new double[dimensions];
	coefficientsDFT[x] = 0; coefficientsDFT[y] = 0;

	for (int i = 0; i < arraySize-1; ++i) {
		double realPart = pointsArray[i][y] * std::cos(2 * PI * k * i / static_cast<double>(arraySize));
		double imPart = pointsArray[i][y] * std::sin(2 * PI * k * i / static_cast<double>(arraySize));

		coefficientsDFT[real] += realPart;
		coefficientsDFT[im] += imPart;
	}
	//std::cout << "current k " << k << std::endl;

	coefficientsDFT[real] *= (2.0 / (static_cast<double>(arraySize)));
	coefficientsDFT[im] *= (2.0 / (static_cast<double>(arraySize)));

	return coefficientsDFT;
}

double** DFTSamplesParallel(double** pointsArray, const int arraySize, const int numberOfElementsPerCPU, const int procRank){
	double** fourierCoefficientsPerCPU = new double* [numberOfElementsPerCPU];
	const int offset = procRank * numberOfElementsPerCPU;

	for (int i = 0; i < numberOfElementsPerCPU; ++i) {
		fourierCoefficientsPerCPU[i] = DFTSampleParallel(pointsArray, arraySize, offset+i);
	}


	std::cout << "CPU: " << procRank << " ended DFT calculations\n";

	return fourierCoefficientsPerCPU;
}

double* Linearize2DTo1D(double** pointsArray, const int dimensions, const int arraySize){
	double* linearizedArray = new double[(arraySize * dimensions)];

	for (int i = 0; i < arraySize; ++i) {
		for (int j = 0; j < dimensions; j++) {
			linearizedArray[i * dimensions + j] = pointsArray[i][j];
		}
	}

	return linearizedArray;
}

double IDFTSampleParallel(double* fourierCoefficients, const int approximationAccuracy,
						  const int k, const int numberOfPoints){
	double yCoordinate = 0.0;
	double summ = 0.0;

	for (int i = 1; i <= approximationAccuracy; ++i) {
		summ = fourierCoefficients[i*dimensions + x] * std::cos(i * 2.0 * PI * k / static_cast<double>(numberOfPoints))
			+ fourierCoefficients[i*dimensions + y] * std::sin(i * 2.0 * PI * k / static_cast<double>(numberOfPoints));
		yCoordinate += summ;
	}
	yCoordinate += fourierCoefficients[0] / 2;

	return yCoordinate;
}

double* IDFTSamplesParallel(double* fourierCoefficients, const int samplesSize, const int approximationAccuracy,
							const double leftBound, const double rightBound, const int numberOfPoints, const int pointsPerCPU, const int procRank){
	if (approximationAccuracy > samplesSize) {
		throw std::length_error("Approximation accuracy must not be bigger than sample size!\n");
		return nullptr;
	}
	double* approximatedPointsArray = new double[pointsPerCPU*2];

	double xAxisOffset = (std::abs(rightBound - leftBound)) / static_cast<double>(numberOfPoints);
	double elementsArrayOffset = procRank * static_cast<double>(pointsPerCPU);
	double xCoordinate = leftBound + elementsArrayOffset*xAxisOffset;

	for (int i = 0; i < pointsPerCPU; ++i) {
			approximatedPointsArray[i * dimensions + x] = xCoordinate;
			approximatedPointsArray[i * dimensions + y] = IDFTSampleParallel(fourierCoefficients, approximationAccuracy,
																			 elementsArrayOffset + i, numberOfPoints);
			xCoordinate += xAxisOffset;
	}
		
	std::cout << "CPU: " << procRank << " succesfully approximated points!" << std::endl;

	return approximatedPointsArray;
}
