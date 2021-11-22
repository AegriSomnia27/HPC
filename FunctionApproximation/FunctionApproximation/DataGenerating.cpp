#include "dataGenerating.h"
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>


const int dimensions = 2;
const int x = 0;
const int y = 1;


inline double Function(double x) {
	// 3x^3 + 7x^2 + 5x + 1.1
	//return 3 * pow(x, 3) + 7 * pow(x, 2) + 5 * x + 1.1;
	return (std::sin(5*x)+std::cos(x));
}

double** GenerateData(const int arraySize, const double leftBound, const double rightBound,
					  const bool isNoisy, const double mean, const double sigma){
	std::random_device randomDevice;
	std::mt19937 gen(randomDevice());
	std::normal_distribution<double> noise(mean, sigma);

	double** pointsArray = new double* [arraySize];

	double offset = (std::abs(rightBound - leftBound)) / static_cast<double>(arraySize);
	double point = leftBound;

	for (int i = 0; i < arraySize; ++i) {
		pointsArray[i] = new double[dimensions];

		pointsArray[i][x] = point;
		if (isNoisy) {
			pointsArray[i][y] = Function(point) + noise(gen);
		}
		else {
			pointsArray[i][y] = Function(point);
		}

		point += offset;
	}

	return pointsArray;
}

void DeleteArray2D(double** pointsArray, const int arraySize) {
	for (int i = 0; i < arraySize; ++i) {
		delete[] pointsArray[i];
	}
	delete[] pointsArray;
}

void MakeCSV(double** pointsArray, const int arraySize, const char* fileName, const char* fileHeader){
	// create a CSV file and open it
	std::ofstream fileCSV;
	fileCSV.open(fileName);

	// raise an error if the file cannot be oppened
	if (!fileCSV.is_open()) {
		std::cerr << "File cannot be opened";
		return;
	}

	// write a header
	fileCSV << fileHeader;

	//write data
	for (int i = 0; i < arraySize; ++i) {
		fileCSV << pointsArray[i][x] << ',' << pointsArray[i][y] << '\n';
	}

	fileCSV.close();
}
