#include "data_generation.h"
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>


const int dimensions = 2;
const int x = 0;
const int y = 1;


inline double Function(double x) {
	// non-periodic function
	// 3x^3 + 7x^2 + 5x + 1.1
	//return 3 * pow(x, 3) + 7 * pow(x, 2) + 5 * x + 1.1;
	
	// piecewise function
	if (x < -5) {
		return -2;
	}
	else if (x >= -5 && x < 5)
		return 2;
	else
		return -2;

	// periodic function
	//sin(5x) + cos(x)
	//return (std::sin(5*x)+std::cos(x));
}

double** GenerateData(const int arraySize, const double leftBound, const double rightBound,
					  const bool isNoisy, const double mean, const double sigma){
	// initializing RNG if we need to add some noise to a function
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

	pointsArray = nullptr;
}

void MakeCSV_2D(double** pointsArray, const int arraySize, const char* fileName, const char* fileHeader){
	// this function is used for generating a csv file based on an array of points
	// in my case it was used for data visualization with python matplotlib 

	// creates a CSV file and opens it
	std::ofstream fileCSV;
	fileCSV.open(fileName);

	// raises an error if the file cannot be oppened
	if (!fileCSV.is_open()) {
		std::cerr << "File cannot be opened";
		return;
	}

	// writes a header
	fileCSV << fileHeader;

	//writes data
	for (int i = 0; i < arraySize; ++i) {
		fileCSV << pointsArray[i][x] << ',' << pointsArray[i][y] << '\n';
	}

	fileCSV.close();
}

void MakeCSV_1D(double* pointsArray, const int arraySize, const char* fileName, const char* fileHeader){
	// this function is used for generating a csv file based on an array of points
	// in my case it was used for data visualization with python matplotlib 

	// creates a CSV file and opens it
	std::ofstream fileCSV;
	fileCSV.open(fileName);

	// raises an error if the file cannot be oppened
	if (!fileCSV.is_open()) {
		std::cerr << "File cannot be opened";
		return;
	}

	// writes a header
	fileCSV << fileHeader;

	//writes data
	for (int i = 0; i < arraySize; ++i) {
		for (int j = 0; j < dimensions; j++) {
			fileCSV << pointsArray[i * dimensions + j];
			if (j == 0) {
				fileCSV << ",";
			}
		}
		fileCSV << std::endl;
	}

	fileCSV.close();
}
