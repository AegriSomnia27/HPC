#include "data_generation.h"
#include "fourier_series_serial.h"

#include "mpi.h"
#include <iostream>
#include <typeinfo>

int main(int argc, char** argv) {
	const int arraySize = 200;
	const int approximatedArraySize = 100;
	const int approximationAccuracy = 50;
	const double leftBound = -10.0;
	const double rightBound = 10.0;

	// Serial parth of the program
	// Generating an array of {x,y}
	double** pointsArray = GenerateData(arraySize, leftBound, rightBound, false);
	MakeCSV(pointsArray, arraySize, "generated_points.csv");

	// Generating an array of DFT coefficients {real, im}
	double** fourierSeriesArray = DFTSamples(pointsArray, arraySize, arraySize);
	MakeCSV(fourierSeriesArray, arraySize, "fourier_coefficients_serial.csv", "real,im\n");

	// Generating a points array {x,y} of the approximated func
	double** approximatedPoints = IDFTSamples(fourierSeriesArray, arraySize, approximationAccuracy, leftBound, rightBound, approximatedArraySize);
	MakeCSV(approximatedPoints, approximatedArraySize, "approximated_function_points_serial.csv");

	// Free allocated memory
	DeleteArray2D(fourierSeriesArray, arraySize);
	DeleteArray2D(approximatedPoints, approximatedArraySize);

	// TODO: add parallel algorithm support using OpenMPI 
	int procNum, procRank;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Status status;
	
	MPI_Finalize();

	return 0;
}