#include "data_generation.h"
#include "fourier_series_serial.h"
#include "fourier_series_parallel.h"

#include "mpi.h"

#include <iostream>
#include <typeinfo>
#include <string>
#include <cmath>

int main(int argc, char** argv) {
	const int arraySize = 200;
	const int approximatedArraySize = 100;
	const int approximationAccuracy = 50;
	const double leftBound = -10.0;
	const double rightBound = 10.0;

	// Serial parth of the program
	// Generating an array of {x,y}
	double** pointsArray = GenerateData(arraySize, leftBound, rightBound, false);
	MakeCSV_2D(pointsArray, arraySize, "generated_points.csv");

	// Generating an array of DFT coefficients {real, im}
	double** fourierSeriesArray = DFTSamples(pointsArray, arraySize, arraySize);
	MakeCSV_2D(fourierSeriesArray, arraySize, "fourier_coefficients_serial.csv", "real,im\n");

	// Generating a points array {x,y} of the approximated func
	double** approximatedPoints = IDFTSamples(fourierSeriesArray, arraySize, approximationAccuracy, leftBound, rightBound, approximatedArraySize);
	MakeCSV_2D(approximatedPoints, approximatedArraySize, "approximated_function_points_serial.csv");

	// Free allocated memory
	DeleteArray2D(fourierSeriesArray, arraySize);
	DeleteArray2D(approximatedPoints, approximatedArraySize);

	// TODO: add parallel algorithm support using OpenMPI 
	int procNum, procRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Status status;


	// How many elements a single CPU must handle
	const int numberOfElementsPerCPU = arraySize / procNum;
	const int elementsRemaind = arraySize - numberOfElementsPerCPU;
	std::cout << "processing unit " << procRank << " number of elements " << numberOfElementsPerCPU << std::endl;
	

	// DFT of the data array
	double** procFourierSeries = DFTSamplesParallel(pointsArray, arraySize, numberOfElementsPerCPU, procRank);
	double* linearizedArray = Linearize2DTo1D(procFourierSeries, 2, numberOfElementsPerCPU);

	// Gather all data into a single array
	//double* DFTCoefficientsArray = new double[(static_cast<double>(numberOfElementsPerCPU)*procNum*2)];
	double* DFTCoefficientsArray = new double[numberOfElementsPerCPU*procNum*2];
	MPI_Allgather(linearizedArray, numberOfElementsPerCPU*2, MPI_DOUBLE, DFTCoefficientsArray, numberOfElementsPerCPU*2, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (procRank==0) {
		// Create one copy of CSV 
		MakeCSV_1D(DFTCoefficientsArray, (numberOfElementsPerCPU * procNum), "fourier_coefficients_parallel.csv", "real,im\n");
	}

	// IDFT part
	const int pointsPerCPU = approximatedArraySize / procNum;
	
	double* procApproximatedPointsArray = IDFTSamplesParallel(DFTCoefficientsArray, numberOfElementsPerCPU*procNum*2, approximationAccuracy,
		leftBound, rightBound, approximatedArraySize, pointsPerCPU, procRank);

	// Gather the data and make CSV

	if (procRank == 0) {
		double* approximatedPointsParallel = new double[pointsPerCPU * procNum*2];
		MPI_Gather(procApproximatedPointsArray, pointsPerCPU*2, MPI_DOUBLE, approximatedPointsParallel, pointsPerCPU*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MakeCSV_1D(approximatedPointsParallel, (pointsPerCPU*procNum), "approximated_function_points_parallel.csv");
		delete[] approximatedPointsParallel;
	}

	DeleteArray2D(procFourierSeries, numberOfElementsPerCPU);
	delete[] linearizedArray;
	delete[] DFTCoefficientsArray;

	MPI_Finalize();

	return 0;
}