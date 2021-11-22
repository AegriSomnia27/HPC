#include "DataGenerating.h"
#include "FourierSeries.h"

int main() {
	const int arraySize = 200;
	const int approximatedArraySize = 100;
	const int approximationAccuracy = 50;
	const double leftBound = -10.0;
	const double rightBound = 10.0;

	// Generating an array of {x,y}
	double** pointsArray = GenerateData(arraySize, leftBound, rightBound, false);
	MakeCSV(pointsArray, arraySize, "generated_points.csv");

	// Generating an array of DFT coefficients {real, im}
	double** fourierSeriesArray = DFTSamples(pointsArray, arraySize);
	MakeCSV(fourierSeriesArray, arraySize, "fourier_coefficients.csv", "real,im\n");

	// Generating a points array {x,y} of the approximated func
	double** approximatedPoints = IDFTSamples(fourierSeriesArray, arraySize, approximationAccuracy, leftBound, rightBound, approximatedArraySize);
	MakeCSV(approximatedPoints, approximatedArraySize, "approximated_function_points.csv");

	// Free allocated memory
	DeleteArray2D(pointsArray, arraySize);
	DeleteArray2D(fourierSeriesArray, arraySize);
	DeleteArray2D(approximatedPoints, approximatedArraySize);

	return 0;
}