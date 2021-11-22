#include "DataGenerating.h"
#include "FourierSeries.h"

int main() {
	int arraySize = 10;
	int testSize = 4;

	//double** testArray = new double* [testSize];
	//for (int i = 0; i < testSize; ++i) {
	//	testArray[i] = new double[2];
	//	testArray[i][1] = 0;
	//}

	//testArray[0][0] = 1.0;
	//testArray[1][0] = -0.0;
	//testArray[2][0] = -1.0;
	//testArray[3][0] = 0.0;

	double** pointsArray = GenerateData(arraySize, -10.0, 10.0, false);
	MakeCSV(pointsArray, arraySize, "first_try.csv");

	double** fourierSeriesArray = DFTSamples(pointsArray, arraySize);
	//double** fourierSeriesArray = DFTSamples(testArray, 4);

	MakeCSV(fourierSeriesArray, arraySize, "testFourier.csv", "real,im\n");

	return 0;
}