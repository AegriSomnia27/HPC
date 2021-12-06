#ifndef DATA_GENERATION_H
#define DATA_GENERATION_H


inline double Function(const double x); 
double** GenerateData(const int arraySize, const double leftBound=-10.0, const double rightBound=10.0, 
					  const bool isNoisy = false, const double mean = 0.0, const double sigma = 1.0);

void DeleteArray2D(double** pointsArray, const int arraySize);

void MakeCSV_2D(double** pointsArray, const int arraySize, const char* fileName="test_file.csv", 
			 const char* fileHeader="x,y\n");

void MakeCSV_1D(double* pointsArray, const int arraySize, const char* fileName = "test_file.csv",
				const char* fileHeader = "x,y\n");


#endif // !DATA_GENERATION_H