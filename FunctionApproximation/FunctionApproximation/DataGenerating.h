#ifndef DATA_GENERATING_H
#define DATA_GENERATING_H

inline double Function(const double x);
double** GenerateData(const int arraySize, const double leftBound=-10.0, const double rightBound=10.0, 
					  const bool isNoisy = true, const double mean = 0.0, const double sigma = 5.0);
void DeleteArray2D(double** pointsArray, const int arraySize);

void MakeCSV(double** pointsArray, const int arraySize, const char* fileName="test_file.csv", 
			 const char* fileHeader="x,y\n");

#endif // !DATA_GENERATING_H
