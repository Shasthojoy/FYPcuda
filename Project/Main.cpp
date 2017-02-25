#include "mat.h"
#include <iostream>
#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <vector> /* For STL */
#include <algorithm>
#include "readmat.h"
#include <valarray>
#include "opencv2/highgui/highgui.hpp"
#include <ctime>
using namespace std;
using namespace cv;
int main()
{
	const clock_t begin_time = clock();
	int i;
	int m;
	int elements;
	int numberofdimension;
	double* arraypint;
	const char *file = "Bracelet.mat";
	const size_t* dimepointer;
	Mat out;
	
	try{
		readmat thismat(file);
		numberofdimension = thismat.getnumbrofdimensions();
		dimepointer = thismat.dimensionpointer();



		cout << "The dimensions are  ";
		for (m = 0; m < numberofdimension; m++)
		{
			cout << *(dimepointer + m);
			cout << "  ";
		}

		// Dimensions

		size_t X, Y, U, V;
		X = *(dimepointer);
		Y = *(dimepointer + 1);
		U = *(dimepointer + 2);
		V = *(dimepointer + 3);
		// Dimensions end

	
	
		

	
	
	




		
	}
	catch (exception& e)
	{
		cout << "ERROR:"<< e.what() << endl;
	}


	cin >> i;
	return 0;

}

