#include "mat.h"
#include <iostream>
#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <vector> /* For STL */
#include "readmat.h"


using namespace std;

int main()
{
	int i;
	int m;
	int elements;
	int numberofdimension;
	double* arraypint;
	const char *file = "Example.mat";
	const size_t* dimepointer;

	readmat thismat(file);
	numberofdimension = thismat.getnumbrofdimensions();
	dimepointer = thismat.dimensionpointer();
	cout << "The dimensions are  ";
	for (m = 0; m < numberofdimension; m++)
	{
		cout << *(dimepointer + m);
		cout << "  ";
	}
	cout << endl;
	cout << "This is the number of dimensions  ";
	cout << numberofdimension << endl;
	arraypint = thismat.getarraypointer();
	elements = thismat.numberofelements();
	cout << "The array elements ";
	cout << elements;
	cout << endl;
	
	for (int n = 0; n < elements; n++)
	{
	cout << *(arraypint + n);
	cout << " ";
	}

	cin >> i;
	return 0;

}

