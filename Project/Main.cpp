#include "mat.h"
#include <iostream>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h> 
#include <vector> 
#include <algorithm>
#include "readmat.h"
#include <math.h>
#include <valarray>
#include "opencv2/highgui/highgui.hpp"
#include <ctime>
#include "omp.h"

using namespace std;
using namespace cv;

typedef struct {
	size_t X;
	size_t Y;
	size_t U;
	size_t V;


} DataIn;

inline float getelement(DataIn data, int x, int v, int u, int y, uint8_t * __restrict dataElements)
{
	int index = data.X*data.Y*data.U*(v)+data.X*data.Y*(u)+data.X*(y)+x;
	return (float)dataElements[index];
}


int main()
{
	const clock_t begin_time = clock();
	int i;
	double start, end, start1;

	start1 = omp_get_wtime();
	int elements;
	int numberofdimension;

	const char *file = "BraceleU.mat";
	const size_t* dimepointer;
	Mat out;


	readmat thismat(file);
	numberofdimension = thismat.getnumbrofdimensions();
	dimepointer = thismat.dimensionpointer();



	cout << "The dimensions are  ";
	for (int k = 0; k < numberofdimension; k++)
	{
		cout << *(dimepointer + k);
		cout << "  ";
	}
	elements = thismat.numberofelements();

	// Dimensions
	DataIn data;
	data.X = *(dimepointer);
	data.Y = *(dimepointer + 1);
	data.U = *(dimepointer + 2);
	data.V = *(dimepointer + 3);

	uint8_t* dataelements;
	dataelements = thismat.getarraypointer();
	cout << "The dimensions are  \n";

	int size;
	size = data.X*data.Y*data.U*data.V;
	int align = 64;
	uint8_t*  __restrict ReOrderedImage = (uint8_t*)_mm_malloc(size, align);
	uint8_t*  __restrict IntegerShiftedImage = (uint8_t*)_mm_malloc(size, align);
	float* __restrict FinalImage = new float[data.U*data.V];
	start = omp_get_wtime();
	omp_set_dynamic(0);
	omp_set_num_threads(4);
	int X_shift, Y_shift;
	int m;
	m = 0;
#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		for (int u = id; u < data.U; u += num_threads)
		{
			for (int v = 0; v < data.V; v++)
			{
				for (int y = 0; y < data.Y; y++)
				{
					for (int x = 0; x < data.X; x++)
					{
						*(ReOrderedImage + v + (data.V*u) + (data.V*data.U*x) + (data.V*data.U*data.X*y)) = getelement(data, x, v, u, y, dataelements);
					}
				}
			}
		}
	}
	start1 = omp_get_wtime();
	//if (m != 0)
	{
#pragma omp parallel
		{
			int X_shift, Y_shift;
			int id = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
			int RowStartSrc, RowShiftSrc, ImageStart, ImageShift, RowStartDst, RowShiftDst;

			for (int x = id; x < data.X; x+=num_threads)
			{
				X_shift = x*m;
				for (int y = 0; y < data.Y; ++y)
				{
					Y_shift = y*m;
					ImageStart = (data.V*data.U*x) + (data.V*data.U*data.X*y);
					for (int u = 0; u < data.U; ++u)
					{

						RowStartDst = (data.V*u) + (data.V*data.U*x) + (data.V*data.U*data.X*y);
						RowShiftDst = RowStartDst + data.V - Y_shift;
						if (u >= X_shift)
						{
							RowStartSrc = (data.V*(u - X_shift)) + (data.V*data.U*x) + (data.V*data.U*data.X*y);
						}
						else
						{
							RowStartSrc = (data.V*(data.U - (X_shift - u))) + (data.V*data.U*x) + (data.V*data.U*data.X*y);
						}
						RowShiftSrc = RowStartSrc + data.V - Y_shift;

						memcpy(IntegerShiftedImage + RowStartDst, ReOrderedImage + RowShiftSrc, sizeof(uint8_t)*Y_shift);

						memcpy(IntegerShiftedImage + RowStartDst + Y_shift, ReOrderedImage + RowStartSrc, sizeof(uint8_t)*(data.V - Y_shift));
					}

				}
			}
		}
	}


	end = omp_get_wtime();
	int h;
	float Temp;
	for (int u = 0; u < data.U; u++)
	{
		for (int v = 0; v < data.V; v++)
		{
			Temp = 0;
			for (int y = 0; y < data.Y; y++)
			{
				for (int x = 0; x < data.X; x++)
				{
					Temp = Temp + (float)*(IntegerShiftedImage + v + (data.V*u) + (data.V*data.U*x) + (data.V*data.U*data.X*y));
				}
			}
			*(FinalImage + v + (data.V*u)) = Temp/(255*(data.X*data.Y));
		}
	}
	
	out = Mat(data.U, data.V, CV_32FC1, FinalImage); //create an image
	namedWindow("Refocused Image", CV_WINDOW_AUTOSIZE);
	//resizeWindow("Refocused Image", 600, 600);
	imshow("Refocused Image", out); //display the image which is stored in the 'img' in the "MyWindow" window

	waitKey(0);  //wait infinite time for a keyress

	destroyWindow("Refocused Image"); //destroy the window with the name, "MyWindow"*/
	cout << "The end of array \n FINISHED\n";
	for (int k = 0; k < 100; k++)
	{
		cout << unsigned(*(IntegerShiftedImage + k + (data.V*data.U * 0) + (data.V*data.U*data.X * 1)));
		cout << " ";
	}


	cout << "\nSource data\n";
	for (int k = 0; k < 100; k++)
	{
		cout << unsigned(*(ReOrderedImage + k + (data.V*data.U * 0) + (data.V*data.U*data.X * 1)));
		cout << " ";
	}
	for (int k = 0; k < 5; k++)







		cout << endl;
	cout << "\nThis is the number of dimensions  ";
	cout << numberofdimension << endl;

	cout << "Data preparation time " << (end - start);
	cout << "\nInteger shift time " << (end - start1);
	cin >> h;




	//cin >> i;
	return 0;

}

