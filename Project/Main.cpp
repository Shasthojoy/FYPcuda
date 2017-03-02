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
#include <errno.h>  
#include <immintrin.h>
using namespace std;
using namespace cv;

typedef struct {
	size_t X;
	size_t Y;
	size_t U;
	size_t V;
	size_t KernalNoOfFilters;
	size_t KernalFilterLength;
	
} DataIn;

inline uint8_t getelement(DataIn data, int x, int v, int u, int y, uint8_t * __restrict dataElements)
{
	int index = data.X*data.Y*data.U*(v)+data.X*data.Y*(u)+data.X*(y)+x;
	return dataElements[index];
}


int main()
{
	const clock_t begin_time = clock();
	int i;
	double start, end, start1;
	errno_t err;
	start1 = omp_get_wtime();
	int elements;
	int numberofdimension;

	const char *fileImage = "BraceleU.mat";
	const size_t* dimepointer;
	Mat out;
	const char *Filter = "filter.mat";

	readmat* thismat = new readmat(fileImage);
	readmat filter(Filter);



	numberofdimension = filter.getnumbrofdimensions();
	dimepointer = thismat->dimensionpointer();
	const size_t * FilterDimensions;
	FilterDimensions = filter.dimensionpointer();
	float* FilterCoefficients = (float*)filter.getarraypointer();
	size_t FilterSizes[7] = { 1, 49, 11, 31, 17, 19, 7 };
	cout << "The filters are  ";
	for (int k = 0; k < 7; k++)
	{
		cout << *(FilterCoefficients + k + 49);
		cout << "  ";
	}
	elements = thismat->numberofelements();

	// Dimensions
	DataIn data;
	data.KernalFilterLength = *(FilterDimensions);
	data.KernalNoOfFilters = *(FilterDimensions + 1);
	data.X = *(dimepointer);
	data.Y = *(dimepointer + 1);
	data.U = *(dimepointer + 2);
	data.V = *(dimepointer + 3);
	int size;
	size = (data.X*data.Y*data.U*data.V * 4);
	int align = 32;
	//float*  __restrict IntegerShiftedImage = (float*)malloc(size);;
	_set_errno(0);
	float*  __restrict IntegerShiftedImage = (float*)_aligned_malloc(size, align);
	if (IntegerShiftedImage == NULL)
	{
		_get_errno(&err);
		cout << "\nThe aligned memory allocation failed\n" << err;
	}
	uint8_t* dataelements;
	dataelements = (uint8_t*)thismat->getarraypointer();
	cout << "The dimensions are  \n";




	float*  __restrict ReOrderedImage = new float[data.X*data.Y*data.U*data.V];
	//float*  __restrict IntegerShiftedImage = (float*)_mm_malloc(size, align);
	//float*  __restrict IntegerShiftedImage = new float[data.X*data.Y*data.U*data.V];
	//delete thismat;

	float* __restrict FinalImage = new float[data.U*data.V];
	start = omp_get_wtime();
	omp_set_dynamic(0);
	omp_set_num_threads(4);
	int X_shift, Y_shift;
	int m;
	m = 1;
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
	int k = 0;
	delete thismat;

	start1 = omp_get_wtime();
	//if (m != 0)
	{
#pragma omp parallel
		{
			int X_shift, Y_shift;
			int id = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
			int RowStartSrc, RowShiftSrc, ImageStart, ImageShift, RowStartDst, RowShiftDst;

			for (int x = id; x < data.X; x += num_threads)
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

						memcpy(IntegerShiftedImage + RowStartDst, ReOrderedImage + RowShiftSrc, sizeof(float)*Y_shift);

						memcpy(IntegerShiftedImage + RowStartDst + Y_shift, ReOrderedImage + RowStartSrc, sizeof(float)*(data.V - Y_shift));
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
			*(FinalImage + v + (data.V*u)) = Temp / (255 * (data.X*data.Y));
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

void ConvolutionAVX(DataIn Data, float* FilterKernal, int* FilterSize, float *out, float *in,int Delay)
{
	__m256 filterCoef;
	__m256 Accumilation[6], InputDataLeft, InputDataRight,SymmetricAdd;
	int PointerFilterCoefficient;
	int offset, LoadPtr, RightPtr;
	for (int u = 0; u < Data.U; u++)
	{
		for (int v = 0; v < Data.V; v++)
		{
			offset = Data.V*u + v;
			for (int k = 0; k < Data.KernalNoOfFilters - 1; k++)
			{
				Accumilation[k] = _mm256_setzero_ps();
				PointerFilterCoefficient = k + 1;
				for (int l = 0; l < ((*(FilterSize + PointerFilterCoefficient) - 1) / 2) + 1; l++)
				{
					if ((v >= (*(FilterSize + PointerFilterCoefficient) - 1) / 2) && (v <= ((Data.V - 1) - (*(FilterSize + PointerFilterCoefficient) - 1) / 2)))
					{
						LoadPtr = ((*(FilterSize + PointerFilterCoefficient) - 1) / 2) - l;
						filterCoef = _mm256_broadcast_ss(FilterKernal + ((Data.KernalFilterLength*PointerFilterCoefficient) + 1));
						InputDataLeft = _mm256_loadu_ps((in + offset - LoadPtr));

						if (l == ((*(FilterSize + PointerFilterCoefficient) - 1) / 2))
						{
							Accumilation[k] = _mm256_mul_ps(InputDataLeft, filterCoef);
						}
						else
						{
							InputDataRight = _mm256_loadu_ps((in + offset + LoadPtr));
							SymmetricAdd = _mm256_add_ps(InputDataLeft, InputDataRight);
							Accumilation[k] = _mm256_mul_ps(SymmetricAdd, filterCoef);

						}
					}
					else
					{

					}
					
					
				}

					
			}

		}
	}


}
