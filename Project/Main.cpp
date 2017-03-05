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
#include <math.h>
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
void ConvolutionAVX(DataIn Data, float* FilterKernal, size_t* FilterSize, float *out, float *in, float* Delay);
void ColumnConvilutionIntegerShiftAVX(DataIn Data, float ShiftAmount, float* __restrict In, float* __restrict Out, size_t* FilterSize, float* __restrict FilterKernal);
int main()
{
	const clock_t begin_time = clock();
	int i;
	double start, end1, start1, start2, end2, end3;
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
	size_t FilterSizes[7] = { 1, 49, 23, 31, 17, 19, 7 };
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
	//float*  __restrict IntegerShiftedImage = (float*)_aligned_malloc(size, align);
	float*  __restrict IntegerShiftedImage = new float[data.X*data.Y*data.U*data.V];
	if (IntegerShiftedImage == NULL)
	{
		_get_errno(&err);
		cout << "\nThe aligned memory allocation failed\n" << err;
	}
	uint8_t* dataelements;
	dataelements = (uint8_t*)thismat->getarraypointer();
	cout << "The dimensions are  \n" << data.V << endl;




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
	int k = 0;
	delete thismat;
	end3 = omp_get_wtime();
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
	end1 = omp_get_wtime();
	// Fractional Shift CPU
	float* FractionalShifted = new float[data.U*data.V*data.Y];
	float* Delay = new float[1];
	*Delay = 0.4;

	start2 = omp_get_wtime();
	//ConvolutionAVX(data, FilterCoefficients, FilterSizes, FractionalShifted, IntegerShiftedImage, Delay);
	end2 = omp_get_wtime();
	// void ConvolutionAVX(DataIn Data, float* FilterKernal, size_t* FilterSize, float *out, float *in, float Delay)
	//(DataIn Data, float ShiftAmount, float* __restrict In, float* __restrict Out, int* FilterSize, float* __restrict FilterKernal)
	ColumnConvilutionIntegerShiftAVX(data, 1.4, ReOrderedImage, FractionalShifted, FilterSizes, FilterCoefficients);
	cout << "\n The data is good\n";
	for (int g = 0; g < 50; g++)
	{
		cout << *(FractionalShifted + g) << " ";
	}


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

	//out = Mat(data.U, data.V,CV_32FC1, FinalImage); //create an image
	//namedWindow("Refocused Image", CV_WINDOW_AUTOSIZE);
	//resizeWindow("Refocused Image", 600, 600);
//	imshow("Refocused Image", out); //display the image which is stored in the 'img' in the "MyWindow" window

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

	cout << "Data preparation time " << (end3 - start);
	cout << "\nInteger shift time " << (end1 - start1);
	cout << "\nFractional shift for time " << (end2 - start2);
	cin >> h;




	//cin >> i;
	return 0;

}

void ConvolutionAVX(DataIn Data, float* FilterKernal, size_t* FilterSize, float *out, float *in, float* Delay)
{

#pragma omp parallel
	{
		int X_shift, Y_shift;
		int id = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		__m256 filterCoef, DelayReg, DelayRegPower;
		__m256 Accumilation[7], InputDataLeft, InputDataRight, SymmetricAdd;
		__m256 ShiftEdge, ShiftEdgeReverse;
		int PointerFilterCoefficient;
		int offset, LoadPtr, RightPtr;
		float *EdgeRight, *EdgeLeft;
		EdgeLeft = new float[16];
		EdgeRight = new float[16];
		int LoadTemp;
		float FilterCoefCheck;
		float FilterNoCheck;
		DelayReg = _mm256_broadcast_ss(Delay);
		for (int u = id; u < Data.U; u += num_threads)
		{
			//To manage the edges
			ShiftEdge = _mm256_loadu_ps(in + (Data.V*u));
			ShiftEdgeReverse = _mm256_permute_ps(ShiftEdge, 27);
			ShiftEdgeReverse = _mm256_permute2f128_ps(ShiftEdgeReverse, ShiftEdgeReverse, 1);
			_mm256_storeu_ps(EdgeLeft, ShiftEdgeReverse);
			_mm256_storeu_ps(EdgeLeft + 8, ShiftEdge);
			ShiftEdge = _mm256_loadu_ps(in + (Data.V*(u + 1)) - 8);
			ShiftEdgeReverse = _mm256_permute_ps(ShiftEdge, 27);
			ShiftEdgeReverse = _mm256_permute2f128_ps(ShiftEdgeReverse, ShiftEdgeReverse, 1);
			_mm256_storeu_ps(EdgeRight, ShiftEdge);
			_mm256_storeu_ps(EdgeRight + 8, ShiftEdgeReverse);

			for (int v = 0; v < Data.V; v += 8)
			{
				offset = Data.V*u + v;
				Accumilation[0] = _mm256_loadu_ps((in + offset));
				DelayRegPower = _mm256_broadcast_ss(Delay);
				for (int k = 1; k < Data.KernalNoOfFilters; k++)
				{
					Accumilation[k] = _mm256_setzero_ps();
					PointerFilterCoefficient = k;
					FilterNoCheck = k % 2;
					for (int l = 0; l < ((*(FilterSize + PointerFilterCoefficient) - 1) / 2) + 1; l++)
					{
						//if ((v >= (*(FilterSize + PointerFilterCoefficient) - 1) / 2) && (v <= ((Data.V - 1) - (*(FilterSize + PointerFilterCoefficient) - 1) / 2)))

						LoadPtr = ((*(FilterSize + PointerFilterCoefficient) - 1) / 2) - l;
						filterCoef = _mm256_broadcast_ss(FilterKernal + ((Data.KernalFilterLength*PointerFilterCoefficient) + l));
						FilterCoefCheck = *(FilterKernal + ((Data.KernalFilterLength*PointerFilterCoefficient) + l));

						if ((v >= LoadPtr))// && )
						{
							InputDataLeft = _mm256_loadu_ps((in + offset - LoadPtr));
						}
						else
						{
							LoadTemp = LoadPtr - v;

							if (LoadTemp <= 7)
								InputDataLeft = _mm256_loadu_ps(EdgeLeft + 8 - LoadTemp);
							else
							{
								LoadTemp = LoadTemp - 7;
								InputDataLeft = _mm256_loadu_ps((in + (Data.V*u)) + LoadTemp - 1);
								InputDataLeft = _mm256_permute_ps(InputDataLeft, 27);
								InputDataLeft = _mm256_permute2f128_ps(InputDataLeft, InputDataLeft, 1);

							}

						}

						if (l == ((*(FilterSize + PointerFilterCoefficient) - 1) / 2))
						{
							Accumilation[k] = _mm256_add_ps(_mm256_mul_ps(InputDataLeft, filterCoef), Accumilation[k]);
						}
						else
						{
							if (((v + 7) <= ((Data.V - 1) - LoadPtr)))
							{
								InputDataRight = _mm256_loadu_ps((in + offset + LoadPtr));
							}
							else
							{
								LoadTemp = v - (Data.V - 1) + LoadPtr;
								if (LoadTemp <= 0)
								{
									InputDataRight = _mm256_loadu_ps(EdgeRight + 7 + LoadTemp);
								}
								else
								{
									InputDataRight = _mm256_loadu_ps(in + (Data.V*(u + 1) - 8) - (LoadTemp - 1));
									InputDataRight = _mm256_permute_ps(InputDataRight, 27);
									InputDataRight = _mm256_permute2f128_ps(InputDataRight, InputDataRight, 1);
								}
							}
							if (FilterNoCheck == 1)
							{
								SymmetricAdd = _mm256_sub_ps(InputDataRight, InputDataLeft);
							}
							else
							{
								SymmetricAdd = _mm256_add_ps(InputDataRight, InputDataLeft);
							}
							Accumilation[k] = _mm256_add_ps(_mm256_mul_ps(SymmetricAdd, filterCoef), Accumilation[k]);

						}


					}
					Accumilation[k] = _mm256_mul_ps(DelayRegPower, Accumilation[k]);
					DelayRegPower = _mm256_mul_ps(DelayRegPower, DelayReg);
					Accumilation[0] = _mm256_add_ps(Accumilation[k], Accumilation[0]);
				}
				_mm256_storeu_ps(out + offset, Accumilation[0]);
			}
		}

	}
}
inline float* PickRow(size_t u, size_t v, float* ImageStart, int IntegerShift, int U, int V)
{
	if (u >= IntegerShift)
	{
		return ImageStart + V*(u - IntegerShift) + v;
	}
	else
	{
		return ImageStart + (V*(U - (IntegerShift - u))) + v;
	}
}
void ColumnConvilutionIntegerShiftAVX(DataIn Data, float ShiftAmount, float* __restrict In, float* __restrict Out, size_t* FilterSize, float* __restrict FilterKernal)
{
	
	int ImageIntergerShift;
	float ImageShiftAmount;
	float* ImageDelay= new float[1];
	float* LightFieldSizeInverse = new float[1];
	*LightFieldSizeInverse = 1 / 17;
	cout << "tttttt " << *LightFieldSizeInverse << endl;
	float* ImageStart;
	__m256 filterCoef, DelayReg, DelayRegPower;
	__m256 Accumilation[7], InputDataLeft, InputDataRight, SymmetricAdd;
	__m256 ColumnSumAccumilator,Divider;
	Divider = _mm256_broadcast_ss(LightFieldSizeInverse);
	int MaxFilterOrder = 25;
	//Check if Filter is Odd or Even
	int FilterNoCheck;
	//size of each filter
	int ThisFilterSize;
	//Calculate Row for the Convolution
	int LoadLeft;
	int LoadRight;
	int ConvolutionFilterIndex;
	/*Here in this code y iterates through image columns and x iterates through image Rows.*/
	
	{
		for (int y = 0; y < Data.Y; y++)
		{
			for (int v = 0; v < Data.V; v += 8)
			{
				for (int u = 0; u < Data.U; u++)
				{
					ColumnSumAccumilator = _mm256_setzero_ps();
					for (int x = 0; x < Data.X; x++)
					{
						ImageShiftAmount = ShiftAmount*x;
						ImageIntergerShift = round(ImageShiftAmount);
						*ImageDelay = ImageShiftAmount - ImageIntergerShift;
						ImageStart = In + (Data.V*Data.U*x) + (Data.V*Data.U*Data.X*y);
						DelayReg = _mm256_broadcast_ss(ImageDelay);
						DelayRegPower = _mm256_broadcast_ss(ImageDelay);
						Accumilation[0] = _mm256_loadu_ps(PickRow(u, v, ImageStart, ImageIntergerShift, Data.U, Data.V));
						if (*ImageDelay != 0)
						{
							for (int l = 0; l < MaxFilterOrder; l++)
							{

								LoadLeft = u - l;
								LoadRight = u + l;
								if ((LoadLeft >= 0))
								{
									InputDataLeft = _mm256_loadu_ps(PickRow(LoadLeft, v, ImageStart, ImageIntergerShift, Data.U, Data.V));
								}
								else
								{
									LoadLeft = u - LoadLeft - 1;
									InputDataLeft = _mm256_loadu_ps(PickRow(LoadLeft, v, ImageStart, ImageIntergerShift, Data.U, Data.V));
								}
								if (l != 0)
								{
									if (LoadRight <= Data.U - 1)
									{
										InputDataRight = _mm256_loadu_ps(PickRow(LoadRight, v, ImageStart, ImageIntergerShift, Data.U, Data.V));
									}
									else
									{
										LoadRight = Data.U - (LoadRight - Data.U - 1);
										InputDataRight = _mm256_loadu_ps(PickRow(LoadRight, v, ImageStart, ImageIntergerShift, Data.U, Data.V));
									}
								}
								for (int k = 1; k < Data.KernalNoOfFilters; k++)
								{
									if (l == 0)
									{
										Accumilation[k] = _mm256_setzero_ps();
									}

									ThisFilterSize = (*(FilterSize + k) - 1) / 2;
									if (l <= ThisFilterSize)
									{
										filterCoef = _mm256_broadcast_ss(FilterKernal + ((Data.KernalFilterLength*k) + ThisFilterSize - l));
										FilterNoCheck = k % 2;
										if (l != 0)
										{
											if (FilterNoCheck == 1)
											{
												SymmetricAdd = _mm256_sub_ps(InputDataRight, InputDataLeft);
											}
											else
											{
												SymmetricAdd = _mm256_add_ps(InputDataRight, InputDataLeft);
											}
											Accumilation[k] = _mm256_add_ps(_mm256_mul_ps(SymmetricAdd, filterCoef), Accumilation[k]);
										}
										else
										{
											Accumilation[k] = _mm256_add_ps(_mm256_mul_ps(InputDataLeft, filterCoef), Accumilation[k]);
										}
										

									}
									if (l == MaxFilterOrder - 1)
									{
										Accumilation[k] = _mm256_mul_ps(DelayRegPower, Accumilation[k]);
										DelayRegPower = _mm256_mul_ps(DelayRegPower, DelayReg);
										Accumilation[0] = _mm256_add_ps(Accumilation[k], Accumilation[0]);
									}
								}
							}
						}
						ColumnSumAccumilator = _mm256_add_ps(ColumnSumAccumilator, Accumilation[0]);
					}
					ColumnSumAccumilator = _mm256_mul_ps(Divider, ColumnSumAccumilator);
					_mm256_storeu_ps(Out + v + (Data.V*u) + (Data.V*Data.U*y), ColumnSumAccumilator);
					if ((v + (Data.V*u) + (Data.V*Data.U*y)) < 51)
					{
						cout << "///" << *(Out + v + (Data.V*u) + (Data.V*Data.U*y));
					}
				}
			}
		}
	}
}

