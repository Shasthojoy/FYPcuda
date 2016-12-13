#include "mat.h"
#include <iostream>
#include <stdio.h>
#include <string.h> /* For strcmp() */
#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <vector> /* For STL */
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
	//double* elements;

} DataIn;

inline float getelement(DataIn data, int x, int v, int u, int y, double * __restrict dataElements)
{
	int index = data.X*data.Y*data.U*(v)+data.X*data.Y*(u)+data.X*(y)+x;
	return (float)dataElements[index];
}


int main()
{
	const clock_t begin_time = clock();
	int i;
	double start, end,start1;
		
	start1 = omp_get_wtime();
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
		for (int k = 0; k < numberofdimension; k++)
		{
			cout << *(dimepointer + k);
			cout << "  ";
		}
		elements = thismat.numberofelements();
		// Dimensions
		DataIn data;
		size_t X, Y, U, V;
		data.X = *(dimepointer);
		data.Y = *(dimepointer + 1);
		data.U = *(dimepointer + 2);
		data.V = *(dimepointer + 3);
		double* dataelements;
		dataelements = thismat.getarraypointer();
		
		float*  __restrict image = new float[data.U*data.V];
		int m = 0;
		
		int X_shift, Y_shift;
		omp_set_dynamic(0);
		omp_set_num_threads(4);
		start = omp_get_wtime();
#pragma omp parallel
		{
			float temp_sum;
			int id = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
			std::cout << "thread is " << id << endl;

			for (int v = 0; v < data.V; v += num_threads)
			{
				for (int u = 0; u < data.U; u++)
				{
					temp_sum = 0;

					for (int y = 0; y < data.Y; y++)
					{
						Y_shift = y*m;

						for (int x = 0; x < data.X; x++)
						{
							
							X_shift = x*m;
							if (v < Y_shift && u < X_shift){
								temp_sum += getelement(data, x, (data.V - (Y_shift - (v + id))), (data.U - (X_shift - u)), y, dataelements);


							}


							else if (v >= Y_shift && u < X_shift){
								temp_sum += getelement(data, x, ((v + id) - (Y_shift)), (data.U - (X_shift - u)), y, dataelements);

							}

							else if (v < Y_shift &&  u >= X_shift){
								temp_sum += getelement(data, x, (data.V - (Y_shift - (v + id))), (u - (X_shift)), y, dataelements);

							}

							else if (v >= Y_shift && u >= X_shift){
								temp_sum += getelement(data, x, ((v + id) - (Y_shift)), (u - (X_shift)), y, dataelements);

							}


						}

					}
					//#pragma omp critical
					{
						*(image + (u*data.V) + (v + id)) = temp_sum / (data.X*data.Y * 255);
					}

				}
			}
		}
		cout << "data is" << endl;
		for (int t = 0; t < 15; t++)
			cout << image[t] << " ";
		
		cout << "The end of array \n FINISHED\n";


		 end = omp_get_wtime();
		std::cout << "\n The wall time algo: " << (end - start)<<endl;
		std::cout << "\n The wall time: " << (end - start1) << endl;
		out = Mat(data.U, data.V, CV_32FC1, image); //create an image
		//cout << "out = " << endl << " " << out << endl << endl;
		if (out.empty()) //check whether the image is loaded or not
		{
			cout << "Error : Image cannot be loaded..!!" << endl;
			//system("pause"); //wait for a key press
			return -1;
		}
		std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;



		//namedWindow("Refocused Image", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
		imshow("Refocused Image", out); //display the image which is stored in the 'img' in the "MyWindow" window

		waitKey(0);  //wait infinite time for a keyress

		destroyWindow("Refocused Image"); //destroy the window with the name, "MyWindow"*/




		cout << endl;
		cout << "This is the number of dimensions  ";
		cout << numberofdimension << endl;

		cout << "The array elements ";
		cout << elements;
		cout << endl;

		for (int n = 0; n < elements; n++)
		{
			//cout << *(arraypint + n);
			//cout << " ";
		}
	}
	catch (exception& e)
	{
		cout << "ERROR:"<< e.what() << endl;
	}


	//cin >> i;
	return 0;

}

