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

		vector<vector<vector <double*>>> image4d;
		image4d.resize(Y);
		for (int y = 0; y < Y; y++)
		{
			image4d[y].resize(V);
			for (int v = 0; v < V; v++)
				image4d[y][v].resize(X*U);
		}

		//get values matlab
		arraypint = thismat.getarraypointer();
		elements = thismat.numberofelements();

		// Assigning values
		int index;

		for (int v = 0; v < V; v++)
		{
			for (int y = 0; y < Y; y++)
			{
				for (int u = 0; u < U; u++)
				{
					// index start
					//index_end = index_start + (X - 1);
					for (int x = 0; x < X; x++)
					{
						index = X*Y*U*(v)+X*Y*(u)+X*(y)+x;
						image4d[y][v][x*U + u] = index + arraypint;


					}
				}
			}
		}
		

		// preprocessing data for shift the columns
		vector <signed int> shiftcolvector;
		int x, u;

		//vector<vector <int*>> column_shift_data;
		//column_shift_data.resize(Y);
		valarray<valarray<double**>> column_shift_data;
		column_shift_data.resize(Y);

		int alpha = 1;
		int m_shift = alpha - 1;
		int temp;

		for (int y = 0; y < Y; y++)
		{
			temp = round(y*m_shift);
			shiftcolvector.push_back(temp);
			//cout << "data shift " << round(y*m_shift);
			//cout << "data shift real " << shiftcolvector[y];
		}

		//shift the pointers
		for (int y = 0; y < Y; y++)
		{
			column_shift_data[y].resize(V);
			for (int v = 0; v < V; v++)
			{
				column_shift_data[y][v] = image4d[y][v].data();
			}
			// shift the pointers
			//rotate(column_shift_data[y].begin(), column_shift_data[y].begin() - shiftcolvector[y], column_shift_data[y].end());
			column_shift_data[y] = column_shift_data[y].cshift(-shiftcolvector[y]);
		}
		//take the sum
		vector<vector <vector<double>>> sum_data;
		sum_data.resize(X);
		for (int x = 0; x < X; x++)
		{
			sum_data[x].resize(U);
			for (int u = 0; u < U; u++)
			{
				sum_data[x][u].resize(V);
			}

		}

		for (int xu = 0; xu < X*U; xu++)
		{
			for (int v = 0; v < V; v++)
			{
				x = xu / U;
				u = xu % U;
				sum_data[x][u][v] = 0;
				for (int y = 0; y < Y; y++)
				{
					sum_data[x][u][v] = **(column_shift_data[y][v] + xu) + sum_data[x][u][v];
				}
			}
		}
		//print the sum
		/*
		cout << "\nThis is the SUM column start\n";
		for (int x = 0; x < X; x++)
		{
		for (int u = 0; u < U; u++)
		{
		for (int v = 0; v < V; v++)

		{
		cout << " " << sum_data[x][u][v];
		}
		cout << "\n";
		}
		cout << "\n";
		}
		cout << "\nThis is the SUM column\n";
		*/
		// data is in format to shift and take sum along y axis. create 2d pointer array and then shift the pointers and take sum to a 2d array.
		vector<signed int>shiftrowvector;

		for (int x = 0; x < X; x++){

			temp = round(x*m_shift);
			shiftrowvector.push_back(temp);

		}

		//shift along u dimension

		valarray<valarray<double*>>row_shift_data;

		row_shift_data.resize(X);

		for (int x = 0; x < X; x++)
		{
			row_shift_data[x].resize(U);
			for (int u = 0; u < U; u++)
			{
				row_shift_data[x][u] = sum_data[x][u].data();//2d array of pointers pointing to the rows
				//cout << "\nthis is pointer\n" << *(row_shift_data[x][u]);
			}

			row_shift_data[x] = row_shift_data[x].cshift(-shiftrowvector[x]);//shift the images along u


		}

		//create the image matrix

		float* image = new float[U*V];
		//for (int i = 0; i < U; ++i)
		//image[i] = new double[V];

		//sum along the third dimension in sum_data

		for (int u = 0; u < U; u++)
		{
			for (int v = 0; v < V; v++)
			{

				*(image + (u*V) + v) = 0;
				for (int x = 0; x < X; x++)
				{
					*(image + (u*V) + v) = *(row_shift_data[x][u] + v) + *(image + (u*V) + v);// sum along the rows
				}
				*(image + (u*V) + v) = *(image + (u*V) + v) / (X*Y * 255);
			}
		}

		//print the image
		/*
		cout << "\n";

		for (int u = 0; u < U; u++)
		{
		for (int v = 0; v < V; v++)
		{
		cout << " " << image[(u*V) + v];
		}
		cout << "\n";
		}

		*/


		cout << "The end of array \n FINISHED\n";



		out = Mat(U, V, CV_32FC1, image); //create an image
		//cout << "out = " << endl << " " << out << endl << endl;
		if (out.empty()) //check whether the image is loaded or not
		{
			cout << "Error : Image cannot be loaded..!!" << endl;
			//system("pause"); //wait for a key press
			return -1;
		}
		std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;



		namedWindow("Refocused Image", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
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


	cin >> i;
	return 0;

}

