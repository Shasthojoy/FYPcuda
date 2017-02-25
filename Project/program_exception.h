#ifndef PROGRAM_EXCEPTION_H
#define PROGRAM_EXCEPTION_H
#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
using namespace std;



class no_mat_exception : public exception
{
	virtual const char* what() const throw()
	{

		return "Cannot read the mat file";

	}

};

class device_exception
{
public:
	device_exception();

	void throw_cuda_error(cudaError_t code)
	{
		if (code != cudaSuccess)
		throw thrust::system_error(code, thrust::cuda_category());
	}
};




#endif