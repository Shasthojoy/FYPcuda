#ifndef CUDAACCELERATOR_H
#define CUDAACCELERATOR_H
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include"readmat.h"
#include "program_exception.h"




typedef struct {
	size_t X;
	size_t Y;
	size_t U;
	size_t V;

	int alpha;
} Data_Struct_In;

class CudaAccelerator
{
private:
	size_t X, Y, U, V; //Dimensions of the Light Field
	size_t size_of_elements; // Size of the of light field
	float * Device_Input; // Device input data pointer
	float * HostPinnedImage; // Host pinned memory
	device_exception Exception_Device; // Cuda exception object
	float AlphaValue;
public:
	CudaAccelerator(const Data_Struct_In ImageMetaData, float alpha);
	~CudaAccelerator();
	void FractionalShift();
	


};



#endif