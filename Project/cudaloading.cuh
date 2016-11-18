#ifndef CUDALOADING_H
#define CUDALOADING_H
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include"readmat.h"

typedef struct {
	size_t X;
	size_t Y;
	size_t U;
	size_t V;
	double* elements;
	int no_of_elements;
	int alpha;
} DataIn;

class cudaloading
{
private:
	//readmat thisMat;
	const size_t* dimensionpointer;
	DataIn thisData;
	double* device_elements;
	cudaError_t cudaStatus;

public:
	cudaloading(const readmat& MAT, int alpha);
	void cudaprepareInt();
	void cudacopy();


};



#endif