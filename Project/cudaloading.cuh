#ifndef CUDALOADING_H
#define CUDALOADING_H
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include"readmat.h"
#include "program_exception.h"


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
	cudaError_t cudaStatus;
	size_t size_of_elements;
	
	

public:
	cudaloading(const readmat& MAT, int alpha);
	~cudaloading();
	void cudaprepareInt();
	void cudacopy();


};



#endif