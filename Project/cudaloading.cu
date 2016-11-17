#include"cudaloading.cuh"
#include "device_launch_parameters.h"
#include "cuda_runtime.h"

cudaloading::cudaloading(const readmat& MAT)
{
	arraypointer = MAT.getarraypointer;
	dimensionpointer = MAT.dimensionpointer;

}




