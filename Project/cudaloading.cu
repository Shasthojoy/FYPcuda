#include"cudaloading.cuh"
#include "device_launch_parameters.h"
#include "cuda_runtime.h"

cudaloading::cudaloading(const readmat& MAT,int alpha)
{
	thisData.elements = MAT.getarraypointer();
	dimensionpointer = MAT.dimensionpointer();
	thisData.X = *(dimensionpointer);
	thisData.Y = *(dimensionpointer + 1);
	thisData.U = *(dimensionpointer + 2);
	thisData.V = *(dimensionpointer + 3);
	thisData.alpha = alpha;



}

void cudaloading::cudacopy()
{

}


