#include"cudaloading.cuh"


cudaloading::cudaloading(const readmat& MAT,int alpha)
{
	
	dimensionpointer = MAT.dimensionpointer();
	thisData.X = *(dimensionpointer);
	thisData.Y = *(dimensionpointer + 1);
	thisData.U = *(dimensionpointer + 2);
	thisData.V = *(dimensionpointer + 3);
	thisData.alpha = alpha;
	// allocation of memory in device for data
	size_of_elements = thisData.X*thisData.Y*thisData.U*thisData.V;
	
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		throw device_operation_exception("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	}
	
	cudaStatus = cudaMalloc((void**)&thisData.elements, size_of_elements * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		throw device_operation_exception("Device Memory allocation failed");
	}

	cudaStatus = cudaMemcpy(thisData.elements, MAT.getarraypointer(), size_of_elements * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		throw device_operation_exception("Device Memory copy to device failed");
	}



}

void cudaloading::cudacopy()
{

}


