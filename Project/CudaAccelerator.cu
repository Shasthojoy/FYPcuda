#include"CudaAccelerator.cuh"


CudaAccelerator::CudaAccelerator(Data_Struct_In ImageMetaData, float Alpha)
{

		// Aquiring LightField Meta data and Other information
		size_of_elements = ImageMetaData.X*ImageMetaData.Y*ImageMetaData.U*ImageMetaData.V;
		AlphaValue = Alpha;
		X = ImageMetaData.X;
		Y = ImageMetaData.Y;
		U = ImageMetaData.U;
		V = ImageMetaData.V;
		Exception_Device.throw_cuda_error(cudaSetDevice(0));
		//Allocate memory In device
		Exception_Device.throw_cuda_error(cudaMalloc((void**)&Device_Input, size_of_elements * sizeof(float)));
		// Allocate host pinned memory
		Exception_Device.throw_cuda_error(cudaMallocHost((void**)&HostPinnedImage, size_of_elements * sizeof(float)));


}

void CudaAccelerator::FractionalShift()
{

}


