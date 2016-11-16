#ifndef CUDALOADING_H
#define CUDALOADING_H
#include"readmat.h"


class cudaloading
{
private:
	//readmat thisMat;
	int* arraypointer;
	const size_t* dimensionpointer;

public:
	cudaloading(readmat& MAT);
	void cudaprepareInt();


};



#endif