#ifndef CUDALOADING_H
#define CUDALOADING_H
#include"readmat.h"
typedef struct {
	size_t X;
	size_t Y;
	size_t U;
	size_t V;
	double* elements;
} DataIn;

class cudaloading
{
private:
	//readmat thisMat;
	
	DataIn thisData;

public:
	cudaloading(const readmat& MAT);
	void cudaprepareInt();
	void cudacopy();


};



#endif