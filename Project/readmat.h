#ifndef READMAT_H
#define READMAT_H
#include "mat.h"
#include "matrix.h"
#include "mex.h"
#include "program_exception.h"


class readmat
{
private:
	const size_t *dimarray;
	const char **dir;
	MATFile *pmat;
	int	ndir;
	mxArray *painfo,*pa;
	const char *file;
	int minute;
	int second;
	const char *name;
	bool isdouble;
	void getmxpointer();
	no_mat_exception noMAT;
public:
	
	//with default value
	readmat(const char *f);
	//	get number of dimensions
	int getnumbrofdimensions();
	// get pointer to array
	double* getarraypointer();
	// get pointer to each dimension size
	const size_t* dimensionpointer();
	//number of elements
	int numberofelements();
	

};

#endif