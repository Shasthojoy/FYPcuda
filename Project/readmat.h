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
	const char *name1;
	bool isdouble;
	no_mat_exception noMAT;
public:
	
	//with default value
	readmat(const char *f);
	//	get number of dimensions
	int getnumbrofdimensions() const;
	// get pointer to array
	double*  getarraypointer() const;
	// get pointer to each dimension size
	const size_t* dimensionpointer() const;
	//number of elements
	int numberofelements() const;
	

};

#endif