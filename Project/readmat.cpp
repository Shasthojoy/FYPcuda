#include "readmat.h"
#include <iostream>
#include "mat.h"
#include "matrix.h"
#include "mex.h"
using namespace std;

// set the file name
readmat::readmat(const char *f)
{
	file = f;
	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		throw noMAT;
	}
	else
	{
		dir = (const char **)matGetDir(pmat, &ndir);
		if (dir == NULL) {
			printf("Error reading directory of file %s\n", file);

		}
		else if (ndir > 1)
		{
			cout << "The number of variables are larger than 1" << endl;
		}
		else
		{
			mxFree(dir);
			matClose(pmat);
			pmat = matOpen(file, "r");
			if (pmat == NULL) {
				throw noMAT;
			}
			else
			{
				painfo = matGetNextVariableInfo(pmat, &name);
				matClose(pmat);
			}
			pmat = matOpen(file, "r");
			if (pmat == NULL) {
				throw noMAT;
			}
			else
			{
				pa = matGetNextVariable(pmat, &name1);
				matClose(pmat);
			}
		}

	}

}



int readmat::getnumbrofdimensions() const
{

	return mxGetNumberOfDimensions(painfo);
}

double * readmat::getarraypointer() const
{
		return mxGetPr(pa);
}

const size_t*  readmat::dimensionpointer() const
{
	return mxGetDimensions(painfo);
}

int readmat::numberofelements() const
{
	return mxGetNumberOfElements(painfo);
}
