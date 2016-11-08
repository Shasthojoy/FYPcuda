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
}

void readmat::getmxpointer()
{
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
				for (int q = 0; q < 1; q++)
				{
					painfo = matGetNextVariableInfo(pmat, &name);

				}
				matClose(pmat);
			}

		}

	}


}

int readmat::getnumbrofdimensions()
{

	getmxpointer();
	return mxGetNumberOfDimensions(painfo);



}

double * readmat::getarraypointer()
{
	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		throw noMAT;
	}
	else
	{
		for (int q = 0; q < 1; q++)
		{

			pa = matGetNextVariable(pmat, &name);
		}
		return mxGetPr(pa);
	}


}

const size_t*  readmat::dimensionpointer()
{

	getmxpointer();
	dimarray = mxGetDimensions(painfo);
	return dimarray;

}

int readmat::numberofelements()
{
	getmxpointer();
	return mxGetNumberOfElements(painfo);

}
