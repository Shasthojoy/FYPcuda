#ifndef PROGRAM_EXCEPTION_H
#define PROGRAM_EXCEPTION_H
#include <exception>
using namespace std;
class no_mat_exception : public exception
{
	virtual const char* what() const throw()
	{

		return "Cannot read the mat file";

	}

};


#endif