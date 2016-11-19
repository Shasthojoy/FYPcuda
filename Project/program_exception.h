#ifndef PROGRAM_EXCEPTION_H
#define PROGRAM_EXCEPTION_H
#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;
class no_mat_exception : public exception
{
	virtual const char* what() const throw()
	{

		return "Cannot read the mat file";

	}

};

class device_operation_exception
{
private:
	string errorMessage;
public:
	device_operation_exception(string Error_message)
	{
		errorMessage = Error_message;
	}
	string what() const throw()
	{

		return errorMessage;

	}

};




#endif