MATFile *pmat;
	int i;
	int	  ndir;
	const char **dir;
	const char *file = "pqfile.mat";
	mxArray *pa;
	const char *name;
	double *in1pr;
	cout << "Hello World!" << endl;
	pmat = matOpen("test.mat", "r");
	if (pmat == NULL) {
		printf("Error opening file\n");

	}
	dir = (const char **)matGetDir(pmat, &ndir);
	if (dir == NULL) {
		printf("Error reading directory of file %s\n", file);
		return(1);
	}
	else {
		printf("Directory of %s:\n", file);
		for (i = 0; i < ndir; i++)
			printf("%s\n", dir[i]);
	}

	mxFree(dir);
	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n", file);
		return(1);
	}
	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
		return(1);
	}
	printf("\nExamining the header for each variable:\n");
	for (i = 0; i < ndir; i++) {
		pa = matGetNextVariableInfo(pmat, &name);
		if (pa == NULL) {
			printf("Error reading in file %s\n", file);
			return(1);
		}
		/* Diagnose header pa */
		printf("According to its header, array %s has %d dimensions\n",
			name, mxGetNumberOfDimensions(pa));
		cout << "Value of the dimensions are: ";
		cout << *(mxGetDimensions(pa)+3) << endl;
		if (mxIsFromGlobalWS(pa))
			printf("  and was a global variable when saved\n");
		else
			printf("  and was a local variable when saved\n");
		mxDestroyArray(pa);
	}
	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n", file);
		return(1);
	}
	pmat = matOpen(file, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
		return(1);
	}

	printf("\nReading in the actual array contents:\n");
	for (i = 0; i<ndir; i++) {
		pa = matGetNextVariable(pmat, &name);
		in1pr = mxGetPr(pa);
		cout << "Value of end variable: ";
		cout << *(in1pr+2) << endl;
		if (pa == NULL) {
			printf("Error reading in file %s\n", file);
			return(1);
		}
		/*
		* Diagnose array pa
		*/
		printf("According to its contents, array %s has %d dimensions\n",
			name, mxGetNumberOfDimensions(pa));
		if (mxIsFromGlobalWS(pa))
			printf("  and was a global variable when saved\n");
		else
			printf("  and was a local variable when saved\n");
		//mxDestroyArray(pa);
	}

	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n", file);
		return(1);
	}
	cout << "Value of var variable: ";
	cout << &pa << endl;
	//mexPrintf("in1=%f\n", *in1pr);
	cin >> i;
	return 0;