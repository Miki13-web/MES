#include "SystemEquations.h"
SystemEquations::SystemEquations(int size) : nN(size) {
    HG = new double* [nN];
    Pg = new double[nN]();
	t = new double[nN]();

	Cg = new double* [nN];

    //zerowanie HG i Pg
    for (int i = 0; i < nN; i++) {
        Pg[i] = 0.0;
        HG[i] = new double[nN]();
		Cg[i] = new double[nN]();
    }
   
}

// Destruktor
SystemEquations::~SystemEquations() {
    for (int i = 0; i < nN; i++) {
        delete[] HG[i];
		delete[] Cg[i];
    }
    delete[] HG;
    delete[] Pg;
	delete[] t;
	delete[] Cg;
}