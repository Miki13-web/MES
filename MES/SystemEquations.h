#pragma once

struct SystemEquations {
    double** HG;
    int nN;
    double* Pg;
    double* t;
    double** Cg;

    SystemEquations(int size);
    ~SystemEquations();
};