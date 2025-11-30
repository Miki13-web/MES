#pragma once
#include <iostream>

struct node
{
    double x, y;
    bool BC = false;
};

struct Jakobian {
    double J[2][2];
    double J1[2][2];
    double detJ;
};

struct element
{
    int ID[4];
    Jakobian* Jaco; //wskaünik na tab jakobianow
    double H[4][4];
    //inicjuje Hbc zerami
    double Hbc[4][4] = { 0.0 };
    double P[4] = { 0.0 };
    double C[4][4] = { 0.0 };
};
