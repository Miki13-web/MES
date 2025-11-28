#pragma once
#include "Structures.h"

struct grid
{
    int nN;
    int nE;
    node* nodes;
    element* elements;
    int* BC;
    int BC_count;

    // Warto w przysz³oœci dodaæ tu destruktor zwalniaj¹cy pamiêæ nodes/elements
    // grid() { ... }
    // ~grid() { ... }
};