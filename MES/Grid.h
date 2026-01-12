#pragma once
#include "Structures.h"
#include <vector>

using namespace std;

struct grid
{
    int nN;
    int nE;
    node* nodes;
    element* elements;
    int* BC;
    int BC_count;
    vector<int> inductionNodes;
    int probeNodeID;

    grid() {
        nN = 0;
        nE = 0;
        nodes = nullptr;
        elements = nullptr;
        BC = nullptr;
        BC_count = 0;
        probeNodeID = -1;
    }

    ~grid() { 
        if (nodes != nullptr) {
            delete[] nodes;
            nodes = nullptr; 
        }

        if (elements != nullptr) {
            for (int i = 0; i < nE; i++) {
                if (elements[i].Jaco != nullptr) {
                    delete[] elements[i].Jaco;
                }
            }
            delete[] elements;
            elements = nullptr;
        }

        if (BC != nullptr) {
            delete[] BC;
            BC = nullptr;
        }
    }
};