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

    grid() {
        nN = 0;
        nE = 0;
        nodes = nullptr;
        elements = nullptr;
        BC = nullptr;
        BC_count = 0;
    }

    ~grid() { 
        // Sprz¹tanie wêz³ów
        if (nodes != nullptr) {
            delete[] nodes;
            nodes = nullptr; // Dobra praktyka
        }

        // Sprz¹tanie elementów i ich Jakobianów
        if (elements != nullptr) {
            // Najpierw musimy usun¹æ tablice Jaco wewn¹trz ka¿dego elementu
            // (bo element nie ma w³asnego destruktora)
            for (int i = 0; i < nE; i++) {
                if (elements[i].Jaco != nullptr) {
                    delete[] elements[i].Jaco;
                }
            }
            // Dopiero teraz usuwamy tablicê elementów
            delete[] elements;
            elements = nullptr;
        }

        // Sprz¹tanie tablicy BC
        if (BC != nullptr) {
            delete[] BC;
            BC = nullptr;
        }
    }
};