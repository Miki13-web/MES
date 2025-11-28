#pragma once

struct elemUniv {
    int npc;
    double** dKsi;
    double** dEta;

    double* eta_pc;
    double* ksi_pc;
    double* wagi;

    double** fN;

    struct Surface {
        int npc_edge;
        double** N;
        double* wagi;

		//konstruktor i destruktor
        Surface(int n);
        ~Surface();

    };

    // Tablica 4 wskaŸników na powierzchnie (Dó³, Prawa, Góra, Lewa)
    Surface* surfaces[4];

    elemUniv(int npc);
    ~elemUniv();
};