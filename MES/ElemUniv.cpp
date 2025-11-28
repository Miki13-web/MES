#include "ElemUniv.h"
#include <iostream>
#include <cmath>

// Implementacja konstuktory Surface
elemUniv::Surface::Surface(int n) : npc_edge(n) {
    wagi = new double[npc_edge];
    N = new double* [npc_edge];
    for (int i = 0; i < npc_edge; i++) {
        N[i] = new double[4];
    }
}

elemUniv::Surface::~Surface() {
    for (int i = 0; i < npc_edge; i++) {
        delete[] N[i];
    }
    delete[] N;
    delete[] wagi;
}

// Implementacja głównego elemUniv
elemUniv::elemUniv(int npc) : npc(npc) {
	int npc2 = npc * npc; // całkowita liczba punktów całkowania w 2D
    dKsi = new double* [npc2];
    dEta = new double* [npc2];
    // funkcje kształtu 
    fN = new double* [npc2];

    for (int i = 0; i < npc2; i++) {
        dKsi[i] = new double[4];
        dEta[i] = new double[4];
		fN[i] = new double[4];
    }
    eta_pc = new double[npc2];
    ksi_pc = new double[npc2];
    wagi = new double[npc2];

	int npc_edge = npc; // punkty całkowania na krawędzi
    double* points_1d = new double[npc_edge];
    double* weights_1d = new double[npc_edge];

    //ustalanie punktów całkowania i wag w zależności od liczby punktów całkowania
    switch (npc) {
    case 2: { // Schemat 2x2
        double p = 1.0 / sqrt(3.0); // wspolrzedna PC 
        ksi_pc[0] = -p; eta_pc[0] = -p; // PC 1
        ksi_pc[1] = p; eta_pc[1] = -p; // PC 2 
        ksi_pc[2] = -p; eta_pc[2] = p; // PC 3 
        ksi_pc[3] = p; eta_pc[3] = p; // PC 4

        //punkty na krawedziach dla macierzy Hbc
        points_1d[0] = -p; weights_1d[0] = 1.0;
        points_1d[1] = p; weights_1d[1] = 1.0;

        // dla macierzy funkcji kształtu N

        for (int i = 0; i < 4; ++i) wagi[i] = 1.0; // wagi
        break;
    }
    case 3: { // Schemat 3x3
        double p1 = sqrt(3.0 / 5.0);
        double p2 = 0.0;
        double w1 = 5.0 / 9.0;
        double w2 = 8.0 / 9.0;

        ksi_pc[0] = -p1; eta_pc[0] = -p1; wagi[0] = w1 * w1;
        ksi_pc[1] = p2;  eta_pc[1] = -p1; wagi[1] = w2 * w1;
        ksi_pc[2] = p1;  eta_pc[2] = -p1; wagi[2] = w1 * w1;
        ksi_pc[3] = -p1; eta_pc[3] = p2;  wagi[3] = w1 * w2;
        ksi_pc[4] = p2;  eta_pc[4] = p2;  wagi[4] = w2 * w2;
        ksi_pc[5] = p1;  eta_pc[5] = p2;  wagi[5] = w1 * w2;
        ksi_pc[6] = -p1; eta_pc[6] = p1;  wagi[6] = w1 * w1;
        ksi_pc[7] = p2;  eta_pc[7] = p1;  wagi[7] = w2 * w1;
        ksi_pc[8] = p1;  eta_pc[8] = p1;  wagi[8] = w1 * w1;

        //punkty na krawedziach dla macierzy Hbc
        points_1d[0] = -p1; weights_1d[0] = w1;
        points_1d[1] = p2; weights_1d[1] = w2;
        points_1d[2] = p1; weights_1d[2] = w1;

        break;
    }
    default:
        std::cerr << "Nieobsługiwana liczba punktów całkowania: " << npc << std::endl;
        // Domyślnie 4, aby uniknąć awarii
        npc2 = 4;
        double p = 1.0 / sqrt(3.0);
        ksi_pc[0] = -p; eta_pc[0] = -p;
        ksi_pc[1] = p; eta_pc[1] = -p;
        ksi_pc[2] = p; eta_pc[2] = p;
        ksi_pc[3] = -p; eta_pc[3] = p;
        for (int i = 0; i < 4; ++i) wagi[i] = 1.0;
        break;
    }

	// Obliczanie wartości pochodnych w każdym punkcie całkowania i funkcji kształtu dla macierzy C
    for (int i = 0; i < npc2; i++) {
        double eta = eta_pc[i];
        double ksi = ksi_pc[i];

        // Pochodne po KSI 
        dKsi[i][0] = -0.25 * (1.0 - eta);
        dKsi[i][1] = 0.25 * (1.0 - eta);
        dKsi[i][2] = 0.25 * (1.0 + eta);
        dKsi[i][3] = -0.25 * (1.0 + eta);

        // Pochodne po ETA (η)
        dEta[i][0] = -0.25 * (1.0 - ksi);
        dEta[i][1] = -0.25 * (1.0 + ksi);
        dEta[i][2] = 0.25 * (1.0 + ksi);
        dEta[i][3] = 0.25 * (1.0 - ksi);

		// Funkcje kształtu N do macierzy C
		fN[i][0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
		fN[i][1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
		fN[i][2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
		fN[i][3] = 0.25 * (1.0 - ksi) * (1.0 + eta);

    }

    //tworzenie 4 powierzchni
    for (int i = 0; i < 4; i++) {
        surfaces[i] = new Surface(npc_edge);
    }

    // Pętla po 4 krawędziach
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < npc_edge; j++) {
            double ksi = 0.0, eta = 0.0;

            // Mapowanie punktu 1D na odpowiednią krawędź 2D
            if (i == 0) { // Dół (Bottom): eta = -1, ksi zmienne
                ksi = points_1d[j];
                eta = -1.0;
            }
            else if (i == 1) { // Prawa (Right): ksi = 1, eta zmienne
                ksi = 1.0;
                eta = points_1d[j];
            }
            else if (i == 2) { // Góra (Top): eta = 1, ksi zmienne (często odwrócone, ale dla skalara N nie ma znaczenia)
                ksi = points_1d[npc_edge - 1 - j]; // Opcjonalnie odwrócona kolejność
                eta = 1.0;
            }
            else if (i == 3) { // Lewa (Left): ksi = -1, eta zmienne
                ksi = -1.0;
                eta = points_1d[npc_edge - 1 - j];
            }

            // Zapisanie wagi dla punktu powierzchni
            surfaces[i]->wagi[j] = weights_1d[j];

            // Obliczenie funkcji kształtu N w tym punkcie brzegowym
            // Wzory: N1=0.25(1-k)(1-e), N2=0.25(1+k)(1-e), N3=0.25(1+k)(1+e), N4=0.25(1-k)(1+e)
            surfaces[i]->N[j][0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
            surfaces[i]->N[j][1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
            surfaces[i]->N[j][2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
            surfaces[i]->N[j][3] = 0.25 * (1.0 - ksi) * (1.0 + eta);
        }
    }
    delete[] points_1d;
    delete[] weights_1d;
}

elemUniv::~elemUniv() {
    // Sprzątanie wnętrza
    for (int i = 0; i < npc*npc; i++) {
        delete[] dKsi[i];
        delete[] dEta[i];
		delete[] fN[i];
    }
    delete[] dKsi;
    delete[] dEta;
    delete[] eta_pc;
    delete[] ksi_pc;
    delete[] wagi;
	delete[] fN;

    // Sprzątanie powierzchni
    for (int i = 0; i < 4; i++) {
        delete surfaces[i];
    }
}