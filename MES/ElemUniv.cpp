#include "ElemUniv.h"
#include <iostream>
#include <cmath>

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

elemUniv::elemUniv(int npc) : npc(npc) {
	int npc2 = npc * npc; // liczba punktów całkowania w 2D
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

	int npc_edge = npc; // liczba punktow całkowania na krawędzi
    double* points_1d = new double[npc_edge];
    double* weights_1d = new double[npc_edge];

    //ustalanie punktów całkowania i wag w zależności od liczby punktów całkowania
    switch (npc) {
    case 2: {
        double p = 1.0 / sqrt(3.0); 
        ksi_pc[0] = -p; eta_pc[0] = -p; 
        ksi_pc[1] = p; eta_pc[1] = -p;  
        ksi_pc[2] = -p; eta_pc[2] = p;  
        ksi_pc[3] = p; eta_pc[3] = p;

        //punkty na krawedziach dla Hbc
        points_1d[0] = -p; weights_1d[0] = 1.0;
        points_1d[1] = p; weights_1d[1] = 1.0;


        for (int i = 0; i < 4; ++i) wagi[i] = 1.0; // wagi
        break;
    }
    case 3: {
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

        //punkty na krawedziach dla Hbc
        points_1d[0] = -p1; weights_1d[0] = w1;
        points_1d[1] = p2; weights_1d[1] = w2;
        points_1d[2] = p1; weights_1d[2] = w1;

        break;
    }
    case 4: {
		double p1 = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
		double p2 = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
		double w1 = (18.0 - sqrt(30.0)) / 36.0;
		double w2 = (18.0 + sqrt(30.0)) / 36.0;

		ksi_pc[0] = -p1; eta_pc[0] = -p1; wagi[0] = w1 * w1;
		ksi_pc[1] = -p2; eta_pc[1] = -p1; wagi[1] = w2 * w1;
		ksi_pc[2] = p2;  eta_pc[2] = -p1; wagi[2] = w2 * w1;
		ksi_pc[3] = p1;  eta_pc[3] = -p1; wagi[3] = w1 * w1;
		ksi_pc[4] = -p1; eta_pc[4] = -p2; wagi[4] = w1 * w2;
		ksi_pc[5] = -p2; eta_pc[5] = -p2; wagi[5] = w2 * w2;
		ksi_pc[6] = p2;  eta_pc[6] = -p2; wagi[6] = w2 * w2;
		ksi_pc[7] = p1;  eta_pc[7] = -p2; wagi[7] = w1 * w2;
		ksi_pc[8] = -p1; eta_pc[8] = p2;  wagi[8] = w1 * w2;
		ksi_pc[9] = -p2; eta_pc[9] = p2;  wagi[9] = w2 * w2;
		ksi_pc[10] = p2; eta_pc[10] = p2; wagi[10] = w2 * w2;
		ksi_pc[11] = p1; eta_pc[11] = p2; wagi[11] = w1 * w2;
		ksi_pc[12] = -p1; eta_pc[12] = p1; wagi[12] = w1 * w1;
		ksi_pc[13] = -p2; eta_pc[13] = p1; wagi[13] = w2 * w1;
		ksi_pc[14] = p2;  eta_pc[14] = p1; wagi[14] = w2 * w1;
		ksi_pc[15] = p1;  eta_pc[15] = p1; wagi[15] = w1 * w1;

        //punkty krawedzi dla Hbc
		points_1d[0] = -p1; weights_1d[0] = w1;
		points_1d[1] = -p2; weights_1d[1] = w2;
		points_1d[2] = p2;  weights_1d[2] = w2;
		points_1d[3] = p1;  weights_1d[3] = w1;

        break;
    }
    default:
        std::cerr << "Nieobsługiwana liczba punktów całkowania: " << npc << std::endl;
        npc2 = 4;
        double p = 1.0 / sqrt(3.0);
        ksi_pc[0] = -p; eta_pc[0] = -p;
        ksi_pc[1] = p; eta_pc[1] = -p;
        ksi_pc[2] = p; eta_pc[2] = p;
        ksi_pc[3] = -p; eta_pc[3] = p;
        for (int i = 0; i < 4; ++i) wagi[i] = 1.0;
        break;
    }

	// wartości pochodnych w każdym punkcie całkowania i funkcji kształtu dla macierzy C
    for (int i = 0; i < npc2; i++) {
        double eta = eta_pc[i];
        double ksi = ksi_pc[i];
 
        dKsi[i][0] = -0.25 * (1.0 - eta);
        dKsi[i][1] = 0.25 * (1.0 - eta);
        dKsi[i][2] = 0.25 * (1.0 + eta);
        dKsi[i][3] = -0.25 * (1.0 + eta);

        dEta[i][0] = -0.25 * (1.0 - ksi);
        dEta[i][1] = -0.25 * (1.0 + ksi);
        dEta[i][2] = 0.25 * (1.0 + ksi);
        dEta[i][3] = 0.25 * (1.0 - ksi);

		fN[i][0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
		fN[i][1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
		fN[i][2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
		fN[i][3] = 0.25 * (1.0 - ksi) * (1.0 + eta);

    }

    // 4 powierzchnie
    for (int i = 0; i < 4; i++) {
        surfaces[i] = new Surface(npc_edge);
    }

    // pętla po krawędziach
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < npc_edge; j++) {
            double ksi = 0.0, eta = 0.0;

            if (i == 0) { 
                ksi = points_1d[j];
                eta = -1.0;
            }
            else if (i == 1) { 
                ksi = 1.0;
                eta = points_1d[j];
            }
            else if (i == 2) { 
                ksi = points_1d[npc_edge - 1 - j]; 
                eta = 1.0;
            }
            else if (i == 3) {
                ksi = -1.0;
                eta = points_1d[npc_edge - 1 - j];
            }

            // waga dla punktu powierzchni
            surfaces[i]->wagi[j] = weights_1d[j];

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

    for (int i = 0; i < 4; i++) {
        delete surfaces[i];
    }
}