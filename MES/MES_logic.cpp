#include "MES_logic.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Oblicza macierze lokalne (H, Hbc, P) dla pojedynczego elementu
void calculateElementMatrices(int i, grid& gri, GlobalData& gData, elemUniv& elemU) {
    
    int nr_elementu = i;

    //Zerowanie macierzy lokalnych H i Hbc oraz C
    for (int m = 0; m < 4; m++) {
        //zerowanie wektora P
        gri.elements[i].P[m] = 0.0;
        for (int n = 0; n < 4; n++) {
            gri.elements[i].H[m][n] = 0.0;
            gri.elements[i].Hbc[m][n] = 0.0;
			gri.elements[i].C[m][n] = 0.0;
        }
    }

    // Pêtla po punktach ca³kowania (Jacobian, H)
    //  Pêtla po wszystkich punktach ca³kowania (PC) dla danego elementu
    for (int j = 0; j < (gData.npc*gData.npc); j++) {
        //cout << "\n  Punkt Calkowania (PC) " << j + 1 << " (ksi=" << elemU.ksi_pc[j] << ", eta=" << elemU.eta_pc[j] << "):" << endl;

        int numer_wagi = j;

        node nodes_xy[4];
        for (int k = 0; k < 4; k++) {
            int nodeId = gri.elements[i].ID[k] - 1;
            nodes_xy[k] = gri.nodes[nodeId];
        }

        // lokalne pochodne dla bie¿¹cego PC
        double* dN_dKsi = elemU.dKsi[j];
        double* dN_dEta = elemU.dEta[j];

		//lokalne wartoœci funkcji kszta³tu N dla macierzy C
		double* N_local = elemU.fN[j];

        // komponenty Jakobianu 
        double dxdksi = 0.0, dydksi = 0.0, dxdeta = 0.0, dydeta = 0.0;


        for (int k = 0; k < 4; k++) {
            dxdksi += dN_dKsi[k] * nodes_xy[k].x;
            dydksi += dN_dKsi[k] * nodes_xy[k].y;
            dxdeta += dN_dEta[k] * nodes_xy[k].x;
            dydeta += dN_dEta[k] * nodes_xy[k].y;
        }


        Jakobian& jaco = gri.elements[i].Jaco[j];
        jaco.J[0][0] = dxdksi;
        jaco.J[0][1] = dydksi;
        jaco.J[1][0] = dxdeta;
        jaco.J[1][1] = dydeta;

        // wyznacznik
        jaco.detJ = (jaco.J[0][0] * jaco.J[1][1]) - (jaco.J[0][1] * jaco.J[1][0]);

        // macierz odwrotna J^-1
        double invDetJ = 1.0 / jaco.detJ;
        jaco.J1[0][0] = jaco.J[1][1] * invDetJ;
        jaco.J1[0][1] = -jaco.J[0][1] * invDetJ;
        jaco.J1[1][0] = -jaco.J[1][0] * invDetJ;
        jaco.J1[1][1] = jaco.J[0][0] * invDetJ;

        // wyniki dla tego PC
        //cout << "    Jakobian J = [ " << setw(12) << jaco.J[0][0] << ", " << setw(12) << jaco.J[0][1] << " ]" << endl;
        //cout << "                 [ " << setw(12) << jaco.J[1][0] << ", " << setw(12) << jaco.J[1][1] << " ]" << endl;
        //cout << "    det(J)     = " << jaco.detJ << endl;
        //cout << "    J^-1       = [ " << setw(12) << jaco.J1[0][0] << ", " << setw(12) << jaco.J1[0][1] << " ]" << endl;
        //cout << "                 [ " << setw(12) << jaco.J1[1][0] << ", " << setw(12) << jaco.J1[1][1] << " ]" << endl;

        // Obliczamy dla kazdego punktu pochodne globalne
        double dNdx[4], dNdy[4];//tablica na pochodne globalne dla punktu
        for (int k = 0; k < 4; k++) { // Pêtla po 4 funkcjach kszta³tu (N1..N4)
            dNdx[k] = jaco.J1[0][0] * dN_dKsi[k] + jaco.J1[0][1] * dN_dEta[k];
            dNdy[k] = jaco.J1[1][0] * dN_dKsi[k] + jaco.J1[1][1] * dN_dEta[k];
        }

        // macierz H dla punktu
        //macierz C dla punktu
        double H_pc[4][4];
		double C_pc[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H_pc[i][j] = (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j]) * gData.Conductivity * jaco.detJ;
				//macierz C
				C_pc[i][j] = N_local[i] * N_local[j] * gData.Density * gData.SpecificHeat * jaco.detJ;
                //dodaje sobie odrazu do macierzy dla elementu
                gri.elements[nr_elementu].H[i][j] += H_pc[i][j] * elemU.wagi[numer_wagi];
				gri.elements[nr_elementu].C[i][j] += C_pc[i][j] * elemU.wagi[numer_wagi];
            }
        }

		
        
        // wypisanie macierzy H dla punktu
        //mozna zakomentowac
       /* cout << "\n    Macierz H dla tego punktu calkowania:" << endl;
        for (int m = 0; m < 4; m++) {
            cout << "      [ ";
            for (int n = 0; n < 4; n++) {
                cout << setw(12) << fixed << setprecision(6) << H_pc[m][n] << (n == 3 ? " ]" : ",");
            }
            cout << endl << endl;
        }*/
		//wypisanie macierzy C dla punktu
        /*cout << "\n    Macierz C dla tego punktu calkowania:" << endl;
        for (int m = 0; m < 4; m++) {
            cout << "      [ ";
            for (int n = 0; n < 4; n++) {
                cout << setw(12) << fixed << setprecision(6) << C_pc[m][n] << (n == 3 ? " ]" : ",");
            }
            cout << endl << endl;
		}*/
    }
    
    //Pêtla po bokach (Hbc)
    // ==============================================
    // OBLICZANIE MACIERZY H_BC (KONWEKCJA) i dodatkowo tutaj wektor P
    // ==============================================


    // Bok 0 (Dó³): wêz³y 0-1
    // Bok 1 (Prawa): wêz³y 1-2
    // Bok 2 (Góra): wêz³y 2-3
    // Bok 3 (Lewa): wêz³y 3-0
    int edge_nodes[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

    // Pêtla po 4 bokach elementu
    for (int edge = 0; edge < 4; edge++) {

        //lokalne indeksy wêz³ów
        int local_n1 = edge_nodes[edge][0];
        int local_n2 = edge_nodes[edge][1];

        //globalne ID
        int id1 = gri.elements[i].ID[local_n1];
        int id2 = gri.elements[i].ID[local_n2];

        // spr BC
        // flag bool BC ze struktury node
        if (gri.nodes[id1 - 1].BC && gri.nodes[id2 - 1].BC) {

            // d³ugoœæ boku (potrzebne do Jakobianu 1D
            node n1 = gri.nodes[id1 - 1];
            node n2 = gri.nodes[id2 - 1];
            double length = sqrt(pow(n2.x - n1.x, 2) + pow(n2.y - n1.y, 2));

            double detJ_surf = length / 2.0; // Jakobian 

            // Pêtla po punktach ca³kowania na powierczhni
            int npc_edge = gData.npc;

            // cout << "  -> Bok " << edge << " jest brzegowy (BC). Dlugosc=" << length << endl;

            for (int p = 0; p < npc_edge; p++) {
                double waga = elemU.surfaces[edge]->wagi[p];
                double* N = elemU.surfaces[edge]->N[p]; // Wartoœci funkcji kszta³tu na krawêdzi

                // sumowanie do macierzy Hbc elementu
                // Wzór: alfa * {N} * {N}^T * waga * detJ_surf
                for (int a = 0; a < 4; a++) {
                    //wektor p sumowanie
                    gri.elements[i].P[a] += gData.Alfa * N[a] * gData.Tot * waga * detJ_surf;
                    for (int b = 0; b < 4; b++) {
                        double val = gData.Alfa * N[a] * N[b] * waga * detJ_surf;
                        // Dodajemy do Hbc
                        gri.elements[i].Hbc[a][b] += val;
                    }
                }
            }
        }
    }
}

// Agreguje macierze elementu do macierzy globalnej
void aggregateToGlobal(int i, grid& gri, SystemEquations& sysEq) {
    for (int m = 0; m < 4; m++) {
        int globalRow = gri.elements[i].ID[m] - 1;

		//agregacja wektora P do Pg
        sysEq.Pg[globalRow] += gri.elements[i].P[m];

		// Agregacja macierzy H oraz Hbc do macierzy globalnej HG
        for (int n = 0; n < 4; n++) {
            int globalCol = gri.elements[i].ID[n] - 1;

            // Sumujemy H (przewodzenie) oraz Hbc (konwekcja)
            double value = gri.elements[i].H[m][n] + gri.elements[i].Hbc[m][n];
            sysEq.HG[globalRow][globalCol] += value;

			//Sumujemy macierz C do Cg
			sysEq.Cg[globalRow][globalCol] += gri.elements[i].C[m][n];
        }
    }
}

//============================================================================================================================================
//  Glowan Petla po wszystkich elementach w siatce
//============================================================================================================================================

void runCalculations(grid& gri, GlobalData& gData, elemUniv& elemU, SystemEquations& sysEq) {
    cout << "\n======\nRozpoczynam obliczenia dla " << gri.nE << " elementow...\n======" << endl;

    //dla kazdego elementu
    for (int i = 0; i < gri.nE; i++) {
		//Wypisanie informacji o elemencie
       /* cout << "==================================" << endl;
        cout << "--- Element " << i + 1 << " [ID: "
            << gri.elements[i].ID[0] << "," << gri.elements[i].ID[1] << ","
            << gri.elements[i].ID[2] << "," << gri.elements[i].ID[3] << "]" << " ---" << endl;
        cout << "==================================" << endl;*/
        
        //oblicz lokalne macierze
        calculateElementMatrices(i, gri, gData, elemU);

        //agreguj do globalnej
        aggregateToGlobal(i, gri, sysEq);

        // Opcjonalnie: wypisywanie kontrolne co np. 10 element lub dla wszystkich
        // cout << "Element " << i+1 << " policzony." << endl;
    }

    cout << "Agregacja zakonczona." << endl;
}

void solveEquation(SystemEquations& sysEq) {
    int n = sysEq.nN;
    //double** A = sysEq.HG; // Skrót do macierzy H globalnej
    double* b = sysEq.t;   // Skrót do wektora wyników (temperatur)

    // 1. TWORZENIE KOPII MACIERZY ABY NIE NISZCZYÆ ORYGINA£U
    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        for (int j = 0; j < n; j++) {
            A[i][j] = sysEq.HG[i][j]; // Kopiujemy wartoœci z globalnej HG
        }
    }
    // =============================================================
    // KROK 1: Przygotowanie prawej strony równania
    // Równanie MES: [H]{t} + {P} = 0  ===>  [H]{t} = -{P}
    // Wektor 'b' (czyli t) na pocz¹tku s³u¿y jako prawa strona równania
    // =============================================================
    for (int i = 0; i < n; i++) {
        b[i] = sysEq.Pg[i]; // Przepisujemy Pg ze zmian¹ znaku na minus
    }

    // =============================================================
    // KROK 2: Eliminacja Gaussa z wyborem elementu g³ównego (Pivoting)
    // =============================================================
    for (int i = 0; i < n - 1; i++) {
        
        // --- Pivoting (szukanie najwiêkszego elementu w kolumnie) ---
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Jeœli znaleziono lepszy element g³ówny, zamieñ wiersze
        if (maxRow != i) {
            std::swap(b[i], b[maxRow]); // Zamiana w wektorze wyników
            for (int k = i; k < n; k++) {
                std::swap(A[i][k], A[maxRow][k]); // Zamiana wierszy w macierzy
            }
        }

        // Zabezpieczenie przed dzieleniem przez zero (macierz osobliwa)
        if (std::abs(A[i][i]) < 1e-12) {
            cerr << "Blad krytyczny: Macierz ukladu jest osobliwa (dzielenie przez zero)!" << endl;
            for (int x = 0; x < n; x++) delete[] A[x];
            delete[] A;
            return;
        }

        // --- W³aœciwa Eliminacja ---
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];

            // Zerowanie elementów pod przek¹tn¹
            // Zaczynamy pêtlê od j=i, bo elementy przed 'i' s¹ ju¿ wyzerowane
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            // Ta sama operacja na wektorze wyrazów wolnych
            b[k] -= factor * b[i];
        }
    }

    // =============================================================
    // KROK 3: Postêpowanie odwrotne (Back Substitution)
    // Obliczamy niewiadome od koñca (od do³u macierzy trójk¹tnej)
    // =============================================================
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        // Sumujemy iloczyny obliczonych ju¿ temperatur i wspó³czynników w wierszu
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * b[j];
        }
        // Wyliczamy t[i]
        b[i] = (b[i] - sum) / A[i][i];
    }
}

void adjustHg(SystemEquations& sysEq, double dt) {
    for (int m = 0; m < sysEq.nN; m++) {
        // dodanie do macierzy Hg czasu
        for (int n = 0; n < sysEq.nN; n++) {
            // Sumujemy Hg oraz C/dt
            sysEq.HG[m][n] += (sysEq.Cg[m][n] / dt);
        }
    }
}

void adjustTime(SystemEquations& sysEq, double dt, double* t0, vector<double>& oldPg) {
    vector<double> C(sysEq.nN, 0.0);
    for (int m = 0; m < sysEq.nN; m++) {
        for (int n = 0; n < sysEq.nN; n++) {
            C[m] += (sysEq.Cg[m][n] / dt) * t0[n];
        }
    }

    for(int i =0; i < sysEq.nN; i++) {
        sysEq.Pg[i] = oldPg[i] + C[i];
	}
}