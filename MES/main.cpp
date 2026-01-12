#include <iostream>
#include <iomanip>
#include <vector>
#include "GlobalData.h"
#include "Structures.h"
#include "ElemUniv.h"
#include "SystemEquations.h"
#include "Grid.h"
#include "LoadData.h"
#include "MES_logic.h" 
#include "GmshExporter.h"

using namespace std;

void przewrocRybe(grid& gri, SystemEquations& sysEq) {
    // Zakresy Y (Wysokoœæ ryby)
    double y_bottom = 0.006;
    double y_top = 0.036;
    double y_center = (y_bottom + y_top) / 2.0;

    // Zakresy X (Szerokoœæ ryby) - TO JEST NOWOŒÆ
    double x_fish_start = 0.01; // Ryba zaczyna siê od 10mm (0.01m)
    double x_fish_end = 0.06;   // Ryba koñczy siê na 60mm (0.06m)

    double epsilon = 0.0001;

    for (int i = 0; i < gri.nN; i++) {
        double y1 = gri.nodes[i].y;
        double x1 = gri.nodes[i].x;

        if (y1 >= y_bottom && y1 < y_center - epsilon &&
            x1 >= x_fish_start - epsilon && x1 <= x_fish_end + epsilon) {

            double target_y = y_top - (y1 - y_bottom);

            for (int j = 0; j < gri.nN; j++) {
                if (abs(gri.nodes[j].x - x1) < epsilon &&
                    abs(gri.nodes[j].y - target_y) < epsilon) {

                    double temp_dol = sysEq.t[i];
                    double temp_gora = sysEq.t[j];

                    sysEq.t[i] = temp_gora; 
                    sysEq.t[j] = temp_dol; 

                    break;
                }
            }
        }
    }
    cout << ">>> RYBA PRZEWROCONA! ZIMNA STRONA JEST TERAZ NA DNE. <<<" << endl;
}

int main() {
    GlobalData gData;
    grid gri1;

	string filename;

    //wybor siatki
    //siatka sandacz
    filename = "siatka_6x6_v2.txt";

    //filename = "Test1_4_4.txt";
    //filename="Test3_31_31_kwadrat.txt";
   // filename="Test2_4_4_MixGrid.txt";

    cout << "\nPodaj ile chcesz punktow calkowania: ";
    cin >> gData.npc;

    readData(filename, gData, gri1);

    elemUniv elemU(gData.npc);
    SystemEquations sysEq(gData.nN);

    //zaby bylo jak w rzeczywistosci to klade na rozgrzna patelnie wiec daje na poczatke inna temperatuire patelni
    double tempPatelni = 60.0;
    double tempSandacz = 15.0;
    double tempOtoczenia = 21.5;

    int probeNodeID = gri1.probeNodeID;
 
    for (int i = 0; i < gData.nN; i++) {
        double y = gri1.nodes[i].y;
        double x = gri1.nodes[i].x;

        // geometria 6x6 cm, dno < 6mm
        if (y <= 0.006) {
            sysEq.t[i] = tempPatelni;
        }
        else {
            if (x >= 0.01 && y <= 0.036) {
                sysEq.t[i] = tempSandacz;
            }
            else {
                sysEq.t[i] = tempOtoczenia;
            }
        }
    }
    runCalculations(gri1, gData, elemU, sysEq);

	//z eksperymentu z patelnia i woda liczê strumieñ ciepla i dodaje do globalnego wektora obci¹¿en Pg
    
    double HeatFlux = 26500;
    double dx = 0.001;

    for (int nodeID : gri1.inductionNodes) {
        sysEq.Pg[nodeID] += HeatFlux * dx;
    }

	//============================================================================================================================================
	// Rozwi¹zanie uk³adu równañ w pêtli po kroku czasowym
	//============================================================================================================================================

    //eksporter do pliku dla stworzenia symulacji
    GmshExporter exporter("wynik_sandacz.msh");
    exporter.exportMesh(gri1);                 
    exporter.exportSolution(gri1, sysEq, 0.0, 0); 


	cout << "\nRozwiazywanie ukladu rownan..." << endl;
	
	adjustHg(sysEq, gData.SimulationStepTime);

    vector<double> oldPg(gData.nN);
    for (int i = 0; i < gData.nN; i++) {
        oldPg[i] = sysEq.Pg[i];
    }

    double t_min, t_max, t_srodek_sandacza;
    bool czyJuzPrzewrocono = false;
    int step = 1;

    for (double i = 0; i < gData.SimulationTime; i += gData.SimulationStepTime) {
        cout << "Czas symulacji: " << i + gData.SimulationStepTime << " s" << "\t";
		double* t = sysEq.t; //poprzednie temperatury do kolejnego kroku czasowego
        adjustTime(sysEq, gData.SimulationStepTime, t, oldPg);
		solveEquation(sysEq);

        //zapis do pliku
        exporter.exportSolution(gri1, sysEq, i + gData.SimulationStepTime, step);
        step++;

        //wyniki temperatur w wêz³ach
        //cout << "\nTemperatury w wezlach po rozwiazaniu ukladu rownan:" << endl;
		t_min = sysEq.t[0];
		t_max = sysEq.t[0];
        for (int i = 0; i < gData.nN; i++) {
			if (sysEq.t[i] > t_max) t_max = sysEq.t[i];
			if (sysEq.t[i] < t_min) t_min = sysEq.t[i];
           // cout << setprecision(6) << "t[" << i << "] = " << sysEq.t[i] << endl;
        }

        //pobranie temp ze œrodka
        t_srodek_sandacza = sysEq.t[probeNodeID];

        if (!czyJuzPrzewrocono && t_srodek_sandacza >= 36.0) {

            przewrocRybe(gri1, sysEq); 

            czyJuzPrzewrocono = true; 
        }
        
        //wypisanie tylko min i max
        cout << setprecision(6) << "T srodek: " << t_srodek_sandacza << " C" << "\t";
		cout << setprecision(6) << "T min: " << t_min << " C" << "\t";
        cout << setprecision(6) << "T max: " << t_max << " C" << endl;
       
        if (t_srodek_sandacza >= 57.0) {
            cout << "\n>>> SANDACZ GOTOWY! Osiagnieto 57 stopni w srodku! <<<" << endl;
            break;
        }
    }

    return 0;
}