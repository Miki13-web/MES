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

using namespace std;

int main() {
    GlobalData gData;
    grid gri1;

	string filename;

    //wybor siatki

    //filename = "Test1_4_4.txt";
    //filename="Test3_31_31_kwadrat.txt";
    filename="Test2_4_4_MixGrid.txt";

    cout << "\nPodaj ile chcesz punktow calkowania: ";
    cin >> gData.npc;

    //wczytanie danych
    readData(filename, gData, gri1);

	// inicjalizacja struktur
    elemUniv elemU(gData.npc);
    SystemEquations sysEq(gData.nN);

	// sprawdzenie dla jakiego kroku czasowego temperatura nie spednie poni¿ej 100 C
    gData.SimulationStepTime = 55;
    // uruchamianie obliczen
    runCalculations(gri1, gData, elemU, sysEq);

	//============================================================================================================================================
    // Wypisywanie wynikow do sprawdzenia
	//============================================================================================================================================
	// wypisanie wyników macierz Cg
	
    /*cout << "\nMacierz C Globalna Cg:" << endl;
    for (int i = 0; i < gData.nN; i++) {
        for (int j = 0; j < gData.nN; j++) {
            cout << setprecision(3) << fixed << sysEq.Cg[i][j] << "  ";
        }
        cout << endl;
	}
	cout << endl;*/

    // wpisanie wyników macierz HG
    // Mo¿esz to te¿ przenieœæ do funkcji np. printResults(sysEq) w Logic.cpp
    /*cout << "\nMacierz Globalna HG:" << endl;
    for (int i = 0; i < gData.nN; i++) {
        for (int j = 0; j < gData.nN; j++) {
            cout << setprecision(3) << fixed << sysEq.HG[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;*/

    //wypisanie Pg
   /* for(int i=0; i<gData.nN; i++) {
        cout << "Pg[" << i << "] = " << sysEq.Pg[i] << endl;
	}
	cout << endl;*/

	//============================================================================================================================================
	// Rozwi¹zanie uk³adu równañ w pêtli po kroku czasowym
	//============================================================================================================================================

    for(int i =0; i<gData.nN; i++) {
        sysEq.t[i] = gData.InitialTemp; //ustawienie poczatkowej temperatury
	}

	cout << "\nRozwiazywanie ukladu rownan..." << endl;

	//gData.SimulationStepTime = 75; //ustawienie kroku czasowego symulacji
	adjustHg(sysEq, gData.SimulationStepTime); //modyfikacja macierzy HG o czas

    vector<double> oldPg(gData.nN);
    for (int i = 0; i < gData.nN; i++) {
        oldPg[i] = sysEq.Pg[i];
    }

    double t_min, t_max;

    for (double i = 0; i < gData.SimulationTime; i += gData.SimulationStepTime) {
        cout << "Czas symulacji: " << i + gData.SimulationStepTime << " s" << "\t";
		double* t = sysEq.t; //poprzednie temperatury do kolejnego kroku czasowego
        adjustTime(sysEq, gData.SimulationStepTime, t, oldPg);
		solveEquation(sysEq);

        // wypisanie wyników temperatury w wêz³ach
        //cout << "\nTemperatury w wezlach po rozwiazaniu ukladu rownan:" << endl;
		t_min = sysEq.t[0];
		t_max = sysEq.t[0];
        for (int i = 0; i < gData.nN; i++) {
			if (sysEq.t[i] > t_max) t_max = sysEq.t[i];
			if (sysEq.t[i] < t_min) t_min = sysEq.t[i];
           // cout << setprecision(6) << "t[" << i << "] = " << sysEq.t[i] << endl;
        }
        
        //wypisanie tylko min i max
		cout << setprecision(6) << "T min: " << t_min << " C" << "\t";
        cout << setprecision(6) << "T max: " << t_max << " C" << endl;
       
    }

    return 0;
}