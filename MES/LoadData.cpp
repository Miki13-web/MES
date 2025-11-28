#include "LoadData.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

ifstream file("Test1_4_4.txt");
//ifstream file("Test1_4_4_inny.txt");
//ifstream file("Test3_31_31_kwadrat.txt");
//ifstream file("Test2_4_4_MixGrid.txt");

/*string filename;
cout << "Podaj nazwê pliku (np. Test2_4_4_MixGrid.txt): ";
cin >> filename;

ifstream file(filename);
*/

void readData(string filename, GlobalData& gData, grid& gri) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        exit(1); // Zakoñcz program, jeœli brak pliku
    }

    string text, line;

    // Wczytywanie zmiennych globalnych
    while (file >> text) {
        if (text == "SimulationTime") file >> gData.SimulationTime;
        else if (text == "SimulationStepTime") file >> gData.SimulationStepTime;
        else if (text == "Conductivity") file >> gData.Conductivity;
        else if (text == "Alfa") file >> gData.Alfa;
        else if (text == "Tot") file >> gData.Tot;
        else if (text == "InitialTemp") file >> gData.InitialTemp;
        else if (text == "Density") file >> gData.Density;
        else if (text == "SpecificHeat") file >> gData.SpecificHeat;
        else if (text == "Nodes") {
            file >> text >> gData.nN;
            gri.nN = gData.nN; // Przypisanie do siatki
        }
        else if (text == "Elements") {
            file >> text >> gData.nE;
            gri.nE = gData.nE; // Przypisanie do siatki
        }
        else if (text == "*Node") break; 
    }

    // Alokacja pamiêci dla wêz³ów i elementów
    // Robimy to tutaj, bo dopiero teraz znamy nN i nE
    gri.nodes = new node[gri.nN];
    gri.elements = new element[gri.nE];
    gri.BC = new int[gri.nN]; // Alokujemy z zapasem na max liczbê wêz³ów
    gri.BC_count = 0;

    // Alokacja pamiêci dla Jakobianów wewn¹trz elementów
    for (int i = 0; i < gri.nE; i++) {
        gri.elements[i].Jaco = new Jakobian[gData.npc*gData.npc];
    }

    // Wczytywanie wêz³ów 
    getline(file, line);
    for (int i = 0; i < gri.nN; i++) {
        getline(file, line);
        if (line.empty()) { i--; continue; }
        stringstream ss(line);
        int id;
        char comma;
        ss >> id >> comma >> gri.nodes[i].x >> comma >> gri.nodes[i].y;
    }

    // Przeskoczenie do sekcji *Element
    while (file >> text) if (text == "*Element,") break;
    string type;
    file >> type; 

    //Wczytywanie elementow
    getline(file, line);
    for (int i = 0; i < gri.nE; i++) {
        getline(file, line);
        if (line.empty()) { i--; continue; }
        stringstream ss(line);
        int id;
        char comma;
        // Format: ID, n1, n2, n3, n4
        ss >> id >> comma >> gri.elements[i].ID[0] >> comma >> gri.elements[i].ID[1]
            >> comma >> gri.elements[i].ID[2] >> comma >> gri.elements[i].ID[3];
    }

    // Wczytywanie BC
    while (file >> text) if (text == "*BC") break;

    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        int id;
        char comma;
        while (ss >> id) {
            gri.BC[gri.BC_count++] = id;
            // Ustawienie flagi w wêŸle
            gri.nodes[id - 1].BC = true;
            if (ss.peek() == ',') ss >> comma;
        }
    }

    file.close();
}