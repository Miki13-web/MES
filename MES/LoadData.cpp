#include "LoadData.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


void readData(string filename, GlobalData& gData, grid& gri) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        exit(1); 
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
            gri.nN = gData.nN; 
        }
        else if (text == "Elements") {
            file >> text >> gData.nE;
            gri.nE = gData.nE; 
        }
        else if (text == "*Node") break; 
    }

    // alokacja pamieci dla wez³ow i elementow
    gri.nodes = new node[gri.nN];
    gri.elements = new element[gri.nE];
    gri.BC = new int[gri.nN]; 
    gri.BC_count = 0;

    // pamiec dla Jakobianów wewn¹trz elementow
    for (int i = 0; i < gri.nE; i++) {
        gri.elements[i].Jaco = new Jakobian[gData.npc*gData.npc];
    }

    getline(file, line);
    for (int i = 0; i < gri.nN; i++) {
        getline(file, line);
        if (line.empty()) { i--; continue; }
        stringstream ss(line);
        int id;
        char comma;
        ss >> id >> comma >> gri.nodes[i].x >> comma >> gri.nodes[i].y;
    }

    while (file >> text) if (text == "*Element,") break;
    string type;
    file >> type; 

    getline(file, line);
    for (int i = 0; i < gri.nE; i++) {
        getline(file, line);
        if (line.empty()) { i--; continue; }
        stringstream ss(line);
        int id;
        char comma;
        ss >> id >> comma >> gri.elements[i].ID[0] >> comma >> gri.elements[i].ID[1]
            >> comma >> gri.elements[i].ID[2] >> comma >> gri.elements[i].ID[3];
    }

    while (file >> text) if (text == "*BC") break;

    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        int id;
        char comma;
        while (ss >> id) {
            gri.BC[gri.BC_count++] = id;
            gri.nodes[id - 1].BC = true;
            if (ss.peek() == ',') ss >> comma;
        }
    }

    file.close();
}