#include "LoadData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

struct MaterialDef {
    double conductivity, density, specificHeat;
};

void readData(string filename, GlobalData& gData, grid& gri) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Nie mozna otworzyc pliku: " << filename << endl;
        exit(1); 
    }

    //wlasnosci materialow
    vector<MaterialDef> materials = {
        {0.026, 1.2,    1005.0}, // ID 0: Powietrze
        {22.5,  7740.0, 440.0},  // ID 1: Stal (Yato 21/0)
        {220.0, 2700.0, 910.0},  // ID 2: Aluminium
        { 0.56, 1060.0, 3340.0 },  // ID 3: Sandacz
        {1.05,  2500.0, 840.0}  // ID 4: SZKLANA PRZYKRYWKA PATELNI
    };

   //woda do testow {0.6,   1000.0, 4186.0},
    string text, line;

    // Wczytywanie zmiennych globalnych
    while (file >> text) {
        if (text == "SimulationTime") file >> gData.SimulationTime;
        else if (text == "SimulationStepTime") file >> gData.SimulationStepTime;
        else if (text == "Alfa") file >> gData.Alfa;
        else if (text == "Tot") file >> gData.Tot;
        else if (text == "InitialTemp") file >> gData.InitialTemp;
        else if (text == "Nodes") {
            file >> text >> gData.nN;
            gri.nN = gData.nN; 
        }
        else if (text == "Elements") {
            file >> text >> gData.nE;
            gri.nE = gData.nE; 
        }
        else if (text == "*Node") break; 
        //to bylo kiedys potrzebne
        // else if (text == "Conductivity") file >> gData.Conductivity;
        //else if (text == "Density") file >> gData.Density;
        //else if (text == "SpecificHeat") file >> gData.SpecificHeat;
    }

    // alokacja pamieci dla wez³ow i elementow
    gri.nodes = new node[gri.nN];
    gri.elements = new element[gri.nE];
    gri.BC = new int[gri.nN]; 
    gri.BC_count = 0;

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
        int id, matID;
        char comma;
        ss >> id >> comma >> gri.elements[i].ID[0] >> comma >> gri.elements[i].ID[1]
            >> comma >> gri.elements[i].ID[2] >> comma >> gri.elements[i].ID[3] >> comma >> matID;

        if (matID < 0 || matID >= materials.size()) matID = 0;

        gri.elements[i].matID = matID;
        gri.elements[i].conductivity = materials[matID].conductivity;
        gri.elements[i].density = materials[matID].density;
        gri.elements[i].specificHeat = materials[matID].specificHeat;
    }

    while (file >> text) {
        //konwekcja
        if (text == "*BC") {
            while (file.peek() != '*' && !file.eof()) {
                string l;
                getline(file, l);
                if (l.empty()) continue;

                stringstream ss(l);
                int id;
                char comma;
                while (ss >> id) {
                    if (gri.BC_count < gri.nN) {
                        gri.BC[gri.BC_count++] = id;
                        gri.nodes[id - 1].BC = true;
                    }
                    if (ss.peek() == ',') ss >> comma;
                }
            }
        }
        // indukcja patelni
        else if (text == "*BC_Induction") {
            while (file.peek() != '*' && !file.eof()) {
                string l;
                getline(file, l);
                if (l.empty()) continue;

                stringstream ss(l);
                int id;
                char comma;
                while (ss >> id) {
                    gri.inductionNodes.push_back(id - 1);
                    if (ss.peek() == ',') ss >> comma;
                }
            }
        }
        else if (text == "*Probe_Node" || text == "*Probe_Node\n") {
            int id;
            file >> id;
            gri.probeNodeID = id - 1;
        }
    }

    file.close();
}