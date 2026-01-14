#include "GmshExporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

GmshExporter::GmshExporter(string fname) {
    if (fname.find(".msh") == string::npos) {
        fname += ".msh";
    }
    filename = fname;

    ofstream file(filename, ios::trunc);
    file.close();
}

void GmshExporter::exportMesh(const grid& gri) {
    ofstream file(filename, ios::app); 

    if (!file.is_open()) {
        cerr << "BLAD: Nie mozna otworzyc pliku do zapisu GMSH!" << endl;
        return;
    }

    // 1. Nag³ówek formatu MSH 2.2
    file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    // 2. Wêz³y ($Nodes)
    file << "$Nodes\n";
    file << gri.nN << "\n"; // Liczba wêz³ów
    for (int i = 0; i < gri.nN; i++) {
        // Format: ID x y z
        // ID w Gmsh musi byæ od 1.  wiêc ID to i+1.
        file << (i + 1) << " " << gri.nodes[i].x << " " << gri.nodes[i].y << " 0.0\n";
    }
    file << "$EndNodes\n";

    // 3. Elementy ($Elements)
    file << "$Elements\n";
    file << gri.nE << "\n"; // Liczba elementów
    for (int i = 0; i < gri.nE; i++) {
        // Format: ID Typ(3=Quad4) LiczbaTagow(2) TagFizyczny TagGeom n1 n2 n3 n4

        file << (i + 1) << " 3 2 " << (gri.elements[i].matID + 1) << " 0 ";

        file << gri.elements[i].ID[0] << " "
            << gri.elements[i].ID[1] << " "
            << gri.elements[i].ID[2] << " "
            << gri.elements[i].ID[3] << "\n";
    }
    file << "$EndElements\n";

    file.close();
    cout << "GMSH: Zapisano geometrie siatki do " << filename << endl;
}

void GmshExporter::exportSolution(const grid& gri, const SystemEquations& sysEq, double currentTime, int stepNumber) {
    ofstream file(filename, ios::app); // Tryb dopisywania

    if (!file.is_open()) return;

    // Sekcja $NodeData - przechowuje wyniki w wêz³ach
    file << "$NodeData\n";
    file << "1\n";              // Liczba tagów string
    file << "\"Temperature\"\n"; // Nazwa pola
    file << "1\n";              // Liczba tagów real
    file << currentTime << "\n"; // Aktualny czas symulacji
    file << "3\n";              // Liczba tagów int
    file << stepNumber << "\n"; // Numer kroku
    file << "1\n";              // Liczba sk³adowych (1 = skalar temperatury)
    file << gri.nN << "\n";     // Liczba wêz³ów

    // Dane: ID_Wêz³a Wartoœæ
    for (int i = 0; i < gri.nN; i++) {
        file << (i + 1) << " " << sysEq.t[i] << "\n";
    }

    file << "$EndNodeData\n";
    file.close();
}