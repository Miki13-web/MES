#pragma once
#include <string>
#include <vector>
#include "Grid.h"
#include "SystemEquations.h"

class GmshExporter {
private:
    std::string filename;

public:
    GmshExporter(std::string fname);

    // Zapisuje geometriê (Wêz³y i Elementy) - raz
    void exportMesh(const grid& gri);

    // Zapisuje wynik temperatury dla danej chwili czasu - w pêtli
    void exportSolution(const grid& gri, const SystemEquations& sysEq, double time, int stepNumber);
};
