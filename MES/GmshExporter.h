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

    // zapisuje siatkê - raz
    void exportMesh(const grid& gri);

    // zapisuje wynik temperatury dla danej chwili czasu - w pêtli
    void exportSolution(const grid& gri, const SystemEquations& sysEq, double time, int stepNumber);
};
