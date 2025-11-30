#pragma once
#include "Structures.h"
#include "GlobalData.h"
#include "ElemUniv.h"
#include "SystemEquations.h"
#include "Grid.h"
#include <vector>

void runCalculations(grid& gri1, GlobalData& gData, elemUniv& elemU, SystemEquations& sysEq);

void solveEquation(SystemEquations& sysEq);

void adjustTime(SystemEquations& sysEq, double dt, double* t, std::vector<double>& oldPg);
void adjustHg(SystemEquations& sysEq, double dt);
