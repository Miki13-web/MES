#pragma once
#include <string>
#include "GlobalData.h"
#include "Grid.h"

// Funkcja przyjmuje nazwê pliku oraz referencje do struktur, które ma wype³niæ
void readData(std::string filename, GlobalData& gData, grid& gri);