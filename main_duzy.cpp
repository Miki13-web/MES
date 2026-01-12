#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath> 

using namespace std;

struct node
{
	double x, y;
	bool BC;
};

struct Jakobian {
	double J[2][2];
	double J1[2][2];
	double detJ;
};

struct element
{
	int ID[4];
	Jakobian* Jaco; //wskaŸnik na tab jakobianow
	double H[4][4];
	//inicjuje Hbc zerami
	double Hbc[4][4] = { 0.0 };
	double P[4] = { 0.0 };
};

struct elemUniv {
	int npc;
	double** dKsi;
	double** dEta;

	double* eta_pc;
	double* ksi_pc;
	double* wagi;

	struct Surface {
		int npc_edge;
		double** N;
		double* wagi;

		Surface(int n) : npc_edge(n) {
			wagi = new double[npc_edge];
			N = new double* [npc_edge];
			for (int i = 0; i < npc_edge; i++) {
				N[i] = new double[4];
			}
		}

		~Surface() {
			for (int i = 0; i < npc_edge; i++) {
				delete[] N[i];
			}
			delete[] N;
			delete[] wagi;
		}

	};

	// Tablica 4 wskaŸników na powierzchnie (Dó³, Prawa, Góra, Lewa)
	Surface* surfaces[4];

	elemUniv(int npc) : npc(npc) {
		dKsi = new double* [npc];
		dEta = new double* [npc];
		for (int i = 0; i < npc; i++) {
			dKsi[i] = new double[4];
			dEta[i] = new double[4];
		}
		eta_pc = new double[npc];
		ksi_pc = new double[npc];
		wagi = new double[npc];

		// Punkty i wagi 1D (pomocnicze)
		int npc_edge = static_cast<int>(sqrt(npc));

		double* points_1d = new double[npc_edge];
		double* weights_1d = new double[npc_edge];

		//ustalanie punktów ca³kowania i wag w zale¿noœci od liczby punktów ca³kowania
		switch (npc) {
		case 4: { // Schemat 2x2
			double p = 1.0 / sqrt(3.0); // wspolrzedna PC 
			ksi_pc[0] = -p; eta_pc[0] = -p; // PC 1
			ksi_pc[1] = p; eta_pc[1] = -p; // PC 2 
			ksi_pc[2] = -p; eta_pc[2] = p; // PC 3 
			ksi_pc[3] = p; eta_pc[3] = p; // PC 4

			//punkty na krawedziach dla macierzy Hbc
			points_1d[0] = -p; weights_1d[0] = 1.0;
			points_1d[1] = p; weights_1d[1] = 1.0;

			// dla macierzy funkcji kszta³tu N

			for (int i = 0; i < 4; ++i) wagi[i] = 1.0; // wagi
			break;
		}
		case 9: { // Schemat 3x3
			double p1 = sqrt(3.0 / 5.0);
			double p2 = 0.0;
			double w1 = 5.0 / 9.0;
			double w2 = 8.0 / 9.0;

			ksi_pc[0] = -p1; eta_pc[0] = -p1; wagi[0] = w1 * w1;
			ksi_pc[1] = p2;  eta_pc[1] = -p1; wagi[1] = w2 * w1;
			ksi_pc[2] = p1;  eta_pc[2] = -p1; wagi[2] = w1 * w1;
			ksi_pc[3] = -p1; eta_pc[3] = p2;  wagi[3] = w1 * w2;
			ksi_pc[4] = p2;  eta_pc[4] = p2;  wagi[4] = w2 * w2;
			ksi_pc[5] = p1;  eta_pc[5] = p2;  wagi[5] = w1 * w2;
			ksi_pc[6] = -p1; eta_pc[6] = p1;  wagi[6] = w1 * w1;
			ksi_pc[7] = p2;  eta_pc[7] = p1;  wagi[7] = w2 * w1;
			ksi_pc[8] = p1;  eta_pc[8] = p1;  wagi[8] = w1 * w1;

			//punkty na krawedziach dla macierzy Hbc
			points_1d[0] = -p1; weights_1d[0] = w1;
			points_1d[1] = p2; weights_1d[1] = w2;
			points_1d[2] = p1; weights_1d[2] = w1;

			break;
		}
		default:
			cerr << "Nieobs³ugiwana liczba punktów ca³kowania: " << npc << endl;
			// Domyœlnie 4, aby unikn¹æ awarii
			npc = 4;
			double p = 1.0 / sqrt(3.0);
			ksi_pc[0] = -p; eta_pc[0] = -p;
			ksi_pc[1] = p; eta_pc[1] = -p;
			ksi_pc[2] = p; eta_pc[2] = p;
			ksi_pc[3] = -p; eta_pc[3] = p;
			for (int i = 0; i < 4; ++i) wagi[i] = 1.0;
			break;
		}

		// Obliczanie wartoœci pochodnych w ka¿dym punkcie ca³kowania
		for (int i = 0; i < npc; i++) {
			double eta = eta_pc[i];
			double ksi = ksi_pc[i];

			// Pochodne po KSI 
			dKsi[i][0] = -0.25 * (1.0 - eta);
			dKsi[i][1] = 0.25 * (1.0 - eta);
			dKsi[i][2] = 0.25 * (1.0 + eta);
			dKsi[i][3] = -0.25 * (1.0 + eta);

			// Pochodne po ETA (?)
			dEta[i][0] = -0.25 * (1.0 - ksi);
			dEta[i][1] = -0.25 * (1.0 + ksi);
			dEta[i][2] = 0.25 * (1.0 + ksi);
			dEta[i][3] = 0.25 * (1.0 - ksi);
		}

		//tworzenie 4 powierzchni
		for (int i = 0; i < 4; i++) {
			surfaces[i] = new Surface(npc_edge);
		}

		// Pêtla po 4 krawêdziach
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < npc_edge; j++) {
				double ksi = 0.0, eta = 0.0;

				// Mapowanie punktu 1D na odpowiedni¹ krawêdŸ 2D
				if (i == 0) { // Dó³ (Bottom): eta = -1, ksi zmienne
					ksi = points_1d[j];
					eta = -1.0;
				}
				else if (i == 1) { // Prawa (Right): ksi = 1, eta zmienne
					ksi = 1.0;
					eta = points_1d[j];
				}
				else if (i == 2) { // Góra (Top): eta = 1, ksi zmienne (czêsto odwrócone, ale dla skalara N nie ma znaczenia)
					ksi = points_1d[npc_edge - 1 - j]; // Opcjonalnie odwrócona kolejnoœæ
					eta = 1.0;
				}
				else if (i == 3) { // Lewa (Left): ksi = -1, eta zmienne
					ksi = -1.0;
					eta = points_1d[npc_edge - 1 - j];
				}

				// Zapisanie wagi dla punktu powierzchni
				surfaces[i]->wagi[j] = weights_1d[j];

				// Obliczenie funkcji kszta³tu N w tym punkcie brzegowym
				// Wzory: N1=0.25(1-k)(1-e), N2=0.25(1+k)(1-e), N3=0.25(1+k)(1+e), N4=0.25(1-k)(1+e)
				surfaces[i]->N[j][0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
				surfaces[i]->N[j][1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
				surfaces[i]->N[j][2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
				surfaces[i]->N[j][3] = 0.25 * (1.0 - ksi) * (1.0 + eta);
			}
		}
		delete[] points_1d;
		delete[] weights_1d;
	}

	~elemUniv() {
		// Sprz¹tanie wnêtrza
		for (int i = 0; i < npc; i++) {
			delete[] dKsi[i];
			delete[] dEta[i];
		}
		delete[] dKsi;
		delete[] dEta;
		delete[] eta_pc;
		delete[] ksi_pc;
		delete[] wagi;

		// Sprz¹tanie powierzchni
		for (int i = 0; i < 4; i++) {
			delete surfaces[i];
		}
	}
};

struct grid
{
	int nN;
	int nE;
	node* nodes;
	element* elements;
	int* BC;
	int BC_count;
};

struct GlobalData {
	double SimulationTime;
	double SimulationStepTime;
	double Conductivity;
	double Alfa;
	double Tot;
	double InitialTemp;
	double Density;
	double SpecificHeat;
	int nN;
	int nE;
	int npc;
};

// tu jest ta macierz globalna Hg która mielismy zagregwac
struct SystemEquations {
	double** HG;
	int nN;
	double* Pg;

	SystemEquations(int size) : nN(size) {
		HG = new double* [nN];
		Pg = new double[nN]();
		for (int i = 0; i < nN; i++) {
			HG[i] = new double[nN]();
		}
	}

	// Destruktor
	~SystemEquations() {
		for (int i = 0; i < nN; i++) {
			delete[] HG[i];
		}
		delete[] HG;
		delete[] Pg;
	}
};

int main() {
	GlobalData gData;
	grid gri1;
	int npc;
	ifstream file("Test1_4_4.txt");
	//ifstream file("Test1_4_4_inny.txt");
	//ifstream file("Test3_31_31_kwadrat.txt");
	//ifstream file("Test2_4_4_MixGrid.txt");

    /*string filename;
    cout << "Podaj nazwê pliku (np. Test2_4_4_MixGrid.txt): ";
    cin >> filename;

    ifstream file(filename);
*/

	cout << "\nPodaj ile chcesz punktow calkowania: "; // 4 lub 9 punktów
	cin >> npc;
	gData.npc = npc;

	if (!file.is_open()) {
		cerr << "Unable to open the file!" << endl;
		return 1;
	}

	string text, line;
	while (file >> text) {
		if (text == "SimulationTime") file >> gData.SimulationTime;
		else if (text == "SimulationStepTime") file >> gData.SimulationStepTime;
		else if (text == "Conductivity") file >> gData.Conductivity;
		else if (text == "Alfa") file >> gData.Alfa;
		else if (text == "Tot") file >> gData.Tot;
		else if (text == "InitialTemp") file >> gData.InitialTemp;
		else if (text == "Density") file >> gData.Density;
		else if (text == "SpecificHeat") file >> gData.SpecificHeat;
		else if (text == "Nodes") { file >> text >> gData.nN; }
		else if (text == "Elements") { file >> text >> gData.nE; }
		else if (text == "*Node") break;
	}

	gri1.nN = gData.nN;
	gri1.nE = gData.nE;
	gri1.nodes = new node[gri1.nN];
	gri1.elements = new element[gri1.nE];
	gri1.BC = new int[gri1.nN];
	gri1.BC_count = 0;

	// pamiec dla jakobianow
	for (int i = 0; i < gri1.nE; i++) {
		// Dla ka¿dego elementu pamiêæ na jego Jakobiany (po 1 na PC)
		gri1.elements[i].Jaco = new Jakobian[gData.npc];
	}

	getline(file, line);
	for (int i = 0; i < gri1.nN; i++) {
		getline(file, line);
		if (line.empty()) { i--; continue; }
		stringstream ss(line);
		int id;
		char comma;
		ss >> id >> comma >> gri1.nodes[i].x >> comma >> gri1.nodes[i].y;
	}

	while (file >> text) if (text == "*Element,") break;
	string type;
	file >> type;

	getline(file, line);
	for (int i = 0; i < gri1.nE; i++) {
		getline(file, line);
		if (line.empty()) { i--; continue; }
		stringstream ss(line);
		int id;
		char comma;
		ss >> id >> comma >> gri1.elements[i].ID[0] >> comma >> gri1.elements[i].ID[1]
			>> comma >> gri1.elements[i].ID[2] >> comma >> gri1.elements[i].ID[3];
	}

	while (file >> text) if (text == "*BC") break;
	while (getline(file, line)) {
		if (line.empty()) continue;
		stringstream ss(line);
		int id;
		char comma;
		while (ss >> id) {
			gri1.BC[gri1.BC_count++] = id;
			// zaznaczenie w strukturze wezla
			gri1.nodes[id - 1].BC = true;
			if (ss.peek() == ',') ss >> comma;
		}
	}

	file.close();

	cout << "=== DANE GLOBALNE ===" << endl;
	cout << "SimulationTime: " << gData.SimulationTime << endl;
	cout << "SimulationStepTime: " << gData.SimulationStepTime << endl;
	cout << "Conductivity: " << gData.Conductivity << endl;
	cout << "Alfa: " << gData.Alfa << endl;
	cout << "Tot: " << gData.Tot << endl;
	cout << "InitialTemp: " << gData.InitialTemp << endl;
	cout << "Density: " << gData.Density << endl;
	cout << "SpecificHeat: " << gData.SpecificHeat << endl;
	cout << "Liczba wezlow (nN): " << gData.nN << endl;
	cout << "Liczba elementow (nE): " << gData.nE << endl;

	cout << "\n======\nWezly\n======" << endl;
	for (int i = 0; i < gri1.nN; i++)
		cout << i + 1 << ": (" << setprecision(10) << fixed << gri1.nodes[i].x << ", " << gri1.nodes[i].y << ")\n";

	cout << "\n======\nElementy\n======" << endl;
	for (int i = 0; i < gri1.nE; i++) {
		cout << i + 1 << ": ";
		for (int j = 0; j < 4; j++) cout << setprecision(10) << fixed << gri1.elements[i].ID[j] << " ";
		cout << endl;
	}

	cout << "\n======\nWarunki brzegowe (BC)\n======" << endl;
	for (int i = 0; i < gri1.BC_count; i++)
		cout << setprecision(10) << fixed << gri1.BC[i] << " ";
	cout << endl;

	// element uniwersalny
	elemUniv elemU(gData.npc);

	//macierz HG
	SystemEquations sysEq(gData.nN);

	cout << "\nKolejnosc punktow calkowania\n";
	for (int i = 0; i < npc; i++) {
		cout << "Punkt " << i + 1 << ": Ksi = " << elemU.ksi_pc[i] << ", Eta = " << elemU.eta_pc[i] << endl;
	}


	cout << "\nPochodne dN/dKsi" << endl;
	for (int i = 0; i < npc; i++) {
		cout << elemU.dKsi[i][0] << " " << elemU.dKsi[i][1] << " " << elemU.dKsi[i][2] << " " << elemU.dKsi[i][3] << " " << endl;
	}

	cout << "\nPochodne dN/dEta" << endl;
	for (int i = 0; i < npc; i++) {
		cout << elemU.dEta[i][0] << " " << elemU.dEta[i][1] << " " << elemU.dEta[i][2] << " " << elemU.dEta[i][3] << " " << endl;
	}

	cout << "\n======\nObliczenia Jakobianu dla Elementow\n======" << endl;

	//============================================================================================================================================
	//  Pêtla po wszystkich elementach w siatce
	//============================================================================================================================================
	for (int i = 0; i < gri1.nE; i++) {
		cout << "==================================" << endl;
		cout << "--- Element " << i + 1 << " [ID: "
			<< gri1.elements[i].ID[0] << "," << gri1.elements[i].ID[1] << ","
			<< gri1.elements[i].ID[2] << "," << gri1.elements[i].ID[3] << "]" << " ---" << endl;
		cout << "==================================" << endl;

		int nr_elementu = i;

		//zerowanie macierzy H dla elementu i Hbc
		for (int m = 0; m < 4; m++) {
			//zerowanie P lokalny
			gri1.elements[i].P[m] = 0.0;
			for (int n = 0; n < 4; n++) {
				gri1.elements[i].H[m][n] = 0.0;
				gri1.elements[i].Hbc[m][n] = 0.0;
			}
		}

		//  Pêtla po wszystkich punktach ca³kowania (PC) dla danego elementu
		for (int j = 0; j < gData.npc; j++) {
			cout << "\n  Punkt Calkowania (PC) " << j + 1 << " (ksi=" << elemU.ksi_pc[j] << ", eta=" << elemU.eta_pc[j] << "):" << endl;

			int numer_wagi = j;

			node nodes_xy[4];
			for (int k = 0; k < 4; k++) {
				int nodeId = gri1.elements[i].ID[k] - 1;
				nodes_xy[k] = gri1.nodes[nodeId];
			}

			// lokalne pochodne dla bie¿¹cego PC
			double* dN_dKsi = elemU.dKsi[j];
			double* dN_dEta = elemU.dEta[j];

			// komponenty Jakobianu 
			double dxdksi = 0.0, dydksi = 0.0, dxdeta = 0.0, dydeta = 0.0;

			for (int k = 0; k < 4; k++) {
				dxdksi += dN_dKsi[k] * nodes_xy[k].x;
				dydksi += dN_dKsi[k] * nodes_xy[k].y;
				dxdeta += dN_dEta[k] * nodes_xy[k].x;
				dydeta += dN_dEta[k] * nodes_xy[k].y;
			}


			Jakobian& jaco = gri1.elements[i].Jaco[j];
			jaco.J[0][0] = dxdksi;
			jaco.J[0][1] = dydksi;
			jaco.J[1][0] = dxdeta;
			jaco.J[1][1] = dydeta;

			// wyznacznik
			jaco.detJ = (jaco.J[0][0] * jaco.J[1][1]) - (jaco.J[0][1] * jaco.J[1][0]);

			// macierz odwrotna J^-1
			double invDetJ = 1.0 / jaco.detJ;
			jaco.J1[0][0] = jaco.J[1][1] * invDetJ;
			jaco.J1[0][1] = -jaco.J[0][1] * invDetJ;
			jaco.J1[1][0] = -jaco.J[1][0] * invDetJ;
			jaco.J1[1][1] = jaco.J[0][0] * invDetJ;

			// wyniki dla tego PC
			cout << "    Jakobian J = [ " << setw(12) << jaco.J[0][0] << ", " << setw(12) << jaco.J[0][1] << " ]" << endl;
			cout << "                 [ " << setw(12) << jaco.J[1][0] << ", " << setw(12) << jaco.J[1][1] << " ]" << endl;
			cout << "    det(J)     = " << jaco.detJ << endl;
			cout << "    J^-1       = [ " << setw(12) << jaco.J1[0][0] << ", " << setw(12) << jaco.J1[0][1] << " ]" << endl;
			cout << "                 [ " << setw(12) << jaco.J1[1][0] << ", " << setw(12) << jaco.J1[1][1] << " ]" << endl;

			            /*
			            // Wypisywanie pochodnych lokalnych (po ksi i eta)
			            cout << "    Pochodne lokalne {dN/d(ksi)} [dKsi]:" << endl;
			            cout << "      [ ";
			            for (int k = 0; k < 4; k++) {
			                // dN_dKsi to tablica double* dN_dKsi = elemU.dKsi[j];
			                cout << setw(10) << fixed << setprecision(6) << dN_dKsi[k] << (k == 3 ? " ]" : ",");
			            }
			            cout << endl;

			            cout << "    Pochodne lokalne {dN/d(eta)} [dEta]:" << endl;
			            cout << "      [ ";
			            for (int k = 0; k < 4; k++) {
			                // dN_dEta to tablica double* dN_dEta = elemU.dEta[j];
			                cout << setw(10) << fixed << setprecision(6) << dN_dEta[k] << (k == 3 ? " ]" : ",");
			            }
			            cout << endl;
			            */

						            /*
						             // pochodne globalne (dN/dx, dN/dy)
						             cout << "    Pochodne globalne {dN/dx} i {dN/dy}:" << endl;
						             for (int k = 0; k < 4; k++) { // Pêtla po 4 funkcjach kszta³tu (N1..N4)
						                 double dNdx = jaco.J1[0][0] * dN_dKsi[k] + jaco.J1[0][1] * dN_dEta[k];
						                 double dNdy = jaco.J1[1][0] * dN_dKsi[k] + jaco.J1[1][1] * dN_dEta[k];
						                 cout << "      N" << k + 1 << ": (dN/dx=" << setw(10) << dNdx << ", dN/dy=" << setw(10) << dNdy << ")" << endl;
						             }
						             cout << endl;
						             */

									 // Obliczamy dla kazdego punktu pochodne globalne
			double dNdx[4], dNdy[4];//tablica na pochodne globalne dla punktu
			for (int k = 0; k < 4; k++) { // Pêtla po 4 funkcjach kszta³tu (N1..N4)
				dNdx[k] = jaco.J1[0][0] * dN_dKsi[k] + jaco.J1[0][1] * dN_dEta[k];
				dNdy[k] = jaco.J1[1][0] * dN_dKsi[k] + jaco.J1[1][1] * dN_dEta[k];
			}

			// macierz H dla punktu
			double H_pc[4][4];
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H_pc[i][j] = (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j]) * gData.Conductivity * jaco.detJ;
					//dodaje sobie odrazu do macierzy dla elementu
					gri1.elements[nr_elementu].H[i][j] += H_pc[i][j] * elemU.wagi[numer_wagi];
				}
			}

			// wypisanie macierzy H dla punktu
			cout << "\n    Macierz H dla tego punktu calkowania:" << endl;
			for (int m = 0; m < 4; m++) {
				cout << "      [ ";
				for (int n = 0; n < 4; n++) {
					cout << setw(12) << fixed << setprecision(6) << H_pc[m][n] << (n == 3 ? " ]" : ",");
				}
				cout << endl << endl;
			}
		}
		// ==============================================
		// OBLICZANIE MACIERZY H_BC (KONWEKCJA)
		// ==============================================


		// Bok 0 (Dó³): wêz³y 0-1
		// Bok 1 (Prawa): wêz³y 1-2
		// Bok 2 (Góra): wêz³y 2-3
		// Bok 3 (Lewa): wêz³y 3-0
		int edge_nodes[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

		// Pêtla po 4 bokach elementu
		for (int edge = 0; edge < 4; edge++) {

			// Pobierz lokalne indeksy wêz³ów tworz¹cych ten bok (0, 1, 2 lub 3)
			int local_n1 = edge_nodes[edge][0];
			int local_n2 = edge_nodes[edge][1];

			// Pobierz ich globalne ID (z pliku, np. 1, 5, 16...)
			int id1 = gri1.elements[i].ID[local_n1];
			int id2 = gri1.elements[i].ID[local_n2];

			// SPRAWDZENIE WARUNKU BRZEGOWEGO:
			// flag bool BC ze struktury node
			if (gri1.nodes[id1 - 1].BC && gri1.nodes[id2 - 1].BC) {

				// d³ugoœæ boku (potrzebne do Jakobianu 1D
				node n1 = gri1.nodes[id1 - 1];
				node n2 = gri1.nodes[id2 - 1];
				double length = sqrt(pow(n2.x - n1.x, 2) + pow(n2.y - n1.y, 2));

				double detJ_surf = length / 2.0; // Jakobian 

				// Pêtla po punktach ca³kowania NA POWIERZCHNI
				// U¿ywamy struktury surfaces z elemUniv
				int npc_edge = sqrt(gData.npc);

				// cout << "  -> Bok " << edge << " jest brzegowy (BC). Dlugosc=" << length << endl;

				for (int p = 0; p < npc_edge; p++) {
					double waga = elemU.surfaces[edge]->wagi[p];
					double* N = elemU.surfaces[edge]->N[p]; // Wartoœci funkcji kszta³tu na krawêdzi

					// sumowanie do macierzy Hbc elementu
					// Wzór: alfa * {N} * {N}^T * waga * detJ_surf
					for (int a = 0; a < 4; a++) {
						for (int b = 0; b < 4; b++) {
							double val = gData.Alfa * N[a] * N[b] * waga * detJ_surf;
							// ZMIANA: Dodajemy do Hbc, a nie do H
							gri1.elements[i].Hbc[a][b] += val;
							// Dodajemy bezpoœrednio do macierzy H elementu
							//gri1.elements[i].H[a][b] += val;
						}
					}
				}
			}
		}


		// Wypisanie skumulowanej macierzy H dla elementu
		cout << "\n################################################" << endl;
		cout << "  Skumulowana macierz H dla elementu " << i + 1 << ":" << endl;
		cout << "################################################" << endl;
		for (int m = 0; m < 4; m++) {
			cout << "      [ ";
			for (int n = 0; n < 4; n++) {
				cout << setw(12) << fixed << setprecision(6) << gri1.elements[i].H[m][n] << (n == 3 ? " ]" : ",");
			}
			cout << endl << endl;
		}

		// wypisanie hbc
		cout << "\n################################################" << endl;
		cout << "  Macierz Hbc dla elementu " << i + 1 << ":" << endl;
		cout << "################################################" << endl;
		for (int m = 0; m < 4; m++) {
			cout << "      [ ";
			for (int n = 0; n < 4; n++) {
				cout << setw(12) << fixed << setprecision(6) << gri1.elements[i].Hbc[m][n] << (n == 3 ? " ]" : ",");
			}
			cout << endl << endl;
		}

		// agregacja do macierzy HG
		for (int m = 0; m < 4; m++) {
			int globalRow = gri1.elements[i].ID[m] - 1; // indeks globalny w macierzy HG musi byc -1 bo w siatce id wezlow zaczyna sie od 1 a tablice indeksuje sie od 0
			for (int n = 0; n < 4; n++) {
				int globalCol = gri1.elements[i].ID[n] - 1; // indeks globalny w macierzy HG z dodaniem Hbc
				sysEq.HG[globalRow][globalCol] += (gri1.elements[i].H[m][n]);
			}
		}


	}

	// wypisanie macierzy hg do sprawdzenia
	cout << "\n=================================\nMacierz globalna HG po agregacji wszystkich elementow\n=================================" << endl;
	for (int i = 0; i < gData.nN; i++) {
		for (int j = 0; j < gData.nN; j++) {
			cout << setw(12) << fixed << setprecision(6) << sysEq.HG[i][j] << "  ;  ";
		}
		cout << endl;
	}

	//zwolnienie pamieci
	delete[] gri1.nodes;

	// Najpierw tablice Jakobianów w ka¿dym elemencie
	for (int i = 0; i < gri1.nE; i++) {
		delete[] gri1.elements[i].Jaco;
	}
	// Dopiero potem tablica elementów
	delete[] gri1.elements;

	delete[] gri1.BC;

	return 0;

}

mam tutaj policzyæ macierz hbc dla ka¿dego elementu ale mi nie dzia³a i jest taka sama niezaleznie od elementu