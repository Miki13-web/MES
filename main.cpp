#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <functional>

using namespace std;

struct tableGauss
{
	const double* x;
	const double* w;
	int n;
};

class GaussIntegration {
private:
	// Tablice z punktami i wagami dla n = 2, 3, 4
	const double x2[2] = { -0.5773502691896257, 0.5773502691896257 };
	const double w2[2] = { 1.0, 1.0 };
	const double x3[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
	const double w3[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556 };
	const double x4[4] = { -0.8611363115940526, -0.3399810435848563,
						  0.3399810435848563, 0.8611363115940526 };
	const double w4[4] = { 0.3478548451374539, 0.6521451548625461,
						  0.6521451548625461, 0.3478548451374539 };
	tableGauss getPointsNumber(int n) const {
		if (n == 2) return { x2, w2, n };
		if (n == 3) return { x3, w3, n };
		if (n == 4) return { x4, w4, n };

		cerr << "Unsupported number of points: " << n << endl;
		return { nullptr, nullptr, 0 };
	}

public:
    // Wariant 1: Ca³kowanie 1D
    
    double integrate(int n, std::function<double(double)> f) const
    {
        tableGauss data = getPointsNumber(n);
        if (data.n == 0) return 0.0;

        double suma = 0.0;

        for (int i = 0; i < n; ++i) {
            suma += data.w[i] * f(data.x[i]);
        }

        return suma;
    }

    // Wariant 2: Ca³kowanie 2D 
    double integrate(int n, std::function<double(double, double)> f) const
    {
        tableGauss data = getPointsNumber(n);
        if (data.n == 0) return 0.0;

        double suma = 0.0;

        // Zagnie¿d¿one pêtle (n x n)
        for (int i = 0; i < n; ++i) { // Wêz³y dla X
            for (int j = 0; j < n; ++j) { // Wêz³y dla Y

                // Iloczyn wag
                double waga_iloczynowa = data.w[i] * data.w[j];

                // Sumowanie: w_i * w_j * f(data.x[i], data.x[j])
                suma += waga_iloczynowa * f(data.x[i], data.x[j]);
            }
        }
        return suma;
    }
};

	
double fun(double x) {
	return 5 * pow(x, 2) + 3 * x + 6;
}

double fun2D(double x, double y) {
	return 5 * pow(x, 2)*pow(y,2) + 3 * x*y + 6;
}

int main() {
    GaussIntegration integration;
    int n, dim;

    cout << "Podaj wymiar calkowania (1D lub 2D): ";
    cin >> dim;

    cout << "Podaj ilosc punktow (n = 2, 3 lub 4): ";
    cin >> n;

    cout << setprecision(10);

    if (dim == 1) {

        double result = integration.integrate(n, fun);
        cout << "\nWynik calkowania 1D na [" << -1 << ", " << 1 << "] wynosi: " << result << endl;

    }
    else if (dim == 2) {

        double result = integration.integrate(n, fun2D);
        cout << "\nWynik calkowania 2D na prostokacie [" << -1 << ", " << 1 << "] x [" << -1 << ", " << 1 << "] wynosi: " << result << endl;

    }
    else {
        cerr << "Nieobslugiwany wymiar calkowania.\n";
    }
	return 0;
}
