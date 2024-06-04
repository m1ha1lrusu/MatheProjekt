#include <iostream>
#include <valarray>
#include "CMyVektor.h"
#include "CMyMatrix.h"

using namespace std;

// Für Aufg. 2 (aus Test: Aufg. 1)
double funktionTest(CMyVektor x) {
    return 4 * x[0] * x[0] + 5 * x[1] * x[2] * x[2];
}

// Für Aufg. 3
double funktion1(CMyVektor x) {
    return sin(x[0] * x[1]) + sin(x[0]) + cos(x[1]);
}

// Für Aufg. 3
double funktion2(CMyVektor x) {
    return -(2 * x[0] * x[0] - 2 * x[0] * x[1] + x[1] * x[1] + x[2] * x[2] - 2 * x[0] - 4 * x[2]);
}

// Funktion f : R^4 -> R^3
CMyVektor funktion3(CMyVektor x) {
    CMyVektor result(3);
    result[0] = x[0] * x[1] * exp(x[2]);
    result[1] = x[1] * x[2] * x[3];
    result[2] = x[3];
    return result;
}

CMyVektor funktion_mumie(CMyVektor x) {
    CMyVektor result(2);
    result[0] = 3 * x[0] * x[0] * x[1];
    result[1] = 5 * x[2] * x[2] + 4;
    return result;
}

CMyVektor funktion4(CMyVektor x) {
    CMyVektor result(2);
    result[0] = x[0] * x[0] * x[0] * x[1] * x[1] * x[1] - 2 * x[1];
    result[1] = x[0] - 2;
    return result;
}

int main() {

    /// Praktikum 1
    /*
    ///--------------------------- 3 ------------------------------
    // Funktion 1:
    CMyVektor start1(2);
    start1[0] = 0.2;
    start1[1] = -2.1;
    cout << "\n\nStarte Gradientenverfahren fuer Funktion 1:" << endl;
    CMyVektor max1 = gradient_max(start1, funktion1);
    CMyVektor grad1 = gradient(start1, funktion1);

    // Funktion 2:
    CMyVektor start2(3);
    start2[0] = 0.0;
    start2[1] = 0.0;
    start2[2] = 0.0;
    cout << "\n\nStarte Gradientenverfahren fuer Funktion 2:" << endl;
    CMyVektor max2 = gradient_max(start2, funktion2, 0.1);
    CMyVektor grad2 = gradient(start2, funktion2);
    */

    /// Praktikum 2
    /*
    // Testen der Jacobi-Matrix
    CMyVektor x(4);
    x[0] = 1;
    x[1] = 2;
    x[2] = 0;
    x[3] = 3;

    CMyMatrix J = jacobi(x, funktion3);
    cout << "Jacobi-Matrix J = " << endl;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << J.at(i, j) << "; ";
        }
        cout << endl;
        */

    /*
    // Testen der Jacobi-Matrix mit Funktion aus Praktikum-Test
    CMyVektor x(3);
    x[0] = 3.0;
    x[1] = 1.0;
    x[2] = 2.0;

    CMyMatrix J = jacobi(x, funktion_mumie);
    cout << "\nJacobi-Matrix J = " << endl;

    for (int i = 0; i < 2; i++) {
        cout << "\t";
        for (int j = 0; j < 3; j++) {
            cout << J.at(i, j) << "\t ";
        }
        cout << endl;
    }
    */

    // Startwert
    CMyVektor start(2);
    start[0] = 1.0;
    start[1] = 1.0;

    // Newton-Verfahren
    CMyVektor nullstelle = newton(start, funktion4);

    return 0;
}