#include <iostream>
#include <iomanip>
#include <valarray>
#include "CMyVektor.h"
#include "CMyMatrix.h"
#include "C_DGLSolver.h"

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

// Funktion für das DGL-System
CMyVektor dglSystem(CMyVektor y, double x) {
    CMyVektor dydx(2);
    dydx[0] = 2 * y[1] - x * y[0];
    dydx[1] = y[0] * y[1] - 2 * x * x * x;
    return dydx;
}

// Funktion für die DGL dritter Ordnung
double dglDritterOrdnung(CMyVektor y, double x) {
    return 2 * x * y[2] + 2 * y[0] * y[1];
}

double dglDritterOrdnung2(CMyVektor y, double x) {
    return 2 * x * y[1] * y[2] + 2 * y[0] * y[0] * y[1];
}

// Exakte Lösung der DGL dritter Ordnung
double exaktLoesung(double x) {
    return 1. / x;
}

void berechneAbweichung(double xStart, double xEnd, int steps, CMyVektor yStart, double (*exaktLoesung)(double),
                        C_DGLSolver &solver, const std::string &method) {
    double h = (xEnd - xStart) / steps;
    double x = xStart;
    CMyVektor y = yStart;

    for (int i = 0; i < steps; ++i) {
        if (method == "Euler") {
            y = solver.euler(y, x, h);
        } else if (method == "Heun") {
            y = solver.heun(y, x, h);
        }
        x += h;
    }

    double y_numerisch = y[0];
    double y_exakt = exaktLoesung(x);
    double abweichung = y_numerisch - y_exakt;

    std::cout << "Abweichung bei " << method << " bei " << steps << " Schritten: " << abweichung
              << std::endl;
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
            cout << J(i, j) << "; ";
        }
        cout << endl;


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

/*
        // Startwert
        CMyVektor start(2);
        start[0] = 1.0;
        start[1] = 1.0;

        // Newton-Verfahren
        CMyVektor nullstelle = newton(start, funktion4);
*/

    ///Praktikum 3
    // Test des DGL-Systems
    CMyVektor yStartSystem(2);
    yStartSystem[0] = 0; // y1(0) = 0
    yStartSystem[1] = 1; // y2(0) = 1
    C_DGLSolver solverSystem(dglSystem);

    //std::cout << "\nTest des DGL-Systems mit Euler-Verfahren:" << std::endl;
    //CMyVektor result = solverSystem.solveEuler(1, 2, 100, yStartSystem);
    //std::cout << result[0] - 0.5;

    //std::cout << "\nTest des DGL-Systems mit Heun-Verfahren:" << std::endl;
    //solverSystem.solveHeun(0, 2, 100, yStartSystem);


    // Test der DGL dritter Ordnung
    CMyVektor yStartDritterOrdnung(3);
    yStartDritterOrdnung[0] = 1; // y(1)
    yStartDritterOrdnung[1] = -1; // y'(1)
    yStartDritterOrdnung[2] = 2; // y''(1)
    C_DGLSolver solverDritterOrdnung(dglDritterOrdnung2);

    std::cout << "\nTest der DGL dritter Ordnung mit Euler-Verfahren:" << std::endl;
    CMyVektor result = solverDritterOrdnung.solveEuler(1, 2, 10000, yStartDritterOrdnung);
    std::cout << result[0] - 0.5 << std::endl;

    //std::cout << "\nTest der DGL dritter Ordnung mit Heun-Verfahren:" << std::endl;
    //solverDritterOrdnung.solveHeun(1, 2, 100, yStartDritterOrdnung);


    /*
    // Abweichungen
    berechneAbweichung(1, 2, 10, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Euler");
    berechneAbweichung(1, 2, 10, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Heun");


    berechneAbweichung(1, 2, 100, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Euler");
    berechneAbweichung(1, 2, 100, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Heun");

    berechneAbweichung(1, 2, 1000, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Euler");
    berechneAbweichung(1, 2, 1000, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Heun");

    berechneAbweichung(1, 2, 10000, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Euler");
    berechneAbweichung(1, 2, 10000, yStartDritterOrdnung, exaktLoesung, solverDritterOrdnung, "Heun");
    */
    return 0;
}
