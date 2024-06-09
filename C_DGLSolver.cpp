#include "C_DGLSolver.h"
#include <iostream>

// Konstruktor für ein DGL-System
C_DGLSolver::C_DGLSolver(CMyVektor (*f)(CMyVektor y, double x)) : f_DGL_System(f), istSystem(true) {}

// Konstruktor für eine DGL höherer Ordnung
C_DGLSolver::C_DGLSolver(double (*f)(CMyVektor y, double x)) : f_DGL_nterOrdnung(f), istSystem(false) {}

// Ableitung von Vektor an der stelle x.
CMyVektor C_DGLSolver::ableitungen(CMyVektor y, double x) {
    if (istSystem) {
        return f_DGL_System(y, x);
    } else {
        CMyVektor result(y.getDimension());

        // Umwandeln (außer der letzten Komponente)
        for (int i = 0; i < y.getDimension() - 1; ++i) {
            result[i] = y[i + 1];
        }
        result[y.getDimension() - 1] = f_DGL_nterOrdnung(y, x);
        return result;
    }
}

// Berechnet die nächste Iteration der Lösung mit dem Euler-Verfahren.
CMyVektor C_DGLSolver::euler(CMyVektor y, double x, double h) {
    return y + h * ableitungen(y, x);
}

// Berechnet die nächste Iteration der Lösung mit dem Heun-Verfahren.
CMyVektor C_DGLSolver::heun(CMyVektor y, double x, double h) {
    CMyVektor k1 = ableitungen(y, x);
    CMyVektor k2 = ableitungen(y + h * k1, x + h); // Vorhersage für den nächsten Schritt, basierend auf k1.
    return y + (h / 2) * (k1 + k2);
}

// Methode zur Lösung der DGL mit dem Euler-Verfahren
CMyVektor C_DGLSolver::solveEuler(double xStart, double xEnd, int steps, CMyVektor yStart) {
    double h = (xEnd - xStart) / steps;
    double x = xStart;
    CMyVektor y = yStart;

    std::cout << "Euler-Verfahren:" << std::endl;
    std::cout << "h = " << h << std::endl;

    for (int i = 0; i < steps; ++i) {
        CMyVektor dydx = ableitungen(y, x);

        std::cout << "\nSchritt " << i << ":" << std::endl;
        std::cout << "\t x = " << x << std::endl;
        std::cout << "\t y = " << y << std::endl;
        std::cout << "\t y' = " << dydx << std::endl;

        y = euler(y, x, h);
        x += h;
    }
    return y;
}

// Methode zur Lösung der DGL mit dem Heun-Verfahren
void C_DGLSolver::solveHeun(double xStart, double xEnd, int steps, CMyVektor yStart) {
    double h = (xEnd - xStart) / steps;
    double x = xStart;
    CMyVektor y = yStart;

    std::cout << "Heun-Verfahren:" << std::endl;
    std::cout << "h = " << h << std::endl;

    for (int i = 0; i <= steps; ++i) {
        CMyVektor k1 = ableitungen(y, x);
        CMyVektor y_Test = y + h * k1;
        CMyVektor k2 = ableitungen(y_Test, x + h);
        CMyVektor y_Mittel = 1. / 2 * (k1 + k2);

        std::cout << "\nSchritt " << i << ":" << std::endl;
        std::cout << "\t x = " << x << std::endl;
        std::cout << "\t y = " << y << std::endl;
        std::cout << "\t y'_orig = " << k1 << std::endl;

        std::cout << "\n\t y_Test = " << y_Test << std::endl;
        std::cout << "\t y'_Test = " << k2 << std::endl;

        std::cout << "\n\t y'_Mittel = " << y_Mittel << std::endl;

        y = heun(y, x, h);
        x += h;
    }
}