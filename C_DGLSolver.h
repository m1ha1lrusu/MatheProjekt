#ifndef MATHEPROJEKT_C_DGLSOLVER_H
#define MATHEPROJEKT_C_DGLSOLVER_H

#include "CMyVektor.h"

class C_DGLSolver {
private:
    /* Funktionspointer für DGL-Systeme.
     * Verweist auf eine Funktion, die die rechte Seite eines DGL-Systems erster Ordnung definiert.
     */
    CMyVektor (*f_DGL_System)(CMyVektor y, double x);

    /* Funktionspointer für DGL höherer Ordnung.
     * Verweist auf eine Funktion, die die rechte Seite einer DGL höherer Ordnung definiert.
     */
    double (*f_DGL_nterOrdnung)(CMyVektor y, double x);

    // Flag zur Unterscheidung der DGL-Typen
    bool istSystem;

    CMyVektor ableitungen(CMyVektor y, double x);

public:
    // Konstruktor für ein DGL-System
    C_DGLSolver(CMyVektor (*f)(CMyVektor y, double x));

    // Konstruktor für eine DGL höherer Ordnung
    C_DGLSolver(double (*f)(CMyVektor y, double x));

    // Methode zur Lösung des DGL-Systems oder der DGL höherer Ordnung mit dem Euler-Verfahren
    CMyVektor euler(CMyVektor y, double x, double h);

    // Methode zur Lösung des DGL-Systems oder der DGL höherer Ordnung mit dem Heun-Verfahren
    CMyVektor heun(CMyVektor y, double x, double h);

    // Methode zur Lösung der DGL mit dem Euler-Verfahren über mehrere Schritte
    CMyVektor solveEuler(double xStart, double xEnd, int steps, CMyVektor yStart);

    // Methode zur Lösung der DGL mit dem Heun-Verfahren über mehrere Schritte
    void solveHeun(double xStart, double xEnd, int steps, CMyVektor yStart);
};

#endif //MATHEPROJEKT_C_DGLSOLVER_H
