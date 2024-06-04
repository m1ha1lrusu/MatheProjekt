#include "CMyVektor.h"
#include "CMyMatrix.h"
#include <iostream>

CMyVektor::CMyVektor(int size, double initial_value) : data(size, initial_value) {}

double &CMyVektor::operator[](int index) {
    return data[index];
}

const double &CMyVektor::operator[](int index) const {
    return data[index];
}

CMyVektor operator+(const CMyVektor &a, const CMyVektor &b) {
    int size = a.getDimension();
    CMyVektor result(size);
    for (int i = 0; i < size; i++) {
        result[i] = a[i] + b[i]; // Komponenten der beiden Vektoren addieren.
    }
    return result;
}

CMyVektor operator-(const CMyVektor &b, const CMyVektor &a) {
    int size = a.getDimension();
    CMyVektor result(size);
    for (int i = 0; i < size; i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

CMyVektor operator*(double lambda, const CMyVektor &a) {
    CMyVektor result(a.getDimension());
    for (int i = 0; i < a.getDimension(); i++) {
        result[i] = lambda * a[i]; // Jede Komponente von a mit lambda multiplizieren.
    }
    return result;
}

CMyVektor operator*(const CMyVektor &a, double lambda) {
    return lambda * a;
}

/// Funktion zur Berechnung des Gradienten einer Funktion an einer gegebenen Stelle
/*
 * Parameter:
 * x - Der Punkt, an dem der Gradient berechnet werden soll.
 * funktion - Ein Funktionszeiger auf die Funktion, deren Gradient berechnet wird.
 *
 * Rückgabe:
 * grad - Ein CMyVektor, der den Gradienten der Funktion am Punkt x darstellt.
 */
CMyVektor gradient(CMyVektor x, double (*funktion)(CMyVektor x)) {
    const double h = 1e-8;  // Schrittweite für die numerische Differentiation
    CMyVektor grad(x.getDimension()); // Vektor für den Gradienten

    double fx = funktion(x);        // f(x1, ..., xn). Funktionswert am ursprünglichen Punkt x.


    // Zentraldifferenz anstatt Vorwärtsdifferenz
    // Über jede Dimension des Punktes x durchlaufen, um jede partielle Ableitung zu berechnen.
    for (int i = 0; i < x.getDimension(); i++) {
        CMyVektor xh1 = x;  // Den aktuellen Vektor kopieren, um die i-te Dimension zu modifizieren.
        xh1[i] += h;  // h zur i-ten Komponente hinzufügen.
        double fxh1 = funktion(xh1);   // f(x1, ..., xi + h, ..., xn). Funktionswert an der modifizierten Stelle xh1.

        CMyVektor xh2 = x;  // Den aktuellen Vektor kopieren, um die i-te Dimension zu modifizieren.
        xh2[i] -= h;  // h von der i-ten Komponente abziehen.
        double fxh2 = funktion(xh2);   // f(x1, ..., xi - h, ..., xn). Funktionswert an der modifizierten Stelle xh2.

        grad[i] = (fxh1 - fxh2) / (2 * h);  // Die i-te Komponente des Gradienten.
    }

    return grad;
}


/// Funktion für das Gradientenverfahren mit Schrittweitensteuerung
CMyVektor gradient_max(CMyVektor x, double (*funktion)(CMyVektor), double lambda, int max_steps, double grad_grenze) {
    CMyVektor current_x = x; // Startpunkt initialisieren.
    double current_value = funktion(current_x); // Funktionswert am Startpunkt.
    int step = 0;

    while (step < max_steps) {

        // Den Gradienten am aktuellen Punkt berechnen.
        CMyVektor grad = gradient(current_x, funktion);
        if (grad.length() < grad_grenze) {
            std::cout << "\nEnde wegen ||grad f(x)||<1e-5 bei" << std::endl;
            std::cout << "\tx = " << current_x << "\n";
            std::cout << "\tlambda = " << lambda << "\n";
            std::cout << "\tf(x) = " << current_value << "\n";
            std::cout << "\tgrad f(x) = " << grad << "\n";
            std::cout << "\t||grad f(x)|| = " << grad.length() << "\n";
            break;
        }

        // Den nächsten Punkt berechnen durch Addition des Gradienten skaliert mit lambda.
        CMyVektor next_x = current_x + lambda * grad;
        // Funktionswert am neuen Punkt.
        double new_value = funktion(next_x);

        // Ausgabe der aktuellen Informationen
        std::cout << "\nSchritt " << step << ":\n";
        std::cout << "\tx = " << current_x << "\n";
        std::cout << "\tlambda = " << lambda << "\n";
        std::cout << "\tf(x) = " << current_value << "\n";
        std::cout << "\tgrad f(x) = " << grad << "\n";
        std::cout << "\t||grad f(x)|| = " << grad.length() << "\n";
        std::cout << "\n\tx_neu = " << next_x << "\n";
        std::cout << "\tf(x_neu) = " << new_value << "\n";

        if (new_value <= current_value) {
            // 2. Schrittweite solange halbieren, bis man eine Stelle mit größerem Funktionswert gefunden hat
            while (new_value <= current_value) {
                lambda /= 2;
                std::cout << "\n\thalbiere Schrittweite, lambda = " << lambda;
                next_x = current_x + lambda * grad;
                new_value = funktion(next_x);
                std::cout << "\n\tx_neu = " << next_x << "\n";
                std::cout << "\tf(x_neu) = " << new_value << "\n";
            }
            current_x = next_x; // Den aktuellen Punkt aktualisieren.
            current_value = new_value;
        }

        // 3. Verdoppelte Schrittweite testen
        double lambda_test = lambda * 2;
        CMyVektor x_test = current_x + lambda_test * grad;
        double test_value = funktion(x_test);

        std::cout << "\n\tTest mit doppelter Schrittweite, lambda = " << lambda_test << "\n";
        std::cout << "\tx_test = " << x_test << "\n";
        std::cout << "\tf(x_test) = " << test_value << "\n";

        if (test_value > new_value) {
            std::cout << "\tverdopple Schrittweite!\n";
            lambda = lambda_test;
            current_x = x_test;
            current_value = test_value;
        } else {
            std::cout << "\tbehalte alte Schrittweite!\n";
            current_x = next_x;
            current_value = new_value;
        }

        step++;

        if (step == max_steps) {
            std::cout << "\nEnde wegen Schrittanzahl = 25 bei" << std::endl;
            std::cout << "\tx = " << current_x << "\n";
            std::cout << "\tlambda = " << lambda << "\n";
            std::cout << "\tf(x) = " << current_value << "\n";
            std::cout << "\tgrad f(x) = " << grad << "\n";
            std::cout << "\t||grad f(x)|| = " << grad.length() << "\n";
            break;
        };
    }

    return current_x;
}