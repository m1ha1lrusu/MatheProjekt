#include "CMyMatrix.h"
#include <stdexcept>
#include <iostream>

// Konstruktor
CMyMatrix::CMyMatrix(int rows, int cols) : matrix(rows, std::vector<double>(cols)) {}

// Zugriff auf Elemente der Matrix
double &CMyMatrix::operator()(int row, int col) {
    return matrix[row][col];
}

const double &CMyMatrix::operator()(int row, int col) const {
    return matrix[row][col];
}

// Berechnung der Inversen einer 2x2-Matrix
CMyMatrix CMyMatrix::invers() const {
    if (matrix.size() != 2 || matrix[0].size() != 2) {
        throw std::runtime_error("Inversion nur fuer 2x2-Matrizen moeglich.");
    }

    double a = matrix[0][0];
    double b = matrix[0][1];
    double c = matrix[1][0];
    double d = matrix[1][1];

    double det = a * d - b * c;
    if (det == 0) {
        throw std::runtime_error("Determinatnte = 0. Matrix kann nicht inversiert werden.");
    }

    CMyMatrix inv(2, 2);
    inv(0, 0) = d / det;
    inv(0, 1) = -b / det;
    inv(1, 0) = -c / det;
    inv(1, 1) = a / det;

    return inv;
}

// Matrix-Vektor-Multiplikation
CMyVektor operator*(const CMyMatrix &A, const CMyVektor &x) {
    int rows = A.getRows();
    int cols = A.getCols();

    if (cols != x.getDimension()) {
        throw std::runtime_error("Matrix- und Vektordimension sind nicht gleich. Multiplikation unmoeglich.");
    }

    CMyVektor result(rows);
    for (int i = 0; i < rows; ++i) {
        result[i] = 0;  // Initialisiere das Ergebnis-Element
        for (int j = 0; j < cols; ++j) {
            result[i] += A(i, j) * x[j];
        }
    }

    return result;
}

CMyMatrix jacobi(CMyVektor x, CMyVektor (*funktion)(CMyVektor)) {
    const double h = 1e-4; // Schrittweite für Differenzquotienten
    int m = x.getDimension(); // Dimension des Eingangsvektors
    int n = funktion(x).getDimension(); // Dimension des Ausgangsvektors

    CMyMatrix J(n, m); // Jacobi-Matrix

    for (int j = 0; j < m; ++j) {
        CMyVektor xh1 = x;
        xh1[j] += h;
        CMyVektor fxh1 = funktion(xh1);

        CMyVektor xh2 = x;
        xh2[j] -= h;
        CMyVektor fxh2 = funktion(xh2);

        for (int i = 0; i < n; ++i) {
            J(i, j) = (fxh1[i] - fxh2[i]) / (2 * h);
        }
    }

    return J;
}

int CMyMatrix::getRows() const {
    return matrix.size();
}

int CMyMatrix::getCols() const {
    return matrix[0].size();
}

const std::vector<double> &CMyMatrix::getRow(int row) const {
    return matrix[row];
}

// Formatierung der Matrix für die Ausgabe
std::ostream &operator<<(std::ostream &os, const CMyMatrix &matrix) {
    for (int i = 0; i < matrix.getRows(); ++i) {
        os << "\t\t";
        for (int j = 0; j < matrix.getCols(); ++j) {
            os << matrix(i, j);
            if (j < matrix.getCols() - 1) {
                os << ", ";
            }
        }
        os << std::endl;
    }
    return os;
}

CMyVektor newton(CMyVektor x, CMyVektor (*funktion)(CMyVektor), int max_steps, double epsilon) {
    for (int i = 0; i < max_steps; ++i) {
        CMyVektor fx = funktion(x);
        if (fx.length() < epsilon) {
            std::cout << "Ende wegen ||f(x)||<" << epsilon << " bei " << std::endl;
            std::cout << "\t x = " << x << std::endl;
            std::cout << "\t f(x) = " << fx << std::endl;
            std::cout << "\t ||f(x)|| = " << fx.length() << std::endl;
            return x;
        }

        CMyMatrix J = jacobi(x, funktion);
        CMyMatrix J_inv = J.invers();
        CMyVektor dx = -1 * (J_inv * fx);

        std::cout << "\nSchritt " << i << ":" << std::endl;
        std::cout << "\t x = " << x << std::endl;
        std::cout << "\t f(x) = " << fx << std::endl;
        std::cout << "\t f'(x) = " << std::endl << J;
        std::cout << "\t (f'(x))^(-1) = " << std::endl << J_inv;
        std::cout << "\t dx = " << dx << std::endl;
        std::cout << "\t ||f(x)|| = " << fx.length() << std::endl;

        x = x + dx;
    }
    std::cout << "Maximale Schrittanzahl erreicht." << std::endl;
    return x;
}