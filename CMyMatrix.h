#ifndef MATHEPROJEKT_CMYMATRIX_H
#define MATHEPROJEKT_CMYMATRIX_H


#include "CMyVektor.h"
#include <vector>


class CMyMatrix {
private:
    std::vector<std::vector<double>> matrix;

public:
    CMyMatrix(int rows, int cols);

    double& operator()(int row, int col);
    const double& operator()(int row, int col) const;

    CMyMatrix invers() const;

    friend CMyVektor operator*(const CMyMatrix &A, const CMyVektor &x);

    int getRows() const;
    int getCols() const;
    const std::vector<double>& getRow(int row) const;
};

CMyMatrix jacobi(CMyVektor x, CMyVektor (*funktion)(CMyVektor));
CMyVektor newton(CMyVektor x, CMyVektor (*funktion)(CMyVektor), int max_steps = 50, double epsilon = 1e-5);

#endif //MATHEPROJEKT_CMYMATRIX_H
