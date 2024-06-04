#ifndef MATHEPROJEKT_CMYVEKTOR_H
#define MATHEPROJEKT_CMYVEKTOR_H

#include <vector>
#include <cmath>
#include <ostream>

class CMyVektor {
private:
    std::vector<double> data;

public:
    CMyVektor(int size, double initial_value = 0.0);

    double &operator[](int index);

    const double &operator[](int index) const;

    void setComponent(int index, double value) {
        if (index >= 0 && index < data.size()) {
            data[index] = value;
        }
    }

    double getComponent(int index) const {
        if (index >= 0 && index < data.size()) {
            return data[index];
        }
        return 0;
    }

    int getDimension() const {
        return data.size();
    }

    double length() const {
        double sum = 0;
        for (double val: data) {
            sum += val * val;
        }
        return sqrt(sum);
    }

    // Überladung des << Operators
    friend std::ostream &operator<<(std::ostream &os, const CMyVektor &v) {
        os << "(";
        for (int i = 0; i < v.getDimension(); ++i) {
            os << v[i];
            if (i != v.getDimension() - 1) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

CMyVektor operator+(const CMyVektor &a, const CMyVektor &b);

CMyVektor operator-(const CMyVektor &b, const CMyVektor &a);

CMyVektor operator*(double lambda, const CMyVektor &a);

CMyVektor operator*(const CMyVektor &a, double lambda);

CMyVektor gradient(CMyVektor x, double (*funktion)(CMyVektor x));

CMyVektor gradient_max(CMyVektor x, double (*funktion)(CMyVektor), double lambda = 1.0,
                       int max_steps = 25, double grad_grenze = 1e-5);


#endif //MATHEPROJEKT_CMYVEKTOR_H
