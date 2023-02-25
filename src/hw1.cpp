#include "hw1.h"
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;
using std::cin;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::setw;

Matrix algebra::zeros(size_t n, size_t m) 
{
    vector<vector<double>> Matrix (n);
    for (vector<double> &rowVector : Matrix)
    {
        rowVector = vector<double>(m, 0);
    }
    return Matrix;
}

Matrix algebra::ones(size_t n, size_t m) 
{
    vector<vector<double>> Matrix (n);
    for (vector<double> &rowVector : Matrix)
    {
        rowVector = vector<double>(m, 1);
    }
    return Matrix;
}

Matrix algebra::random(size_t n, size_t m, double min, double max) 
{
    if (min > max) {
        throw std::logic_error("min can't be greater than max!");
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    Matrix matrix;
    for (int i = 0; i < n; i++) {
        vector<double> rowVec ;
        for (int j = 0; j < m; j++) {
            rowVec.push_back(dis(gen));
        }
        matrix.push_back(rowVec);
    }
    return matrix;
}

void algebra::show(const Matrix& matrix) 
{
    for (const vector<double> &rowVec : matrix) {
        for (double element : rowVec) {
            cout << setw(7) << fixed << setprecision(3) << element << "";
        }
        cout << endl;
    }
}

Matrix algebra::multiply(const Matrix& matrix, double c)
{
    Matrix mulMatrix;
    for (const vector<double>& rowVec : matrix) {
        vector<double> mulMatrixRowVec;
        for (double element : rowVec)
        {
            mulMatrixRowVec.push_back(element * c);
        }
        mulMatrix.push_back(mulMatrixRowVec);
    }
    return mulMatrix;
}

