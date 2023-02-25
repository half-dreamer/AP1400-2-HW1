#include "hw1.h"
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;
using std::cin;
using std::cout;
using std::endl;
using std::fixed;
using std::logic_error;
using std::setprecision;
using std::setw;
using std::swap;
using std::transform;

// user-defined functions
int algebra::rowNum(const Matrix& matrix) 
{
    if (matrix.empty()) {
        return 0;
    }
    return matrix.size();
}

int algebra::colNum(const Matrix& matrix) 
{
    if (matrix.empty()) {
        return 0;
    }
    return matrix[0].size();
}

bool algebra::isEmpty(const Matrix& matrix) 
{
    return matrix.empty();
}

double& algebra::element(Matrix& matrix, int row, int col) 
{
    return matrix[row][col];
}

double algebra::dotProduct(const Matrix& matrix1, const Matrix& matrix2, 
                    int row, int col) 
{
    double dotProduct = 0;
    for (int i = 0; i < colNum(matrix1); ++i) {
        dotProduct += matrix1[row][i] * matrix2[i][col];
    }
    return dotProduct;
}
// to implemented functions
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

Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2) 
{
    if (colNum(matrix1) != rowNum(matrix2)) {
        throw logic_error("these two matrix can't be multiplied!");
    }
    Matrix muledMatrix = zeros(rowNum(matrix1), colNum(matrix2));
    for (int row = 0; row < rowNum(matrix1); ++row) {
        for (int col = 0; col < colNum(matrix2); ++col) {
            muledMatrix[row][col] = dotProduct(matrix1, matrix2, row, col);
        }
    }
    return muledMatrix;
}

Matrix algebra::sum(const Matrix& matrix, double c)
{
    Matrix sumedMatrix;
    for (vector<double> rowVec : matrix) {
        vector<double> sumedVec;
        for (double element : rowVec) {
            sumedVec.push_back(element + c);
        }
        sumedMatrix.push_back(sumedVec);
    }
    return sumedMatrix;
}

Matrix algebra::sum(const Matrix &matrix1, const Matrix &matrix2)
{
    // validate input
    if (rowNum(matrix1) != rowNum(matrix2) || colNum(matrix1) != colNum(matrix2)) {
        throw logic_error("matrixs with two different dimension can't be sumed!");
    }
    Matrix sumedMatrix;
    for (int i = 0; i < matrix1.size(); ++i) {
        vector<double> newRowVec(colNum(matrix1));
        std::transform(matrix1[i].begin(), matrix1[i].end(),
                       matrix2[i].begin(), newRowVec.begin(), std::plus<double>());
        sumedMatrix.push_back(newRowVec);
    }
    return sumedMatrix;
}

Matrix algebra::transpose(const Matrix &matrix)
{
    Matrix transedMatrix = zeros(colNum(matrix), rowNum(matrix));
    for (int i = 0; i < rowNum(matrix); ++i) {
        for (int j = 0; j < colNum(matrix); j++) {
            transedMatrix[j][i] = matrix[i][j];
        }
    }
    return transedMatrix;
}

Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m) 
{
    Matrix minoredMatrix;
    int rowCount = 0;
    for (vector<double> rowVec : matrix)
    {
        int colCount = 0;
        vector<double> newRowVec;
        if (rowCount == n) {
            ++rowCount;
            continue;
        } else {
            for (double element : rowVec) {
                if (colCount == m) {
                    ++colCount;
                    continue;
                } else {
                    newRowVec.push_back(element);
                    ++colCount;
                }
            }
        }
        minoredMatrix.push_back(newRowVec);
        ++rowCount;
    }
    return minoredMatrix;
}

double algebra::determinant(const Matrix &matrix)
{
    //validate input matrix 
    if (colNum(matrix) != rowNum(matrix)) {
        throw logic_error("input matrix is not a square!");
    }
    // corner case
    if (isEmpty(matrix)) {
        return 1;
    }
    // base case
    if (colNum(matrix) == 1) {
        return matrix[0][0];
    }
    double determinant = 0;
    int flag = 1;
    for (int i = 0; i < colNum(matrix); ++i) {
        determinant += flag * matrix[0][i] * algebra::determinant(minor(matrix, 0, i));
        flag = -flag;
    }
    return determinant;
}

Matrix algebra::inverse(const Matrix& matrix)
{
    throw "I forget how to calculate the inverse matrix🤡";
}

Matrix algebra::ero_swap(const Matrix &matrix, size_t r1, size_t r2) 
{
    if (r1 > rowNum(matrix) - 1 || r2 > rowNum(matrix) - 1) {
        throw logic_error("toSwap row number is larger than matrix size!");
    }
    Matrix swapedMatrix(matrix);
    swap(swapedMatrix[r1], swapedMatrix[r2]);
    return swapedMatrix;
}

Matrix algebra::ero_multiply(const Matrix &matrix, size_t r, double c)
{
    Matrix muledMatrix(matrix);
    transform(muledMatrix[r].begin(), muledMatrix[r].end(), muledMatrix[r].begin(),
              std::bind(std::multiplies<double>(), std::placeholders::_1, c));
    // this unary function can also use [c](auto &a){return a*c;}
    return muledMatrix;
}

Matrix algebra::ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2)
{
    Matrix sumedMatrix(matrix);
    transform(sumedMatrix[r1].begin(), sumedMatrix[r1].end(), sumedMatrix[r2].begin(), sumedMatrix[r2].begin(),
              [c](double ele1, double ele2)
              { return c * ele1 + ele2;});
    return sumedMatrix;
}
Matrix algebra::upper_triangular(const Matrix& matrix)
{
    Matrix step1 = ero_sum(matrix, )
}