//  Matrix.cpp
//  Matrix
#include "Matrix.hpp"
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdint.h>
#include <array>
#include <iterator>
using namespace std;

Matrix::Matrix(){
    rowSize = 0;
    colSize = 0;
}
Matrix::Matrix(const std::vector<std::vector<double>>& data){
    rowSize = (int)data.size();
    colSize = (int)data[0].size();
    matrix = data;
}
Status Matrix::add(const Matrix& other){
    std::vector<std::vector<double>> otherMatrix = other.matrix;
    int otherRowSize = (int)otherMatrix.size();
    int otherColSize = (int)otherMatrix[0].size();
    if((rowSize == colSize) && (otherRowSize == otherColSize)){
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix.size(); j++){
                matrix[i][j] += otherMatrix[i][j];
            }
        }
        return NoError;
    }else
        return DimensionError;
}
Status Matrix::subtract(const Matrix& other){
    std::vector<std::vector<double>> otherMatrix = other.matrix;
    int otherRowSize = (int)otherMatrix.size();
    int otherColSize = (int)otherMatrix[0].size();
    if((rowSize == colSize) && (otherRowSize == otherColSize)){
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix.size(); j++){
                matrix[i][j] -= otherMatrix[i][j];
            }
        }
        return NoError;
    }else
        return DimensionError;
}
Status Matrix::multiply(const Matrix& other){
    std::vector<std::vector<double>> otherMatrix = other.matrix;
    int otherRowSize = (int)otherMatrix.size();
    int otherColSize = (int)otherMatrix[0].size();
    std::vector<std::vector<double>> resultMatrix;
    resultMatrix = matrix;
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            resultMatrix[i][j] = 0;
        }
    }
    if((rowSize == otherColSize)&&(colSize == otherRowSize)){
        for (int i = 0; i < rowSize; i++){
            for (int j = 0; j < otherColSize; j++){
                for (int k = 0; k < otherRowSize; k++){
                    resultMatrix[i][j] += matrix[i][k] * otherMatrix[k][j];
                }
            }
        }
        matrix = resultMatrix;
        return NoError;
    }else
        return DimensionError;
}
Status Matrix::divide(double scalar){
    if(scalar != 0){
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix.size(); j++){
                matrix[i][j] /= scalar;
            }
        }
        return NoError;
    }else
        return DivideByZeroError;
}
void Matrix::multiply(double scalar){
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            matrix[i][j] *= scalar;
        }
    }
}
double calcDeterminant(std::vector<std::vector<double>> matrix){
    int determinant = 0;
    determinant = ( matrix[0][0] * ( (matrix[1][1]*matrix[2][2])-(matrix[1][2]*matrix[2][1]) ) - (matrix[0][1] * ((matrix[1][0]*matrix[2][2])-(matrix[1][2]*matrix[2][0]))) + (matrix[0][2] * ((matrix[1][0]*matrix[2][1])-(matrix[1][1]*matrix[2][0]))) );
    return determinant;
}
void Matrix::transpose(){
    Matrix transposed;
    transposed = matrix;
    transposed.zero();
    transposed.rowSize = colSize;
    transposed.colSize = rowSize;
    for (int i = 0; i < rowSize; i++){
        for (int j = 0; j < colSize; j++){
            transposed.matrix[j][i] = matrix[i][j];
        }
    }
    rowSize = transposed.rowSize;
    colSize = transposed.colSize;
    for (int i = 0; i < rowSize; i++){
        for (int j = 0; j < colSize; j++){
            matrix[i][j] = transposed.matrix[i][j];
        }
    }
}
void Matrix::zero(){
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.size(); j++){
            matrix[i][j] = 0;
        }
    }
}
void Matrix::show() const{
    for(int i = 0; i < rowSize; i++){
        for(int j = 0; j < colSize; j++){
            cout << left << setw(3) << matrix[i][j];
        }
        cout << endl;
    }
}
double Matrix::getDeterminant(){
    int determinant = 0;
    if((rowSize == 2)&&(colSize == 2)){
        determinant = ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
        return determinant;
    }else if((rowSize == 3)&&(colSize == 3)){
        determinant = ( matrix[0][0] * ( (matrix[1][1]*matrix[2][2])-(matrix[1][2]*matrix[2][1]) ) - (matrix[0][1] * ((matrix[1][0]*matrix[2][2])-(matrix[1][2]*matrix[2][0]))) + (matrix[0][2] * ((matrix[1][0]*matrix[2][1])-(matrix[1][1]*matrix[2][0]))) );
        return determinant;
    }else if((rowSize==4)&&(colSize==4)){
        std::vector<std::vector<double>> sub33matrix1 {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        std::vector<std::vector<double>> sub33matrix2 {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        std::vector<std::vector<double>> sub33matrix3 {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        std::vector<std::vector<double>> sub33matrix4 {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        sub33matrix1[0][0] = matrix[1][1];
        sub33matrix1[0][1] = matrix[1][2];
        sub33matrix1[0][2] = matrix[1][3];
        sub33matrix1[1][0] = matrix[2][1];
        sub33matrix1[1][1] = matrix[2][2];
        sub33matrix1[1][2] = matrix[2][3];
        sub33matrix1[2][0] = matrix[3][1];
        sub33matrix1[2][1] = matrix[3][2];
        sub33matrix1[2][2] = matrix[3][3];
        double det1 = calcDeterminant(sub33matrix1);
        sub33matrix2[0][0] = matrix[1][0];
        sub33matrix2[0][1] = matrix[1][2];
        sub33matrix2[0][2] = matrix[1][3];
        sub33matrix2[1][0] = matrix[2][0];
        sub33matrix2[1][1] = matrix[2][2];
        sub33matrix2[1][2] = matrix[2][3];
        sub33matrix2[2][0] = matrix[3][0];
        sub33matrix2[2][1] = matrix[3][2];
        sub33matrix2[2][2] = matrix[3][3];
        double det2 = calcDeterminant(sub33matrix2);
        sub33matrix3[0][0] = matrix[1][0];
        sub33matrix3[0][1] = matrix[1][1];
        sub33matrix3[0][2] = matrix[1][3];
        sub33matrix3[1][0] = matrix[2][0];
        sub33matrix3[1][1] = matrix[2][1];
        sub33matrix3[1][2] = matrix[2][3];
        sub33matrix3[2][0] = matrix[3][0];
        sub33matrix3[2][1] = matrix[3][1];
        sub33matrix3[2][2] = matrix[3][3];
        double det3 = calcDeterminant(sub33matrix3);
        sub33matrix4[0][0] = matrix[1][0];
        sub33matrix4[0][1] = matrix[1][1];
        sub33matrix4[0][2] = matrix[1][2];
        sub33matrix4[1][0] = matrix[2][0];
        sub33matrix4[1][1] = matrix[2][1];
        sub33matrix4[1][2] = matrix[2][2];
        sub33matrix4[2][0] = matrix[3][0];
        sub33matrix4[2][1] = matrix[3][1];
        sub33matrix4[2][2] = matrix[3][2];
        double det4 = calcDeterminant(sub33matrix4);
        double determinant = (matrix[0][0] * det1) - (matrix[0][1] * det2) + (matrix[0][2] * det3) - (matrix[0][3] * det4);
        return determinant;
    }else
        return NAN;
}
double Matrix::getAt(int row, int column) const{
    return matrix[row][column];
}
int Matrix::getRowSize() const{
    return rowSize;
}
int Matrix::getColSize() const{
    return colSize;
}
bool Matrix::isSquare() const{
    if((rowSize==colSize)&&(rowSize>0)&&(colSize>0))
        return true;
    else
        return false;
}
bool Matrix::hasSameDimensionAs(const Matrix& other){
    std::vector<std::vector<double>> otherMatrix = other.matrix;
    int otherRowSize = (int)otherMatrix.size();
    int otherColSize = (int)otherMatrix[0].size();
    if((rowSize == otherRowSize) && (colSize == otherColSize)){
        return true;
    }else
        return false;
}
Matrix Matrix::getMinor(int row, int column){
    std::vector<std::vector<double>> minor;
    minor = matrix;
    if(minor.size() > row){
        minor.erase(minor.begin() + row);
    }
    for(int i = 0; i < minor.size(); ++i){
        if(minor[i].size() > column){
            minor[i].erase(minor[i].begin() + column);
        }
    }
    matrix = minor;
    return matrix;
}
double HelpingFunctions::dotProduct(const std::vector<double>& x, const std::vector<double>& y){
    std::vector<double> vectorX = x;
    std::vector<double> vectorY = y;
    int resultofDP = 0;
    int vectorxSize = (int)vectorX.size();
    int vectorySize = (int)vectorY.size();
    if(vectorxSize == vectorySize){
        for (int i = 0; i < vectorxSize; i++){
            resultofDP += vectorX[i] * vectorY[i] ;
        }
        return resultofDP;
    }else
        return 0;
}
