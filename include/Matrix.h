#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
using namespace std;
#include <vector>
using std::vector;

class Matrix {
  private:
    int ncol;
    int nrow;
    int ntot;
    vector<vector<double>> values;

  public:
    Matrix();
    Matrix(int);
    Matrix(int, int);
    Matrix(const Matrix &); // copying constructor
    Matrix(vector<vector<double>>);
    ~Matrix();

    int get_nrow() const;
    int get_ncol() const;
    double get_vij(int, int) const;
    void print(string) const;
    void print() const;
    void print(int, int, string) const;
    double &operator()(const int &, const int &);

    void fill();
    void fill(double);
    Matrix transpose();
    Matrix inverse();
    Matrix eigenSystem(double *);
    void Eye();
    void Random();
    void Resize(int, int);

    // Scalar operations
    Matrix operator+(double);
    Matrix operator-(double);
    Matrix operator*(double);
    Matrix operator/(double);

    // Matrix Operatios
    Matrix operator+(Matrix &);
    Matrix operator-(Matrix &);
    Matrix operator*(Matrix &);
    Matrix &operator=(const Matrix &);
};

#endif