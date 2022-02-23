#include "Matrix.h"
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <lapacke.h>

using namespace std;

Matrix::Matrix() {
    ncol = 1;
    nrow = 1;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
};

Matrix::Matrix(int n) {
    ncol = n;
    nrow = n;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
}

Matrix::Matrix(int n, int m) {
    nrow = n;
    ncol = m;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
}
Matrix::Matrix(const Matrix &matB) {
    nrow = matB.nrow;
    ncol = matB.ncol;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            values[i][j] = matB.values[i][j];
        }
    }
}
Matrix::Matrix(vector<vector<double>> list) {
    nrow = list.size();
    int nmax = -1;
    for (int i = 0; i < nrow; i++) {
        if (int(list[i].size()) > nmax)
            nmax = int(list[i].size());
    }
    ncol = nmax;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < int(list[i].size()); j++) {
            values[i][j] = list[i][j];
        }
    }
}

Matrix::~Matrix() {}

int Matrix::get_ncol() const { return this->ncol; }

int Matrix::get_nrow() const { return this->nrow; }

double Matrix::get_vij(int i, int j) const { return this->values[i][j]; }

double &Matrix::operator()(const int &i, const int &j) {
    return this->values[i][j];
}

void Matrix::print(string st) const {
    cout << "  Matrix: " << st << endl;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            cout << setw(10) << setprecision(6) << fixed << values[i][j];
        }
        cout << endl;
    }
}

void Matrix::print() const {
    cout << "  Matrix: " << endl;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            cout << setw(10) << setprecision(6) << fixed << values[i][j];
        }
        cout << endl;
    }
}
void Matrix::print(int _row, int _col, string name) const {
    cout << "  Matrix: " << name << endl;
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _col; j++) {
            cout << setw(18) << setprecision(6) << fixed << scientific
                 << values[i][j];
        }
        cout << endl;
    }
}
void Matrix::fill() {
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            values[i][j] = 10 * (i + 1) + (j + 1);
        }
    }
}

void Matrix::fill(double k) {
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            values[i][j] = k;
        }
    }
}

void Matrix::Random() {
    srand((unsigned int)time(NULL));
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            values[i][j] = (double(rand()) / double(RAND_MAX));
        }
    }
}

Matrix Matrix::transpose() {
    Matrix matT(ncol, nrow);

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            matT(j, i) = this->values[i][j];
        }
    }
    return matT;
}

Matrix Matrix::inverse() {
    Matrix mInv(ncol, nrow);
    int m = nrow;
    int n = ncol;
    int nn = ntot;
    int error;
    int *pivot = new int[nrow];
    double *mat = new double[ntot];
    double *works = new double[ntot];

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            mat[i * ncol + j] = values[i][j];
        }
    }

    dgetrf_(&m, &n, mat, &n, pivot, &error);
    // cout << " dgetrf error : " << error << " , should be zero" << endl;
    if (error != 0) {
        cout << "  [ERROR] dgetrf FAILS to factorize! " << endl;
    }

    dgetri_(&n, mat, &n, pivot, works, &nn, &error);
    // cout << " dgetri error : " << error << " , should be zero" << endl;
    if (error != 0) {
        cout << "  [ERROR] dgetri FAILS to compute the inverse! " << endl;
    }

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            mInv(i, j) = mat[i * ncol + j];
        }
    }

    delete[] mat;

    delete[] pivot;

    delete[] works;

    return mInv;
}

Matrix Matrix::eigenSystem(double *eval) {

    if (ncol != nrow) {
        cout << " The matrix is not squared!" << endl;
        return *this;
    }
    Matrix mevec(ncol, ncol);
    const char jobz = 'V'; // Eigenvalues and eigenvectors
    const char uplo = 'U'; // Upper triangular

    int n = ncol;
    int lwork = ntot;
    int error;
    double *mat = new double[ntot];
    double *works = new double[lwork];
    if (!eval) {
        eval = new double[ncol];
    }

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            mat[i * ncol + j] = values[i][j];
        }
    }

    dsyev_(&jobz, &uplo, &n, mat, &n, eval, works, &lwork, &error);
    if (error != 0) {
        cout << "  [ERROR] dsyev FAILS to diagonalize!" << endl;
    }

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            mevec(i, j) = mat[i * ncol + j];
        }
    }

    delete[] mat;

    delete[] works;

    return mevec;
}

void Matrix::Eye() {
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (i == j)
                this->values[i][j] = 1.;
        }
    }
}

void Matrix::Resize(int irow, int icol) {
    ncol = irow;
    nrow = icol;
    ntot = ncol * nrow;

    values.resize(nrow);
    for (int i = 0; i < nrow; i++) {
        values[i].resize(ncol, 0.);
    }
}
// Scalar Operations
Matrix Matrix::operator+(double val) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = this->values[i][j] + val;
        }
    }
    return res;
}
Matrix Matrix::operator-(double val) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = this->values[i][j] - val;
        }
    }
    return res;
}

Matrix Matrix::operator*(double k) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = k * (this->values[i][j]);
        }
    }
    return res;
}

Matrix Matrix::operator/(double k) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = (this->values[i][j]) / k;
        }
    }
    return res;
}

// Matrix Operation
Matrix Matrix::operator+(Matrix &B) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = this->values[i][j] + B(i, j);
        }
    }
    return res;
}
Matrix Matrix::operator-(Matrix &B) {
    Matrix res(nrow, ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            res(i, j) = this->values[i][j] - B(i, j);
        }
    }
    return res;
}

Matrix Matrix::operator*(Matrix &B) {
    Matrix res(nrow, B.get_ncol());

    if (ncol == B.get_nrow()) {
        int i, j, k;
        double tmp;
        for (i = 0; i < nrow; i++) {
            for (j = 0; j < B.get_ncol(); j++) {
                tmp = 0.;
                for (k = 0; k < ncol; k++) {
                    tmp += values[i][k] * B(k, j);
                }
                res(i, j) = tmp;
            }
        }
    } else {
        cout << "Error the matrix do not have the correct size!" << endl;
    }
    return res;
}

Matrix &Matrix::operator=(const Matrix &matA) {
    // Matrix res(matA.nrow, matA.ncol);
    this->nrow = matA.nrow;
    this->ncol = matA.ncol;
    this->ntot = matA.ntot;
    this->values.resize(matA.nrow);
    for (int i = 0; i < matA.nrow; i++) {
        this->values[i].resize(matA.ncol, 0.);
    }

    for (int i = 0; i < matA.nrow; i++) {
        for (int j = 0; j < matA.ncol; j++) {
            this->values[i][j] = matA.values[i][j];
        }
    }
    return *this;
}