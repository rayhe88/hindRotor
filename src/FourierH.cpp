#include "FourierH.h"
#include "Matrix.h"
#include "screen.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;

Hamiltonian::Hamiltonian(int n, double red, vector<double> cA,
                         vector<double> cB) {
    nx = n;
    if (n % 2 == 0) {
        nx++;
        cout << "The value of n must be Odd, n is incremented in 1" << endl;
    }
    nk = (int)(nx - 1) / 2;
    reduceM = red;
    deltax = (xf - xi) / nx;
    deltak = (2. * M_PI) / (nx * deltax);

    setCoefA(cA);
    setCoefB(cB);

    setKinetic();
    setPotential();

    printI();
    cout << " Hamilton size : " << nx << endl;
    printI();
    cout << "  Kinetic size : " << nk << endl;
    printI();
    cout << "   delta   x   : " << deltax << endl;
    printI();
    cout << "   delta   k   : " << deltak << endl;
    printI();
    cout << "   Data in pot : " << ndat << endl;
}

void Hamiltonian::setCoefA(vector<double> cA) {
    coefA.clear();
    for (int i = 0; i < (int)cA.size(); i++) {
        coefA.push_back(cA[i]);
    }
    if (ndat == 0) {
        ndat = (int)cA.size();
    }
}

void Hamiltonian::setCoefB(vector<double> cB) {
    coefB.clear();
    for (int i = 0; i < (int)cB.size(); i++) {
        coefB.push_back(cB[i]);
    }
    if (ndat == 0) {
        ndat = (int)cB.size();
    }
}

double Hamiltonian::Potential(double theta) {

    double vt = (double)0.f;

    for (int i = 0; i < ndat; i++) {
        vt += (coefA[i] * (1. - cos(theta * (i + 1))) +
               coefB[i] * sin(theta * (i + 1)));
    }

    return vt;
}

void Hamiltonian::setKinetic() {
    // This is the original value, the very original
    double factor = 6509.657026642;
    // To obtaing this value, I used the next transformation:
    // I(Da A^2) (1822.889888817 me / 1 Da) (1.8897259886 a0 / 1 A)^2
    // I(Da A^2) (6509.657026642 me a0^2 / Da A^2) = I (me a0^2)
    // The value of reduceM needs to be in Da A^2.

    kinetic.clear();

    for (int l = 0; l <= nk; l++) {
        kinetic.push_back((1. / (2. * factor * reduceM)) *
                          pow(l * deltak, 2.0));
    }
}

void Hamiltonian::setPotential() {
    vpot.clear();
    grid.clear();
    double theta;
    double vtheta;

    for (int i = 0; i < nx; i++) {
        theta = xi + i * deltax;
        vtheta = Potential(theta);

        grid.push_back(theta);
        vpot.push_back(vtheta);
    }
}

Matrix Hamiltonian::SolveFGH() {
    int i, j, l;
    double tmp;
    double hij;

    hmat.Resize(nx, nx);

    for (i = 0; i < nx; i++) {
        for (j = i; j < nx; j++) {
            hij = (double)0.;

            tmp = 2. * M_PI * (i - j) / (double)nx;

            for (l = 1; l <= nk; l++) {
                hij += cos(l * tmp) * kinetic[l];
            }

            hij *= (2.f / (double)nx);
            if (i == j)
                hij += vpot[i];
            hmat(i, j) = hij;
            hmat(j, i) = hij;
        }
    }

    double *val = new double[nx];

    evec = hmat.eigenSystem(val);

    normalizeWF(1.);

    eval.clear();

    for (i = 0; i < nx; i++)
        eval.push_back(val[i]);

    delete[] val;

    return evec;
}

void Hamiltonian::normalizeWF(double factor) {
    double cnorm;
    double psi;

    for (int i = 0; i < nx; i++) {
        double sumPsi2 = 0.;
        for (int j = 0; j < nx; j++) {
            psi = evec(i, j);
            sumPsi2 += psi * psi;
        }
        cnorm = sqrt(factor / (deltax * sumPsi2));
        for (int j = 0; j < nx; j++) {
            psi = evec(i, j);
            psi *= cnorm;
            evec(i, j) = psi;
        }
    }
}

vector<double> &Hamiltonian::getEigenVal() { return eval; }

void Hamiltonian::printPot(string init) {
    string ext("Pot.dat");
    string name = init + ext;

    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(8) << fixed << scientific << grid[i]
             << setw(16) << setprecision(8) << fixed << scientific << vpot[i]
             << endl;
    }
    fout.close();
}

void Hamiltonian::printWaveFunction(string init) {
    string ext("Wf.dat");
    string name = init + ext;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);
    double wf;
    // Converts the Energy in Hartree to kJ mol^-1
    static double Eh2kJmol = 2625.5;
    Matrix evecT = evec.transpose();
    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(8) << fixed << scientific << grid[i];
        for (int j = 0; j < nx; j++) {
            wf = evecT(i, j);
            fout << setw(16) << setprecision(8) << fixed << scientific
                 << wf + (eval[j] * Eh2kJmol);
        }
        fout << endl;
    }
    fout.close();
}

void Hamiltonian::printDensity(string init) {
    string ext("Rho.dat");
    string name = init + ext;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);
    double wf;
    // Converts the Energy in Hartree to kJ mol^-1
    static double Eh2kJmol = 2625.5;
    Matrix evecT = evec.transpose();
    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(8) << fixed << scientific << grid[i];
        for (int j = 0; j < nx; j++) {
            wf = evecT(i, j);
            fout << setw(16) << setprecision(8) << fixed << scientific
                 << (wf * wf) + (eval[j] * Eh2kJmol);
        }
        fout << endl;
    }
    fout.close();
}

void Hamiltonian::printEnergyPlot(string init) {
    string ext("EnergyPlot.dat");
    string name = init + ext;
    ofstream fout;
    // Converts the Energy in Hartree to kJ mol^-1
    static double Eh2kJmol = 2625.5;
    fout.open(name.c_str(), fstream::out);

    fout << setw(16) << setprecision(8) << fixed << scientific << grid[0];
    for (int j = 0; j < nx; j++) {
        fout << setw(16) << setprecision(8) << fixed << scientific
             << eval[j] * 2625.5;
    }
    fout << endl;
    fout << setw(16) << setprecision(8) << fixed << scientific << grid[nx - 1];
    for (int j = 0; j < nx; j++) {
        fout << setw(16) << setprecision(8) << fixed << scientific
             << eval[j] * Eh2kJmol;
    }
    fout << endl;

    fout.close();
}

void Hamiltonian::printEnergy(string init) {
    string ext("Ene.dat");
    string name = init + ext;
    ofstream fout;
    // Converts the Energy in Hartree to kJ mol^-1
    static double Eh2kJmol = 2625.5;
    fout.open(name.c_str(), fstream::out);

    fout << setw(10) << "Eigenvalue" << setw(20) << "Energy (kJ/mol) " << endl;
    for (int j = 0; j < nx; j++) {
        fout << setw(10) << fixed << j << setw(20) << setprecision(8) << fixed
             << eval[j] * Eh2kJmol << endl;
    }

    fout.close();
}

void Hamiltonian::printDistribution(string init, double temp) {
    string ext("DistEneT.dat");
    string name = init + ext;
    // Converts the Energy in Hartree to kJ mol^-1
    const double Eh2kJmol = 2625.5;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    double BoltzmannK = 0.008314462618153242; // kB in kJ mol^-1 K^-1
    double val, sum;

    sum = 0;
    for (int i = 0; i < nx; i++) {
        val = exp(-(eval[i] * Eh2kJmol) / (BoltzmannK * temp));
        sum += val;
        fout << setw(16) << setprecision(6) << fixed << scientific
             << (double)i / (nx - 1) << setw(16) << setprecision(6) << fixed
             << scientific << val << setw(16) << setprecision(6) << fixed
             << scientific << sum << endl;
    }
    fout.close();
}

void Hamiltonian::printDistribution(string init, double Ti, double Tf, int nT) {
    string ext("DistEne.dat");
    string name = init + ext;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    double dT = (Tf - Ti) / (double)nT - 1;
    const double BoltzmannK = 0.008314462618153242; // kB in kJ mol^-1 K^-1
    // Converts the Energy in Hartree to kJ mol^-1
    const double Eh2kJmol = 2625.5;

    double val, sum;
    double temp;

    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(6) << fixed << scientific
             << (double)i / (nx - 1);
        for (int j = 0; j < nT; j++) {
            temp = Ti + j * dT;
            sum = 0;
            for (int k = 0; k < i; k++) {
                val = -eval[k] * Eh2kJmol / (BoltzmannK * temp);
                val = exp(val);
                sum += val;
            }

            fout << setw(16) << setprecision(6) << fixed << scientific << sum;
        }
        fout << endl;
    }
    fout.close();
}
