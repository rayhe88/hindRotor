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
    // This is the original value that I used... I don't remember how I get the
    // number...
    //// double factor = 0.000153618; // h^2/( 4 Pi^2 Ired)
    // New calculations take in mind the value of reduceM = reduce Inertia
    // Moment (Da A^2) Inside the Schrodinger equation I need to use atomi
    // units. 1 Da = 1822.888 486 209 me 1 A = 0.529 177 249 a0
    double factor = 1. / 510.460839399;
    // The value 510.4608393999 is the factor to convert (Da A^2) -> (me a0^2)
    // reduceM needs to enter in Da A^2 units.
    kinetic.clear();

    for (int l = 0; l <= nk; l++) {
        kinetic.push_back((factor / (2. * reduceM)) * pow(l * deltak, 2));
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
    Matrix evecT = evec.transpose();
    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(8) << fixed << scientific << grid[i];
        for (int j = 0; j < nx; j++) {
            wf = evecT(i, j);
            fout << setw(16) << setprecision(8) << fixed << scientific
                 << wf + (eval[j] * 2625.5);
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
    Matrix evecT = evec.transpose();
    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(8) << fixed << scientific << grid[i];
        for (int j = 0; j < nx; j++) {
            wf = evecT(i, j);
            fout << setw(16) << setprecision(8) << fixed << scientific
                 << (wf * wf) + (eval[j] * 2626.6);
        }
        fout << endl;
    }
    fout.close();
}

void Hamiltonian::printEnergyPlot(string init) {
    string ext("EnergyPlot.dat");
    string name = init + ext;
    ofstream fout;
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
             << eval[j] * 2625.5;
    }
    fout << endl;

    fout.close();
}

void Hamiltonian::printEnergy(string init) {
    string ext("Ene.dat");
    string name = init + ext;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    fout << setw(10) << "Eigenvalue" << setw(20) << "Energy (kJ/mol) " << endl;
    for (int j = 0; j < nx; j++) {
        fout << setw(10) << fixed << j << setw(20) << setprecision(8) << fixed
             << eval[j] * 2625.5 << endl;
    }

    fout.close();
}

void Hamiltonian::printDistribution(string init, double temp) {
    string ext("DistEneT.dat");
    string name = init + ext;
    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    double BoltzmannK = 0.008314462618153242;
    double val, sum;

    sum = 0;
    for (int i = 0; i < nx; i++) {
        val = exp(-(eval[i] * 2625.5) / (BoltzmannK * temp));
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
    double BoltzmannK = 0.008314462618153242;
    double val, sum;
    double temp;

    for (int i = 0; i < nx; i++) {
        fout << setw(16) << setprecision(6) << fixed << scientific
             << (double)i / (nx - 1);
        for (int j = 0; j < nT; j++) {
            temp = Ti + j * dT;
            sum = 0;
            for (int k = 0; k < i; k++) {
                val = -eval[k] * 2625.5;
                val /= (BoltzmannK * temp);
                val = exp(val);
                sum += val;
            }

            fout << setw(16) << setprecision(6) << fixed << scientific << sum;
        }
        fout << endl;
    }
    fout.close();
}