#include "HPot.h"
#include "Matrix.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;

HPot::HPot() {
    ntot = 0;
    ndat = 0;
}

void HPot::loadPotential(const string &name) {
    ntot = 0;
    double theta, energy;
    double emin = 1.E2;
    ifstream finp;

    finp.open(name.c_str(), std::ios::in);

    if (!finp.good()) {
        cout << " The file [" << name << "] can't be opened!" << endl;
    }

    finp.seekg(finp.beg);
    while (!finp.eof()) {
        finp >> theta >> energy;
        hval.push_back(theta);
        epot.push_back(energy);
        if (emin > energy)
            emin = energy;
        ntot++;
    }
    ntot--;
    hval.pop_back();
    epot.pop_back();

    finp.close();

    ndat = int(ntot / 2) - 1;

    for (int i = 0; i < ntot; i++) {
        epot[i] = epot[i] - emin;
        hval[i] = (hval[i] * M_PI) / 180.;
    }
}

void HPot::getCoeffs() {
    Matrix matA(ntot, 2 * ndat);

    int nj = 1;
    for (int j = 0; j < 2 * ndat; j++) {
        if (j % 2 == 0) {
            for (int i = 0; i < ntot; i++) {
                matA(i, j) = 1 - cos(nj * hval[i]);
            }
        } else {
            for (int i = 0; i < ntot; i++) {
                matA(i, j) = sin(nj * hval[i]);
            }
        }
        if (j % 2 == 1)
            nj++;
    }

    Matrix matAt = matA.transpose();
    Matrix matQ = matAt * matA;
    Matrix matQi = matQ.inverse();

    Matrix tmp = matQi * matAt;
    Matrix energy(ntot, 1);

    for (int i = 0; i < ntot; i++)
        energy(i, 0) = epot[i];

    Matrix Coeff = tmp * energy;

    for (int i = 0; i < ndat; i++) {
        coeffA.push_back(Coeff(2 * i, 0));
    }
    for (int i = 0; i < ndat; i++) {
        coeffB.push_back(Coeff(2 * i + 1, 0));
    }
}

void HPot::checkCoeffs(){
    cout << " Check the conditions for Coefficients" << endl;
    double cond1 = 0., cond2 = 0.;
    for (int k = 0; k < ndat; k++) {
        cond1 += coeffB[k] * (k);
        cond2 += coeffA[k] * (k);
    }

    cout << " First  condition: " << setw(12)
         << setprecision(6) << fixed << cond1 << endl;
    cout << " Second condition: " << setw(12)
         << setprecision(6) << fixed << cond2 << endl;

}

void HPot::printCoeffA() {
    cout << " COEFFICIENTS  A " << endl;
    for (int i = 0; i < int(coeffA.size()); i++)
        cout << setw(4) << i << setw(12) << setprecision(6) << fixed
             << coeffA[i] << endl;
}

void HPot::printCoeffB() {
    cout << " COEFFICIENTS  B " << endl;
    for (int i = 0; i < int(coeffB.size()); i++)
        cout << setw(4) << i << setw(12) << setprecision(6) << fixed
             << coeffB[i] << endl;
}

vector<double> &HPot::getCoeffA() { return coeffA; }
vector<double> &HPot::getCoeffB() { return coeffB; }

double HPot::potentialV(double theta) {
    double v = 0.;
    for (int k = 0; k < ndat; k++) {
        v += (coeffA[k] * (1. - cos(theta * (k + 1))) +
              coeffB[k] * sin(theta * (k + 1)));
    }

    return v;
}

double HPot::der1V(double theta) {
    double v = 0.;
    for (int k = 0; k < ndat; k++) {
        v += (coeffA[k] * (k + 1) * sin(theta * (k + 1)) +
              coeffB[k] * (k + 1) * cos(theta * (k + 1)));
    }

    return v;
}

double HPot::der2V(double theta) {
    double v = 0.;
    for (int k = 0; k < ndat; k++) {
        double k2 = (k + 1) * (k + 1);
        v += (coeffA[k] * k2 * cos(theta * (k + 1)) -
              coeffB[k] * k2 * sin(theta * (k + 1)));
    }

    return v;
}

void HPot::printVder(int n, const string &init) {
    string ext("der.dat");
    string name = init + ext;

    double delta, theta, pot, der1, der2;

    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    delta = 2. * M_PI / double(n - 1);

    fout << "# File with the Potential, first and second derivative" << endl;
    fout << "#  Theta(degree)   V(theta)   V'(theta)   V''(theta)" << endl;

    for (int i = 0; i < n; i++) {
        theta = i * delta;
        pot = potentialV(theta);
        der1 = der1V(theta);
        der2 = der2V(theta);
        theta *= (180. / M_PI);
        fout << setw(12) << setprecision(8) << fixed << theta << setw(12)
             << setprecision(8) << fixed << pot;
        fout << setw(12) << setprecision(8) << fixed << der1;
        fout << setw(12) << setprecision(8) << fixed << der2 << endl;
    }
    fout.close();
}

void HPot::printV(int n, const string &init) {
    string ext("pot.dat");
    string name = init + ext;

    double delta, theta, pot;

    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    delta = 2. * M_PI / double(n - 1);

    for (int i = 0; i < n; i++) {
        theta = i * delta;
        pot = potentialV(theta);
        theta *= (180. / M_PI);
        fout << setw(12) << setprecision(8) << fixed << theta << setw(12)
             << setprecision(8) << fixed << pot << endl;
    }
    fout.close();
}

double HPot::getMaxV() {
    int n = 1000;
    double delta, theta, pot, max;
    max = -(double)n;

    delta = 2. * M_PI / double(n - 1);
    for (int i = 0; i < n; i++) {
        theta = i * delta;
        pot = potentialV(theta);
        if (pot > max) {
            max = pot;
        }
    }

    return max;
}
