#ifndef _FOURIER_H_H_
#define _FOURIER_H_H_
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;
#include <vector>
using std::vector;
#include "Matrix.h"

class Hamiltonian {
  private:
    int nx;
    int nk;
    int ndat = 0;
    double deltax, deltak;
    double reduceM;
    double xi = 0.;
    double xf = 2 * M_PI;
    vector<double> coefA;
    vector<double> coefB;
    vector<double> eval;
    vector<double> grid;
    vector<double> vpot;
    vector<double> kinetic;
    Matrix hmat;
    Matrix evec;

  public:
    Hamiltonian(int nx, double red, vector<double> coefA, vector<double> coefB);

    void setCoefA(vector<double>);
    void setCoefB(vector<double>);
    void setPotential();
    void setKinetic();
    void printPot(string name);
    void printWaveFunction(string name);
    void printDensity(string name);
    void printEnergy(string name);
    void printEnergyPlot(string name);
    void normalizeWF(double);
    void printDistribution(string name, double, double, int);
    void printDistribution(string, double);

    vector<double> &getEigenVal();

    Matrix SolveFGH();

    double Potential(double theta);
};

#endif
