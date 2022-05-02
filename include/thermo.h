#ifndef _THERMO_CH_H_
#define _THERMO_CH_H_
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;
#include <vector>
using std::vector;
#include "Matrix.h"

class ThermoCh {
  private:
    int ndat = 0;
    int nrel = 0;
    int sigma = 1;
    double redI;
    double kBoltzmann;
    vector<double> energy;

  public:
    ThermoCh(int, int, double, vector<double>);
    void change_energy(double);

    void getPropertiesAtT(double, double, int);

    double getFunctionQhr(double temp);

    double getEntropy_hr(double q, double temp);

    double getEnergy_hr(double temp);

    double getCapacity_hr(double temp);

    double getFunctionQfr(double temp);

    double getEntropy_fr(double q);

    double getEnergy_fr(double temp);

    double getCapacity_fr();

    void getPropertiesAtRange(double, double, int, int, string);
};

#endif
