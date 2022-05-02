#ifndef _HPOT_H_
#define _HPOT_H_
#include <fstream>
#include <iostream>
using namespace std;
#include <vector>
using std::vector;

class HPot {
  private:
    int ntot;
    int ndat;
    vector<double> epot;
    vector<double> hval;
    vector<double> coeffA;
    vector<double> coeffB;

  public:
    HPot();
    void loadPotential(const string &);
    void getCoeffs();
    void printCoeffA();
    void printCoeffB();
    void printV(int, const string &);
    vector<double> &getCoeffA();
    vector<double> &getCoeffB();
    double potentialV(double t);
    double getMaxV();
};

#endif
