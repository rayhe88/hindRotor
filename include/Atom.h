#ifndef _ATOM_H_
#define _ATOM_H_

#include "Rvector.h"
#include <iostream>
#include <string>
using std::string;

using namespace std;

class Atom {
  public:
    Rvector coor;

  private:
    int zatom = -1;
    double mass = 0.1;
    string symb = "00";
    string setSymbol(int);

    int setAtomicNumberfromSymbol(string);
    double setAtomicMass(int);

  public:
    //  Rvector coor;
    Atom(int, double, double, double);
    Atom(int, double *);
    Atom(int, Rvector);
    ~Atom();
    //***********************************************
    Atom(string, double, double, double);
    Atom(string, double *);
    Atom(string, Rvector);
    Atom(const char *, double, double, double);
    Atom(const char *, double *);
    Atom(const char *, Rvector);
    //***********************************************
    Atom &operator=(const Atom &at);
    //***********************************************
    double getMass();
    string getSymbol();
    Rvector getCoors();
    double get_x();
    double get_y();
    double get_z();
    //***********************************************
    friend ostream &operator<<(ostream &o, const Atom &);
};

#endif
