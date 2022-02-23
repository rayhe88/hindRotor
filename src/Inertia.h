#ifndef _INERTIA_H_
#define _INERTIA_H_

#include "Atom.h"
#include "Molecule.h"

#include <iostream>
using namespace std;

class Inertia {
  private:
    double momentI21;
    double momentI22;
    double momentI23;
    Rvector Cross(Rvector, Rvector);
    double DtoAxis(Rvector, Rvector, Rvector);

  public:
    Inertia(Molecule top1, Molecule top2, int at1, int at2);

    void setRed21(Molecule top1, Molecule top2, int at1, int at2);
    void setRed22(Molecule top1, Molecule top2, int at1, int at2);
    void setRed23(Molecule top1, Molecule top2, int at1, int at2);

    double getRed21();
    double getRed22();
    double getRed23();

    double getIred(int);
};

#endif
