#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include "Atom.h"
#include <fstream>
#include <iostream>
#include <string>
using std::string;
#include <vector>
using std::vector;

using namespace std;

class Molecule {
  private:
    int natoms;
    bool isInAtomicUnits;
    bool isCentered;
    vector<Atom> atoms;
    double totalMass;
    Rvector centreMass;

    void readWFNfile(const string &);
    void readWFXfile(const string &);
    void readXYZfile(const string &);

  public:
    Molecule();
    ~Molecule();
    void addAtom(Atom);
    double getTotalMass();
    Rvector getCentreMass();
    void Translate2CentreMass();
    void ConvertAngstrom2Bohr();
    void status();
    void update();

    int getTotalAtoms();
    string getSymbolAtomi(int i);
    double getMassAtomi(int i);
    double getXcoorAtomi(int i);
    double getYcoorAtomi(int i);
    double getZcoorAtomi(int i);
    double getXcoorCM();
    double getYcoorCM();
    double getZcoorCM();

    Atom getAtom(int i);
    Rvector getCoors(int i);

    void loadMolecule(const string &);
};
#endif
