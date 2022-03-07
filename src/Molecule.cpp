#include "Molecule.h"
#include "Rvector.h"
#include "screen.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
using namespace std;

inline bool checkExtention(string const &name, string const &ending) {
    if (ending.size() > name.size())
        return false;
    return std::equal(ending.rbegin(), ending.rend(), name.rbegin());
}

inline bool is_number(const string &s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

Molecule::Molecule() {
    natoms = 0;
    isInAtomicUnits = false;
    isCentered = false;
    totalMass = 0.;
}

Molecule::~Molecule() {
    atoms.clear();
    centreMass.~Rvector();
}

void Molecule::addAtom(Atom a) {
    atoms.push_back(a);
    natoms++;
}

double Molecule::getTotalMass() {
    double mass = 0.;
    for (int i = 0; i < natoms; i++) {
        mass += atoms[i].getMass();
    }

    return mass;
}

Rvector Molecule::getCentreMass() {
    double mass = 0.;
    Rvector centreMass(0., 0., 0.);

    for (int i = 0; i < natoms; i++) {
        centreMass += (atoms[i].coor) * (atoms[i].getMass());
        mass += atoms[i].getMass();
    }

    centreMass /= mass;
    return centreMass;
}

void Molecule::Translate2CentreMass() {
    Rvector cm = getCentreMass();

    if (!isCentered) {
        for (int i = 0; i < natoms; i++) {
            atoms[i].coor -= cm;
        }
        isCentered = true;
    }
}

void Molecule::ConvertAngstrom2Bohr() {
    if (isInAtomicUnits) {
        for (int i = 0; i < natoms; i++) {
            atoms[i].coor *= 0.529177248994098;
        }
        isInAtomicUnits = false;
    }
}

void Molecule::update() {
    cout << " === Update the molecular data ===" << endl;
    ConvertAngstrom2Bohr();
    totalMass = getTotalMass();
    centreMass = getCentreMass();
    Translate2CentreMass();
    cout << " =================================" << endl;
}

void Molecule::status() {
    cout << " Number of total Atoms  : " << natoms << endl;
    cout << " Is it in atomic units? : " << boolalpha << isInAtomicUnits
         << endl;
    cout << " Is it centered?        : " << boolalpha << isCentered << endl;
    cout << " Total mass             : " << totalMass << " Da" << endl;
    cout << " Mass center            : " << centreMass << endl;
    for (int i = 0; i < natoms; i++) {
        cout << " Atom " << setw(4) << i + 1 << ": " << atoms[i] << endl;
    }
}

int Molecule::getTotalAtoms() { return natoms; }

string Molecule::getSymbolAtomi(int i) { return atoms[i].getSymbol(); }

double Molecule::getXcoorAtomi(int i) { return atoms[i].get_x(); }

double Molecule::getYcoorAtomi(int i) { return atoms[i].get_y(); }

double Molecule::getZcoorAtomi(int i) { return atoms[i].get_z(); }

double Molecule::getMassAtomi(int i) { return atoms[i].getMass(); }

double Molecule::getXcoorCM() { return centreMass.get_x(); }

double Molecule::getYcoorCM() { return centreMass.get_y(); }

double Molecule::getZcoorCM() { return centreMass.get_z(); }

Atom Molecule::getAtom(int idx) { return atoms[idx]; }

Rvector Molecule::getCoors(int idx) { return atoms[idx].getCoors(); }

void Molecule::loadMolecule(const string &name) {

    if (checkExtention(name, string(".xyz"))) {
        readXYZfile(name);
        isInAtomicUnits = false;
    }
    if (checkExtention(name, string(".XYZ"))) {
        readXYZfile(name);
        isInAtomicUnits = false;
    }

    if (checkExtention(name, string(".wfn"))) {
        readWFNfile(name);
        isInAtomicUnits = true;
    }
    if (checkExtention(name, string(".WFN"))) {
        readWFNfile(name);
        isInAtomicUnits = true;
    }

    if (checkExtention(name, string(".wfx"))) {
        readWFXfile(name);
        isInAtomicUnits = true;
    }
    if (checkExtention(name, string(".WFX"))) {
        readWFXfile(name);
        isInAtomicUnits = true;
    }

    if ((int)atoms.size() == 0 || natoms == 0) {
        printE();
    }

    ConvertAngstrom2Bohr();
    totalMass = getTotalMass();
    centreMass = getCentreMass();
    Translate2CentreMass();
}

void Molecule::readWFNfile(const string &name) {
    int nnuclei;
    int nprim;
    int nmolorb;
    char tmp[4];
    int itmp;
    double xt, yt, zt, q;

    string line;
    ifstream finp;

    finp.open(name.c_str(), std::ios::in);
    if (!finp.good()) {
        cout << " The file [" << name << "] can't be opened!" << endl;
    }
    finp.seekg(finp.beg);

    getline(finp, line); // the title of file is readed
    getline(finp, line); // read the line with information of nnuclei, prim;

    sscanf(line.c_str(), "GAUSSIAN %d MOL ORBITALS %d PRIMITIVES %d NUCLEI",
           &nmolorb, &nprim, &nnuclei);
    for (int i = 0; i < nnuclei; i++) {

        getline(finp, line); // read the line with information of centres
        sscanf(line.c_str(), "%s %d (CENTRE %d) %lf %lf %lf CHARGE = %lf", tmp,
               &itmp, &itmp, &xt, &yt, &zt, &q);

        Atom atm(tmp, xt, yt, zt);

        addAtom(atm);
    }

    finp.close();
}

void Molecule::readXYZfile(const string &name) {
    int nnuclei;
    char tmp[10];
    string s;
    double xt, yt, zt;

    string line;
    ifstream finp;

    finp.open(name.c_str(), std::ios::in);
    if (!finp.good()) {
        cout << " The file [" << name << "] can't be opened!" << endl;
    }
    finp.seekg(finp.beg);

    getline(finp, line); // read the line with the number of nuclei
    sscanf(line.c_str(), "%d", &nnuclei);
    getline(finp, line); // read the comment line;
    for (int i = 0; i < nnuclei; i++) {
        getline(finp, line); // read the line with information of centres
        sscanf(line.c_str(), "%s %lf %lf %lf", tmp, &xt, &yt, &zt);

        if (!is_number(tmp)) {
            s = string(tmp);
            s.erase(std::remove_if(s.begin(), s.end(),
                                   [](char ch) { return std::isdigit(ch); }),
                    s.end());
            addAtom(Atom(s, xt, yt, zt));
        } else {
            addAtom(Atom(atoi(tmp), xt, yt, zt));
        }
    }

    finp.close();
}

void Molecule::readWFXfile(const string &name) {
    int nnuclei;
    int tmp;
    double xt, yt, zt;

    vector<int> atomicNumbers;
    vector<Rvector> atomicCoordinates;

    string line;
    ifstream finp;

    finp.open(name.c_str(), std::ios::in);
    if (!finp.good()) {
        cout << " The file [" << name << "] can't be opened!" << endl;
    }
    finp.seekg(finp.beg);

    while (!finp.eof()) {
        getline(finp, line); // read the line with the number of nuclei
        if (line.find(string("<Number of Nuclei>")) != std::string::npos) {
            getline(finp, line);
            sscanf(line.c_str(), "%d", &nnuclei);
        }

        if (line.find(string("<Nuclear Cartesian Coordinates>")) !=
            std::string::npos) {
            for (int i = 0; i < nnuclei; i++) {
                getline(finp, line);
                sscanf(line.c_str(), "%lf %lf %lf", &xt, &yt, &zt);
                atomicCoordinates.push_back(Rvector(xt, yt, zt));
            }
        }

        if (line.find(string("<Atomic Numbers>")) != std::string::npos) {
            for (int i = 0; i < nnuclei; i++) {
                getline(finp, line);
                sscanf(line.c_str(), "%d", &tmp);
                atomicNumbers.push_back(tmp);
            }
        }
    }

    for (int i = 0; i < nnuclei; i++) {
        addAtom(Atom(atomicNumbers[i], atomicCoordinates[i]));
    }

    finp.close();
}