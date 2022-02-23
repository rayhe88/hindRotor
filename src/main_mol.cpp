#include "Atom.h"
#include "Inertia.h"
#include "Matrix.h"
#include "Molecule.h"
#include "Rvector.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

int main() {
    /*
        double x = 1., y = 2., z = 3.;
        double a = 4., b = 5., c = 6.;

        Rvector vecA(x, y, z);
        Rvector vecB(a, b, c);

        Rvector vec0;
        Rvector res;

        res = vecA + vecB;

        cout << " Vector   A   : " << vecA << endl;
        cout << " Vector   B   : " << vecB << endl;
        cout << " Vector A + B : " << res << endl;

        res = vecA - vecB;

        cout << " Vector   A   : " << vecA << endl;
        cout << " Vector   B   : " << vecB << endl;
        cout << " Vector A + B : " << res << endl;

        res = vecA * 2;

        cout << " Vector   A   : " << vecA << endl;
        cout << " Vector  A 2  : " << res << endl;

        res = vecA / 3;

        cout << " Vector   A   : " << vecA << endl;
        cout << " Vector  A/3  : " << res << endl;

        double dprod = vecA.dot(vecB);
        cout << " VectorA . VectorB : " << dprod << endl;
        dprod = vecB.dot(vecA);
        cout << " VectorB . VectorA : " << dprod << endl;

        dprod = vecA.norm();
        cout << " Norma del vector A: " << dprod << endl;

        vecA.normalize();
        cout << " normalizado VectorA : " << vecA << "\n" << endl;

        cout << "=====================================================" << endl;

        Atom atom1(1, 1., 2., 3.);

        cout << "atom1" << atom1 << endl;

        Rvector coord(0., -1., -2.);

        Atom atom2(2, coord);
        cout << "atom2" << atom2 << endl;

        double r[3] = {10., 20., 30.};
        Atom atom3(3, r);
        cout << "atom3" << atom3 << endl;

        cout << "=====================================================" << endl;

        Atom atom4("H", 1., 2., 3.);

        cout << "atom4" << atom4 << endl;

        Atom atom5("He", coord);
        cout << "atom5" << atom5 << endl;

        Atom atom6("Li", r);
        cout << "atom6" << atom6 << endl;

        cout << "=====================================================" << endl;
        string st7 = "Be";
        string st8 = "B";
        string st9 = "C";
        Atom atom7(st7, 1., 2., 3.);

        cout << "atom7" << atom7 << endl;

        Atom atom8(st8, coord);
        cout << "atom8" << atom8 << endl;

        Atom atom9(st9, r);
        cout << "atom9" << atom9 << endl;

        cout << "=====================================================" << endl;
        Rvector dist = atom3.coor - atom2.coor;

        cout << "Distancia entre atom3 y atom2 : " << dist.norm() << " u.a" <<
       endl; for (int idx = 0; idx < 119; idx++) { Rvector r(0, 0, 0); Atom
       atom(idx, r); cout << "Atom " << setw(4) << idx << ": " << atom << endl;
        }

        Molecule mol;
        Rvector r1(0., 0., 0.);
        Rvector r2(2.59808, 1.5, 0.);
        Rvector r3(0., 3., 0.);
        Rvector r4(0.866025, 1.5, 2.44949);

        mol.addAtom(Atom("H", r1));
        mol.addAtom(Atom(1, r2));
        mol.addAtom(Atom(1, r3));
        mol.addAtom(Atom(1, r4));
        mol.status();
        mol.update();
        mol.status();

        ofstream myfile;
        myfile.open("Tetra1.xyz");

        myfile << mol.getTotalAtoms() + 1 << endl;
        myfile << "Testfile" << endl;

        myfile << setw(4) << fixed << "CM" << setw(10) << fixed <<
       setprecision(6)
               << mol.getXcoorCM() << setw(10) << fixed << setprecision(6)
               << mol.getYcoorCM() << setw(10) << fixed << setprecision(6)
               << mol.getZcoorCM() << endl;

        for (int i = 0; i < mol.getTotalAtoms(); i++) {
            myfile << setw(4) << fixed << mol.getSymbolAtomi(i) << setw(10) <<
       fixed
                   << setprecision(6) << mol.getXcoorAtomi(i) << setw(10) <<
       fixed
                   << setprecision(6) << mol.getYcoorAtomi(i) << setw(10) <<
       fixed
                   << setprecision(6) << mol.getZcoorAtomi(i) << endl;
        }

        myfile.close();
    */

    /* cout << "Archivo con ext xyz" << endl;
     mol.loadMolecule(string("hola.xyz"));

     cout << "Archivo con ext XYZ" << endl;
     mol.loadMolecule(string("hola.XYZ"));

     cout << "Archivo con ext wfn" << endl;
     mol.loadMolecule(string("hola.wfn"));

     cout << "Archivo con ext WFN" << endl;
     mol.loadMolecule(string("hola.WFN"));

     cout << "Archivo con ext wfx" << endl;
     mol.loadMolecule(string("hola.wfx"));

     cout << "Archivo con ext WFX" << endl;
     mol.loadMolecule(string("hola.WFX"));*/
    /*
        Molecule mol1;
        mol1.loadMolecule(string("examples/methanol-opt.wfn"));
        mol1.update();
        mol1.status();
        for (int i = 0; i < mol1.getTotalAtoms(); i++) {
            cout << setw(4) << fixed << mol1.getSymbolAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol1.getXcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol1.getYcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol1.getZcoorAtomi(i) << endl;
        }

        Molecule mol2;
        mol2.loadMolecule(string("examples/waterDimer.xyz"));
        for (int i = 0; i < mol2.getTotalAtoms(); i++) {
            cout << setw(4) << fixed << mol2.getSymbolAtomi(i) << setw(10)
   << fixed
                 << setprecision(6) << mol2.getXcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol2.getYcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol2.getZcoorAtomi(i) << endl;
        }

        Molecule mol3;
        mol3.loadMolecule(string("examples/waterDimer.wfx"));
        for (int i = 0; i < mol3.getTotalAtoms(); i++) {
            cout << setw(4) << fixed << mol3.getSymbolAtomi(i) << setw(10)
   << fixed
                 << setprecision(6) << mol3.getXcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol3.getYcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << mol3.getZcoorAtomi(i) << endl;
        }
     */
    /*
        Molecule top1, top2;

        top1.addAtom(mol1.getAtom(0));
        top1.addAtom(mol1.getAtom(1));
        top1.addAtom(mol1.getAtom(2));
        top1.addAtom(mol1.getAtom(3));

        top2.addAtom(mol1.getAtom(4));
        top2.addAtom(mol1.getAtom(5));
        cout << " =========   TOP 1    ===========" << endl;
        for (int i = 0; i < top1.getTotalAtoms(); i++) {
            cout << setw(4) << fixed << top1.getSymbolAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top1.getXcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top1.getYcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top1.getZcoorAtomi(i) << endl;
        }
        cout << " =========   TOP 2    ===========" << endl;
        for (int i = 0; i < top2.getTotalAtoms(); i++) {
            cout << setw(4) << fixed << top2.getSymbolAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top2.getXcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top2.getYcoorAtomi(i) << setw(10) <<
       fixed
                 << setprecision(6) << top2.getZcoorAtomi(i) << endl;
        }
        cout << " ================================" << endl;

        Inertia moments(top1, top2, 0, 0);

        cout << " Moments of Inertia" << endl;
        cout << " I(2,1) : " << setw(10) << setprecision(7) << fixed
             << moments.getRed21() << endl;
        cout << " I(2,2) : " << setw(10) << setprecision(7) << fixed
             << moments.getRed22() << endl;
        cout << " I(2,3) : " << setw(10) << setprecision(7) << fixed
             << moments.getRed23() << endl;
        */
    /*
    Matrix matA(3, 4);
    matA.fill();
    matA.print(" A ");

    Matrix matT = matA.transpose();

    matT.print(" At ");

    matA.fill(3);

    matA.print(" A ");

    // Matrix matB = matA;
    Matrix matB = matA;

    matB.print(" B ");

    Matrix matC;

    matC = matB;
    matC.print(" C ");

    Matrix matI(4);
    matI.Eye();

    matI.print(" Identidad ");

    Matrix suma = matA + matB;

    suma.print(" suma ");

    Matrix prod = suma * matI;

    prod.print(" Producto matA . mat I");*/
    cout << "=================================================================="
         << endl;
    Molecule mol3;
    mol3.loadMolecule(string("examples/waterDimer.wfx"));
    mol3.status();
    cout << "=================================================================="
         << endl;

    Molecule mol2;

    mol2.loadMolecule(string("examples/waterDimer2.xyz"));
    mol2.status();
    cout << "=================================================================="
         << endl;
    return (EXIT_SUCCESS);
}
