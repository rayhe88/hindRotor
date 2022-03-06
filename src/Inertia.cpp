#include "Inertia.h"
#include "Molecule.h"
#include "Rvector.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double Inertia::getIred(int typeI) {
    double redI;
    switch (typeI) {
    case 1:
        redI = getRed21();
        break;
    case 2:
        redI = getRed22();
        break;
    case 3:
        redI = getRed23();
        break;
    default:
        redI = getRed21();
    }
    return redI;
}

Inertia::Inertia(Molecule top1, Molecule top2, int at1, int at2) {
    setRed21(top1, top2, at1, at2);
    setRed22(top1, top2, at1, at2);
    setRed23(top1, top2, at1, at2);
}

Rvector Inertia::Cross(Rvector u, Rvector v) {

    double ti, tj, tk;

    ti = u[1] * v[2] - u[2] * v[1];
    tj = u[2] * v[0] - u[0] * v[2];
    tk = u[0] * v[1] - u[1] * v[0];

    return Rvector(ti, tj, tk);
}

double Inertia::DtoAxis(Rvector p, Rvector q, Rvector r) {
    Rvector v1;
    Rvector v2;
    Rvector v3;

    v1 = p - q;
    v2 = q - r;

    v3 = Cross(v1, v2);
    v3 /= v2.norm();

    return v3.norm();
}

void Inertia::setRed21(Molecule top1, Molecule top2, int at1, int at2) {
    double reduced;
    // Internal moment of inertia I(2,1)
    // The moment of inertia of the rotating group is computed about the axis
    // containing the twisting bond.
    double it1 = (double)0.;
    double it2 = (double)0.;

    for (int i = 0; i < top1.getTotalAtoms(); i++) {
        it1 += (top1.getMassAtomi(i) *
                pow(DtoAxis(top1.getCoors(i), top1.getCoors(at1),
                            top2.getCoors(at2)),
                    2.));
    }

    for (int i = 0; i < top2.getTotalAtoms(); i++) {
        it2 += (top2.getMassAtomi(i) *
                pow(DtoAxis(top2.getCoors(i), top1.getCoors(at1),
                            top2.getCoors(at2)),
                    2.));
    }

    reduced = (it1 * it2) / (it1 + it2);

    momentI21 = reduced;
}

void Inertia::setRed22(Molecule top1, Molecule top2, int at1, int at2) {
    // Internal moment of inertia I(2,2)
    // The moment of inertia of the rotating group is computed about the axis
    // parallel to the bound but passing through the center of mass of rotating
    // group.
    Rvector vecu = top1.getCoors(at1) - top2.getCoors(at2);
    Rvector vecp;
    vecu.normalize();
    vecp = top1.getCentreMass() + vecu * 0.2;

    double reduced;
    double it1 = 0., it2 = 0.;

    for (int i = 0; i < top1.getTotalAtoms(); i++) {
        it1 += (top1.getMassAtomi(i) *
                pow(DtoAxis(top1.getCoors(i), top1.getCentreMass(), vecp), 2.));
    }

    for (int i = 0; i < top2.getTotalAtoms(); i++) {
        it2 += (top2.getMassAtomi(i) *
                pow(DtoAxis(top2.getCoors(i), top1.getCentreMass(), vecp), 2.));
    }

    reduced = (it1 * it2) / (it1 + it2);

    momentI22 = reduced;
}

void Inertia::setRed23(Molecule top1, Molecule top2, int at1, int at2) {
    // the moment of inertia of the rotating group is computed about the axis
    // passing through the centers-of-mass of both the rotating group and the
    // remainder of the molecule

    double reduced;
    double it1 = 0., it2 = 0.;

    for (int i = 0; i < top1.getTotalAtoms(); i++) {
        it1 += (top1.getMassAtomi(i) *
                pow(DtoAxis(top1.getCoors(i), top1.getCentreMass(),
                            top2.getCentreMass()),
                    2.));
    }

    for (int i = 0; i < top2.getTotalAtoms(); i++) {
        it2 += (top2.getMassAtomi(i) *
                pow(DtoAxis(top2.getCoors(i), top1.getCentreMass(),
                            top2.getCentreMass()),
                    2.));
    }

    reduced = (it1 * it2) / (it1 + it2);

    momentI23 = reduced;
}

double Inertia::getRed21() { return momentI21; }

double Inertia::getRed22() { return momentI22; }

double Inertia::getRed23() { return momentI23; }
