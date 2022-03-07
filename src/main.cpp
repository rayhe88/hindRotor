#include "FourierH.h"
#include "HPot.h"
#include "Inertia.h"
#include "Molecule.h"
#include "runCommands.h"
#include "thermo.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {

    runCommands run(argc, argv);

    HPot vpotentialT;
    Molecule mol;
    Molecule top1, top2;
    double redI, vmax;

    mol.loadMolecule(run.geometry_name);

    for (int i = 0; i < (int)run.top1.size(); i++)
        top1.addAtom(mol.getAtom(run.top1[i]));
    for (int i = 0; i < (int)run.top2.size(); i++)
        top2.addAtom(mol.getAtom(run.top2[i]));

    Inertia moments(top1, top2, run.atom1Bond, run.atom2Bond);
    // the Inertia moments are in Da A^2 units
    redI = moments.getIred(run.typeI);

    /*double redI21 = 0.671762;
    double redI22 = 0.674461;
    double redI23 = 0.626877;

    cout << moments.getRed21() << "    " << redI21 << endl;
    cout << moments.getRed22() << "    " << redI22 << endl;
    cout << moments.getRed23() << "    " << redI23 << endl;*/

    vpotentialT.loadPotential(run.potential_name);

    vpotentialT.getCoeffs();

    vmax = vpotentialT.getMaxV();

    Hamiltonian fourier(run.hsize, redI, vpotentialT.getCoeffA(),
                        vpotentialT.getCoeffB());

    Matrix evec = fourier.SolveFGH();

    vector<double> eval = fourier.getEigenVal();

    fourier.printPot(run.output_name);

    fourier.normalizeWF(1.0);

    fourier.printWaveFunction(run.output_name);
    fourier.printDensity(run.output_name);
    fourier.printEnergy(run.output_name);
    fourier.printEnergyPlot(run.output_name);

    ThermoCh chemistry((int)eval.size(), run.sigma, redI, eval);

    if (run.singleT)
        chemistry.getPropertiesAtT(vmax, run.temp, run.typeI);

    if (run.rangeT)
        chemistry.getPropertiesAtRange(run.tempi, run.tempf, run.ntemp,
                                       run.typeI, run.output_name);

    return (EXIT_SUCCESS);
}
