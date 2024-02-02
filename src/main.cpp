#include "FourierH.h"
#include "HPot.h"
#include "Inertia.h"
#include "Molecule.h"
#include "runCommands.h"
#include "thermo.h"

#include "Matrix.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

void testEigen(){
    Matrix testA;
    Matrix evec;
    testA.Resize(5,5);
    testA(0,0) =  1.96; testA(0,1) = -6.49; testA(0,2) = -0.47; testA(0,3) = -7.20; testA(0,4) = -0.65;
    testA(1,0) = -6.49; testA(1,1) =  3.80; testA(1,2) = -6.39; testA(1,3) =  1.50; testA(1,4) = -6.34;
    testA(2,0) = -0.47; testA(2,1) = -6.39; testA(2,2) =  4.17; testA(2,3) = -1.51; testA(2,4) =  2.67;
    testA(3,0) = -7.20; testA(3,1) =  1.50; testA(3,2) = -1.51; testA(3,3) =  5.70; testA(3,4) =  1.80;
    testA(4,0) = -0.65; testA(4,1) = -6.34; testA(4,2) =  2.67; testA(4,3) =  1.80; testA(4,4) = -7.10;

    double *val = new double[5];
    evec = testA.eigenSystem(val);

    cout << "Eigenvalues" << endl;
    for(int i=0; i < 5; i++)
        cout << setw(6) << setprecision(2) << fixed << val[i] << endl;

    evec.print("Eigenvectors");

    delete [] val;

    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {

    //testEigen();

    runCommands run(argc, argv);

    HPot vpotentialT;
    Molecule mol;
    Molecule top1, top2;
    double redI, vmax;

    mol.loadMolecule(run.geometry_name);

    if (run.assignMoment == false) {

        for (int i = 0; i < (int)run.top1.size(); i++)
            top1.addAtom(mol.getAtom(run.top1[i]));
        for (int i = 0; i < (int)run.top2.size(); i++)
            top2.addAtom(mol.getAtom(run.top2[i]));

        Inertia moments(top1, top2, run.atom1Bond, run.atom2Bond);
        // the Inertia moments are in Da A^2 units
        redI = moments.getIred(run.typeI);
    } else {
        redI = run.moment;
    }

    /*double redI21 = 0.671762;
    double redI22 = 0.674461;
    double redI23 = 0.626877;

    cout << moments.getRed21() << "    " << redI21 << endl;
    cout << moments.getRed22() << "    " << redI22 << endl;
    cout << moments.getRed23() << "    " << redI23 << endl;*/

    vpotentialT.loadPotential(run.potential_name);

    vpotentialT.getCoeffs();
    vpotentialT.checkCoeffs();

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

    //    chemistry.getPropertiesAtT_range(vmax, run.temp, run.typeI);

    if (run.singleT)
        chemistry.getPropertiesAtT(vmax, run.temp, run.typeI);

    if (run.rangeT)
        chemistry.getPropertiesAtRange(run.tempi, run.tempf, run.ntemp,
                                       run.typeI, run.output_name);

    return (EXIT_SUCCESS);
}
