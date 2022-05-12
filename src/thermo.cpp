#include "thermo.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

ThermoCh::ThermoCh(int n, int sig, double red, vector<double> ener) {
    double unitfactor;
    ndat = (int)ener.size();
    nrel = n;
    redI = red;
    sigma = sig;

    energy.clear();
    for (int i = 0; i < ndat; i++) {
        energy.push_back(ener[i]);
    }

    unitfactor = 2625.5;               // Converts Hartree to kJ mol^-1
    kBoltzmann = 0.008314462618153242; // kB in kJ mol^-1 K^-1

    change_energy(unitfactor);
}

void ThermoCh::change_energy(double factor) {
    for (int i = 0; i < ndat; i++) {
        energy[i] *= factor;
    }
}

/****************************************************************
 *       Hindered rotor  Properties
 ***************************************************************/
double ThermoCh::getFunctionQhr(double temp) {
    double q = 0.f;
    for (int j = 0; j < nrel; j++) {
        q += exp(-energy[j] / (kBoltzmann * temp));
    }

    return q / sigma;
}

double ThermoCh::getEntropy_hr(double temp) {
    double s;
    double num = 0.f, den = 0.f;
    double q;
    for (int j = 0; j < nrel; j++) {
        num += energy[j] * exp(-energy[j] / (kBoltzmann * temp));
        den += exp(-energy[j] / (kBoltzmann * temp));
    }
    q = den / sigma;
    s = kBoltzmann * log(q) + (1. / temp) * (num / den);

    return s;
}

double ThermoCh::getEnergy_hr(double temp) {
    double num = 0.f, den = 0.f;
    for (int j = 0; j < nrel; j++) {
        num += energy[j] * exp(-energy[j] / (kBoltzmann * temp));
        den += exp(-energy[j] / (kBoltzmann * temp));
    }
    return (num / den);
}

double ThermoCh::getCapacity_hr(double temp) {
    double d0q = 0.f, d1q = 0.f, d2q = 0.f;
    double valexp, cv;
    for (int j = 0; j < nrel; j++) {
        valexp = exp(-energy[j] / (kBoltzmann * temp));

        d0q += valexp;
        d1q += energy[j] * valexp;
        d2q += pow(energy[j], 2.f) * valexp;
    }

    cv = (d0q * d2q - d1q * d1q) / (d0q * d0q);

    cv /= (kBoltzmann * temp * temp);

    return cv;
}

/****************************************************************
 *       Free rotor  Properties
 ***************************************************************/
double ThermoCh::getFunctionQfr(double temp) {
    double val;
    val = sqrt(redI * temp * 0.129527);

    return val / sigma;
}

double ThermoCh::getEntropy_fr(double temp) {
    double q = sqrt(redI * temp * 0.129527) / sigma;
    return kBoltzmann * (log(q) + 0.5);
}

double ThermoCh::getEnergy_fr(double temp) { return kBoltzmann * temp / 2.; }

double ThermoCh::getCapacity_fr() { return kBoltzmann / 2.; }

void ThermoCh::getPropertiesAtRange(double ti, double tf, int nstep, int typeI,
                                    string init) {
    string ext("chem.dat");
    string name = init + ext;
    ofstream fout;
    double deltaT = (tf - ti) / (double)(nstep - 1);
    double valT;
    fout.open(name.c_str(), fstream::out);

    fout << "#";
    for (int i = 0; i < 145; i++)
        fout << "=";
    fout << endl;
    fout << "# Print the values of thermo chemical properties" << endl;
    fout << "#     Temperature:  [" << ti << " : " << tf << "] K" << endl;
    fout << "#      I_red(2," << typeI << ")   : " << setw(16)
         << setprecision(8) << fixed << redI << "  Da A^2" << endl;
    fout << "#";
    for (int i = 0; i < 145; i++)
        fout << "=";
    fout << endl;
    fout << "#" << setw(13) << "Temp" << setw(16) << "Q_fr" << setw(16)
         << "Q_hr" << setw(16) << "S_fr" << setw(16) << "S_hr" << setw(16)
         << "E_fr" << setw(16) << "E_hr" << setw(16) << "Cv_fr" << setw(16)
         << "Cv_hr" << endl;
    fout << "#";
    for (int i = 0; i < 145; i++)
        fout << "=";
    fout << endl;

    fout << " # Cv are in J mol-1 K-1" << endl;

    for (int i = 0; i < nstep; i++) {
        valT = ti + i * deltaT;
        double q_hr = getFunctionQhr(valT);
        double q_fr = getFunctionQfr(valT);

        double entropy_hr = getEntropy_hr(valT);
        double entropy_fr = getEntropy_fr(valT);

        double energy_hr = getEnergy_hr(valT);
        double energy_fr = getEnergy_fr(valT);

        double capacity_hr = getCapacity_hr(valT);
        double capacity_fr = getCapacity_fr();
        fout << setw(16) << setprecision(8) << fixed << valT << setw(16)
             << setprecision(8) << fixed << q_fr << setw(16) << setprecision(8)
             << fixed << q_hr << setw(16) << setprecision(8) << fixed
             << entropy_fr << setw(16) << setprecision(8) << fixed << entropy_hr
             << setw(16) << setprecision(8) << fixed << energy_fr << setw(16)
             << setprecision(8) << fixed << energy_hr << setw(16)
             << setprecision(8) << fixed << capacity_fr * 1000.0 << setw(16)
             << setprecision(8) << fixed << capacity_hr * 1000.0 << endl;
    }
    fout.close();
}

void ThermoCh::getPropertiesAtT(double maxV, double valT, int typeI) {
    // Converts the Energy in Hartree to kJ mol^-1
    const double Eh2kJmol = 2625.5;
    cout << "================================================================="
         << endl;
    cout << " Print the values of thermo chemical properties" << endl;
    cout << "      Temperature  : " << setw(16) << setprecision(8) << fixed
         << valT << "  K" << endl;
    cout << "          kB T     : " << setw(16) << setprecision(8) << fixed
         << valT * kBoltzmann << "  kJ / mol" << endl;
    cout << " Highest potential : " << setw(16) << setprecision(8) << fixed
         << maxV * Eh2kJmol << "  kJ / mol" << endl;
    cout << "                     " << setw(16) << setprecision(8) << fixed
         << maxV * Eh2kJmol / kBoltzmann << "  K" << endl;
    cout << "  Symmetry number  : " << setw(16) << setprecision(8) << fixed
         << sigma << endl;
    cout << "      I_red(2," << typeI << ")   : " << setw(16) << setprecision(8)
         << fixed << redI << "  Da A^2" << endl;
    cout << "================================================================="
         << endl;

    double q_hr = getFunctionQhr(valT);
    double q_fr = getFunctionQfr(valT);

    double entropy_hr = getEntropy_hr(valT);
    double entropy_fr = getEntropy_fr(valT);

    double energy_hr = getEnergy_hr(valT);
    double energy_fr = getEnergy_fr(valT);

    double capacity_hr = getCapacity_hr(valT);
    double capacity_fr = getCapacity_fr();

    cout << "                  Free Rotor    Hindered Rotor" << endl;
    cout << "================================================================="
         << endl;

    cout << "    Qr :    " << setw(16) << setprecision(8) << fixed << q_fr
         << setw(16) << setprecision(8) << fixed << q_hr << endl;

    cout << "    S  :    " << setw(16) << setprecision(8) << fixed << entropy_fr
         << setw(16) << setprecision(8) << fixed << entropy_hr
         << "   kJ / (mol K )" << endl;
    cout << "    E  :    " << setw(16) << setprecision(8) << fixed << energy_fr
         << setw(16) << setprecision(8) << fixed << energy_hr << "   kJ / mol "
         << endl;
    cout << "    Cv :    " << setw(16) << setprecision(8) << fixed
         << capacity_fr * 1000.0 << setw(16) << setprecision(8) << fixed
         << capacity_hr * 1000.0 << "   J / (mol K )" << endl;
    cout << "================================================================="
         << endl;
}

void ThermoCh::getPropertiesAtT_range(double maxV, double valT, int typeI) {
    // Converts the Energy in Hartree to kJ mol^-1
    const double Eh2kJmol = 2625.5;
    string name("DataChem2.dat");
    ofstream fout;
    fout.open(name.c_str(), fstream::out);

    fout << "#================================================================="
         << endl;
    fout << "# Print the values of thermo chemical properties" << endl;
    fout << "#      Temperature  : " << setw(16) << setprecision(8) << fixed
         << valT << "  K" << endl;
    fout << "#          kB T     : " << setw(16) << setprecision(8) << fixed
         << valT * kBoltzmann << "  kJ / mol" << endl;
    fout << "# Highest potential : " << setw(16) << setprecision(8) << fixed
         << maxV * Eh2kJmol << "  kJ / mol" << endl;
    fout << "#  Symmetry number  : " << setw(16) << setprecision(8) << fixed
         << sigma << endl;
    fout << "#      I_red(2," << typeI << ")   : " << setw(16)
         << setprecision(8) << fixed << redI << "  Da A^2" << endl;
    fout << "#================================================================="
         << endl;

    fout << "#" << setw(15) << "ndat" << setw(16) << "Qtor" << setw(16) << "S"
         << setw(16) << "E" << setw(16) << "Cv" << endl;

    int nmax = nrel;

    for (int idata = 5; idata < nmax; idata += 5) {
        nrel = idata;
        double q_hr = getFunctionQhr(valT);
        double entropy_hr = getEntropy_hr(valT);
        double energy_hr = getEnergy_hr(valT);
        double capacity_hr = getCapacity_hr(valT);

        fout << setw(16) << setprecision(8) << fixed << nrel;

        fout << setw(16) << setprecision(8) << fixed << q_hr;

        fout << setw(16) << setprecision(8) << fixed << entropy_hr;

        fout << setw(16) << setprecision(8) << fixed << energy_hr;

        fout << setw(16) << setprecision(8) << fixed << capacity_hr << endl;
    }

    fout.close();
}