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

    unitfactor = 2625.5;               // kJ / mol
    kBoltzmann = 0.008314462618153242; // kB in kJ / (mol K)

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
    double q = 0;
    for (int i = 0; i < nrel; i++) {
        q += exp(-energy[i] / (kBoltzmann * temp));
    }

    return q / sigma;
}

double ThermoCh::getEntropy_hr(double q, double temp) {
    double s;
    double num = 0., den = 0.;
    for (int j = 0; j < nrel; j++) {
        num += energy[j] * exp(-energy[j] / (kBoltzmann * temp));
        den += exp(-energy[j] / (kBoltzmann * temp));
    }
    s = kBoltzmann * log(q) + (num / (temp * den));

    return s;
}

double ThermoCh::getEnergy_hr(double temp) {
    double num = 0., den = 0.;
    for (int j = 0; j < nrel; j++) {
        num += energy[j] * exp(-energy[j] / (kBoltzmann * temp));
        den += exp(-energy[j] / (kBoltzmann * temp));
    }
    return (num / den);
}

double ThermoCh::getCapacity_hr(double temp) {
    double t0 = 0., t1 = 0., t2 = 0., t3 = 0.;
    double valexp;
    for (int j = 0; j < nrel; j++) {
        valexp = exp(-energy[j] / (kBoltzmann * temp));

        t0 += valexp;
        t1 += energy[j] * energy[j] * valexp / (kBoltzmann * temp * temp);
        t2 += energy[j] * valexp / (kBoltzmann * temp * temp);
        t3 += energy[j] * valexp;
    }

    return ((t0 * t1) - (t2 * t3)) / (t0 * t0);
}

/****************************************************************
 *       Free rotor  Properties
 ***************************************************************/
double ThermoCh::getFunctionQfr(double temp) {
    double val;
    val = sqrt(redI * temp * 0.129527);

    return val / sigma;
}

double ThermoCh::getEntropy_fr(double q) { return kBoltzmann * (log(q) + 0.5); }

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
    cout << "      I_red(2," << typeI << ")   : " << setw(16) << setprecision(8)
         << fixed << redI << "  Da A^2" << endl;
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

    for (int i = 0; i < nstep; i++) {
        valT = ti + i * deltaT;
        double q_hr = getFunctionQhr(valT);
        double q_fr = getFunctionQfr(valT);

        double entropy_hr = getEntropy_hr(q_hr, valT);
        double entropy_fr = getEntropy_fr(q_fr);

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
             << setprecision(8) << fixed << capacity_fr << setw(16)
             << setprecision(8) << fixed << capacity_hr << endl;
    }
    fout.close();
}

void ThermoCh::getPropertiesAtT(double maxV, double valT, int typeI) {
    double unitfactor = 2625.5; // Transform Hartree to kJ / mol
    double q_hr = getFunctionQhr(valT);
    double q_fr = getFunctionQfr(valT);

    double entropy_hr = getEntropy_hr(q_hr, valT);
    double entropy_fr = getEntropy_fr(q_fr);

    double energy_hr = getEnergy_hr(valT);
    double energy_fr = getEnergy_fr(valT);

    double capacity_hr = getCapacity_hr(valT);
    double capacity_fr = getCapacity_fr();
    cout << "================================================================="
         << endl;
    cout << " Print the values of thermo chemical properties" << endl;
    cout << "      Temperature  : " << setw(16) << setprecision(8) << fixed
         << valT << "  K" << endl;
    cout << "          kB T     : " << setw(16) << setprecision(8) << fixed
         << valT * kBoltzmann << "  kJ / mol" << endl;
    cout << " Highest potential : " << setw(16) << setprecision(8) << fixed
         << maxV * unitfactor << "  kJ / mol" << endl;
    cout << "      I_red(2," << typeI << ")   : " << setw(16) << setprecision(8)
         << fixed << redI << "  Da A^2" << endl;
    cout << "================================================================="
         << endl;
    cout << "                  Free Rotor    Hindered Rotor" << endl;
    cout << "================================================================="
         << endl;

    cout << "    Qr :    " << setw(16) << setprecision(8) << fixed << q_fr
         << setw(16) << setprecision(8) << fixed << q_hr << endl;

    cout << "    S  :    " << setw(16) << setprecision(8) << fixed << entropy_fr
         << setw(16) << setprecision(8) << fixed << entropy_hr << endl;
    cout << "    E  :    " << setw(16) << setprecision(8) << fixed << energy_fr
         << setw(16) << setprecision(8) << fixed << energy_hr << endl;
    cout << "    Cv :    " << setw(16) << setprecision(8) << fixed
         << capacity_fr << setw(16) << setprecision(8) << fixed << capacity_hr
         << endl;
    cout << "================================================================="
         << endl;
}
