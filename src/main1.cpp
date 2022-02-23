#include "Matrix.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

int main() {
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

    prod.print(" Producto suma . mat I");

    Matrix matR(4, 5);
    matR.Random();
    matR.print(" Random ");

    Matrix matRt;
    matRt = matR.transpose();
    matRt.print(" Random transpose");

    Matrix matQ(3, 3);
    matQ.print(" MatQ ");

    matQ(0, 0) = 1.;
    matQ(0, 1) = 2.;
    matQ(0, 2) = 3.;
    matQ(1, 0) = 4.;
    matQ(1, 1) = 5.;
    matQ(1, 2) = 6.;
    matQ(2, 0) = 7.;
    matQ(2, 1) = 8.;
    matQ(2, 2) = 19.;
    matQ.print(" MatQ ");

    Matrix matQi = matQ.inverse();
    matQi.print(" MatQi ");

    Matrix p = matQ * matQi;
    Matrix q = matQi * matQ;

    p.print(" Matrix p = Q . Qi");
    q.print(" Matrix q = Qi . Q");

    Matrix eigA(5, 5);

    eigA(0, 0) = 1.96;
    eigA(0, 1) = -6.49;
    eigA(0, 2) = -0.47;
    eigA(0, 3) = -7.20;
    eigA(0, 4) = -0.65;

    eigA(1, 0) = -6.49;
    eigA(1, 1) = 3.80;
    eigA(1, 2) = -6.39;
    eigA(1, 3) = 1.50;
    eigA(1, 4) = -6.34;

    eigA(2, 0) = -0.47;
    eigA(2, 1) = -6.39;
    eigA(2, 2) = 4.17;
    eigA(2, 3) = -1.51;
    eigA(2, 4) = 2.67;

    eigA(3, 0) = -7.20;
    eigA(3, 1) = 1.50;
    eigA(3, 2) = -1.51;
    eigA(3, 3) = 5.70;
    eigA(3, 4) = 1.80;

    eigA(4, 0) = -0.65;
    eigA(4, 1) = -6.34;
    eigA(4, 2) = 2.67;
    eigA(4, 3) = 1.80;
    eigA(4, 4) = -7.10;

    eigA.print(" EJEMPLO DIAGONALIZACION");
    double *eval = new double[eigA.get_nrow()];

    Matrix evec = eigA.eigenSystem(eval);

    evec.print(" Eigenvectors ");

    for (int i = 0; i < 5; i++) {
        cout << setw(7) << i << setw(12) << setprecision(6) << fixed << eval[i]
             << endl;
    }
    delete[] eval;

    return (EXIT_SUCCESS);
}
