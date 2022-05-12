#include "runCommands.h"
#include "screen.h"
#include "version.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <lapacke.h>

using namespace std;

bool runCommands::exists_file(string name) {
    bool exits = false;
    ifstream file;
    file.open(name.c_str(), std::ios::in);
    exits = file.good();
    if (exits == true)
        file.close();
    return !exits;
}

void runCommands::checkNames() {
    string none = string("NULL");
    bool flagError = false;

    if (!geometry_name.compare(none)) {
        printE();
        cout << "Geometry file name is missing" << endl;
        flagError = true;
    } else {
        if (exists_file(geometry_name)) {
            printE();
            cout << " " << geometry_name << " file not found" << endl;
            flagError = true;
        }
    }

    if (!potential_name.compare(none)) {
        printE();
        cout << "Potential file name is missing" << endl;
        flagError = true;
    } else {
        if (exists_file(potential_name)) {
            printE();
            cout << potential_name << " file not found" << endl;
            flagError = true;
        }
    }

    if (!tops_name.compare(none)) {
        printE();
        cout << "Tops file name is missing" << endl;
        flagError = true;
    } else {
        if (exists_file(tops_name)) {
            printE();
            cout << tops_name << " file not found" << endl;
            flagError = true;
        }
    }

    if (!output_name.compare(none)) {
        output_name = string(DEFAULTNAME);
        printW();
        cout << "Output name is missing" << endl;
        cout << "           " << output_name
             << " will be assigned to the output file" << endl;

    } else {
        if (!exists_file(output_name)) {
            printW();
            cout << " A file with the name '" << tops_name
                 << "' wal already fount, this file will be rewritten" << endl;
        }
    }
    if (flagError == true) {
        printE();
        cout << "End of program: Missing input files" << endl;
        exit(EXIT_FAILURE);
    }
}

void runCommands::getTemperature(char **argv) {
    temp = 298.15;
    bool flagT = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'T') {
            flagT = true;
            if (argv[i + 1] != NULL) {
                temp = atof(argv[i + 1]);
                singleT = true;
                if (argv[i + 1][0] == '-') {
                    printW();
                    cout << " Temperature flag was found but no "
                            "assigned a value."
                         << endl;
                    cout << "           It will evaluate to "
                            " T = 298.15 K"
                         << endl;
                    temp = 298.15;
                }
            }
        }
    }
    if (flagT == false) {
        printW();
        cout << "No temperature assigned" << endl;
        cout << "           The temperature will be taken by default T = "
                "298.15 K"
             << endl;
        singleT = true;
    }
    if (temp < 0) {
        printW();
        cout << " The Temperature needs to be > 0 " << endl;
        cout << "           The temperature will be taken by default T = "
                "298.15 K "
             << endl;
        sigma = 1;
    }
}

void runCommands::getRangeTemperature(char **argv) {
    tempf = 0.;
    tempi = 0.;
    ntemp = 0;
    bool flagError = false;
    bool findFlag = false;
    double swap;

    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'r') {
            findFlag = true;
            if (argv[i + 1] != NULL)
                tempi = atof(argv[i + 1]);
            if (argv[i + 2] != NULL)
                tempf = atof(argv[i + 2]);
            if (argv[i + 3] != NULL)
                ntemp = atoi(argv[i + 3]);
        }
    }

    if (findFlag == false) {
        return;
    }

    if (tempf < tempi) {
        printW();
        cout << "The final temperature cannot be less than initial "
                "temperature"
             << endl;
        printW();
        cout << "The temperature values will be swapped" << endl;
        if (tempi > tempf) {
            swap = tempi;
            tempi = tempf;
            tempf = swap;
        }
        cout << "         T initial : " << setw(12) << setprecision(6) << fixed
             << tempi << endl;
        cout << "         T final   : " << setw(12) << setprecision(6) << fixed
             << tempf << endl;
    }
    if (tempi < 0) {
        printE();
        cout << "The initial temperature cannot be less than zero" << endl;
        cout << "         T initial : " << setw(12) << setprecision(6) << fixed
             << tempi << endl;
        flagError = true;
    }
    if (tempf <= 0) {
        printE();
        cout << "The final temperature cannot be less than or equal "
                "zero"
             << endl;
        cout << "         T final : " << setw(12) << setprecision(6) << fixed
             << tempf << endl;
        flagError = true;
    }

    if (ntemp <= 1) {
        printE();
        cout << "The number of evaluations must be greater than 1" << endl;
        cout << "         evaluations : " << setw(12) << setprecision(6)
             << fixed << ntemp << endl;
        flagError = true;
    }

    if (flagError == false) {
        rangeT = true;

        htemp = (tempf - tempi) / (double)ntemp - 1;
    }
}

string runCommands::getName(const char *flag, char **argv) {

    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == flag[0] && argv[i][1] == flag[1]) {
            if (argv[i + 1] != NULL) {
                if (argv[i + 1][0] != '-') {
                    return (string(argv[i + 1]));
                }
            } else {
                return (string("NULL"));
            }
        }
    }
    return (string("NULL"));
}

void runCommands::getVersion(char **argv) {
    bool flagV = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && (argv[i][1] == 'v' || argv[i][1] == 'V')) {
            flagV = true;
        }
    }
    if (flagV) {
        printI();
        cout << "Version: " << PROJECT_VER << endl;
        printI();
        cout << "Compilation Date: " << __DATE__ << "  " << __TIME__ << endl;
        printI();
        cout << "Git SHA1: " << GIT_SHA1 << endl;
        exit(EXIT_SUCCESS);
    }
}

void runCommands::getLapackVersion(char **argv) {
    bool flagV = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'l') {
            flagV = true;
        }
    }
    if (flagV) {
        printI();
        int major, minor, patch;
        ilaver_(&major, &minor, &patch);

        cout << "Lapack version: " << major << "." << minor << "." << patch
             << endl;
        exit(EXIT_SUCCESS);
    }
}

void runCommands::selectInertia(char **argv) {
    typeI = 1;
    bool flag = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'I') {
            if (argv[i + 1] != NULL) {
                typeI = atoi(argv[i + 1]);
                flag = true;
                if (argv[i + 1][0] == '-') {
                    printW();
                    cout << "Inertia moment flag was found but no "
                            "assigned a value."
                         << endl;
                    cout << "          It will evaluate to "
                            " I(2,1)"
                         << endl;
                    typeI = 1;
                    flag = false;
                }
            } else {
                flag = false;
            }
        }
    }
    if (flag == false) {
        printW();
        cout << "No Inertia moment assigned" << endl;
        cout
            << "           The Inertia moment will be taken by default I (2,1) "
            << endl;
        typeI = 1;
    }
    if (typeI < 1 || typeI > 3) {
        printW();
        cout << "The Inertia moment only can take the values of I(2,n) n = 1, "
                "2 or 3"
             << endl;
        cout
            << "           The Inertia moment will be taken by default I (2,1) "
            << endl;
        typeI = 1;
    }
}

void runCommands::assignInertia(char **argv) {
    assignMoment = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'm') {
            if (argv[i + 1] != NULL) {
                moment = atof(argv[i + 1]);
                assignMoment = true;
            }
        }
    }
}

void runCommands::getSigma(char **argv) {
    sigma = 1;
    bool flag = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 's') {
            flag = true;
            if (argv[i + 1] != NULL) {
                sigma = atoi(argv[i + 1]);
                if (argv[i + 1][0] == '-') {
                    printW();
                    cout << "Symmetry number flag was found but no "
                            "assigned a value."
                         << endl;
                    cout << "          It will evaluate to "
                            " sigma = 1"
                         << endl;
                    sigma = 1;
                }
            }
        }
    }
    if (flag == false) {
        printW();
        cout << "No Symmetry number assigned" << endl;
        cout << "           The Symmetry number will be taken by default sigma "
                "= 1"
             << endl;
        sigma = 1;
    }
    if (sigma < 1) {
        printW();
        cout << "The Symmetry number only can take the values sigma >= 0 "
             << endl;
        cout << "           The Symmetry number will be taken by default sigma "
                "= 1 "
             << endl;
        sigma = 1;
    }
}

void runCommands::getSizeH(char **argv) {
    hsize = 501;
    bool flag = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'n') {
            flag = true;
            if (argv[i + 1] != NULL) {
                hsize = atoi(argv[i + 1]);
                if (argv[i + 1][0] == '-') {
                    printW();
                    cout << "Hamiltonian size flag was found but no "
                            "assigned a value."
                         << endl;
                    cout << "          It will evaluate to "
                            " H size = 501"
                         << endl;
                    hsize = 501;
                }
            }
        }
    }
    if (flag == false) {
        printW();
        cout << "No Hamiltonian size assigned" << endl;
        cout << "           The Hamiltonian size will be taken by default H "
                "size "
                "= 501"
             << endl;
        hsize = 501;
    }
    if (hsize < 1) {
        printW();
        cout << "The Hamiltonian size only can take the values H size > 1 "
             << endl;
        cout << "           The Hamiltonian size will be taken by default H "
                "size "
                "= 501 "
             << endl;
        hsize = 501;
    }
    if (hsize % 2 == 0) {
        hsize++;
        printW();
        cout << "The Hamiltonian size needs to be odd" << endl;
        cout << "           The Hamiltonian size will be increase in one H "
                "size = "
             << hsize << endl;
    }
}

void runCommands::helper(char **argv) {
    bool flagH = false;
    for (int i = 0; i < nargc; i++) {
        if (argv[i][0] == '-' && (argv[i][1] == 'h' || argv[i][1] == 'H')) {
            flagH = true;
        }
    }
    if (flagH) {
        cout << "\n  Usage: " << argv[0] << " FLAGS [VALUES]" << endl;
        cout
            << "  The program computes the rotational partition function for a "
            << endl;
        cout << "  molecular system at a temperature T or a range of "
                "temperatures.\n"
             << endl;
        cout << RED << "  Mandatory flags! " << RST << endl;
        cout
            << "  -g [string]   Geometry file. Formats supported xyz, wfn, wfx."
            << endl;
        cout << "  -p [string]   Potential file. File whit the rotational "
                "potential."
             << endl;
        cout << "                  In the format :" << endl;
        cout << "                   __________________________ " << endl;
        cout << "                  |  theta_0   potential_0   |" << endl;
        cout << "                  |  theta_1   potential_1   |" << endl;
        cout << "                  |  theta_2   potential_2   |" << endl;
        cout << "                  |    ...         ...       |" << endl;
        cout << "                  |  theta_m   potential_m   |" << endl;
        cout << "                  |__________________________|\n" << endl;
        cout << "  -t [string]   File with the TOPS information. TOP A is the "
                "group of"
             << endl;
        cout << "                nA atoms, the second line corresponds to the "
                "indexes of each atom."
             << endl;
        cout << "                TOP B is the group of nB atoms, the fourth "
                "line correspond to the"
             << endl;
        cout << "                indexes of each atom in this group." << endl;
        cout << "                After BOND the atoms that form the axis of "
                "rotation are assigned."
             << endl;
        cout << "                i.e. 5 6" << endl;
        cout << "                  The file's format is:" << endl;
        cout << "                   __________________________ " << endl;
        cout << "                  |>> TOP A  nA              |" << endl;
        cout << "                  |1 2 5 ... nA              |" << endl;
        cout << "                  |>> TOP B  nB              |" << endl;
        cout << "                  |3 4 6 ... nB              |" << endl;
        cout << "                  |>> BOND                   |" << endl;
        cout << "                  |5 6                       |" << endl;
        cout << "                  |__________________________|\n" << endl;
        cout << GRE << "  Optional flags!" << RST << endl;
        cout << "  -o [string]   Output name, generic name for the output. "
                "Default '"
             << DEFAULTNAME << "'." << endl;
        cout << "  -T [float]    Select a temperature in Kelvin. Default T = "
                "298.15 K."
             << endl;
        cout << "  -r [float] [float] [int] " << endl;
        cout << "                Range of temperature." << endl;
        cout << "                 The first value is initial temperature."
             << endl;
        cout << "                 The second value is final temperature."
             << endl;
        cout << "                 The third value is number of steps > 1."
             << endl;
        cout << "  -s  [int]     Symmetry number (sigma)." << endl;
        cout << "  -I  [int]     Indicates the kind of Inertia moment I(2,n). "
                "n can be = 1, 2 or 3."
             << endl;
        cout << "                  - I(2,1) the moment of inertia of the "
                "rotating group is computed"
             << endl;
        cout << "                    about the axis containing the twisting "
                "bond."
             << endl;
        cout << "                  - I(2,2) the moment of inertia of the "
                "rotating group is computed"
             << endl;
        cout
            << "                    about the axis  parallel to the  bound but "
               "passing through the"
            << endl;
        cout << "                    center of mass of rotating group." << endl;
        cout << "                  - I(2,3) the moment of inertia of the "
                "rotating group is computed"
             << endl;
        cout << "                    about the axis  passing  through  the  "
                "centers-of-mass of both"
             << endl;
        cout << "                    rotating groups and the remainder of the "
                "molecule. "
             << endl;
        cout << "                  - Default n = 1, equivalent to I(2,1). "
             << endl;
        cout << "  -n  [int]     Indicates the number  of points  in the "
                "Fourier  Grid  Hamiltonian"
             << endl;
        cout << "                method (H size). Default value is 501 points."
             << endl;
        cout << "  -h, -H        Display help." << endl;
        cout << "  -v, -V        Display version." << endl;
        cout << "  -l            Display LAPACK version." << endl;

        exit(EXIT_SUCCESS);
    }
}

runCommands::runCommands(int argc, char *argv[]) {
    nargc = argc;
    if (argc == 1) {
        printE();
        cout << " The program needs more input parameters" << endl;
        cout << "            Try " << argv[0] << " -h to display the help."
             << endl;
        exit(EXIT_FAILURE);
    }

    helper(argv);

    getVersion(argv);
    getLapackVersion(argv);

    geometry_name = getName("-g", argv);
    potential_name = getName("-p", argv);
    output_name = getName("-o", argv);
    tops_name = getName("-t", argv);

    checkNames();

    printI();
    cout << "Geometry  name : " << geometry_name << endl;
    printI();
    cout << "Potential name : " << potential_name << endl;
    printI();
    cout << "Tops      name : " << tops_name << endl;
    printI();
    cout << "Output    name : " << output_name << endl;

    getTemperature(argv);

    getRangeTemperature(argv);
    printI();
    cout << "Single    Temp :  " << boolalpha << singleT << endl;
    printI();
    cout << "Range     Temp :  " << boolalpha << rangeT << endl;

    readTopsFile();
    selectInertia(argv);
    assignInertia(argv);
    getSigma(argv);
    getSizeH(argv);
}

void runCommands::readTopsFile() {
    ifstream finp;
    string line;
    int nA, nB;
    int at1, at2;
    int tmp;
    bool flagE = false;
    atom1Bond = -1;
    atom2Bond = -1;

    top1.clear();
    top2.clear();

    finp.open(tops_name.c_str(), std::ios::in);

    if (!finp.good()) {
        printE();
        cout << "The file [" << tops_name << "] can't be opened!" << endl;
        exit(EXIT_FAILURE);
    }

    finp.seekg(finp.beg);

    while (!finp.eof()) {
        getline(finp, line);
        if (line.find(string(">> TOP_A")) != string::npos) {
            sscanf(line.c_str(), ">> TOP_A %d", &nA);
            for (int i = 0; i < nA; i++) {
                finp >> tmp;
                top1.push_back(tmp - 1);
            }
        }
        if (line.find(string(">> TOP_B")) != string::npos) {
            sscanf(line.c_str(), ">> TOP_B %d", &nB);
            for (int i = 0; i < nB; i++) {
                finp >> tmp;
                top2.push_back(tmp - 1);
            }
        }
        if (line.find(string(">> BOND")) != string::npos) {
            finp >> at1 >> at2;
            at1--;
            at2--;
        }
    }
    finp.close();

    for (long unsigned int i = 0; i < top1.size(); i++) {
        if (at1 == top1[i])
            atom1Bond = i;
    }
    for (long unsigned int i = 0; i < top2.size(); i++) {
        if (at2 == top2[i])
            atom2Bond = i;
    }

    if (top1.size() == 0) {
        printE();
        cout << "No atoms assigned to TOP 1" << endl;
        flagE = true;
    }
    if (top2.size() == 0) {
        printE();
        cout << "No atoms assigned to TOP 2" << endl;
        flagE = true;
    }
    if (atom1Bond == -1) {
        printE();
        cout << "Atom 1 was not assigned to define the axis of rotation"
             << endl;
        flagE = true;
    }
    if (atom2Bond == -1) {
        printE();
        cout << "Atom 2 was not assigned to define the axis of rotation"
             << endl;
        flagE = true;
    }
    if (flagE == true) {
        exit(EXIT_FAILURE);
    }
}