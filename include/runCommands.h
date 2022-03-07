#ifndef _RUN_COMMANDS_H_
#define _RUN_COMMANDS_H_

#include <vector>
using std::vector;
#include <fstream>
#include <iostream>
#include <string>
using std::string;

#include "version.h"

using namespace std;

#define DEFAULTNAME PROJECT_NAME

class runCommands {
  private:
    string getName(const char *, char **);
    void getTemperature(char **);
    void getRangeTemperature(char **);
    void checkNames();
    bool exists_file(string);
    void helper(char **);
    void getVersion(char **);
    void getLapackVersion(char **);
    void readTopsFile();
    void selectInertia(char **);
    void getSigma(char **);
    void getSizeH(char **);

  public:
    int nargc;
    int typeI;
    int hsize;
    string geometry_name;
    string potential_name;
    string tops_name;
    string output_name;
    int atom1Bond;
    int atom2Bond;
    int sigma;
    vector<int> top1;
    vector<int> top2;
    double temp;

    int ntemp;
    double tempi, tempf, htemp;

    bool singleT = false;
    bool rangeT = false;
    runCommands(int, char **);
};
#endif