#ifndef _RVECTOR_H_
#define _RVECTOR_H_

#include <iostream>
using namespace std;

class Rvector {
  private:
    double x;
    double y;
    double z;

  public:
    Rvector(double rx = 0.0, double ry = 0.0, double rz = 0.0);
    Rvector(double *);

    void set_x(double);
    void set_y(double);
    void set_z(double);

    double get_x(void);
    double get_y(void);
    double get_z(void);

    double operator[](int);

    void operator+=(const Rvector &);
    void operator-=(const Rvector &);
    void operator*=(int k);
    void operator*=(double k);
    void operator/=(double k);
    Rvector operator+(const Rvector &);
    Rvector operator-(const Rvector &);
    Rvector operator*(int k);
    Rvector operator*(double k);
    Rvector operator/(double k);

    double dot(const Rvector &);
    double norm();

    void normalize();

    friend ostream &operator<<(ostream &o, const Rvector &);
};

#endif
