#ifndef _SCREEN_H_
#define _SCREEN_H_
#include <iostream>
//#ifdef COLOR
// #define RED "\033[1;3;91m"
// #define PIN "\033[1;3;95m"
// #define GRE "\033[1;3;92m"
// #define RST "\033[0m"
//#else
#define RED ""
#define GRE ""
#define PIN ""
#define RST ""
//#endif

inline void printW() { std::cout << "[" << PIN << "WARNING" << RST << "]  "; }

inline void printE() { std::cout << "[" << RED << " ERROR " << RST << "]  "; }

inline void printI() { std::cout << "[" << GRE << " INFO  " << RST << "]  "; }

#endif
