#ifndef __USRDEF_H
#define __USRDEF_H

#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <float.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
//#include <direct.h>
#include <iomanip>

const double CERO   = 0.0E+000;
const double CEROP  = 1.0E-006;	// 0+ : small enough??????
const double P5     = 0.5E+000;
const double UNO    = 1.0E+000;
const double PI   = 3.14159265358979323846E+000;
const double PI_2 = 1.57079632679489661923E+000;
const double EGAMMA = 0.57721566490153286060E+000;

// some physical constant
const double ME   = 9.1093897540E-031;	// KG
const double CE   = 1.6021773349E-019;	// C
const double HBAR = 1.0545726663E-034;	// J * SEC
const double VC	  = 2.9979245800E+008;	// M / SEC
const double EPSL = 8.8541878170E-012;	// F / M

// the unit of length is A
// the unit of energy is EV
const double LUNIT = 1.0E-010;
const double HSD2M = (((HBAR * HBAR) / ME * P5) / CE) / LUNIT / LUNIT;
const double HBCME = 4.0 * HBAR / VC / ME / LUNIT / LUNIT;
const double ESEPS = CE / (4.0 * PI) / EPSL / LUNIT;
const double MIUBR  = P5 * CE * HBAR / ME;

using namespace std;

// ==> declare a structure;
struct tCubicSplineData{
	double vm;	// ==> v_m;
	double vm1;	// ==> v_{m+1};

	double a0;
	double a1;
	double a2;
	double a3;
};

#endif /* __USRDEF_H */
