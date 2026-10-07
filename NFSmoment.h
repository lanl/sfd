/*
	Includes and Defines
*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include <complex>
#ifdef OMP
#include <omp.h>
#endif

using namespace std;

/* Sav Data structure */
typedef struct { double x; double y; double z;} Vector;
typedef struct { int type; Vector r; Vector s;
} SaveParticle;

