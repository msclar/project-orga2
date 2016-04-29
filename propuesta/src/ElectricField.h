#ifndef ELECTRICFIELD_H
#define ELECTRICFIELD_H

#include <iostream>
#include "Malla2D.h"
#include "Entorno.h"
#include "Potential.h"

using namespace std;

class ElectricField {

    public:
        Malla2D<long double> mx;     // Malla 2D
        Malla2D<long double> my;     // Malla 2D
        Malla2D<long double> mx_k1;  // Malla 2D
        Malla2D<long double> my_k1;  // Malla 2D
        Malla2D<long double> mx_k2;  // Malla 2D
        Malla2D<long double> my_k2;  // Malla 2D
	//ElectricField(int dimx, int dimy):mx(dimx,dimy),my(dimx,dimy),mx_k1(dimx,dimy),my_k1(dimx,dimy),mx_k2(dimx,dimy),my_k2(dimx,dimy) {} ;
	long double errorConvergencia_x;
	long double errorConvergencia_y;
	Malla2D<long double> wt_x;
	Malla2D<long double> wt_y;
	
	ElectricField(Entorno& e, Potential& Phi, const long double xecx, const long double xecy, const long double xwtx, const long double xwty);
};

// Salida
ostream& operator <<(ostream& o, const ElectricField& s);

#endif // POTENTIAL_H

