#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include "Malla2D.h"
#include "Entorno.h"
#include "Specie.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <map>
using namespace std;

class Specie;
class Potential {

  public:
    Malla2D<long double> m;     // Malla 2D
    Malla2D<long double> m_k1;  // Malla 2D
    Malla2D<long double> m_k2;  // Malla 2D

		long double errorConvergencia;
		long double errorConvergenciaInicial;
		
		Malla2D<long double> wt;
		
		// Optimizacion
		Malla2D<long double> opt_restaEnX; 
		Malla2D<long double> opt_restaEnY;
		Malla2D<long double> opt_sumaPhi;
		Malla2D<long double> opt_sumaPhiConRestaPhi;
	                              
		Potential(Entorno& e, map<string,Specie>& species, const long double xec, const long double xeci, const long double xwt);
		Potential(){};
		
		void optimizarParaCalculosDominioSobreK2(const int i, const int j);
		void assignK2nernstPlankReactionEquation(const int i, const int j, map<string,Specie>& species);
};

// Salida
ostream& operator <<(ostream& o, const Potential& s);

#endif // POTENTIAL_H

