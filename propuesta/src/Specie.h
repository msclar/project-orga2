#ifndef SPECIE_H
#define SPECIE_H

#include <iostream>
#include "Malla2D.h"
#include "Entorno.h"
#include "Potential.h"
#include <map>

using namespace std;

class Potential;
class Specie {

    public:
		string		Name;	// Name
		string		fileName;	// fileName is the same as Name but without symbols + and -
		int		z;	// number of charges carried by ions
		long double	C0;	// Initial concentration
		
		Malla2D<long double> 	m;	// Malla 2D Concentration
		Malla2D<long double> 	m_k1;	// Malla 2D Concentration
		Malla2D<long double> 	m_k2;	// Malla 2D Concentration
		Malla2D<long double> 	D;	// Malla 2D Diffusion Coefficient
		Malla2D<long double>	u;	// Malla 2D Mobility			
		
		long double	errorConvergencia;	// cuando itera
		Malla2D<long double>  wt;
		
		Malla2D<long double> opt_z_z_u;
		Malla2D<long double> opt_z_D;
		Malla2D<long double> opt_abs_z_u;
		Malla2D<long double> opt_4_abs_z_u;
		Malla2D<long double> opt_dt_D_Over_h2;
		Malla2D<long double> opt_z_z_u_dt_Over_h2;
		Malla2D<long double> opt_z_z_u_dt_Over_4h2;
		Malla2D<long double> opt_One_Plus_dt_D_Over_h2;
		Malla2D<long double> opt_restaEnX;                 // se actualiza en cada i,j de k2
		Malla2D<long double> opt_restaEnY;                 // se actualiza en cada i,j de k2
		Malla2D<long double> opt_suma_vecinos;             // se actualiza en cada i,j de k2
		Malla2D<long double> opt_RestaEnX_Plus_RestaEnY;   // se actualiza en cada i,j de k2
		
		long double opt_dt_Over_h2;
		long double opt_dt_By_4_Over_h2;
		long double opt_dt_Over_4h2;
		long double opt_z_Over_absz;
		long double opt_4_absz;
		long double opt_z_z_dt_Over_h2;
		long double opt_z_z_dt_Over_4h2;
			
		Specie(const string xn, const string xfn,const Entorno& e, const int xz, const long double xc0, 
			   const long double xd, const long double xec, const long double xwt)
			:Name(xn),
			fileName(xfn),
			z(xz), C0(xc0), 
			m(C0,e.ii+1, e.jj+1),
			m_k1(C0,e.ii+1, e.jj+1),
			m_k2(C0,e.ii+1, e.jj+1),
			D(xd,e.ii+1, e.jj+1), 
			u(xd*e.F/(e.R*e.T),e.ii+1, e.jj+1),
			errorConvergencia(xec), 
			wt(xwt,e.ii+1, e.jj+1), 
			opt_z_z_u(e.ii+1,e.jj+1),
			opt_z_D(e.ii+1,e.jj+1),
			opt_abs_z_u(e.ii+1,e.jj+1),
			opt_4_abs_z_u(e.ii+1,e.jj+1),
			opt_dt_D_Over_h2(e.ii+1,e.jj+1),
			opt_z_z_u_dt_Over_h2(e.ii+1,e.jj+1),
			opt_z_z_u_dt_Over_4h2(e.ii+1,e.jj+1),
			opt_One_Plus_dt_D_Over_h2(e.ii+1,e.jj+1),
			opt_restaEnX(e.ii+1,e.jj+1),                // se actualiza en cada i,j de k2
			opt_restaEnY(e.ii+1,e.jj+1),                // se actualiza en cada i,j de k2
			opt_suma_vecinos(e.ii+1,e.jj+1),            // se actualiza en cada i,j de k2
			opt_RestaEnX_Plus_RestaEnY(e.ii+1,e.jj+1)   // se actualiza en cada i,j de k2
			{ 
				// Optimizacion
				reoptimizar(e);
			} ;
		Specie(){};
		// Funciones calculos dominio

		// calculo dominio sin corriente
		void optimizarParaCalculosDominioSobreK1FickReactionEquation(const Entorno& e);
		void optimizarParaCalculosDominioSobreK2FickReactionEquation(const int i, const int j);
		void assignK2FickReactionEquation(const int i, const int j, const long double terminoReactivoSINSpecieEnCuestion, const long double terminoReactivoCONSpecieEnCuestion);
		// calculo dominio con corriente
		void optimizarParaCalculosDominioSobreK1NernstPlankReactionEquation(const Entorno& e);
		void optimizarParaCalculosDominioSobreK2NernstPlankReactionEquation(const int i, const int j, Potential& Phi);
		void assignK2electroNeutrality(const int i, const int j, map<string,Specie>& species);
		void assignK2nernstPlankReactionEquation(const int i, const int j, Potential& Phi, const long double terminoReactivoSINSpecieEnCuestion, const long double terminoReactivoCONSpecieEnCuestion);

		bool isConvergenceDiffusionConditionOk(const Entorno& e);
		
		// optimizaciones
		void reoptimizar(const Entorno& e);
		
};

// Salida
ostream& operator <<(ostream& o, const Specie& s);
#endif // SPECIE_H
