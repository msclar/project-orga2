#include "Specie.h"
#include <math.h>

ostream& operator <<(ostream& o, const Specie& s){
         o << "Especie: " << s.Name << " z: " << s.z << " D[0][0]: " << s.D[0][0] << " u[0][0]: " << s.u[0][0] << " C0: " << s.C0 << endl;
         return o;
}



void Specie::optimizarParaCalculosDominioSobreK1NernstPlankReactionEquation(const Entorno& e)
{
  
	optimizarParaCalculosDominioSobreK1FickReactionEquation(e);
	
	for(int i=1; i<e.ii; i++){
		for(int j=1; j<e.jj; j++){
		       opt_z_z_u[i][j] =  opt_z_Over_absz *   u[i][j];

		       opt_z_D[i][j] =   z *  D[i][j];

		       opt_abs_z_u[i][j] = fabs(  z) *   u[i][j];

		       opt_4_abs_z_u[i][j] =  opt_4_absz *    u[i][j];

		       opt_z_z_u_dt_Over_h2[i][j]  = opt_z_z_dt_Over_h2  * u[i][j];

		       opt_z_z_u_dt_Over_4h2[i][j] = opt_z_z_dt_Over_4h2 * u[i][j];
		      
		       
// 		      if (i==99 && j==70) {
// 		      
// 			cout << "opt_RestauEnY[i][j] ("<< i << "," << j << ")"<< opt_RestauEnY[i][j]<< endl;
// 			cout << "u[i][j+1] ("<< i << "," << j+1 << ")"<< u[i][j+1]<< endl;
// 			cout << "u[i][j-1] ("<< i << "," << j-1 << ")"<< u[i][j-1]<< endl;
// 		      }

		       
		}
	}
}

void Specie::optimizarParaCalculosDominioSobreK2NernstPlankReactionEquation(const int i, const int j, Potential& Phi){
  optimizarParaCalculosDominioSobreK2FickReactionEquation(i,j);
  opt_RestaEnX_Plus_RestaEnY[i][j]  =  opt_restaEnX[i][j] * Phi.opt_restaEnX[i][j] +  opt_restaEnY[i][j] * Phi.opt_restaEnY[i][j];
}

void Specie::assignK2electroNeutrality(const int i, const int j, map<string,Specie>& species){
    long double concentration=0.0;
    for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
		// Specie
		if(iMapSpec->second.Name!=Name) {
			concentration+=(iMapSpec->second.z * iMapSpec->second.m_k2[i][j]);
		}
    }   
    m_k2[i][j] = - concentration / z;
}

void Specie::assignK2nernstPlankReactionEquation(const int i, const int j, Potential& Phi, const long double terminoReactivoSINSpecieEnCuestion, const long double terminoReactivoCONSpecieEnCuestion){

	m_k2[i][j]  =	(
		+ opt_dt_D_Over_h2[i][j] * opt_suma_vecinos[i][j]
		+ opt_z_z_u_dt_Over_4h2[i][j] * opt_RestaEnX_Plus_RestaEnY[i][j]
		+ terminoReactivoSINSpecieEnCuestion
		+ m[i][j]
		)
		/
		(
		  opt_One_Plus_dt_D_Over_h2[i][j]
		- opt_z_z_u_dt_Over_h2[i][j] * Phi.opt_sumaPhiConRestaPhi[i][j]
		- terminoReactivoCONSpecieEnCuestion
		)
		;
}

void Specie::optimizarParaCalculosDominioSobreK1FickReactionEquation(const Entorno& e){
	for(int i=1; i<e.ii; i++){
		for(int j=1; j<e.jj; j++){
      opt_dt_D_Over_h2[i][j] =  opt_dt_Over_h2 *  D[i][j];
      opt_One_Plus_dt_D_Over_h2[i][j] = 1.0 + opt_dt_By_4_Over_h2 * D[i][j];
		}
	}
}

void Specie::optimizarParaCalculosDominioSobreK2FickReactionEquation(const int i, const int j){
  opt_restaEnX[i][j]     =  m_k2[i+1][j] -  m_k2[i-1][j];
  opt_restaEnY[i][j]     =  m_k2[i][j+1] -  m_k2[i][j-1];
  opt_suma_vecinos[i][j] =  m_k2[i+1][j] +  m_k2[i-1][j] +  m_k2[i][j+1] +  m_k2[i][j-1];
}

void Specie::assignK2FickReactionEquation(const int i, 
																					const int j, 
																					const long double terminoReactivoSINSpecieEnCuestion, 
																					const long double terminoReactivoCONSpecieEnCuestion){

	m_k2[i][j]  =	(
			+ opt_dt_D_Over_h2[i][j] * opt_suma_vecinos[i][j]
			+ terminoReactivoSINSpecieEnCuestion
			+ m[i][j]
			)
			/
			(
			  opt_One_Plus_dt_D_Over_h2[i][j]
			- terminoReactivoCONSpecieEnCuestion
			)
			;

			
}

bool Specie::isConvergenceDiffusionConditionOk(const Entorno& e){
  bool ret=true;
  int i=0;
  int j=0;
  while(i<=e.ii && ret){
    j=0;
    while(j<=e.jj && ret){
      if((e.dt * D[i][j] / (e.h*e.h)) >= 0.25 ){
	ret = false;
      }
      j++;
    }
    i++;
  }
  if(!ret){
    cout << Name << ": ";
    cout << " e.dt = " << e.dt;
    cout << " D[" << (i-1) << "][" << (j-1) << "] = " << D[i-1][j-1];
    cout << " e.h = " << e.h << " -> (e.dt *  D[" << (i-1) << "][" << (j-1) << "] / (e.h*e.h) = " << (e.dt * D[i-1][j-1] / (e.h*e.h));
    cout << "  !!! Problemas de convergencia!!!"; 
    cout << endl;
  }
  else{
    cout << Name << " -> Ok!!!" << endl; 
  }
  return ret;
};


void Specie::reoptimizar(const Entorno& e){
	opt_dt_Over_h2      =       e.dt / ( e.h * e.h);  
	opt_dt_By_4_Over_h2 = 4.0 * opt_dt_Over_h2;  
	opt_dt_Over_4h2     =       e.dt / ( 4.0 * e.h * e.h);  
	if(z!=0){
	  opt_z_Over_absz = (  z / fabs(  z));
	}
	else{
	  opt_z_Over_absz = 0.0; // en el caso de especies electroneutras
	}
	opt_4_absz = 4.0 * fabs(   z);
	opt_z_z_dt_Over_h2 = opt_z_Over_absz * e.dt / ( e.h * e.h);
	opt_z_z_dt_Over_4h2 = opt_z_z_dt_Over_h2 / 4.0;
}