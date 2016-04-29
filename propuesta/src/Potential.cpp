#include "Potential.h"

struct parametroInicializacionAnodo{
	Entorno* e;
};

int funcionObjetivoInitAnodo (const gsl_vector * vectorIncognitas, void *params, gsl_vector * vectorFunciones)
{
	// Parametros
	Entorno*       e    = ((struct parametroInicializacionAnodo *) params)->e    ;

	// Incognitas
	long double X_Phi = gsl_vector_get (vectorIncognitas, 0);

	// Funciones objetivo
	long double F_Phi;

	//Optimizacion de calculos
	long double F_divided_By_2_R_T = e->F / (2.0 * e->R * e->T);
	long double I_1_exp = F_divided_By_2_R_T * (X_Phi + e->E0eq_I  ) ;
	long double I_1 = e->I_0I  * ( exp (-I_1_exp) -  exp (I_1_exp)) ;

	// ....Ecuacion 20 - I = I1 
	F_Phi=  I_1 - e->I;

	// Guarda las ecuaciones dentro del vector de ecuaciones
	gsl_vector_set (vectorFunciones, 0, F_Phi );
	return GSL_SUCCESS;
}



bool superaToleranciaErrorConvergenciaPhi(const Entorno& e, const Potential& Phi){

	bool res=false;

	for(int i=0;!res && (i<=e.ii);i++){
		for(int j=0;!res && (j<=e.jj);j++){
			res = res || (fabs(Phi.m_k1[i][j] - Phi.m_k2[i][j]) > Phi.errorConvergenciaInicial);
		}
	}

	return res;
}


void muestraErrorPorComponente(const Entorno& e, const Potential& Phi){
	long double maximo = -100;
	for(int i=0;i<=e.ii;i++){
		for(int j=0;j<=e.jj;j++){
			long double candidato   = fabs(  Phi.m_k1[i][j] -   Phi.m_k2[i][j]);
			if (candidato   > maximo) maximo=candidato  ;
		}
	}
	cout << "   Phi: " << maximo;
}

Potential::Potential(Entorno& e, map<string,Specie>& species, const long double xec, const long double xeci, const long double xwt):
	      m(e.ii+1,e.jj+1),
	      m_k1(e.ii+1,e.jj+1),
	      m_k2(e.ii+1,e.jj+1),
	      errorConvergencia(xec),errorConvergenciaInicial(xeci),
	      wt(xwt,e.ii+1,e.jj+1),
	      opt_restaEnX(e.ii+1,e.jj+1),
	      opt_restaEnY(e.ii+1,e.jj+1),
	      opt_sumaPhi(e.ii+1,e.jj+1),
	      opt_sumaPhiConRestaPhi(e.ii+1,e.jj+1)
	{
	////////////////////////////////////////////////////////////////////////////////////////////////
	////// IMPORTANTE: Las especies vienen YA cargadas con el valor C0 en toda la malla       //////
	//////             por lo que solo tengo que asignar aquellos puntos que quiero modificar //////
	////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////
	// Calculo Phi en el anodo (electrodo)    //
	////////////////////////////////////////////

	// Utiliza el solver de GSL
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;

	// TamaÃ±o de vector incognita y ecuaciones
	const size_t num = 1; //Phi

	// Parametros
	// Los primeros los paso vacios, ya que no se usan
	struct parametroInicializacionAnodo p = {&e};

	//Vector de Ecuaciones
	gsl_multiroot_function f = {&funcionObjetivoInitAnodo, num, &p};


	//Vector de Incognitas, valores iniciales
	gsl_vector *incognitasX = gsl_vector_alloc (num);

	// Le asigno -2.0 porque se que con ese valor siempre converge
	gsl_vector_set (incognitasX, 0, -2.0);

	//Seteo del solver
	T = gsl_multiroot_fsolver_dnewton;
	s = gsl_multiroot_fsolver_alloc (T, 1);
	gsl_multiroot_fsolver_set (s, &f, incognitasX);

	// Iteracion de convergencia
	do
	{
		iter++;
		// Itera
		status = gsl_multiroot_fsolver_iterate (s);

		if (status)   // check if solver is stuck
		break;

		status = gsl_multiroot_test_residual (s->f, 1e-9);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	// Guarda el valor obtenido en todo el anodo
	long double Phi_inicial_anodo = gsl_vector_get(s->x,0);

	// Libera memoria
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (incognitasX);

	// Calculo Ey  // OJO ESTOY UTILIZANDO u como contasnte en TODA LA MALLA!!!!
	long double valorEy = e.I / (e.F* (  species["H+"  ].u[0][0] * species["H+"  ].C0 
	                                   + species["OH-" ].u[0][0] * species["OH-" ].C0 
	                                   + species["Cl-" ].u[0][0] * species["Cl-" ].C0 
	                                   + species["Na+" ].u[0][0] * species["Na+" ].C0 
					  ));

	long double Phi_inicial_catodo = Phi_inicial_anodo - valorEy * e.minimaDistanciaAnodoCatodo;


	// Seteo el valor de Phi en toda la malla
	m.sameValue(Phi_inicial_anodo - 0.5 * (valorEy * e.minimaDistanciaAnodoCatodo));
	// Seteo el valor de Phi en los electrodos
	for(int i=1; i<e.ii; i++){
	  for(int j=1; j<e.jj; j++){
		if(e.geometriaTipoNodo[i][j]==bordeAnodo ) m[i][j] = Phi_inicial_anodo;
		if(e.geometriaTipoNodo[i][j]==bordeCatodo) m[i][j] = Phi_inicial_catodo;   
	  }
	}

	m.save("data/initial/PhiInicialAntesDeDifFIN.dat");
	// seteo el valor de phi en el dominio
	// REVISARRRRRRRRRRRRRRRRRRRRRRRRRRRR
// 	for(int i=1; i<e.ii; i++){
// 	    for(int j=1; j<e.jj; j++){
// 		Phi.m[i][j] = Phi_inicial_anodo - valorEy * e.h * j;
// 	    }
// 	}
	
	// Itero hasta converger.
	m_k1=m;
	m_k2=m;

	bool superaToleranciaError;
	long it=0;
	do{
		for(int i=1;i<e.ii;i++){
			for(int j=1;j<e.jj;j++){
				
					if(e.geometriaTipoNodo[i][j]==dominio) m_k2[i][j] = 0.25 * (m_k2[i+1][j] + m_k2[i-1][j] + m_k2[i][j+1] + m_k2[i][j-1]);
				
			}
		}
		for(int i=0;i<=e.ii;i++){
			m_k2[i][0]    = m_k2[i][1];
			m_k2[i][e.jj] = m_k2[i][e.jj-1];
		}
		for(int j=0;j<=e.jj;j++){
			m_k2[0][j]    = m_k2[1][j];
			m_k2[e.ii][j] = m_k2[e.ii-1][j];
		}

		superaToleranciaError = superaToleranciaErrorConvergenciaPhi(e, *this);

// 		if(it % 100 == 0){
// 			cout << "it = " << it << " ";
// 			muestraErrorPorComponente(e, *this);
// 			cout << endl;
// 		}

		m_k1 = m_k2;

		it++;
	} while (superaToleranciaError && (it<e.cantidadItConvergencia));

	
	m=m_k2;
	// Guardo Phi
	m.save("data/initial/PhiInicialAntesDeIterar.dat");
}


void Potential::optimizarParaCalculosDominioSobreK2(const int i, const int j){
  opt_restaEnX[i][j]  =  m_k2[i+1][j] -  m_k2[i-1][j];
  opt_restaEnY[i][j]  =  m_k2[i][j+1] -  m_k2[i][j-1];
  opt_sumaPhi[i][j]    = m_k2[i+1][j] +   m_k2[i][j+1] + m_k2[i-1][j] + m_k2[i][j-1];
  opt_sumaPhiConRestaPhi[i][j] = opt_sumaPhi[i][j] - 4.0 * m_k2[i][j];
}

void Potential::assignK2nernstPlankReactionEquation(const int i, const int j, map<string,Specie>& species){
  
  long double suma_abs_z_u_C_Phi=0.0;
  long double suma_z_D_Cs=0.0;
  long double suma_4_abs_z_u_C=0.0;
    
  for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
		// Specie
    suma_abs_z_u_C_Phi+= iMapSpec->second.opt_abs_z_u[i][j] * (0.25 *  iMapSpec->second.opt_RestaEnX_Plus_RestaEnY[i][j] +  iMapSpec->second.m_k2[i][j] * opt_sumaPhi[i][j] );
		suma_z_D_Cs       += iMapSpec->second.opt_z_D[i][j]  * (  iMapSpec->second.opt_suma_vecinos[i][j] - 4.0 *  iMapSpec->second.m_k2[i][j] );
		suma_4_abs_z_u_C  += iMapSpec->second.opt_4_abs_z_u[i][j]  *  iMapSpec->second.m_k2[i][j];
   }   

    m_k2[i][j] = (  suma_abs_z_u_C_Phi
		  + suma_z_D_Cs
		 )
		 /
		 (
		    suma_4_abs_z_u_C
		 )
		;
		
}