#include <cstdlib>
#include <iostream>
#include "Entorno.h"
#include "Specie.h"
#include "Potential.h"
#include "ElectricField.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "Malla2D.h"
#include <map>
#include <vector>

#include <time.h>

using namespace std;

void calculaNodoDominio(const int i, const int j, const int it, const Entorno& e, map<string,Specie>& species, Potential& Phi);
void calculaNodoDominioSinCorriente(const int i, const int j, const int it, const Entorno& e, map<string,Specie>& species);

void calculaDominioSubIteracion(int it, const Entorno& e, map<string,Specie>& species, Potential& Phi);

void calculaElectrodosSubIteracion(int n, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E);

void calculaBordesSubIteracion(Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeSuperiorSubIteracion(int j, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeInferiorSubIteracion(int j, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeIzquierdoSubIteracion(int i, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeDerechoSubIteracion  (int i, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E);

void actualizacionEspacial(map<string,Specie>& species, Potential& Phi, ElectricField& E);
void actualizacionTemporal(map<string,Specie>& species, Potential& Phi, ElectricField& E);
void actualizacionElectricField(const Entorno& e,const Potential& Phi,ElectricField& E);
void actualizacionOptimizacionCalculos(const Entorno& e,  map<string,Specie>& species);
void actualizacionCorrienteVoltaje(const int n, Entorno& e, const long double diferenciaPotencial);
void actualizacionFase(const int n, Entorno& e, map<string,Specie>& species, const long double diferenciaPotencial);
void reoptimizacionTemporal(Entorno& e, map<string,Specie>& species);

void saveAll(const Entorno& e, const int foto, const long& n, const long& it, map<string,Specie>& species, const Potential& Phi, const ElectricField& E);
void saveAll_k2(const Entorno& e, const long& n, const long& it, map<string,Specie>& species, const Potential& Phi, const ElectricField& E);

void calculaBordeAnodo (int i, int j, int n, Entorno& e, vector<string>& speciePosition, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeAnodoSinCorriente(int i, int j, int n, vector<string>& speciePosition, map<string,Specie>& species);

void calculaBordeCatodo (int i, int j, int n, Entorno& e, vector<string>& speciePosition, map<string,Specie>& species, Potential& Phi, ElectricField& E);
void calculaBordeCatodoSinCorriente(int i, int j, int n, vector<string>& speciePosition, map<string,Specie>& species);

int funcionObjetivoAnodo(const gsl_vector * vectorIncognitas, void *params, gsl_vector * vectorFunciones);

int funcionObjetivoCatodo(const gsl_vector * vectorIncognitas, void *params, gsl_vector * vectorFunciones);

void muestraErrorPorComponente(const Entorno& e, map<string,Specie>& species, const Potential& Phi_k, const ElectricField& E);
bool superaToleranciaErrorConvergencia(const Entorno& e, map<string,Specie>& species, const Potential& Phi, const ElectricField& E);




struct parametrosElectrodo{
	int i;
	int j;
	int n;
	Entorno* e;
	vector<string>* speciePosition;
	map<string,Specie>*  species;
	ElectricField* Ex;
};


int main(int argc, char *argv[])
{
	double begin = time(NULL);

	cout << "==================================================================" << endl;
	cout << "=                       Malla                                    =" << endl;
	cout << "==================================================================" << endl;

	Entorno e(50,50,5,2);
	cout << "tamano       : " << e.h*e.ii << " x " << e.h*e.jj << endl;
  
	cout << "==================================================================" << endl;
	cout << "=                     Especies                                   =" << endl;
	cout << "==================================================================" << endl;

	// long double velocidad = 0.5	; //(1/2 es mÃ¡s lento)
	long double velocidad = 1.0/1.2	; //(1/2 es mÃ¡s lento)
	
  long double  D_H_en_agua  = 9.31e-9;
	long double  D_OH_en_agua = 5.26e-9;
	// long double  D_Cl_en_agua = 2.03e-9 * 1e-1;//2.03e-9 * 1e-1; // Cl = SNARF
	long double  D_Cl_en_agua = 2.03e-9;
	long double  D_Na_en_agua = 1.33e-9;

	map<string,Specie> species;
	// Species que estan disueltas en la Solucion. No agregar la Solucion (H2O) xq si no calcula un monton de cosas de mas.
	//                                  map_id         Name    fileName  e   z   C0              D0          errorConvergencia  wt
	species.insert(pair<string,Specie> ("H+"  , Specie("H+ " , "H"     , e,  1,  1.0e-4       ,  D_H_en_agua*velocidad, 5.0e-7, 1.0 )));
	species.insert(pair<string,Specie> ("OH-" , Specie("OH-" , "OH"    , e, -1,  1.0e-4       ,  D_OH_en_agua*velocidad, 5.0e-7, 1.0 )));
	species.insert(pair<string,Specie> ("Cl-" , Specie("Cl-" , "Cl"    , e, -1,  160.0       ,  D_Cl_en_agua*velocidad, 5.0e-6, 1.0 )));
	species.insert(pair<string,Specie> ("Na+" , Specie("Na+" , "Na"    , e,  1,  160.0      ,  D_Na_en_agua*velocidad, 5.0e-6, 1.0 )));

	cout << species["H+"]    << endl;
	cout << species["OH-"]   << endl;
	cout << species["Cl-"]   << endl;
	cout << species["Na+"]   << endl;
	cout << "H2O - C0 = " << e.H2O_C0 << endl;
	cout << "velocidad: " << velocidad << endl;
	cout << endl;

	// Chequeo que se cumple la relacion: (D * dt / (h*h)) < 1/4
	bool ok = true;
	cout << endl;
	cout << "=======================================================" << endl;
	cout << "==  Chequeo de convergencia: (D * dt / (h*h)) < 1/4  ==" << endl;
	cout << "=======================================================" << endl;
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    ok = ok && iMapSpec->second.isConvergenceDiffusionConditionOk(e);
	}
	cout << endl;
	if(!ok) return 1;
	
	cout << "==================================================================" << endl;
	cout << "=      Campo electrico y potencial                              =" << endl;
	cout << "==================================================================" << endl;
	cout << "... calculando Phi inicial (iterando, aguarde por favor) ... " << endl;
	Potential Phi(e, species,  5.0e-5, 1.0e-6, 1.0);

	//Ahora el ElectricField tiene 2 coordenadas pues E = (-Phi'/dx,-Phi'/dy).
	//Usamos mx para representar la primera coordenada (sentido vertical).
	//USamos my para representar la segunda coordenada (sentido horizontal).
	ElectricField E(e, Phi, 1.0e-3, 1.0e-3, 1.0, 1.0);
	
	cout << "I       inicial        = " << e.I << endl;
	cout << "Phi     inicial anodo  = " <<  Phi.m_k2[ e.ii/2 ][        e.separacionBordeLateral ] << endl;
	cout << "Phi     inicial catodo = " <<  Phi.m_k2[ e.ii/2 ][ e.jj - e.separacionBordeLateral ] << endl;
	cout << "Voltaje inicial        = " << (Phi.m_k2[ e.ii/2 ][        e.separacionBordeLateral ] - Phi.m_k2[ e.ii/2 ][ e.jj-e.separacionBordeLateral ])  << endl;
	cout << "distancia electrodos   = " << e.x[e.jj-e.separacionBordeLateral]-e.x[e.separacionBordeLateral] << endl;
	cout << endl;
	// guardo condiciones iniciales

	cout << "==================================================================" << endl;
	cout << "=      Guardando condiciones iniciales (foto 0)                  =" << endl;
	cout << "==================================================================" << endl;
	{
		long int cero=0;
		saveAll(e, cero, cero, cero, species, Phi, E);
	}
	cout << endl;
	cout << "==================================================================" << endl;

	int foto = 1;
	ofstream voltage  ((e.initialValuesFolder+"voltage.dat"  ).c_str());
	ofstream corriente((e.initialValuesFolder+"corriente.dat").c_str());
	ofstream tiempo   ((e.initialValuesFolder+"tiempo.dat"   ).c_str());
	ofstream deltaT   ((e.initialValuesFolder+"deltaT.dat"   ).c_str());
	ofstream status   ((e.initialValuesFolder+"status.dat"   ).c_str());
	ofstream subit    ((e.initialValuesFolder+"subit.dat"   ).c_str());
	
	// Iteracion temporal
	for(long n=e.n0;(n<=e.nn) && (e.ctrl_deltaTAcumulado <= e.umbralTiempoTotal);n++){
		bool superaToleranciaError;
		long it = 0;

		voltage.flush()   << (Phi.m[e.ii/2][e.separacionBordeLateral] - Phi.m[e.ii/2][e.jj-e.separacionBordeLateral]) << " " ;
		corriente.flush() << e.I << " ";
		deltaT.flush()    << e.dt << " ";
		tiempo.flush()    << e.ctrl_deltaTAcumulado << " ";
		status.flush()    << e.ctrl_faseActual << " ";

    // Optimizacion de calculos porque pudo haber cambiado D, u o dt
    actualizacionOptimizacionCalculos(e,species);
		
		do{
			calculaDominioSubIteracion(it, e, species, Phi);
			calculaElectrodosSubIteracion(n, e, species, Phi, E);
			calculaBordesSubIteracion(e, species, Phi, E);

			superaToleranciaError = superaToleranciaErrorConvergencia(e, species, Phi, E);

			if ((it%1000==0) && (it!=0)) {
				cout 	<< "Subit: " << it
							<< "; Tiempo: " << e.ctrl_deltaTAcumulado
							<< "; Fase: " << e.ctrl_faseActual
							<< "; I: " << e.I
							<< "; V: " << (Phi.m_k2[ e.ii/2 ][ e.separacionBordeLateral ]-Phi.m_k2[ e.ii/2 ][ e.jj-e.separacionBordeLateral ])
							<< "; Phi_anodo = " << Phi.m_k2[ e.ii/2 ][ e.separacionBordeLateral ]
							<< "; Phi_catodo = " << Phi.m_k2[ e.ii/2 ][ e.jj- e.separacionBordeLateral ]
							<< "; n: " << n;

				muestraErrorPorComponente(e, species, Phi, E);
				cout << endl;
				saveAll_k2(e, n, it, species, Phi, E);
			}

			actualizacionEspacial(species, Phi, E);
			it++;

		} while (superaToleranciaError && (it<e.cantidadItConvergencia));
		
		//TODO
		//actualizarTemperatura(e,k,sigma,Phi,T)
		
		subit.flush() << it << " ";
		/*status*/ 
		if(n%1000==0)  cout << "n: " << n << " e.ctrl_deltaTAcumulado: " << e.ctrl_deltaTAcumulado << " it: " << it << endl;

		actualizacionElectricField(e,Phi,E);
		actualizacionTemporal(species, Phi, E);
		
		// Testeo si diverge. En caso afirmativo termino la ejecucion
		if (isnan(species["Cl-"].m[0][0]) || isnan(species["Cl-"].m[e.ii][e.jj])) {
			saveAll(e, foto, n, it, species, Phi, E);
			cout << "SALI POR NAN" << endl;
			n=e.nn+1; // termina la ejecucion;
		}

		long double diferenciaPotencial = Phi.m[e.ii/2][e.separacionBordeLateral] - Phi.m[e.ii/2][e.jj-e.separacionBordeLateral];
		actualizacionCorrienteVoltaje(n,e,diferenciaPotencial);
		actualizacionFase(n,e,species,diferenciaPotencial);

		// cuento el tiempo del pulso aplicado.
		e.ctrl_tiempoActualPulso = e.ctrl_tiempoActualPulso + e.dt;
		e.ctrl_deltaTAcumulado = e.ctrl_deltaTAcumulado + e.dt;
		
		/*status*/
		if ((e.ctrl_deltaTAcumulado/foto) >=e.tiempoFoto)  {
			saveAll(e, foto, n, it, species, Phi, E);
			cout << "Aca salvo con foto=" << foto << endl;
			foto++;
		}

		// saveAll(e, n, n, it, species, Phi, E);
	}

	voltage.close();
	corriente.close();
	tiempo.close();
	deltaT.close();
	status.close();

	double time_elapsed = time(NULL) - begin;
	cout << "Tiempo transcurrido: " << time_elapsed << endl;

	return EXIT_SUCCESS;
}

void saveAll(const Entorno& e, const int foto, const long& n, const long& it, 
			 map<string,Specie>& species, const Potential& Phi, const ElectricField& E)
{
	cout << "foto: " << foto << " n: " << n << " subIt: " << it << endl;
	stringstream prefix;prefix << e.commit << "_foto_" << foto << "_";

	string folder = "data/";

	// save species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m.save(folder + prefix.str()+iMapSpec->second.fileName+".dat");
	    // Difusion guv_membraneCoefficient
	    iMapSpec->second.D.save(folder + prefix.str()+"D_"+iMapSpec->second.fileName+".dat");
	}   
	    
	// Grabo campo elecrtico y potencial
	E.mx.save(folder + prefix.str()+ "Ex.dat");
	E.my.save(folder + prefix.str()+ "Ey.dat");
	Phi.m.save(folder + prefix.str()+ "Phi.dat");
}

void saveAll_k2(const Entorno& e, const long& n, const long& it, map<string,Specie>& species, const Potential& Phi, const ElectricField& E){
	cout << "n: " << n << " subIt: " << it << endl;
	stringstream prefix;prefix << e.commit << "_n_" << n << "_subit_" << it << "_";

	string folder = "data/";

	// save species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2.save(folder + prefix.str()+iMapSpec->second.fileName+".dat");
	    // Difusion guv_membraneCoefficient
	    iMapSpec->second.D.save(folder + prefix.str()+"D_"+iMapSpec->second.fileName+".dat");
	    // Difusion guv_membraneCoefficient
	    iMapSpec->second.u.save(folder + prefix.str()+"u_"+iMapSpec->second.fileName+".dat");
	}   

	// ElectricField & Potencial
	E.mx_k2.save(folder + prefix.str()+ "Ex.dat");
	E.my_k2.save(folder + prefix.str()+ "Ey.dat");
	Phi.m_k2.save(folder + prefix.str()+ "Phi.dat");
}


void actualizacionEspacial(map<string,Specie>& species, Potential& Phi, ElectricField& E){
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k1 = iMapSpec->second.m_k2;
	}   
	  E.mx_k1=       E.mx_k2;
	  E.my_k1=       E.my_k2;
	Phi.m_k1 =      Phi.m_k2;
}

void actualizacionTemporal(map<string,Specie>& species, Potential& Phi, ElectricField& E){

	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m = iMapSpec->second.m_k2;
	}   
	  E.mx=       E.mx_k2;
	  E.my=       E.my_k2;
	Phi.m =      Phi.m_k2;
}

void muestraErrorPorComponente(const Entorno& e, map<string,Specie>& species, const Potential& Phi, const ElectricField& E){
	long double maximo;
	
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
    maximo = -100;
    for(int i=0;i<=e.ii;i++){
	    for(int j=0;j<=e.jj;j++){
		    long double candidato   = fabs(  iMapSpec->second.m_k1[i][j] -   iMapSpec->second.m_k2[i][j]);
		    if (candidato > maximo) maximo=candidato;
	    }
    }
    cout << " Errores Maximos: " << iMapSpec->second.Name << ": " << maximo;
	}   

	maximo = -100;
	for(int i=0;i<=e.ii;i++){
		for(int j=0;j<=e.jj;j++){
			long double candidato   = fabs(  Phi.m_k1[i][j] -   Phi.m_k2[i][j]);
			if (candidato   > maximo) maximo=candidato  ;
		}
	}
	cout << "   Phi: " << maximo;

	//No estoy muy seguro de si esto hace falta. Solo tiene sentido en el anodo y catodo para Ex. Habria que revisarlo.
	maximo = -100;
	for(int i=0;i<=e.ii;i++){
		for(int j=0;j<=e.jj;j++){
			long double candidato   = fabs(  E.mx_k1[i][j] -   E.mx_k2[i][j]);
			if (candidato   > maximo) maximo=candidato  ;
		}
	}
	cout << "   Ex: " << maximo;

	maximo = -100;
	for(int i=0;i<=e.ii;i++){
		for(int j=0;j<=e.jj;j++){
			long double candidato   = fabs(  E.my_k1[i][j] -   E.my_k2[i][j]);
			if (candidato   > maximo) maximo=candidato  ;
		}
	}
	cout << "   Ey: " << maximo;
}

// La idea de esta funcion es revisar de manera inteligente las diferencias entre el valor anterior y el actual
bool superaToleranciaErrorConvergencia(const Entorno& e, map<string,Specie>& species, const Potential& Phi, const ElectricField& E){
	// en cuanto la diferencia supera el maximo permitido deja de buscar.

	bool res=false;

	for(int i=0;!res && (i<=e.ii);i++){
		for(int j=0;!res && (j<=e.jj);j++){
			res = res || (fabs(Phi.m_k1[i][j] - Phi.m_k2[i][j]) > Phi.errorConvergencia);
		}
	}

	for(int i=0;!res && (i<=e.ii);i++){
		for(int j=0;!res && (j<=e.jj);j++){
			res = res || (fabs(E.mx_k1[i][j] - E.mx_k2[i][j])>E.errorConvergencia_x);
		}
	}

	for(int i=0;!res && (i<=e.ii);i++){
		for(int j=0;!res && (j<=e.jj);j++){
			res = res || (fabs(E.my_k1[i][j] - E.my_k2[i][j])>E.errorConvergencia_y);
		}
	}

	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); !res && iMapSpec!=species.end(); ++iMapSpec){
	    for(int i=0;!res && (i<=e.ii);i++){
		    for(int j=0;!res && (j<=e.jj);j++){
			    res = res || (fabs(iMapSpec->second.m_k1[i][j] - iMapSpec->second.m_k2[i][j]) > iMapSpec->second.errorConvergencia);
		    }
	    }
	}   
	

	return res;
}

void calculaBordeCatodo (int i, int j, int n, Entorno& e, vector<string>& speciePosition, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	
	// Utiliza el sover de GSL
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;

	// Tamano de vector ingognita y ecuaciones
	const size_t num = speciePosition.size()+1; 

	// ParÃ¡metros
	struct parametrosElectrodo p = {i, j, n, &e, &speciePosition, &species, &E};

	//Vector de Ecuaciones
	gsl_multiroot_function f = {&funcionObjetivoCatodo, num, &p};

	//Vector de Incognitas, valores iniciales
	gsl_vector *incognitasX = gsl_vector_alloc (num);
	// Nodos Ã¡nodo, Dominio	e Catodo
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
	  	gsl_vector_set (incognitasX, iPos,  species[speciePosition[iPos]].m_k1[i][j] );
	}
	gsl_vector_set (incognitasX, speciePosition.size(), E.my_k1[i][j] );
	//Seteo del solver
	T = gsl_multiroot_fsolver_dnewton;
	s = gsl_multiroot_fsolver_alloc (T, speciePosition.size()+1);
	gsl_multiroot_fsolver_set (s, &f, incognitasX);

	// IteraciÃ³n de convergencia
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

	// Guarda los valores obtenidos en el nodo
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
		species[speciePosition[iPos]].m_k2[i][j] = species[speciePosition[iPos]].wt[i][j]   * gsl_vector_get(s->x,iPos) + (1.0 - species[speciePosition[iPos]].wt[i][j])   *  species[speciePosition[iPos]].m_k1[i][j];
	}	
	E.my_k2[i][j] = E.wt_y[i][j]  * gsl_vector_get(s->x,speciePosition.size()) + (1.0 - E.wt_y[i][j])  * E.my_k1[i][j];
	
	Phi.m_k2[i][j] = -E.my_k2[i][j] * e.h + Phi.m_k2[i][j-1];
	//Le asignamos 0 a Ey por completitud.
	E.mx_k2[i][j] = 0;

	// Libera memoria
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (incognitasX);
}

int funcionObjetivoCatodo (const gsl_vector * vectorIncognitas, void *params, gsl_vector * vectorFunciones)
{
	// Parametros
	int                 i              = ((struct parametrosElectrodo *) params)->i              ;
	int                 j              = ((struct parametrosElectrodo *) params)->j              ;
//	int                 n              = ((struct parametrosElectrodo *) params)->n              ;
	Entorno*            e              = ((struct parametrosElectrodo *) params)->e              ;
	vector<string>*     speciePosition = ((struct parametrosElectrodo *) params)->speciePosition ;
	map<string,Specie>* species        = ((struct parametrosElectrodo *) params)->species        ;

	// valores iniciales incognitas
	Malla1D<long double> X_Specie(speciePosition->size());
	long double X_Ey; 
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
		    X_Specie[iPos]=gsl_vector_get (vectorIncognitas, iPos);
	}	
	X_Ey  = gsl_vector_get (vectorIncognitas, speciePosition->size());

	// Funciones objetivo
	Malla1D<long double> F_Specie(speciePosition->size());
	long double F_Ey;

	// Derivadas
	Malla1D<long double> derivada_X(speciePosition->size());
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	    string specieName = (*speciePosition)[iPos];
	    derivada_X[iPos]=( 3.0 * X_Specie[iPos]   - 4.0 *  (*species)[specieName].m_k2[i][j-1] +  (*species)[specieName].m_k2[i][j-2]) / (2.0*e->h);;
	}	

	long double Ey_By_n = - (X_Ey * e->normalNodoX[i][j]);
	
	// Iguala el flujo a cero para todas las especies
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	  string specieName = (*speciePosition)[iPos];  
	  F_Specie[iPos]= -   (*species)[specieName].D[i][j]         * ( derivada_X[iPos] * e->normalNodoX[i][j])
			  -   (*species)[specieName].opt_z_z_u[i][j] *  X_Specie[iPos]    * Ey_By_n;
	}
	
	// En el caso de OH, le agrego al flujo la corriente
	int iPosOH=0;
	for(iPosOH=0;iPosOH<(int)speciePosition->size() && (*speciePosition)[iPosOH]!="OH-";iPosOH++);
	F_Specie[iPosOH]+=e->I / e->F; 
	
	// Electroneutralidad
	F_Ey  =	0.0;
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	  string specieName = (*speciePosition)[iPos]; 
	  F_Ey+=X_Specie[iPos]* (*species)[specieName].z;
	}
	
	// Guarda las ecuaciones dentro del vector de ecuaciones
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	    gsl_vector_set (vectorFunciones, iPos, F_Specie[iPos]);
	}
	
	gsl_vector_set (vectorFunciones, speciePosition->size(), F_Ey   );
	return GSL_SUCCESS;
}

void calculaElectrodosSubIteracion(int n, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	int jAnodo  = 0    + e.separacionBordeLateral;
	int jCatodo = e.jj - e.separacionBordeLateral;
	// electrodos
	
	// Especies que voy a usar, mantengo un orden para el vector de incognitas y funciones de las funciones objetivo del gsl
	vector<string> speciePosition;
	speciePosition.push_back("H+");
	speciePosition.push_back("OH-");
	speciePosition.push_back("Cl-");
	speciePosition.push_back("Na+");
// 	speciePosition.push_back("LH+");
// 	speciePosition.push_back("LOH-");


	if(e.ctrl_faseActual==4){
	    for(int i=e.separacionBordeSup;i<=(e.ii-e.separacionBordeSup);i++){
 		calculaBordeAnodoSinCorriente  (i, jAnodo , n, speciePosition, species);
 		calculaBordeCatodoSinCorriente (i, jCatodo, n, speciePosition, species);
	    }
	}	
	else{
	    for(int i=e.separacionBordeSup;i<=(e.ii-e.separacionBordeSup);i++){
 		calculaBordeAnodo (i,jAnodo , n, e, speciePosition, species, Phi, E);
		calculaBordeCatodo(i,jCatodo, n, e, speciePosition, species, Phi, E);
	    }
	}	
}

void calculaDominioSubIteracion(int it, const Entorno& e, map<string,Specie>& species, Potential& Phi){
	//El if lo meto afuera (por mas que quede feo el codigo) para que no pregunte por cada nodo en que e.ctrl_faseActual esta.
	if(e.ctrl_faseActual==4){
		for(int i=1;i<e.ii;i++){
			for(int j=1;j<e.jj;j++){
				if(e.geometriaTipoNodo[i][j]==dominio)
					calculaNodoDominioSinCorriente(i, j, it, e, species);
			}
		}
	}else{
		for(int i=1;i<e.ii;i++){
			for(int j=1;j<e.jj;j++){
				if(e.geometriaTipoNodo[i][j]==dominio)
				    calculaNodoDominio(i, j, it, e, species, Phi);
			}
		}
	}
}

void calculaNodoDominio(const int i, const int j, const int it, const Entorno& e, map<string,Specie>& species, Potential& Phi){

	// 1. Optimizacion de calculos
	Phi.optimizarParaCalculosDominioSobreK2(i, j);
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
		// Specie
	  iMapSpec->second.optimizarParaCalculosDominioSobreK2NernstPlankReactionEquation(i, j, Phi);
	}   

	// 2. Calculo de especies
	// Ejemplo pra ver como pasar parametros de Termino reactivo: R_H+ = kwb H2O - kwf H OH
	species["H+"].assignK2nernstPlankReactionEquation(i, j, Phi, e.opt_H2O_By_k_wb_dt, - e.opt_k_wf_dt * species["OH-"].m_k2[i][j]);
	species["OH-"].assignK2nernstPlankReactionEquation(i, j, Phi, e.opt_H2O_By_k_wb_dt, - e.opt_k_wf_dt * species["H+"].m_k2[i][j] );
	species["Cl-"].assignK2nernstPlankReactionEquation(i, j, Phi, 0.0, 0.0);
	species["Na+"].assignK2electroNeutrality(i, j, species);

	//Calculo de Phi
	Phi.assignK2nernstPlankReactionEquation(i, j, species);
}

void calculaNodoDominioSinCorriente(const int i, const int j, const int it, const Entorno& e, map<string,Specie>& species){
  	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
			// Specie
	  	iMapSpec->second.optimizarParaCalculosDominioSobreK2FickReactionEquation(i, j);
		}   

	// 2. Calculo de especies
	
	// Ejemplo para ver como pasar parametros de Termino reactivo: R_H+ = kwb H2O - kwf H OH
	// entonces (1) - kwf OH  * dt   (es (1) xq incluia a H+)
	//          (2)   kwb H2O * dt  
	species["H+"].assignK2FickReactionEquation(i, j, e.opt_H2O_By_k_wb_dt, - e.opt_k_wf_dt * species["OH-"].m_k2[i][j]);
	species["OH-"].assignK2FickReactionEquation(i, j, e.opt_H2O_By_k_wb_dt, - e.opt_k_wf_dt * species["H+"].m_k2[i][j] );
	species["Cl-"].assignK2FickReactionEquation(i, j, 0.0, 0.0);
	species["Na+"].assignK2electroNeutrality(i, j, species);  
	
}

void calculaBordeAnodo (int i, int j, int n, Entorno& e, vector<string>& speciePosition, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	
	// Utiliza el sover de GSL
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;

	// Tamano de vector ingognita y ecuaciones
	const size_t num = 6; //(int)(speciePosition.size()+2);; //  species, Ex, Phi

	// Parametros
	struct parametrosElectrodo p = {i, j, n, &e, &speciePosition, &species, &E};

	//Vector de Ecuaciones
	gsl_multiroot_function f = {&funcionObjetivoAnodo, num, &p};

	//Vector de Incognitas, valores iniciales
	gsl_vector *incognitasX = gsl_vector_alloc (num);
	// valores iniciales incognitas 
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
	  	gsl_vector_set (incognitasX, iPos,  species[speciePosition[iPos]].m_k1[i][j] );
	}
	gsl_vector_set (incognitasX, (int)speciePosition.size()  ,  E.my_k1[i][j]);
	gsl_vector_set (incognitasX, (int)speciePosition.size()+1, Phi.m_k1[i][j]);

	//Seteo del solver
	T = gsl_multiroot_fsolver_dnewton;
// 	s = gsl_multiroot_fsolver_alloc (T, speciePosition.size()+2);
	s = gsl_multiroot_fsolver_alloc (T, 6);
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

	// Guarda los valores obtenidos en el nodo
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
	    	  species[speciePosition[iPos]].m_k2[i][j] = species[speciePosition[iPos]].wt[i][j]   * gsl_vector_get(s->x,iPos) + (1.0 - species[speciePosition[iPos]].wt[i][j])   *  species[speciePosition[iPos]].m_k1[i][j];
	}	
	  E.my_k2[i][j] =   E.wt_y[i][j]  * gsl_vector_get(s->x,speciePosition.size()  ) + (1.0 -   E.wt_y[i][j])  *   E.my_k1[i][j];
	 Phi.m_k2[i][j] =   Phi.wt[i][j]  * gsl_vector_get(s->x,speciePosition.size()+1) + (1.0 -   Phi.wt[i][j])  *  Phi.m_k1[i][j];
	
	
	//Le asignamos 0 a Ey por completitud.
	E.mx_k2[i][j] = 0;

	// Libera memoria
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (incognitasX);
}

int funcionObjetivoAnodo (const gsl_vector * vectorIncognitas, void *params, gsl_vector * vectorFunciones)
{
	// Parametros
	int             i    = ((struct parametrosElectrodo *) params)->i    ;
	int             j    = ((struct parametrosElectrodo *) params)->j    ;
//	int             n    = ((struct parametrosElectrodo *) params)->n    ;
	Entorno*        e    = ((struct parametrosElectrodo *) params)->e    ;
	vector<string>*     speciePosition = ((struct parametrosElectrodo *) params)->speciePosition ;
	map<string,Specie>* species        = ((struct parametrosElectrodo *) params)->species        ;

	// valores iniciales incognitas
	Malla1D<long double> X_Specie(speciePosition->size()+1);
	long double X_Ey; 
	long double X_Phi;
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
		    X_Specie[iPos]=gsl_vector_get (vectorIncognitas, iPos);
	  
	}	
	X_Ey  = gsl_vector_get (vectorIncognitas, speciePosition->size()  );
	X_Phi = gsl_vector_get (vectorIncognitas, speciePosition->size()+1);
	// Funciones objetivo
	Malla1D<long double> F_Specie(speciePosition->size()+1);
	long double F_Ey;
	long double F_Phi;

	// Derivadas
	Malla1D<long double> derivada_X(speciePosition->size());
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	    string specieName = (*speciePosition)[iPos];
	    derivada_X[iPos]=( -11.0 * X_Specie[iPos]
	                       +18.0 * (*species)[specieName].m_k2[i][j+1]
	                       - 9.0 * (*species)[specieName].m_k2[i][j+2]
	                       + 2.0 * (*species)[specieName].m_k2[i][j+3]
			     ) / (6.0*e->h);
	}	

	//Optimizacion de calculos
	// H+ tiene un caso especial. Necesito saber su posicion en el vector 
	int iPosH=0;
	for(iPosH=0;iPosH<(int)speciePosition->size() && (*speciePosition)[iPosH]!="H+";iPosH++);
	
	long double I_1_exp = e->opt_F_divided_By_2_R_T * (X_Phi + e->E0eq_I ) ;
	long double I_1 = e->I_0I  * (                      exp (-I_1_exp) -  (X_Specie[iPosH] / (*species)["H+"].C0)  * exp (I_1_exp)) ;
	
	long double Ey_By_n = - (X_Ey * e->normalNodoX[i][j]);

	// Iguala el flujo a cero para todas las especies
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	  string specieName = (*speciePosition)[iPos];  
	  F_Specie[iPos]= -   (*species)[specieName].D[i][j]         * ( derivada_X[iPos] * e->normalNodoX[i][j])
			  -   (*species)[specieName].opt_z_z_u[i][j] *  X_Specie[iPos]    * Ey_By_n;
	}	

	// En el caso de H+, le agrego al flujo la corriente
	F_Specie[iPosH]+= - I_1 / e->F; 


	// Electroneutralidad
	F_Ey  =	0.0;
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	  string specieName = (*speciePosition)[iPos]; 
	  F_Ey+=X_Specie[iPos]* (*species)[specieName].z;
	}
	
	// Corriente parcial == corriente total
	F_Phi  = e->I - I_1; //- I_2;
	
	// Guarda las ecuaciones dentro del vector de ecuaciones
	// Guarda las ecuaciones dentro del vector de ecuaciones
	for(int iPos=0;iPos<(int)speciePosition->size();iPos++){
	    gsl_vector_set (vectorFunciones, iPos, F_Specie[iPos]);
	}
	
	gsl_vector_set (vectorFunciones, speciePosition->size()  , F_Ey   );
	gsl_vector_set (vectorFunciones, speciePosition->size()+1, F_Phi  );

	return GSL_SUCCESS;
}

void calculaBordeAnodoSinCorriente (int i, int j, int n, vector<string>& speciePosition, map<string,Specie>& species){
	// valores iniciales incognitas 
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
	  	species[speciePosition[iPos]].m_k2[i][j] =  (   4.0 * species[speciePosition[iPos]].m_k2[i][j+1]
	  	                                              -       species[speciePosition[iPos]].m_k2[i][j+2]
	  	                                            ) / 3.0;  
	}
}

void calculaBordeCatodoSinCorriente (int i, int j, int n, vector<string>& speciePosition, map<string,Specie>& species){

	// valores iniciales incognitas 
	for(int iPos=0;iPos<(int)speciePosition.size();iPos++){
	  	species[speciePosition[iPos]].m_k2[i][j] =  (   4.0 * species[speciePosition[iPos]].m_k2[i][j-1]
	  	                                              -       species[speciePosition[iPos]].m_k2[i][j-2]
	  	                                            ) / 3.0;  
	}
}

void calculaBordesSubIteracion(Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E){

	for(int j=1;j<e.jj;j++){
		calculaBordeSuperiorSubIteracion(j, species, Phi, E);

		calculaBordeInferiorSubIteracion(j, e, species, Phi, E);
	}
	
	for(int i=1;i<e.ii;i++){
		calculaBordeIzquierdoSubIteracion(i, species, Phi, E);

		calculaBordeDerechoSubIteracion  (i, e, species, Phi, E);
	}
	

	//Esquinas, en sentido horario
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2[0][0]       = (iMapSpec->second.m_k2[0][1]         + iMapSpec->second.m_k2[1][0]        ) / 2;
	    iMapSpec->second.m_k2[0][e.jj]    = (iMapSpec->second.m_k2[1][e.jj]      + iMapSpec->second.m_k2[0][e.jj-1]   ) / 2;
	    iMapSpec->second.m_k2[e.ii][e.jj] = (iMapSpec->second.m_k2[e.ii-1][e.jj] + iMapSpec->second.m_k2[e.ii][e.jj-1]) / 2;
	    iMapSpec->second.m_k2[e.ii][0]    = (iMapSpec->second.m_k2[e.ii-1][0]    + iMapSpec->second.m_k2[e.ii][1]     ) / 2;
	}   

	Phi.m_k2[0][0] = (Phi.m_k2[0][1] + Phi.m_k2[1][0]) / 2;
	Phi.m_k2[0][e.jj] = (Phi.m_k2[1][e.jj] + Phi.m_k2[0][e.jj-1]) / 2;
	Phi.m_k2[e.ii][e.jj] = (Phi.m_k2[e.ii-1][e.jj] + Phi.m_k2[e.ii][e.jj-1]) / 2;
	Phi.m_k2[e.ii][0] = (Phi.m_k2[e.ii-1][0] + Phi.m_k2[e.ii][1]) / 2;
}


void calculaBordeSuperiorSubIteracion(int j, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2[0][j] = iMapSpec->second.m_k2[1][j];
	}   
	//Necesito tener un valor de Phi en este borde porque luego lo voy a usar para propagar por la malla.
	Phi.m_k2[0][j] =      Phi.m_k2[1][j];
}


void calculaBordeInferiorSubIteracion(int j, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2[e.ii][j] = iMapSpec->second.m_k2[e.ii-1][j];
	}   
	Phi.m_k2[e.ii][j] =      Phi.m_k2[e.ii-1][j];
}

void calculaBordeIzquierdoSubIteracion(int i, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2[i][0] = iMapSpec->second.m_k2[i][1];
	}   
	
	// Necesito tener un valor de Phi en este borde porque luego lo voy a usar para propagar por la malla.
	Phi.m_k2[i][0] =      Phi.m_k2[i][1];
}

void calculaBordeDerechoSubIteracion(int i, Entorno& e, map<string,Specie>& species, Potential& Phi, ElectricField& E){
	// species
	for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	    // Specie
	    iMapSpec->second.m_k2[i][e.jj] = iMapSpec->second.m_k2[i][e.jj-1];
	}   
        Phi.m_k2[i][e.jj] =      Phi.m_k2[i][e.jj-1];
}

void actualizacionElectricField(const Entorno& e,const Potential& Phi,ElectricField& E){
		// Calculo Ex en el dominio y bordes laterales
		for(int i=1;i<e.ii;i++){
			for(int j=0;j<=e.jj;j++){
				if(e.geometriaTipoNodo [i][j]==dominio) E.mx_k2[i][j] = (Phi.m_k2[i-1][j] - Phi.m_k2[i+1][j]) / (2 * e.h);
			}
		}


		// Calculo en el domino Ey
		for(int i=0;i<=e.ii;i++){
			for(int j=1;j<e.jj;j++){
				if(e.geometriaTipoNodo [i][j]==dominio) E.my_k2[i][j] =  (Phi.m_k2[i][j-1] - Phi.m_k2[i][j+1]) / (2 * e.h);
			}
		}
		//Seteo los bordes inferior/superior de Ex/Ey con el valor en (1,j)
		for(int j=0;j<=e.jj;j++){
			E.mx_k2[0][j]    = E.mx_k2[1][j];
			E.mx_k2[e.ii][j] = E.mx_k2[e.ii-1][j];
			E.my_k2[0][j]    = E.my_k2[1][j];
			E.my_k2[e.ii][j] = E.my_k2[e.ii-1][j];
		}
}

void actualizacionOptimizacionCalculos(const Entorno& e,  map<string,Specie>& species){
  // Optimizacion de calculos
  for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
  	// Specie
    if(e.ctrl_faseActual==4) iMapSpec->second.optimizarParaCalculosDominioSobreK1FickReactionEquation(e);
    else                     iMapSpec->second.optimizarParaCalculosDominioSobreK1NernstPlankReactionEquation(e);
  }
}

void actualizacionCorrienteVoltaje(const int n, Entorno& e, const long double diferenciaPotencial){
	/*status*/ 
	if(n % 1000 == 0) cout << "Fase " << e.ctrl_faseActual << " - n = " << n << " - diferenciaPotencial = " << diferenciaPotencial << endl;
	
	// Accion a tomar en cada e.ctrl_faseActual
	switch(e.ctrl_faseActual){
	  case 1:{
		    e.I = e.I + 5;
		    break;
	  }

	  case 2:
		{
	    if (diferenciaPotencial < e.umbralVoltageMax) e.I = e.I + 2.5;
	    if (diferenciaPotencial > e.umbralVoltageMax) e.I = e.I - 2.5;
	    break;
		}

    case 3:
		{
	    e.I = e.I - 5;
	    break;
		}

    case 4:
		{
			// nop
	    break;
		}
	}
}

void actualizacionFase(const int n, Entorno& e, map<string,Specie>& species, const long double diferenciaPotencial){
      // Cambio de e.ctrl_faseActual?
      // e.ctrl_faseActual 1 -> e.ctrl_faseActual 2
      if((e.ctrl_faseActual==1) && (diferenciaPotencial>=e.umbralVoltageMax)){
	      e.ctrl_faseActual = 2;
	      // Lo supongo estable -> acelero el tiempo
	      e.dt = 1e-4;
	      e.tiempoFoto = 5e-3;
	      // reoptimizar x cambio de dt
	      reoptimizacionTemporal(e,species);
	      /*status*/ cout << "PASE FASE " << e.ctrl_faseActual << "  ==> tiempo: " << e.ctrl_deltaTAcumulado << "  - diferenciaPotencial: " << diferenciaPotencial << " - corriente: " << e.I << " - n: " << n << endl;
      }

      // e.ctrl_faseActual 2 -> e.ctrl_faseActual 3
      if((e.ctrl_faseActual==2) && (e.ctrl_tiempoActualPulso>=e.duracionPulso)){
	      e.ctrl_faseActual = 3;
	      // Se vuelve mas inestable -> bajo el tiempo
	      e.dt = 1e-12;
	      e.tiempoFoto = 1e-2;
	      // Reseteo la cuenta del pulso
	      e.ctrl_tiempoActualPulso=0.0;
	      // reoptimizar x cambio de dt
	      reoptimizacionTemporal(e,species);
	      /*status*/ cout << "PASE A FASE " << e.ctrl_faseActual << "  ==> tiempo: " << e.ctrl_deltaTAcumulado << "  - diferenciaPotencial: " << diferenciaPotencial << " - corriente: " << e.I << " - n: " << n << endl;
      }

      // e.ctrl_faseActual 3 -> e.ctrl_faseActual 4
      if((e.ctrl_faseActual==3) && (e.I<=e.umbralCorrienteMin)){
	      e.ctrl_faseActual = 4;
	      // Lo supongo estable -> acelero el tiempo
	      e.dt = 1e-5;
	      e.I = 12.0;
	      e.tiempoFoto = 2e-2;
	      // reoptimizar x cambio de dt
	      reoptimizacionTemporal(e,species);
	      /*status*/ cout << "PASE A FASE " << e.ctrl_faseActual << "  ==> tiempo: " << e.ctrl_deltaTAcumulado << "  - diferenciaPotencial: " << diferenciaPotencial << " - corriente: " << e.I << " - n: " << n << endl;
      }

      // e.ctrl_faseActual 4 -> e.ctrl_faseActual 1
      if((e.ctrl_faseActual==4) && (e.ctrl_tiempoActualPulso>(e.duracionCicloECT - e.duracionPulso))){
	      e.ctrl_faseActual = 1;
	      // Se vuelve mas inestable -> bajo el tiempo
	      e.dt = 1e-12;
	      // Reseteo la cuenta del pulso
	      e.ctrl_tiempoActualPulso=0.0;
	      // reoptimizar x cambio de dt
	      reoptimizacionTemporal(e,species);
	      /*status*/ cout << "PASE A FASE " << e.ctrl_faseActual << "  ==> tiempo: " << e.ctrl_deltaTAcumulado << "  - diferenciaPotencial: " << diferenciaPotencial << " - corriente: " << e.I << " - n: " << n << endl;
      }


}

void reoptimizacionTemporal(Entorno& e, map<string,Specie>& species){
    e.reoptimizar();
    // species
    for (std::map<string,Specie>::iterator iMapSpec=species.begin(); iMapSpec!=species.end(); ++iMapSpec){
	// Specie
	iMapSpec->second.reoptimizar(e);
    }   

}
