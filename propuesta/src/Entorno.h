#ifndef ENTORNO_H
#define ENTORNO_H

#include <iostream>
#include <math.h> 
#include "Malla1D.h"
#include "Malla2D.h"
#include "Tipos.h"

using namespace std;

class Entorno {

    public:

        Entorno();
        Entorno(long xCantNodos, long yCantNodos, int xSeparacionBordeLateral, int xSeparacionBordeSup);

        // Constantes Universales
    	long double	T;	// K
    	long double	F;	// A.s / mol
    	long double	R;	// Kg.m2 / K.mol.s2

        // Malla Espacial
        int ii;				// ultimo nodo en el eje x
        int jj;				// ultimo nodo en el eje y

    	int cantAnodos;			// cantidad de anodos en la malla
    	int cantCatodos;		// cantidad de catodos en la malla
    	int separacionBordeLateral;	// separacion Borde Lateral de la barra de electrodo
    	int separacionBordeSup;		// separacion Borde Superior de la barra de electrodo

    	Malla1D<long double>   x;	// m	posición en el eje x
    	Malla1D<long double>   y;	// m	posicion en el eje y

    	long double h;			// m	factor de discretizacion

    	Malla2D<TipoNodo> geometriaTipoNodo;
    	Malla2D<long double> normalNodoX;
    	long double minimaDistanciaAnodoCatodo;

        // Malla Temporal
        long n0;	// Numero de iteracion temporal por la que empezara a iterar
        long nn;	// Cantidad de iteraciones temporales
        long double dt;	// seg. Cuanto es un una iteracion temporal en segundos

        // Corriente
        long double I;		// A / m2
        long double I_0I;	// A / m2   constante I0 correspondiente a H+  anodo
        long double I_0II;	// A / m2   constante I0 correspondiente a Cl- anodo
        long double I_0III;	// A / m2   constante I0 correspondiente a OH- catodo
        long double E0eq_I;	// V        constante E0 correspondiente a H+  anodo
        long double E0eq_II;	// V        constante E0 correspondiente a Cl- anodo
        long double E0eq_III;	// V        constante E0 correspondiente a OH- catodo

        // Concentracion inicial de la Solucion (H2O)
        long double H2O_C0;
	
        // Termino Reactivo
        long double k_wf;	// m3 / mol.s  constante correspondiente al termino reactivo de H2O forward
        long double k_wb;	// s-1         constante correspondiente al termino reactivo de H2O backward

        long double k_L_LH_f;
        long double k_L_LH_b;
    
        long double k_L_LOH_f;
        long double k_L_LOH_b;
    
        long double k_LH_LHOH_f;
        long double k_LH_LHOH_b;

        long double k_LOH_LHOH_f;
        long double k_LOH_LHOH_b;
        
        // Parametro de corrida

        long   cantidadItConvergencia;		// Cantidad de Iteraciones aceptables para convergencia.
        long double errorConvergenciaInicial;	// Error aceptable para convergencia para la condicion inicial.

        long double wtPhi_Inicial;

        Malla2D<long double> wt_Ex;		// Parametro de relajacion SOR.
        Malla2D<long double> wt_Ey;		// Parametro de relajacion SOR.
        Malla2D<long double> wt_Phi;		// Parametro de relajacion SOR.

        string load;				// Nombre del archivo con los datos de la corrida n
        string commit;				// Nombre del archivo que irá guardando
        string dataFolder;			// Nombre de la carpeta donde se guardan los datos          
        string initialValuesFolder;             // Nombre de la carpeta donde se guardan los datos iniciales

        long double factorMigracionElectrodo;

    	// No se usan porque estamos trabajando a corriente constante
        long double umbralVoltageMax; 
        //long double umbralVoltageMin;
        long double duracionPulso;
        long double duracionCicloECT;

        // variable que llevan el valor actual (las denomino de control)
        long double ctrl_tiempoActualPulso;
        long double ctrl_deltaTAcumulado;
        int ctrl_faseActual;
	
        long double tiempoFoto;
        long double umbralCorrienteMin;
        long double umbralTiempoTotal;
	
        // optimizaciones
        long double opt_F_divided_By_2_R_T;
        long double opt_H2O_By_k_wb_dt;
        long double opt_k_wf_dt;

        void reoptimizar();

    private:
        void spaceGenerator();
        void seteoDeRelajacion();
};

// Salida
ostream& operator <<(ostream& o, const Entorno& s);

#endif // ENTORNO_H