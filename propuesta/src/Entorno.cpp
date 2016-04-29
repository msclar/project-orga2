#include "Entorno.h"
#include "math.h"
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

Entorno::Entorno(){
}

Entorno::Entorno(long xCantNodos, long yCantNodos, int xSeparacionBordeLateral, int xSeparacionBordeSup)
                 :  
                 separacionBordeLateral(xSeparacionBordeLateral),
                 separacionBordeSup(xSeparacionBordeSup),
                 x(xCantNodos),
                 y(yCantNodos),
                 geometriaTipoNodo(xCantNodos,yCantNodos),
                 normalNodoX(xCantNodos,yCantNodos)
{
    // Constantes Universales
    T = 298;          // K
    F = 96485.3415;   // A.s / mol
    R = 8.31;         // kg.m2 / K.mol.s2

    // Malla Espacial
    ii = xCantNodos-1 ; //      Cantidad de nodos en el eje x
    jj = yCantNodos-1 ; //      Cantidad de nodos en el eje y
    h = 2.5e-5; // m

    x[0] = 0.000;
    for(int i=1;i<=ii;i++){
       x[i]= x[i-1] + h;
    }
    
    y[0] = 0.000;
    for(int j=1;j<=jj;j++){
       y[j]= y[j-1] + h;
    }

    spaceGenerator();

    // Malla Temporal
    n0 = 1		 ;//     numero de iteracion temporal por la que empezara a iterar
    nn = 100000000;//		 Cantidad de iteraciones temporales
    dt  = 1e-12; //0.0005;//      seg. Cuanto es un una iteracion temporal en segundos

    // Corriente
    I= 0.250; //1494.0; //30.0;//1000.0;
    //I=0.002/0.0000273197;  // A/m2  -->  2mA
    //I=0.004/0.0000086063;  // A/m2  -->  4mA
    //I=0.008/0.0000011742; 10.0 // A/m2  -->  8mA
    //I=0.010/0.000000051819;// A/m2  --> 10mA

    //I_0I         = 1e-006;// A / m2
    I_0I         = 1e-006;// A / m2
    I_0II        = 10.0  ;// A / m2
    I_0III       = 1.0   ;// A / m2

    E0eq_I       = 0.816   ;// V
    E0eq_II      = 1.407   ;// V
    E0eq_III     = -0.828 ;// V

    // Concentracion de la solucion (H2O)
    H2O_C0 = 55500.0;
    
    // Termino Reactivo
    k_wf  = 1.5e8;        // m3 / mol.s
    k_wb  = 2.7e-5;       // s-1

    // Parametros de corrida
    errorConvergenciaInicial = 1e-15; 	     // Error aceptable para convergencia para condicion inicial
    
    cantidadItConvergencia =20000;	     // Cantidad de Iteraciones aceptables para convergencia.

    // fileName
    load	        = "2D_segundos";
    commit	        = "2D_segundos";
    dataFolder          = "data/";
    initialValuesFolder = "data/initial/";
    
    factorMigracionElectrodo = 1.0 / 1.0; // Deprecated

    // variable que llevan el valor actual (los inicializo en cero)
    ctrl_tiempoActualPulso = 0;
    ctrl_deltaTAcumulado = 0;
    // Fase 1:   Pulso - Rampa ascendente para llegar al diferenciaPotencial
    // Fase 2:   Pulso - Mantener el diferenciaPotencial
    // Fase 3: NoPulso - Rampa descendente para bajar el diferenciaPotencial
    // Fase 4: NoPulso - Mantener diferenciaPotencial
    ctrl_faseActual = 1;
    
    umbralCorrienteMin = 1;
    umbralTiempoTotal = 5.00; // tiempo total de simulacion (seg)
    tiempoFoto = 1.0e-10; //va cambiando segun la fase
    
    // Estos son parametros para corriente pulsatil
    //umbralVoltageMin = 0.1;
    umbralVoltageMax = h*100*(xCantNodos-2*xSeparacionBordeLateral)*40;  // 40 V/cm;
    duracionPulso    = 1e-3;
    duracionCicloECT = 1.0;

    // optimizaciones
    opt_F_divided_By_2_R_T = F / (2.0 * R * T);
    reoptimizar();

    // guarda en un archivo los valores de la malla
    x.save(initialValuesFolder+commit+ "_n_x.dat");
    y.save(initialValuesFolder+commit+ "_n_y.dat");
    geometriaTipoNodo.save(initialValuesFolder+"tipoNodo.dat");
    normalNodoX.save(initialValuesFolder+"normalNodoX.dat");
}

void Entorno::reoptimizar(){
    opt_H2O_By_k_wb_dt = H2O_C0 * k_wb * dt;
    opt_k_wf_dt = k_wf * dt;
}

ostream& operator <<(ostream& o, const Entorno& e){

	o << " Constantes Universales" << endl;
	o << " ----------------------" << endl;
	o << " T :" << e.T << " K"                << endl;
	o << " F :" << e.F << " A.s / mol"        << endl;
	o << " R :" << e.R << " Kg.m2 / K.mol.s2" << endl;
	o << endl << endl;

	o << " Malla Espacial" << endl;
	o << " -------------- " << endl;
	o << "  0 ... ii   : " << e.ii << endl;
	o << "  0 ... jj   : " << e.jj << endl;

	o << " Malla Temporal" << endl;
	o << " --------------" << endl;
	o << " n0        : " << e.n0 << endl;
	o << " nn        : " << e.nn << endl;
	o << " dt        : " << e.dt << " s" << endl;
	o << endl << endl;

	o << " Corriente" << endl;
	o << " ---------" << endl;
	o << " I_0I    : " << e.I_0I     << " A/m2 " << endl;
	o << " I_0II   : " << e.I_0II    << " A/m2 " << endl;
	o << " I_0III  : " << e.I_0III   << " A/m2 " << endl;
	o << " E0eq_I  : " << e.E0eq_I   << " V "    << endl;
	o << " E0eq_II : " << e.E0eq_II  << " V "    << endl;
        o << " E0eq_III: " << e.E0eq_III << " V "    << endl;
	//o << endl << endl;+ I_2 

	o << " Termino Reactivo" << endl;
	o << " ----------------" << endl;
	o << " k_wf    : " << e.k_wf << " m3 / mol.s " << endl;
	o << " k_wb    : " << e.k_wb << "	s-1 "        << endl;
	o << endl << endl;

	o << " Parametros de corrida" << endl;
	o << " ---------------------" << endl;
	o << " cantidadItConvergencia  : " << e.cantidadItConvergencia << endl;
	o << " load   : " << e.load   << endl;
	o << " commit : " << e.commit << endl;

	o << endl << " vector I " << endl;
	o         << " ........." << endl;
	o << e.I << endl;

	o << endl << " vector x " << endl;
	o         << " ........." << endl;
	o << e.x << endl;


	o << "....................................................................................................................."  << endl;

    return o;
}

void Entorno::spaceGenerator(){
    
//    int separacionBordeLateral=25; // en nodos
//    int separacionBordeSup=10; // en nodos
    
    // Al principio defino a toda la matriz como dominio
    geometriaTipoNodo.sameValue(dominio);

    for(int i=separacionBordeSup;i<=ii-separacionBordeSup;i++){
       normalNodoX[i][0 + separacionBordeLateral] =  1.0;
       normalNodoX[i][jj-separacionBordeLateral] =  1.0;
    }

    for(int i=separacionBordeSup;i<=ii-separacionBordeSup;i++){
       geometriaTipoNodo [i][ 0+separacionBordeLateral] = bordeAnodo;
       geometriaTipoNodo [i][jj-separacionBordeLateral] = bordeCatodo;
    }

    // Electrodos redondos
    minimaDistanciaAnodoCatodo= h * abs(jj-separacionBordeLateral-separacionBordeLateral);

    cout << "Distancia entre electrodos: " << minimaDistanciaAnodoCatodo << " m" << endl;
}