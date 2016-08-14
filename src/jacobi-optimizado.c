#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define EPS 0.000001

long long tt[16];
struct timespec startt[16], endt[16];
int measure_type = -1;

#define STARTMEASURE(x)\
if (measure_type != x)\
{\
	clock_gettime(CLOCK_MONOTONIC,&(startt[x]));\
}
	
#define ENDMEASURE(x)\
if (measure_type != x)\
{\
	clock_gettime(CLOCK_MONOTONIC, &(endt[x]));\
	tt[x] += endt[x].tv_sec*1000000000LL + endt[x].tv_nsec - startt[x].tv_sec*1000000000LL - startt[x].tv_nsec;\
}

#ifdef ASM
	extern void jacobiStep (double*, double*, double*, double*, int, int);
	extern void laplaceStep (double*, double*, int, int);
	extern void updateB (double*, double*, double*, double, double, int);
	extern void fillWithConstant (double*, double, int);
	extern void calculateVectorError (double*, double*, double*, double*, int, int);
	extern void createA (double*, double*, double, double, double, int, int);
	extern void calculateTInd (double*, double, double, double*, double*, int, int, double);
	extern double vectorNorm (double*, int);
	extern double distanceBetweenVectors (double*, double*, int);
#endif

int indice (int i, int j, int max_j) {
	return i * max_j + j;
}

void print1DColumna(double matrix[], int cant_elems) {
	int i;
	for (i = 0; i < cant_elems; i++) {
		printf("%f\n", matrix[i]);
	}
	printf("\n");
}

void printDbg2DMatlab(double matrix[], int max_i, int max_j) {
	int i, j;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			if (j) printf(" ");
			printf("%f", matrix[i * max_j + j]);
		}
		printf("\n");
	}
	printf("\n\n");	
}

void print2DMatlab(double matrix[], int max_i, int max_j, FILE* out) {
	int i, j;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			if (j) fprintf(out, " ");
			fprintf(out, "%f", matrix[i * max_j + j]);
		}
		fprintf(out, "\n");
	}
	fprintf(out, "\n\n");	
}

double normaVector(double matriz[], int cant_elems) {
	double result = 0.0;
	
	int i;
	for (i = 0; i < cant_elems; i++) {
		result += matriz[i] * matriz[i];
	}
	// se puede evitar tomar raiz cuadrada si tomamos la convencion
	return sqrt(result);
}

double promedioVecinos (int i, int j, double* phi, int max_j) {
	return (phi[indice(i-1, j, max_j)] + phi[indice(i+1, j, max_j)] + phi[indice(i, j-1, max_j)] + phi[indice(i, j+1, max_j)]) / 4.0;
}

void pasoLaplace(double* res, double* phi, int max_i, int max_j) {
	
	#ifdef MEASURE_TIME
		STARTMEASURE(3);
	#endif
	
	#ifdef ASM
		laplaceStep(res, phi, max_i, max_j);
	#endif
	
	#ifndef ASM
	int i, j;
	// promedio los 4 vecinos (estos tienen los 4 vecinos)
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			res[indice(i , j, max_j)] = promedioVecinos(i, j, phi, max_j);
		}
	}
	// bordes izquierdo y derecho
	for (i = 1; i < max_i - 1; i++) {
		res[indice(i, 0, max_j)] = res[indice(i, 1, max_j)];
		res[indice(i, max_j-1, max_j)] = res[indice(i, max_j-2, max_j)]; 
	}
	
	// bordes inferior y superior
	for (j = 0; j < max_j; j++) {
		res[indice(0, j, max_j)] = res[indice(1, j, max_j)];
		res[indice(max_i-1, j, max_j)] = res[indice(max_i-2, j, max_j)]; 
	}
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(3);
	#endif
	
	
	return;
}

double calcErrorLaplace(double* actPhi, double* nextPhi, int max_i, int max_j) {
	double norma = 0.0;
	
	#ifdef MEASURE_TIME
		STARTMEASURE(4);
	#endif
	
	#ifdef ASM
		norma = distanceBetweenVectors(actPhi, nextPhi, max_i * max_j);
	#endif
	
	#ifndef ASM
	int i;
	for (i = 0; i < max_i * max_j; i++) {
		norma += (actPhi[i] - nextPhi[i]) * (actPhi[i] - nextPhi[i]);
	}
	norma = sqrt(norma);
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(4);
	#endif
	
	return norma;
}

void obtenerLaplace(double* res, 
					int anodoi, 
					int anodoj, 
					double anodov, 
					int catodoi, 
					int catodoj, 
					double catodov, 
					int max_i, 
					int max_j) {
	
	#ifdef MEASURE_TIME
		STARTMEASURE(5);
	#endif
						
	#ifndef ASM
	int i;
	for (i = 0; i < max_i * max_j; i++) {
		res[i] = (anodov + catodov) / 2.0;
	}
	#endif
	
	#ifdef ASM
		fillWithConstant(res, (anodov + catodov) / 2.0, max_i * max_j);
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(5);
	#endif
	
	int posanodo = indice(anodoi, anodoj, max_j);
	int poscatodo = indice(catodoi, catodoj, max_j);
	res[posanodo] = anodov;
	res[poscatodo] = catodov;
	
	double *aux = malloc(max_i * max_j * sizeof(double));
	while (calcErrorLaplace(aux, res, max_i, max_j) > EPS) {
		pasoLaplace((double*) aux, res, max_i, max_j);
		aux[posanodo] = anodov;
		aux[poscatodo] = catodov;
		
		pasoLaplace(res, (double*) aux, max_i, max_j);
		res[posanodo] = anodov;
		res[poscatodo] = catodov;
	}
	free(aux);
	return;
}

double calcularPosicion(int i,
						int j,
						double A[],
						double Tn[],
						double B[],
						int max_i,
						int max_j) {
	int s = indice(i, j, max_j);
	double res = -B[s];
	res += A[indice(0, s, max_i*max_j)] * Tn[indice(i - 1, j, max_j)];
	res += A[indice(1, s, max_i*max_j)] * Tn[indice(i, j - 1, max_j)];
	res += A[indice(2, s, max_i*max_j)] * Tn[indice(i + 1, j, max_j)];
	res += A[indice(3, s, max_i*max_j)] * Tn[indice(i, j + 1, max_j)];
	res += A[indice(4, s, max_i*max_j)] * Tn[indice(i, j, max_j)];
	return res;
}

double calcVectorError(double A[], 
					   double Tn[], 
					   double B[],
					   double k[],
					   double res[],
					   int anodo_x,
					   int anodo_y,
					   int catodo_x,
					   int catodo_y,
					   int max_i,
					   int max_j,
					   double delta_x) {
						   
	#ifdef MEASURE_TIME
		STARTMEASURE(6);
	#endif
	
	#ifdef ASM
		calculateVectorError(A, Tn, B, res, max_i, max_j);
	#endif

	#ifndef ASM
	int i, j;
	// calculo A * Tn - B
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			res[s] = calcularPosicion(i, j, A, Tn, B, max_i, max_j);
		}
	}
	for (i = 0; i < max_i; i++) {
		res[indice(i, 0, max_j)] = 0.0;
		res[indice(i, max_j-1, max_j)] = 0.0;
	}
	
	for (j = 0; j < max_j; j++) {
		res[indice(0, j, max_j)] = 0.0;
		res[indice(max_i-1, j, max_j)] = 0.0;
	}
	#endif	
			   
	#ifdef MEASURE_TIME
		ENDMEASURE(6);
	#endif
	
	double r = k[indice(catodo_x, catodo_y, max_j)] / (delta_x * 10);
	res[indice(catodo_x, catodo_y, max_j)] = 
	r * (
			res[indice(catodo_x+1, catodo_y, max_j)] + 
			res[indice(catodo_x-1, catodo_y, max_j)] + 
			res[indice(catodo_x, catodo_y+1, max_j)] + 
			res[indice(catodo_x, catodo_y-1, max_j)]
		) / (4 * (r - 1));
	r = k[indice(anodo_x, anodo_y, max_j)] / (delta_x * 10);
	res[indice(anodo_x, anodo_y, max_j)] =  
	r * (
			res[indice(anodo_x+1, anodo_y, max_j)] + 
			res[indice(anodo_x-1, anodo_y, max_j)] + 
			res[indice(anodo_x, anodo_y+1, max_j)] + 
			res[indice(anodo_x, anodo_y-1, max_j)]
		) / (4 * (r - 1));

	double norma;
		   
	#ifdef MEASURE_TIME
		STARTMEASURE(7);
	#endif
	
	#ifndef ASM
		norma = normaVector(res, max_i * max_j);
	#endif
	
	#ifdef ASM
		norma = vectorNorm(res, max_i * max_j);
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(7);
	#endif
	return norma;
}

/* CONVENCION para los 5 vecinos
 *   0
 * 1 4 3
 *   2
 */

double calcularSiguienteJacobi (int i,
                                int j,
                                double A[],
                                double Tn[],
                                double B[],
                                int max_i,
                                int max_j)  {
	int s = indice(i, j, max_j);
	double sum = 0.0;
	sum += A[indice(0, s, max_i * max_j)] * Tn[indice(i-1, j, max_j)];
	sum += A[indice(2, s, max_i * max_j)] * Tn[indice(i+1, j, max_j)];
	sum += A[indice(1, s, max_i * max_j)] * Tn[indice(i, j-1, max_j)];
	sum += A[indice(3, s, max_i * max_j)] * Tn[indice(i, j+1, max_j)];
	return (B[s] - sum) / A[indice(4, s, max_i * max_j)];
}

void jacobiStepOptimized(double* Tn_sig, 
				double* Tn,
				double* B,
				double* TInd, 
				double A[],
				double k[],
				int catodo_x, 
				int catodo_y, 
				int anodo_x, 
				int anodo_y,
				int max_i,
				int max_j,
				double delta_x,
				double T_aire) {
	// X_i^(iter+1) = (b_i - sum (j != i) a_ij * X_j^iter	
	#ifdef MEASURE_TIME
		STARTMEASURE(8);
	#endif
	
	#ifndef ASM
	
	int i;
	int j;
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			Tn_sig[indice(i, j, max_j)] = calcularSiguienteJacobi(i, j, A, Tn, B, max_i, max_j);
		}
	}

	// bordes izquierdo y derecho	
	for (i = 0; i < max_i; i++) {
		Tn_sig[indice(i, 0, max_j)] = Tn_sig[indice(i, 1, max_j)];
		Tn_sig[indice(i, max_j-1, max_j)] = Tn_sig[indice(i, max_j-2, max_j)];
	}

	// bordes inferior y superior
	for (j = 1; j < max_j - 1; j++) {
		Tn_sig[indice(0, j, max_j)] = Tn_sig[indice(1, j, max_j)];
		Tn_sig[indice(max_i-1, j, max_j)] = Tn_sig[indice(max_i-2, j, max_j)];
	}
	#endif
	
	#ifdef ASM
		jacobiStep(Tn_sig, Tn, B, A, max_i, max_j);
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(8);
	#endif
	
	// temperatura en electrodos
	double r = k[indice(catodo_x, catodo_y, max_j)] / (delta_x * 10); //  h = 10 W / m² K
	Tn_sig[indice(catodo_x, catodo_y, max_j)] = 
		r * (
			Tn_sig[indice(catodo_x+1, catodo_y, max_j)] + 
			Tn_sig[indice(catodo_x-1, catodo_y, max_j)] + 
			Tn_sig[indice(catodo_x, catodo_y+1, max_j)] + 
			Tn_sig[indice(catodo_x, catodo_y-1, max_j)]
		) / (4 * (r - 1))
		- T_aire / (r - 1); 
	
	r = k[indice(anodo_x, anodo_y, max_j)] / (delta_x * 10);
	Tn_sig[indice(anodo_x, anodo_y, max_j)] = 
		r * (
			Tn_sig[indice(anodo_x+1, anodo_y, max_j)] + 
			Tn_sig[indice(anodo_x-1, anodo_y, max_j)] + 
			Tn_sig[indice(anodo_x, anodo_y+1, max_j)] + 
			Tn_sig[indice(anodo_x, anodo_y-1, max_j)]
		) / (4 * (r - 1))
		- T_aire / (r - 1); 
}

#ifndef ASM
void calculateTInd(double* phi, double delta_x, double delta_y, 
				   double* sigma, double* TInd, int max_i, int max_j,
				   double resto) {
	int i, j;
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {			
			double deriv_phi_x = (phi[indice(i+1, j, max_j)] - phi[indice(i-1, j, max_j)]) / (2 * delta_x);
			double deriv_phi_y = (phi[indice(i, j+1, max_j)] - phi[indice(i, j-1, max_j)]) / (2 * delta_y);
			double gradient_phi_quad = deriv_phi_x * deriv_phi_x + deriv_phi_y * deriv_phi_y;
			TInd[indice(i, j, max_j)] = resto + sigma[indice(i, j, max_j)] * gradient_phi_quad;
		}
	}
}
#endif

void updateB_C (double* B, double* Tn, double* TIndAct, double rho_times_C_rho, double delta_t, int max_ij) {
	int i;
	for (i = 0; i < max_ij; i++) {
		B[i] = - rho_times_C_rho * Tn[i] / delta_t - TIndAct[i];
	}
}

void createA_C (double* A, double* k, double indep_term, double delta_x, double delta_y, int max_i, int max_j) {
	// indep_term = - w_b * C_b * rho_b - rho * C_rho / delta_t
	// ya comenzamos con el assembler: ver createA.asm

	int i, j;
	// creacion de la matriz A
	for (i = 0; i < max_i * max_j * 5; i++) {
		A[i] = 0;
	}
	
	// calculo de A
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			
			double derivada_k_x = (k[indice(i+1, j, max_j)] - k[indice(i-1, j, max_j)]) / (2 * delta_x);
			A[indice(0, s, max_i*max_j)] = - derivada_k_x / (2 * delta_x) + k[indice(i-1, j, max_j)] / (delta_x * delta_x);
			A[indice(2, s, max_i*max_j)] = derivada_k_x / (2 * delta_x) + k[indice(i+1, j, max_j)] / (delta_x * delta_x);

			double derivada_k_y = (k[indice(i, j+1, max_j)] - k[indice(i, j-1, max_j)]) / (2 * delta_y);
			A[indice(1, s, max_i*max_j)] = - derivada_k_y / (2 * delta_y) + k[indice(i, j-1, max_j)] / (delta_y * delta_y);
			A[indice(3, s, max_i*max_j)] = derivada_k_y / (2 * delta_y) + k[indice(i, j+1, max_j)] / (delta_y * delta_y);

			A[indice(4, s, max_i*max_j)] = - 2 * k[s] / (delta_x * delta_x) - 2 * k[s] / (delta_y * delta_y) + indep_term;
		}
	}
}

int main (int argc, char** argv) {
	int max_i, max_j;
	double delta_t, delta_x, delta_y;
	int max_cant_delta_t; // cantidad de tiempos a analizar
	int cant_prints = 1; // cantidad de impresiones parciales a hacer
	
	if (argc == 6) {
		max_i = atoi(argv[1]);
		max_j = atoi(argv[2]);
		delta_t = atof(argv[3]);
		max_cant_delta_t = atoi(argv[4]);
		cant_prints = atoi(argv[5]);
	}
	else if (argc == 1) {
		max_i = 53;
		max_j = 59;
		delta_t = 0.000005;
		max_cant_delta_t = 100000;
	}
	else {
		printf("La cantidad de parametros pasados es incorrecta\n");
		printf("usage: max_i max_j delta_t max_cant_delta_t cant_prints\n");
		return 1;
	}
	
	// inicializacion de variables
	int C_rho = 3680;
	int rho = 1039;
	int q_ddd = 10347;
	int C_b = 3840;
	int rho_b = 1060;
	double w_b = 0.00715;
	double T_a = 310.15;
	double T_aire = 296;

	double catodo_x = 0.025, catodo_y = 0.05, anodo_x = 0.075, anodo_y = 0.05;
	
	delta_x = 0.1 / max_i; delta_y = 0.1 / max_j; // recinto de 0.1cm x 0.1cm

	double *Tn, *Tn_sig, *TInd, *B, *k, *sigma, *phi, *A, *phiZero, *TIndPhiZero;
	
	Tn = malloc(max_i * max_j * sizeof(double));
	TInd = malloc(max_i * max_j * sizeof(double));
	Tn_sig = malloc(max_i * max_j * sizeof(double));
	B = malloc(max_i * max_j * sizeof(double));
	k = malloc(max_i * max_j * sizeof(double));
	sigma = malloc(max_i * max_j * sizeof(double));
	phi = malloc(max_i * max_j * sizeof(double));
	A = malloc(max_i * max_j * 5 * sizeof(double));
	phiZero = malloc(max_i * max_j * sizeof(double));
	TIndPhiZero = malloc(max_i * max_j * sizeof(double));

	int catodo_x_idx = catodo_x / delta_x;
	int catodo_y_idx = catodo_y / delta_y;
	double catodo_v = 0; // voltaje en catodo
	
	int anodo_x_idx = anodo_x / delta_x;
	int anodo_y_idx = anodo_y / delta_y;
	double anodo_v = 1000; // voltaje en anodo
	
	int i, j;
	
	#ifndef MEASURE_TIME
		FILE* phi_file = fopen("phi.out", "w");
		print2DMatlab(phi, max_i, max_j, phi_file);
		fclose(phi_file);
	#endif
	
	// lectura de input
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			int s = indice(i, j, max_j);
			Tn[s] = T_a;
			Tn_sig[s] = T_a;
			TInd[s] = 0;
			k[s] = 0.565;
			sigma[s] = 0.75;
			phiZero[s] = 0;
		}
	}

	#ifdef MEASURE_TIME
		//printf("Ingrese el tipo de medicion\n");
		// scanf("%i", &measure_type);
		// Ingreso el tipo de la medicion, los tipos posibles son los siguientes:
		/*
		0 = TODO
		1 = CREACION DE A
		2 = UPDATE B
		3 = PASO LAPLACE
		4 = ERROR LAPLACE
		5 = FILL WITH CONSTANT LAPLACE
		6 = CALC VECTOR ERROR
		7 = NORM VECTOR
		8 = JACOBI STEP
		9 = CALCULATE TIND
		*/
		for (i = 0; i < 16; i++)
			tt[i] = 0;
		
		STARTMEASURE(0);
	#endif

	obtenerLaplace(phi, anodo_x_idx, anodo_y_idx, anodo_v, catodo_x_idx, catodo_y_idx, catodo_v, max_i, max_j);
	// calculo de TInd
	
	#ifdef MEASURE_TIME
		STARTMEASURE(9);
	
	#endif
	calculateTInd(phi, delta_x, delta_y, sigma, TInd, max_i, max_j, w_b * C_b * rho_b * T_a + q_ddd);
	calculateTInd(phiZero, delta_x, delta_y, sigma, TIndPhiZero, max_i, max_j, w_b * C_b * rho_b * T_a + q_ddd);
	
	#ifdef MEASURE_TIME
		ENDMEASURE(9);
	#endif
	
	#ifdef MEASURE_TIME
		STARTMEASURE(1);
	#endif
	
	#ifdef ASM
		createA (A, k, - w_b * C_b * rho_b - rho * C_rho / delta_t, delta_x, delta_y, max_i, max_j);
	#endif
	
	#ifndef ASM
		createA_C (A, k, - w_b * C_b * rho_b - rho * C_rho / delta_t, delta_x, delta_y, max_i, max_j);
	#endif
	
	#ifdef MEASURE_TIME
		ENDMEASURE(1);
	#endif
	
	double *auxVectorError = malloc(max_i * max_j * sizeof(double));

	// resolucion por Jacobi
	int n; // n es la cantidad de delta_t corridos (el tiempo real es n * delta_t * 2)
	for (n = 1; n <= max_cant_delta_t; n++) {
		
		// Para simular los pulsos electricos
		double* TIndAct = TIndPhiZero;
		if (n % 4000 < 40) 
			TIndAct = TInd;
			
		#ifdef MEASURE_TIME
			STARTMEASURE(2);
		#endif
		
		#ifdef ASM
			updateB(B, Tn, TIndAct, (double) rho * C_rho, 1.0 / delta_t, max_i * max_j);
		#endif
		
		#ifndef ASM
			updateB_C(B, Tn, TIndAct, (double) rho * C_rho, delta_t, max_i * max_j);
		#endif
		
		#ifdef MEASURE_TIME
			ENDMEASURE(2);
		#endif
		
		double errorPrevio = 0, errorActual = calcVectorError(A, Tn, B, k, auxVectorError, anodo_x_idx, anodo_y_idx, catodo_x_idx, catodo_y_idx, max_i, max_j, delta_x);
		while (fabs(errorPrevio - errorActual) > EPS) {
			// siempre empieza estando la info bien en Tn, Tn_sig es auxiliar
			// por eso, pongamos una cantidad de iteraciones *par*
			jacobiStepOptimized(Tn_sig, Tn, B, TIndAct, A, k, catodo_x_idx, catodo_y_idx, anodo_x_idx, anodo_y_idx, max_i, max_j, delta_x, T_aire);
			jacobiStepOptimized(Tn, Tn_sig, B, TIndAct, A, k, catodo_x_idx, catodo_y_idx, anodo_x_idx, anodo_y_idx, max_i, max_j, delta_x, T_aire);

			errorPrevio = errorActual;
			errorActual = calcVectorError(A, Tn, B, k, auxVectorError, anodo_x_idx, anodo_y_idx, catodo_x_idx, catodo_y_idx, max_i, max_j, delta_x);		
		}
		
		#ifndef MEASURE_TIME
			// imprime cant_prints elementos, comenzando de max_cant_delta_t
			// y descenciendo de a saltos de ceil(max_cant_delta_t / cant_prints) 
			if ((max_cant_delta_t - n) % ((max_cant_delta_t + cant_prints) / cant_prints) == 0) {
				char str[32];
				sprintf(str, "output/T%d.out", n);
				FILE* Tn_file = fopen(str, "w");
				print2DMatlab(Tn, max_i, max_j, Tn_file);
				fclose(Tn_file);
			}
		#endif
	}
	
	#ifdef MEASURE_TIME
		ENDMEASURE(0);
		printf("%d, %d, %d", max_i, max_j, max_cant_delta_t); // Tiempo en milisegundos
		
		for (i = 0; i < 10; i++) {
			printf(", %.6f", tt[i]/1000000.0); // Tiempo en milisegundos	
		}
		printf("\n");
	#endif
	
	free(Tn);
	free(TInd);
	free(Tn_sig);
	free(B);
	free(k);
	free(sigma);
	free(phi);
	free(A);
	free(phiZero);
	free(TIndPhiZero);
	return 0;
}
