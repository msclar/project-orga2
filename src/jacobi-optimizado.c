#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 0.000001

#ifdef ASM
	extern void jacobiStep (double*, double*, double*, double*, int, int);
	extern void laplaceStep (double*, double*, int, int);
	extern void updateB (double*, double*, double*, double, double, int);
	extern void fillWithZeros (double*, int);
	extern void calculateVectorError (double*, double*, double*, double*, int, int);
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
		for(j = 0; j < max_j; j++) {
			if(j) printf(" ");
			printf("%f", matrix[i * max_j + j]);
		}
		printf("\n");
	}
	printf("\n\n");	
}

void print2DMatlab(double matrix[], int max_i, int max_j, FILE* out) {
	int i, j;
	for (i = 0; i < max_i; i++) {
		for(j = 0; j < max_j; j++) {
			if(j) fprintf(out, " ");
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

void pasoLaplace(double* res, double* phi, int max_i, int max_j) {
	#ifdef ASM
		laplaceStep(res, phi, max_i, max_j);
	#endif
	
	int i, j;
	#ifndef ASM
	// promedio los 4 vecinos (estos tienen los 4 vecinos)
	for(i = 1; i < max_i - 1; i++) {
		for(j = 1; j < max_j - 1; j++) {
			res[indice(i , j, max_j)] = (phi[indice(i-1, j, max_j)] + phi[indice(i+1, j, max_j)] + phi[indice(i, j-1, max_j)] + phi[indice(i, j+1, max_j)]) / 4.0;
		}
	}
	#endif
	
	// bordes izquierdo y derecho
	for(i = 1; i < max_i - 1; i++) {
		res[indice(i, 0, max_j)] = res[indice(i, 1, max_j)];
		res[indice(i, max_j-1, max_j)] = res[indice(i, max_j-2, max_j)]; 
	}
	
	// bordes inferior y superior
	for(j = 0; j < max_j; j++) {
		res[indice(0, j, max_j)] = res[indice(1, j, max_j)];
		res[indice(max_i-1, j, max_j)] = res[indice(max_i-2, j, max_j)]; 
	}
	
	return;
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
	int i;
	for(i = 0; i < max_i * max_j; i++) {
		res[i] = (anodov + catodov) / 2.0;
	}
	
	int posanodo = indice(anodoi, anodoj, max_j);
	int poscatodo = indice(catodoi, catodoj, max_j);
	res[posanodo] = anodov;
	res[poscatodo] = catodov;
	
	double *aux = malloc(max_i * max_j * sizeof(double));
	int corridas;
	for(corridas = 0; corridas < 750; corridas++) {		
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
	int i, j;
	
	#ifdef ASM
		//fillWithZeros(res, max_i * max_j);
		calculateVectorError(A, Tn, B, res, max_i, max_j);
	#endif

	#ifndef ASM
	// calculo A * Tn - B
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			res[s] = calcularPosicion(i, j, A, Tn, B, max_i, max_j);
		}
	}
	#endif

	for (i = 0; i < max_i; i++) {
		res[indice(i, 0, max_j)] = 0.0;
		res[indice(i, max_j-1, max_j)] = 0.0;
	}
	
	for (j = 0; j < max_j; j++) {
		res[indice(0, j, max_j)] = 0.0;
		res[indice(max_i-1, j, max_j)] = 0.0;
	}

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

	return normaVector(res, max_i * max_j);
}

/* CONVENCION para los 5 vecinos
 *   0
 * 1 4 3
 *   2
 */

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
	
	int i;
	int j;
	#ifndef ASM
	
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			double sum = 0.0;
			sum += A[indice(0, s, max_i*max_j)] * Tn[indice(i-1, j, max_j)];
			sum += A[indice(2, s, max_i*max_j)] * Tn[indice(i+1, j, max_j)];
			sum += A[indice(1, s, max_i*max_j)] * Tn[indice(i, j-1, max_j)];
			sum += A[indice(3, s, max_i*max_j)] * Tn[indice(i, j+1, max_j)];
			Tn_sig[s] = (B[s] - sum) / A[indice(4, s, max_i*max_j)];
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

void calcularTInd(double* phi, double delta_x, double delta_y, 
				  double* sigma, double* TInd, int max_i, int max_j,
				  double w_b, int C_b, int rho_b, int q_ddd, double T_a) {
					  
	int i, j;
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {			
			double deriv_phi_x = (phi[indice(i+1, j, max_j)] - phi[indice(i-1, j, max_j)]) / (2 * delta_x);
			double deriv_phi_y = (phi[indice(i, j+1, max_j)] - phi[indice(i, j-1, max_j)]) / (2 * delta_y);
			double gradient_phi_quad = deriv_phi_x * deriv_phi_x + deriv_phi_y * deriv_phi_y;
			TInd[indice(i, j, max_j)] = w_b * C_b * rho_b * T_a + q_ddd + sigma[indice(i, j, max_j)] * gradient_phi_quad;
		}
	}
}

void updateB_C (double* B, double* Tn, double* TIndAct, double rho_times_C_rho, double delta_t, int max_ij) {
	int i;
	for (i = 0; i < max_ij; i++) {
		B[i] = - rho_times_C_rho * Tn[i] / delta_t - TIndAct[i];
	}
}

int main( int argc, char** argv ) {
	// inicializacion de variables
	int C_rho = 3680;
	int rho = 1039;
	int q_ddd = 10347;
	int C_b = 3840;
	int rho_b = 1060;
	double w_b = 0.00715;
	double T_a = 310.15;
	double T_aire = 296;

	double delta_t = 0.000005, delta_x, delta_y; // numeros al azar


	double *Tn, *Tn_sig, *TInd, *B, *k, *sigma, *phi, *A, *phiZero, *TIndPhiZero;
	
	int max_i = 53, max_j = 59;
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

	double catodo_x = 0.025, catodo_y = 0.05, anodo_x = 0.075, anodo_y = 0.05;
	
	delta_x = 0.1 / max_i; delta_y = 0.1 / max_j; // recinto de 0.1cm x 0.1cm
	
	int catodo_x_idx = catodo_x / delta_x;
	int catodo_y_idx = catodo_y / delta_y;
	double catodo_v = 0;
	int anodo_x_idx = anodo_x / delta_x;
	int anodo_y_idx = anodo_y / delta_y;
	double anodo_v = 1000;
	
	printf("%d, %d, %d, %d\n", catodo_x_idx, catodo_y_idx, anodo_x_idx, anodo_y_idx);
	
	int i, j;
	obtenerLaplace(phi, anodo_x_idx, anodo_y_idx, anodo_v, catodo_x_idx, catodo_y_idx, catodo_v, max_i, max_j);
	FILE* phi_file = fopen("phi.out", "w");
	print2DMatlab(phi, max_i, max_j, phi_file);
	
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

	// calculo de TInd
	calcularTInd(phi, delta_x, delta_y, sigma, TInd, max_i, max_j, w_b, C_b, rho_b, q_ddd, T_a);
	calcularTInd(phiZero, delta_x, delta_y, sigma, TIndPhiZero, max_i, max_j, w_b, C_b, rho_b, q_ddd, T_a);
	
	// creacion de la matrix A
	for (i = 0; i < max_i * max_j * 5; i++) {
		A[i] = 0;
	}
	
	// calculo de A
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			
			double derivada_k_x = (k[indice(i+1, j, max_j)] - k[indice(i-1, j, max_j)]) / (2 * delta_x);
			double derivada_k_y = (k[indice(i, j+1, max_j)] - k[indice(i, j-1, max_j)]) / (2 * delta_y);
			
			A[indice(4, s, max_i*max_j)] = - 2 * k[s] / (delta_x * delta_x) - 2 * k[s] / (delta_y * delta_y) - w_b * C_b * rho_b - rho * C_rho / delta_t;
			A[indice(2, s, max_i*max_j)] = derivada_k_x / (2 * delta_x) + k[indice(i+1, j, max_j)] / (delta_x * delta_x);
			A[indice(0, s, max_i*max_j)] = - derivada_k_x / (2 * delta_x) + k[indice(i-1, j, max_j)] / (delta_x * delta_x);
			A[indice(3, s, max_i*max_j)] = derivada_k_y / (2 * delta_y) + k[indice(i, j+1, max_j)] / (delta_y * delta_y);
			A[indice(1, s, max_i*max_j)] = - derivada_k_y / (2 * delta_y) + k[indice(i, j-1, max_j)] / (delta_y * delta_y);
		}
	}

	double *auxVectorError = malloc(max_i * max_j * sizeof(double));

	// resolucion por Jacobi
	int n;
	for (n = 0; n < 100000; n++) {
		double* TIndAct = TIndPhiZero;
		if(n % 4000 < 40) TIndAct = TInd;
			
		#ifdef ASM
			updateB(B, Tn, TIndAct, (double) rho * C_rho, delta_t, max_i * max_j);
		#endif
		
		#ifndef ASM
			updateB_C(B, Tn, TIndAct, (double) rho * C_rho, delta_t, max_i * max_j);
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
		
	}
	
	char str[32];
	sprintf(str, "output/T%d.out", n);
	FILE* Tn_file = fopen(str, "w");
	print2DMatlab(Tn, max_i, max_j, Tn_file);
	fclose(Tn_file);
	fclose(phi_file);
	
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
