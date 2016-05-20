#include <stdio.h>
#include <math.h>
#define EPS 0.001
// inicializacion de variables
int C_rho = 3680;
int rho = 1039;
int q_ddd = 10347;
int C_b = 3840;
int rho_b = 1060;
double w_b = 0.00715;
double T_a = 310.15;
double T_aire = 296;

int max_i = 50, max_j = 50; // numeros al azar
double delta_t = 0.1, delta_x, delta_y; // numeros al azar


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

void print2DMatlab(double matrix[], int max_i, int max_j, FILE* out) {
	int i, j;
	fprintf(out, "[");
	for (i = 0; i < max_i; i++) {
		fprintf(out, "[");
		for(j = 0; j < max_j; j++) {
			if(j) fprintf(out, " ");
			fprintf(out, "%f", matrix[i * max_j + j]);
		}
		fprintf(out, "];\n");
	}
	fprintf(out, "]\n\n");	
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


void pasoLaplace(double* res, double* phi, int imax, int jmax) {
	int i, j;
	for(i = 1; i < imax - 1; i++) {
		for(j = 1; j < jmax - 1; j++) {
			int posactual = indice(i , j, jmax);
			res[posactual] = (phi[posactual - 1] + phi[posactual + 1] + phi[posactual - jmax] + phi[posactual + jmax]) / 4.0;
		}
	}
	
	for(i = 1; i < imax - 1; i++) {
		res[indice(i, 0, jmax)] = res[indice(i, 1, jmax)];
		res[indice(i, jmax-1, jmax)] = res[indice(i, jmax-2, jmax)]; 
	}
	
	for(j = 0; j < imax; j++) {
		res[indice(0, j, jmax)] = res[indice(1, j, jmax)];
		res[indice(imax-1, j, jmax)] = res[indice(imax-2, j, jmax)]; 
	}
	
	return;
}

void obtenerLaplace(double* res, int anodoi, int anodoj, double anodov, int catodoi, int catodoj, double catodov, int imax, int jmax) {
	int i;
	for(i = 0; i < imax*jmax; i++) {
		res[i] = (anodov + catodov) / 2.0;
	}
	
	int posanodo = indice(anodoi, anodoj, jmax);
	int poscatodo = indice(catodoi, catodoj, jmax);
	res[posanodo] = anodov;
	res[poscatodo] = catodov;
	
	double aux[imax][jmax];
	int corridas;
	for(corridas = 0; corridas < 750; corridas++) {		
		pasoLaplace((double*) aux, res, imax, jmax),
		aux[anodoi][anodoj] = anodov;
		aux[catodoi][catodoj] = catodov;
		
		pasoLaplace(res,(double*) aux, imax, jmax),
		res[posanodo] = anodov;
		res[poscatodo] = catodov;
	}	
	return;
}

double calcularPosicion(int i,
						int j,
						double A[],
						double Tn[],
						double B[]) {
	int s = indice(i, j, max_j);
	double res = -B[s];
	res += A[indice(s, 0, 5)] * Tn[indice(i - 1, j, max_j)];
	res += A[indice(s, 1, 5)] * Tn[indice(i, j - 1, max_j)];
	res += A[indice(s, 2, 5)] * Tn[indice(i + 1, j, max_j)];
	res += A[indice(s, 3, 5)] * Tn[indice(i, j + 1, max_j)];
	res += A[indice(s, 4, 5)] * Tn[indice(i, j, max_j)];
	return res;
}

double calcVectorError(double A[], 
					   double Tn[], 
					   double B[],
					   double k[],
					   int anodo_x,
					   int anodo_y,
					   int catodo_x,
					   int catodo_y) {
	int i, j;
	double res[max_i * max_j];
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			res[indice(i, j, max_j)] = 0;
		}
	}
	// calculo A * Tn - B
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			res[s] = calcularPosicion(i, j, A, Tn, B);
		}
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
 *   1
 * 0 4 2
 *   3
 */

void jacobiStepOptimized(double Tn_sig[], 
				double Tn[],
				double B[],
				double TInd[], 
				double A[],
				double k[],
				int catodo_x, 
				int catodo_y, 
				int anodo_x, 
				int anodo_y) {
	// X_i^(iter+1) = (b_i - sum (j != i) a_ij * X_j^iter	
	
	int i, j;
	// los que tienen los 4 vecinos
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {
			int s = indice(i, j, max_j);
			double sum = 0.0;
			sum += A[indice(s, 0, 5)] * Tn[indice(i-1, j, max_j)];
			sum += A[indice(s, 2, 5)] * Tn[indice(i+1, j, max_j)];
			sum += A[indice(s, 1, 5)] * Tn[indice(i, j-1, max_j)];
			sum += A[indice(s, 3, 5)] * Tn[indice(i, j+1, max_j)];
			Tn_sig[s] = (B[s] - sum) / A[indice(s, 4, 5)];
		}
	}
	
	// bordes inferior y superior
	for (j = 1; j < max_j - 1; j++) {
		Tn_sig[indice(0, j, max_j)] = Tn_sig[indice(1, j, max_j)];
		Tn_sig[indice(max_i-1, j, max_j)] = Tn_sig[indice(max_i-2, j, max_j)];
	}

	// bordes izquierdo y derecho	
	for (i = 0; i < max_i; i++) {
		Tn_sig[indice(i, 0, max_j)] = Tn_sig[indice(i, 1, max_j)];
		Tn_sig[indice(i, max_j-1, max_j)] = Tn_sig[indice(i, max_j-2, max_j)];
	}
	
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

int main( int argc, char** argv ) {
	freopen("test2.out", "w", stdout);
	double Tn[max_i * max_j];
	double Tn_sig[max_i * max_j];
	double TInd[max_i * max_j];
	double B[max_i * max_j];
	double k[max_i * max_j];
	double sigma[max_i * max_j];
	double phi[max_i * max_j];
	double A[max_i * max_j * 5];

	double catodo_x = 1.365, catodo_y = 2, anodo_x = 2.731, anodo_y = 2;
	
	delta_x = 4.0 / max_i; delta_y = 4.0 / max_j; // recinto de 4cm x 4cm
	
	int catodo_x_idx = catodo_x / delta_x;
	int catodo_y_idx = catodo_y / delta_y;
	double catodo_v = 0;
	int anodo_x_idx = anodo_x / delta_x;
	int anodo_y_idx = anodo_y / delta_y;
	double anodo_v = 1500;
	
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
			//phi[s] = 1500 * i * delta_x * pow(M_E, - (delta_x * i - 2) * (delta_x * i - 2) - (delta_y * j - 2) * (delta_y * j - 2)); // estan bien usados los delta?
		}
	}

	// calculo de TInd
	for (i = 1; i < max_i - 1; i++) {
		for (j = 1; j < max_j - 1; j++) {			
			double deriv_phi_x = (phi[indice(i+1, j, max_j)] - phi[indice(i-1, j, max_j)]) / (2 * delta_x);
			double deriv_phi_y = (phi[indice(i, j+1, max_j)] - phi[indice(i, j-1, max_j)]) / (2 * delta_y);
			double gradient_phi_quad = deriv_phi_x * deriv_phi_x + deriv_phi_y * deriv_phi_y;
			TInd[indice(i, j, max_j)] = w_b * C_b * rho_b * T_a + q_ddd + sigma[indice(i, j, max_j)] * gradient_phi_quad;
		}
	}
	
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
			
			A[indice(s, 4, 5)] = - 2 * k[s] / (delta_x * delta_x) - 2 * k[s] / (delta_y * delta_y) - w_b * C_b * rho_b - rho * C_rho / delta_t;
			A[indice(s, 2, 5)] = derivada_k_x / (2 * delta_x) + k[indice(i+1, j, max_j)] / (delta_x * delta_x);
			A[indice(s, 0, 5)] = - derivada_k_x / (2 * delta_x) + k[indice(i-1, j, max_j)] / (delta_x * delta_x);
			A[indice(s, 3, 5)] = derivada_k_y / (2 * delta_y) + k[indice(i, j+1, max_j)] / (delta_y * delta_y);
			A[indice(s, 1, 5)] = - derivada_k_y / (2 * delta_y) + k[indice(i, j-1, max_j)] / (delta_y * delta_y);
		}
	}

	// resolucion por Jacobi
	print2DMatlab(Tn, max_i, max_j, phi_file);
	
	int n;
	for (n = 0; n < 400; n++) {
		for (i = 0; i < max_i * max_j; i++) {
			B[i] = - rho * C_rho * Tn[i] / delta_t - TInd[i];
		}
	
		while (calcVectorError(A, Tn, B, k, anodo_x_idx, anodo_y_idx, catodo_x_idx, catodo_y_idx) > EPS) {
			// siempre empieza estando la info bien en Tn, Tn_sig es auxiliar
			// por eso, pongamos una cantidad de iteraciones *par*
			jacobiStepOptimized(Tn_sig, Tn, B, TInd, A, k, catodo_x_idx, catodo_y_idx, anodo_x_idx, anodo_y_idx);
			jacobiStepOptimized(Tn, Tn_sig, B, TInd, A, k, catodo_x_idx, catodo_y_idx, anodo_x_idx, anodo_y_idx);
			printf("%.6f\n", calcVectorError(A, Tn, B, k, anodo_x_idx, anodo_y_idx, catodo_x_idx, catodo_y_idx));
		}
		
		char str[32];
		sprintf(str, "output/T%d.out", n);
		FILE* Tn_file = fopen(str, "w");
		print2DMatlab(Tn, max_i, max_j, Tn_file);
		fclose(Tn_file);
	}
	fclose(phi_file);
	
	return 0;
}