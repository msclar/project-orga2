#include <stdio.h>
#include <math.h>

// inicializacion de variables
int c_rho = 3680;
int rho = 1039;
int q_ddd = 10347;
int c_b = 3840;
int rho_b = 1060;
double w_b = 0.00715;
double T_a = 37.0;

int max_i = 3, max_j = 3; // numeros al azar
double delta_t = 0.1, delta_x = 0.1, delta_y = 0.1; // numeros al azar
	
int indice (int i, int j, int max_j) {
	return i * max_j + j;
}

void print1D(double matrix[], int cant_elems) {
	int i;
	for (i = 0; i < cant_elems; i++) {
		printf("%f\n", matrix[i]);
	}
	printf("\n");
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

void jacobiStep(double Tn_sig[], 
				double Tn[],
				double B[],
				double TInd[], 
				double A[][max_i * max_j]) {
	// X_i^(iter+1) = (b_i - sum (j != i) a_ij * X_j^iter	
	int i, j, i_a, j_a;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			double sum = 0.0;
		
			int s = indice(i, j, max_j);
			for (i_a = 0; i_a < max_i; i_a++) {
				for (j_a = 0; j_a < max_j; j_a++) {
					sum += A[s][indice(i_a, j_a, max_j)] * Tn[indice(i_a, j_a, max_j)];
				}
			}
			sum -= A[s][s] * Tn[s];
			Tn_sig[s] = (B[s] - sum) / A[s][s];
		}
	}
}

void jacobiStepOptimized(double Tn_sig[], 
				double Tn[],
				double B[],
				double TInd[], 
				double A[][max_i * max_j],
				double k[][max_j],
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
			sum += A[s][indice(i-1, j, max_j)] * Tn[indice(i-1, j, max_j)];
			sum += A[s][indice(i+1, j, max_j)] * Tn[indice(i+1, j, max_j)];
			sum += A[s][indice(i, j-1, max_j)] * Tn[indice(i, j-1, max_j)];
			sum += A[s][indice(i, j+1, max_j)] * Tn[indice(i, j+1, max_j)];
			Tn_sig[s] = (B[s] - sum) / A[s][s];
		}
	}
	
	// bordes inferior y superior
	for (j = 1; j < max_j - 1; j++) {
		Tn_sig[indice(0, j, max_j)] = Tn_sig[indice(1, j, max_j)];
		Tn_sig[indice(max_i-1, j, max_j)] = Tn_sig[indice(max_i-2, j, max_j)];
	}

	// bordes izquierdo y derecho	
	for (i = 1; i < max_i - 1; i++) {
		Tn_sig[indice(i, 0, max_j)] = Tn_sig[indice(i, 1, max_j)];
		Tn_sig[indice(i, max_j-1, max_j)] = Tn_sig[indice(i, max_j-2, max_j)];
	}
		
	// temperatura en electrodos
	double r = k[catodo_x][catodo_y] / (delta_x * 10); //  h = 10 W / m² K
	Tn_sig[indice(catodo_x, catodo_y, max_j)] = 
		r * (
			Tn_sig[indice(catodo_x+1, catodo_y, max_j)] + 
			Tn_sig[indice(catodo_x-1, catodo_y, max_j)] + 
			Tn_sig[indice(catodo_x, catodo_y+1, max_j)] + 
			Tn_sig[indice(catodo_x, catodo_y-1, max_j)] - 4 * T_a
		) / (4 * r - 1); 
	
	r = k[anodo_x][anodo_y] / (delta_x * 10);
	Tn_sig[indice(anodo_x, anodo_y, max_j)] = 
		r * (
			Tn_sig[indice(anodo_x+1, anodo_y, max_j)] + 
			Tn_sig[indice(anodo_x-1, anodo_y, max_j)] + 
			Tn_sig[indice(anodo_x, anodo_y+1, max_j)] + 
			Tn_sig[indice(anodo_x, anodo_y-1, max_j)] - 4 * T_a
		) / (4 * r - 1); 
}

double calcVectorError(double A[][max_i * max_j], 
					   double Tn[], 
					   double B[]) {
	int i, j;
	double res[max_i * max_j];
	
	// calculo A * Tn - B
	for (i = 0; i < max_i * max_j; i++) {
		res[i] = -B[i];
		for (j = 0; j < max_i * max_j; j++) {
			res[i] += A[i][j] * Tn[j];
		}
	}
	return normaVector(res, max_i * max_j);
}

double calcVectorErrorOptimized(double A[][max_i * max_j], 
					   double Tn[], 
					   double B[]) {
	// TODO: code this
	return 1.0;
}

int main( int argc, char** argv ) {
	freopen("test.out", "w", stdout);
	double Tn[max_i * max_j];
	double Tn_sig[max_i * max_j];
	double TInd[max_i * max_j];
	double B[max_i * max_j];
	double k[max_i][max_j];
	double sigma[max_i][max_j];
	double phi[max_i][max_j];
	double A[max_i * max_j][max_i * max_j];
	
	int i, j;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			int s = indice(i, j, max_j);
			Tn[s] = T_a;
			Tn_sig[s] = T_a;
			k[i][j] = 0.565;
			sigma[i][j] = 0.75;
			phi[i][j] = 1500 * i * delta_x * pow(M_E, - (delta_x * i - 2) * (delta_x * i - 2) - (delta_y * j - 2) * (delta_y * j - 2)); // estan bien usados los delta?
			
			double deriv_phi_x = (phi[i+1][j] - phi[i-1][j]) / (2 * delta_x);
			double deriv_phi_y = (phi[i][j+1] - phi[i][j-1]) / (2 * delta_y);
			double gradient_phi_quad = deriv_phi_x * deriv_phi_x + deriv_phi_y * deriv_phi_y;
			TInd[s] = w_b * c_b * rho_b * T_a + q_ddd + sigma[i][j] * gradient_phi_quad;
		}
	}
	
	for (i = 0; i < max_i * max_j; i++) {
		for (j = 0; j < max_i * max_j; j++) {
			A[i][j] = 0;
		}
	}
	
	// creacion de la matriz A
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			// lleno la fila s = indice(i, j) de la matriz A
			int s = indice(i, j, max_j);
			double derivada_k_x = (k[i+1][j] - k[i-1][j]) / (2 * delta_x);
			double derivada_k_y = (k[i][j+1] - k[i][j-1]) / (2 * delta_y);
			
			A[s][indice(i, j, max_j)] = - 2 * k[i][j] / (delta_x * delta_x) - 2 * k[i][j] / (delta_y * delta_y) + w_b * c_b * rho_b - rho * c_rho / delta_t; 
			A[s][indice(i+1, j, max_j)] = derivada_k_x / (2 * delta_x) + k[i+1][j] / (delta_x * delta_x);
			A[s][indice(i-1, j, max_j)] = - derivada_k_x / (2 * delta_x) + k[i-1][j] / (delta_x * delta_x);
			A[s][indice(i, j+1, max_j)] = derivada_k_y / (2 * delta_y) + k[i][j+1] / (delta_y * delta_y);
			A[s][indice(i, j-1, max_j)] = - derivada_k_y / (2 * delta_y) + k[i][j+1] / (delta_y * delta_y);
		}
	}

	// resolucion por Jacobi
	int iter, n;
	for (n = 0; n < 1; n++) {
		for (i = 0; i < max_i * max_j; i++) {
			B[i] = Tn[i] - TInd[i];
		}
	
		for (iter = 0; iter < 6; iter++) {
			// siempre empieza estando la info bien en Tn, Tn_sig es auxiliar
			// por eso, pongamos una cantidad de iteraciones *par*
			if (iter % 2 == 0) {
				jacobiStep(Tn_sig, Tn, B, TInd, A);
				print1D(Tn_sig, max_i * max_j);
			} 
			else {
				jacobiStep(Tn, Tn_sig, B, TInd, A);
				print1D(Tn, max_i * max_j);
			} 
		}
	}
	
	return 0;
}
