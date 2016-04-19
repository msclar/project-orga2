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

int max_i = 20, max_j = 20; // numeros al azar
double delta_t = 0.1, delta_x = 0.1, delta_y = 0.1; // numeros al azar
	
int indice (int i, int j, int max_j) {
	return i * max_j + j;
}

void jacobiStep(double Tn_sig[][max_j], 
				double Tn[][max_j], 
				double TInd[][max_j], 
				double A[][max_i * max_j]) {
	// X_i^(iter+1) = (b_i - sum (j != i) a_ij * X_j^iter	
	int i, j, i_a, j_a;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			double sum = 0.0;
		
			int s = indice(i, j, max_j);
			for (i_a = 0; i_a < max_i; i_a++) {
				for (j_a = 0; j_a < max_j; j_a++) {
					sum += A[s][indice(i_a, j_a, max_j)] * Tn[i_a][j_a];
				}
			}
			sum -= A[s][s] * Tn[i][j];
			
			Tn_sig[i][j] = Tn[i][j] - TInd[i][j] - sum;
		}
	}
}

int main( int argc, char** argv ) {
	double Tn[max_i][max_j];
	double Tn_sig[max_i][max_j];
	double TInd[max_i][max_j];
	double k[max_i][max_j];
	double sigma[max_i][max_j];
	double phi[max_i][max_j];
	double A[max_i * max_j][max_i * max_j];
	
	int i, j;
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			Tn[i][j] = T_a;
			Tn_sig[i][j] = T_a;
			k[i][j] = 0.565;
			sigma[i][j] = 0.75;
			phi[i][j] = 1500 * i * delta_x * pow(M_E, - (delta_x * i - 2) * (delta_x * i - 2) - (delta_y * j - 2) * (delta_y * j - 2)); // estan bien usados los delta?
			
			double deriv_phi_x = (phi[i+1][j] - phi[i-1][j]) / (2 * delta_x);
			double deriv_phi_y = (phi[i][j+1] - phi[i][j-1]) / (2 * delta_y);
			double gradient_phi_quad = deriv_phi_x * deriv_phi_x + deriv_phi_y * deriv_phi_y;
			TInd[i][j] = w_b * c_b * rho_b * T_a + q_ddd + sigma[i][j] * gradient_phi_quad;
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
	// CAMBIAR ITER POR EL CALCULO DE LA CONVERGENCIA
	// QUE IMPRIMIMOS COMO OUTPUT?
	int iter, n;
	for (n = 0; n < 10; n++) {
		if (n % 2 == 0) {
			for (iter = 0; iter < 10; iter++)
				jacobiStep(Tn_sig, Tn, TInd, A);
		} 
		else {
			for (iter = 0; iter < 10; iter++)
				jacobiStep(Tn, Tn_sig, TInd, A);
		}
	}
	
	return 0;
}
