#include <stdio.h>

//   a
// b c d
//   e

double jacobiIterCoeficiente(double Aa, double Ab, double Ac, double Ad, double Ae, double xa, double xb, double bc, double xd, double xe) {
	double Rx = Aa * xa + Ab * xb + Ad * xd + Ae * xe;
	return (bc - Rx) / Ac;
}

void jacobiIter(double* A, double* b, double* x0, double* x1, int iMax, int jMax) {
	int i, j;
	for(i = 1; i < iMax-1; i++) {
		for(j = 1; j < jMax-1; j++) {
			int pos = i*jMax+j;
			x1[i] = jacobiIterCoeficiente(
				A[pos * iMax*jMax + pos - jMax], A[pos * iMax*jMax + pos - 1], A[pos * iMax*jMax + pos],
				A[pos * iMax*jMax + pos + 1], A[pos * iMax*jMax + pos + jMax], 
				x0[pos - jMax], x0[pos - 1], b[pos], x0[pos + 1], x0[pos + jMax]);
		}
	}
	return;
} 

void calcVectorError(double* A, double* x, double* b, double* error, int iMax, int jMax) {
	int i;
	for(i = 1; i < iMax * jMax; i++) {
		for(i = 1; i < iMax * jMax; j++) {
			//HOLA MEL, COMPLETA ESTO (?)
		}
	}
}

int main( int argc, char** argv ) {
	return 0;
}
