#include "ElectricField.h"

ElectricField::ElectricField(Entorno& e, Potential& Phi, const long double xecx, const long double xecy, const long double xwtx, const long double xwty):
      mx(e.ii+1,e.jj+1),
      my(e.ii+1,e.jj+1),
      mx_k1(e.ii+1,e.jj+1),
      my_k1(e.ii+1,e.jj+1),
      mx_k2(e.ii+1,e.jj+1),
      my_k2(e.ii+1,e.jj+1), 
      errorConvergencia_x(xecx), 
      errorConvergencia_y(xecy),
      wt_x(xwtx,e.ii+1,e.jj+1),
      wt_y(xwty,e.ii+1,e.jj+1) 
    {
	//Calculo Ey
	for(int i=1;i<=e.ii;i++){
		for(int j=0;j<=e.jj;j++){
		      mx[i][j] = (Phi.m[i-1][j] - Phi.m[i][j]) / e.h;
		}
	}

	for(int j=0;j<=e.jj;j++){
		mx[0][j] = mx[1][j];
	}

	//Calculo Ey
	for(int j=1;j<=e.jj;j++){
		for(int i=0;i<=e.ii;i++){
			my[i][j] = (Phi.m[i][j-1] - Phi.m[i][j]) / e.h;
		}
	}

	for(int i=0;i<=e.ii;i++){
		my[i][0] = my[i][1];
	}
	
	mx_k1=mx;
	mx_k2=mx;
	
	my_k1=my;
	my_k2=my;
}


