#ifndef MALLA2D_H
#define MALLA2D_H

#include "nr.h"

template <class T>
class Malla2D: public NRMat<T> {
public:
	Malla2D():NRMat<T>() {};
	Malla2D(int n, int m):NRMat<T>(n,m){};			// Zero-based array
	Malla2D(const T &a, int n, int m):NRMat<T>(a,n,m){};	// Initialize to constant value
	Malla2D(const T *a, int n, int m):NRMat<T>(a,n,m){};	// Initialize to array
	Malla2D(const Malla2D &rhs):NRMat<T>(rhs){};		// Copy constructor

	void save(const string name) const;
	void load(const string name);
	void sameValue(const T &a);				//assign a to every element
};


template <class T>
ostream& operator <<(ostream& o, const Malla2D<T> & a)
{

	int nn = a.nrows();
	int mm = a.ncols();

	for (int i=0; i< nn; i++)
		{   for (int j=0; j<mm; j++) o << a[i][j] << " ";
            o << endl;
        }
    o << endl;
	return o;
}

/*
template <class T>
void Malla2D<T>::save(const string name) const{
	ofstream of;
	of.open(name.data());
        cout.precision(15);
	for (int i=0; i< this->nn; i++)
		{ for (int j=0; j<this->mm; j++) of << this->v[i][j] << " ";
		  of << endl;
        }
	of.close();
}
*/
template <class T>
void Malla2D<T>::save(const string name) const{
	ofstream of;
	of.open(name.data());
        //cout.precision(50);
	for (int j=0; j< this->nn; j++){
		for (int i=0; i<this->mm; i++) of << this->v[j][i] << " ";
		of << endl;
        }
	of.close();
}

/*
template <class T>
void Malla2D<T>::load(const string name){
	ifstream arch;
	arch.open(name.data());
	for (int i=0; i< this->nn; i++)
		{ for (int j=0; j<this->mm; j++) arch >> this->v[i][j];
		  //of << endl;
        }
	arch.close();
}
*/
template <class T>
void Malla2D<T>::load(const string name){
	ifstream arch;
	arch.open(name.data());
	for (int j=0; j< this->nn; j++){
		for (int i=0; i<this->mm; i++) arch >> this->v[j][i];
		  //of << endl;
        }
	arch.close();
}


template <class T>
void  Malla2D<T>::sameValue(const T &a)	//assign a to every element
{
	for (int i=0; i< this->nn; i++)
		for (int j=0; j<this->mm; j++)
			this->v[i][j] = a;
}

#endif // MALLA2D_H


