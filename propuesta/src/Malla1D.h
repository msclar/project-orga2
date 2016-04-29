#ifndef MALLA1D_H
#define MALLA1D_H

#include "nr.h"

template <class T>
class Malla1D: public NRVec<T> {
public:
	Malla1D():NRVec<T>() {};
	explicit Malla1D(int n):NRVec<T>(n){};		// Zero-based array
	Malla1D(const T &a, int n):NRVec<T>(a,n){};	//initialize to constant value
	Malla1D(const T *a, int n):NRVec<T>(a,n){};	// Initialize to array
	Malla1D(const Malla1D &rhs):NRVec<T>(rhs){};	// Copy constructor

    void save(const string name) const;
    void load(const string name);
	void sameValue(const T &a); //assign a to every element
};


template <class T>
ostream& operator <<(ostream& o, const Malla1D<T> & v)
{
	int nn = v.size();
    for (int i=0; i<nn; i++) o << v[i] << " ";
	o << endl;
	return o;
}

template <class T>
void Malla1D<T>::save(const string name) const{
	ofstream of;
	of.open(name.data());
	for (int i=0; i<this->nn; i++) of << this->v[i] << " ";
	of.close();
}

template <class T>
void Malla1D<T>::load(const string name){
	ifstream arch;
	arch.open(name.data());
	for (int i=0; i<this->nn; i++) arch >> this->v[i];
	arch.close();
}

template <class T>
void  Malla1D<T>::sameValue(const T &a)	//assign a to every element
{
    for (int i=0; i<this->nn; i++) this->v[i]=a;
}

#endif // MALLA1D_H


