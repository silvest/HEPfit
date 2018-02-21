/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Expanded.h
 * Author: enrico
 *
 * Created on 16 novembre 2017, 12.02
 */

#ifndef EXPANDED_H
#define EXPANDED_H

#include <vector>
#include <iostream>

#include "gslpp_complex.h"
#include "gslpp_vector_double.h"
#include "gslpp_vector_complex.h"
#include "gslpp_matrix_complex.h"
#include "gslpp_matrix_double.h"

typedef unsigned int uint;

using namespace gslpp;

#define isMat(X,M) (std::is_same<X,matrix<M> >::value)
#define isVec(X,V) (std::is_same<X,vector<V> >::value)
#define isSc(X,S)  (std::is_same<X,S>::value)
#define isS(X)   (isSc(X,double) || isSc(X,complex))
#define isV(X)   (isVec(X,double) || isVec(X,complex))
#define isM(X)   (isMat(X,double) || isMat(X,complex))
#define isComp(X) (isSc(X,complex) || isVec(X,complex) || isMat(X,complex))
#define sameType(X,Y) ( (isS(X) && isS(Y)) || (isV(X) && isV(Y)) || (isM(X) && isM(Y)))

template <class T> class Expanded {
public:

    Expanded<T>(std::vector<T> & dinp) {
        ord1 = dinp.size();
        ord2 = 0;
        data1 = dinp;
    }

    Expanded<T>(std::vector<std::vector<T> > & dinp) {
        ord1 = dinp.size();
        ord2 = dinp[0].size();
        data2 = dinp;
    }

    T Zero() const {
        return T();
    }
    
    template <class Q>
    typename std::enable_if<
       isSc(T,double) ||
       (isSc(T,complex)  && isComp(Q)) ||
       (isMat(T,double) && !isS(Q))  ||
       (isMat(T,complex) && (isVec(Q,complex) || isMat(Q,complex))), Expanded<Q> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<Q> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<Q>(res);
        } else {
            std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<Q>(res);
        }
    }

    template <class Q>
    typename std::enable_if<
     (!isSc(T,double) && isSc(Q,double)) ||
     (!isS(T) && !isMat(T,double) && isMat(Q,double)) ||
     ((isMat(T,complex) || isVec(T,complex)) && isSc(Q,complex)) || 
     (isVec(T,complex) && isMat(Q,complex)), Expanded<T> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<T> res(up1, Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(up1, std::vector<T>(up2, Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<T>(res);
        }
    }

    template <class Q>
    typename std::enable_if<
    (isVec(T,double) && isSc(Q,complex)) ||
    (isVec(T,double) && isMat(Q,complex)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<vector<complex> > res(up1,Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<vector<complex> >(res);
        } else {
            std::vector<std::vector<vector<complex> > > res(up1, std::vector<vector<complex> >(up2,Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<vector<complex> >(res);
        }
    }

    template <class Q>
    typename std::enable_if<
    (isSc(T,complex) && isVec(Q,double)) ||
    (isMat(T,complex) && isVec(Q,double)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<vector<complex> > res(up1,z.Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<vector<complex> >(res);
        } else {
            std::vector<std::vector<vector<complex> > > res(up1, std::vector<vector<complex> >(up2,z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<vector<complex> >(res);
        }
    }


    template <class Q>
    typename std::enable_if<
    (isVec(T,double) && isVec(Q,complex)) ||
    (isVec(T,complex) && isVec(Q,double)) ||
    (isVec(T,complex) && isVec(Q,complex)), Expanded<complex> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<complex> res(up1, 0.);
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<complex>(res);
        } else {
            std::vector<std::vector<complex> > res(up1, std::vector<complex>(up2, 0.));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<complex>(res);
        }
    }

    template <class Q>
    typename std::enable_if<
    (isMat(T,double) && isSc(Q,complex)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<matrix<complex> > res(up1, Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<matrix<complex> >(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(up1, std::vector<matrix<complex> >(up2, Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<matrix<complex> >(res);
        }
    }
    template <class Q>
    typename std::enable_if<
    (isSc(T,complex) && isMat(Q,double)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<matrix<complex> > res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<matrix<complex> >(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(up1, std::vector<matrix<complex> >(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<matrix<complex> >(res);
        }
    }

    
    template <class Q>
    typename std::enable_if<
    isVec(T,double) && isVec(Q,double), Expanded<double> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<double> res(up1, 0.);
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<double>(res);
        } else {
            std::vector<std::vector<double> > res(up1, std::vector<double>(up2, 0.));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    for (uint j1 = 0; j1 <= i1; j1++)
                        for (uint j2 = 0; j2 <= i2; j2++)
                            res[i1][i2] += data2[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
            return Expanded<double>(res);
        }
    }

// Operator +

    template <class Q>
    typename std::enable_if<
       sameType(T,Q) && isComp(T), Expanded<T> >::type
    operator+(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<T> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                res[i] = data1[i] + z.getOrd(i);
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(up1, std::vector<T>(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    res[i1][i2] = data2[i1][i2] + z.getOrd(i1, i2);
            return Expanded<T>(res);
        }
    }

    template <class Q>
    typename std::enable_if<
       sameType(T,Q) && !isComp(T), Expanded<Q> >::type
    operator+(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<Q> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                res[i] = data1[i] + z.getOrd(i);
            return Expanded<Q>(res);
        } else {
            std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    res[i1][i2] = data2[i1][i2] + z.getOrd(i1, i2);
            return Expanded<Q>(res);
        }
    }

// Operator -

    template <class Q>
    typename std::enable_if<
       sameType(T,Q) && isComp(T), Expanded<T> >::type
    operator-(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<T> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                res[i] = data1[i] - z.getOrd(i);
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(up1, std::vector<T>(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    res[i1][i2] = data2[i1][i2] - z.getOrd(i1, i2);
            return Expanded<T>(res);
        }
    }

    template <class Q>
    typename std::enable_if<
       sameType(T,Q) && !isComp(T), Expanded<Q> >::type
    operator-(const Expanded<Q>& z) const {
        uint up1=std::min(ord1,z.getN1());
        uint up2=std::min(ord2,z.getN2());
        if (up2 == 0) {
            std::vector<Q> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                res[i] = data1[i] - z.getOrd(i);
            return Expanded<Q>(res);
        } else {
            std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
            for (uint i1 = 0; i1 < up1; i1++)
                for (uint i2 = 0; i2 < up2; i2++)
                    res[i1][i2] = data2[i1][i2] - z.getOrd(i1, i2);
            return Expanded<Q>(res);
        }
    }
    
    Expanded<T>  operator-() const {
        if (ord2 == 0) {
            std::vector<T> res(ord1, Zero());
            for (uint i = 0; i < ord1; i++)
                res[i] = -data1[i];
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                    res[i1][i2] = -data2[i1][i2];
            return Expanded<T>(res);
        }
    }


    
/***************** End Expanded-Expanded *****************/

/***************** Begin Expanded * UnExpanded *****************/
    
// Return T.
    
    template <class Q>
    typename std::enable_if<
     (std::is_same<T,Q>::value  && !isV(T)) ||
     isSc(Q,double) ||
     (!isS(T) && isMat(Q,double)) ||
     ((isMat(T,complex) || isVec(T,complex)) && isS(Q)) || 
     (isVec(T,complex) && isMat(Q,complex)) ||
     (isV(T) && isSc(Q,double)) ||
     (isMat(T,double) && isSc(Q,double)),Expanded<T> >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<T> res(ord1, Zero());
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] * z;
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<T>(res);
        }
    }
    
// Return double

    template <class Q>
    typename std::enable_if<
      isVec(T,double) && isVec(Q,double) ,Expanded<double> >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<double> res(ord1, 0.);
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<double>(res);
        } else {
            std::vector<std::vector<double> > res(ord1, std::vector<double>(ord2, 0.));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<double>(res);
        }
    }

// Return complex

    template <class Q>
    typename std::enable_if<
      (isSc(T,double) && isSc(Q,complex)) ||
      (isV(T) && isV(Q) && !(isVec(T,double) && isVec(Q,double))),Expanded<complex> >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<complex> res(ord1, 0.);
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<complex>(res);
        } else {
            std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2, 0.));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<complex>(res);
        }
    }
// Return vector<double>

    template <class Q>
    typename std::enable_if<
      isVec(Q,double) && (isSc(T,double) || isMat(T,double)), Expanded<vector<double> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<vector<double> > res(ord1, vector<double>(z.size(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<vector<double> >(res);
        } else {
            std::vector<std::vector<vector<double> > > res(ord1, std::vector<vector<double> >(ord2, vector<double>(z.size(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<vector<double> >(res);
        }
    }


// Return vector<complex>

    template <class Q>
    typename std::enable_if<
      (isVec(Q,complex) && (isS(T) || isM(T))) ||
      (isVec(Q,double) && (isSc(T,complex) || isMat(T,complex))),Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<vector<complex> > res(ord1, vector<complex>(z.size(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<vector<complex> >(res);
        } else {
            std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(z.size(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<vector<complex> >(res);
        }
    }

    template <class Q>
    typename std::enable_if<
      isVec(T,double) && (isSc(Q,complex) || isMat(Q,complex)),Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<vector<complex> > res(ord1, vector<complex>(data1[0].size(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<vector<complex> >(res);
        } else {
            std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(data2[0][0].size(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<vector<complex> >(res);
        }
    }

// Return matrix<double>
    
    template <class Q>
    typename std::enable_if<
      (isSc(T,double) && isMat(Q,double)) ,Expanded<matrix<double> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<matrix<double> > res(ord1, matrix<double>(z.size_i(),z.size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<matrix<double> >(res);
        } else {
            std::vector<std::vector<matrix<double> > > res(ord1, std::vector<matrix<double> >(ord2, matrix<double>(z.size_i(),z.size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<matrix<double> >(res);
        }
    }

// Return matrix<complex>
    
    template <class Q>
    typename std::enable_if<
      (isS(T) && isMat(Q,complex)) ||
      (isSc(T,complex) && isMat(Q,double)) ||
      (isMat(T,double) && isMat(Q,complex)) ,Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<matrix<complex> > res(ord1, matrix<complex>(z.size_i(),z.size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<matrix<complex> >(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(z.size_i(),z.size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<matrix<complex> >(res);
        }
    }

    template <class Q>
    typename std::enable_if<
      isMat(T,double) && isSc(Q,complex) ,Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<matrix<complex> > res(ord1, matrix<complex>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<matrix<complex> >(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<matrix<complex> >(res);
        }
    }

/***************** End Expanded * UnExpanded *****************/

/***************** Begin Expanded / UnExpanded-Scalar *****************/

    
    Expanded<T> operator/(const double& z) const {
        if (ord2 == 0) {
            std::vector<T> res(ord1, Zero());
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] / z;
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] / z;
            return Expanded<T>(res);
        }
    }
    
    template <class Q>
    typename std::enable_if<isComp(T) && isSc(Q,complex),Expanded<T> >::type
      operator/(const Q& z) const {
        if (ord2 == 0) {
            std::vector<T> res(ord1, Zero());
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] / z;
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] / z;
            return Expanded<T>(res);
        }
    }
    template <class Q>
    typename std::enable_if<isSc(T,double) && isSc(Q,complex),Expanded<complex > >::type
    operator/(const Q& z) const {
        if (ord2 == 0) {
            std::vector<complex> res(ord1, 0.);
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] / z;
            return Expanded<complex>(res);
        } else {
            std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2, 0.));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] / z;
            return Expanded<complex>(res);
        }
      }
    template <class Q>
    typename std::enable_if<isVec(T,double) && isSc(Q,complex),Expanded<vector<complex> > >::type
      operator/(const Q& z) const {
        if (ord2 == 0) {
            std::vector<vector<complex> > res(ord1, vector<complex>(data1[0].size(),0.));
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] / z;
            return Expanded<vector<complex> >(res);
        } else {
            std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(data2[0][0].size(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] / z;
            return Expanded<vector<complex> >(res);
        }
    }
    template <class Q>
    typename std::enable_if<isMat(T,double) && isSc(Q,complex),Expanded<matrix<complex> > >::type
      operator/(const Q& z) const {
        if (ord2 == 0) {
            std::vector<matrix<complex> > res(ord1, matrix<complex>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                res[i] = data1[i] / z;
            return Expanded<matrix<complex> >(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] / z;
            return Expanded<matrix<complex> >(res);
        }
    }

    /***************** End Expanded / UnExpanded-Scalar  *****************/
    Expanded<double> real() const;

    Expanded<double> abs2() const;

    Expanded<T> conjugate() const;

    Expanded<T> transpose() const;
    
    T getAllOrd() {
        T res();
        if(ord2 == 0) {
            for (uint i=0;i<ord1;i++)
                res += data1[i];
        } else {
            for (uint i=0;i<ord1;i++)
                for (uint j=0;j<ord2;j++)
                    res += data2[i][j];
        }
        return res;
    }

    T getOrd(uint i) const {
        return data1[i];
    }

    T getOrd(uint i, uint j) const {
        return data2[i][j];
    }

    void setOrd(uint i, T value) {
        data1[i] = value;
    }

    void setOrd(uint i, uint j, T value) {
        data2[i][j] = value;
    }

    const uint getN1() const {
        return ord1;
    }

    const uint getN2() const {
        return ord2;
    }
    
    private:
    uint ord1, ord2;
    std::vector<T> data1;
    std::vector<std::vector<T> > data2;
};

// Template specialization

template <> Expanded<complex> Expanded<complex>::conjugate() const {
    uint ord1 = getN1();
    uint ord2 = getN2();
    if (ord2 == 0) {
        std::vector<complex> res(ord1);
        for (uint i = 0; i < ord1; i++)
            res[i] = getOrd(i).conjugate();
        return Expanded<complex>(res);
    } else {
        std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2));
        for (uint i = 0; i < ord1; i++)
            for (uint j = 0; j < ord2; j++)
                res[i][j] = getOrd(i, j).conjugate();
        return Expanded<complex>(res);
    }
}

template <> Expanded<matrix<double> > Expanded<matrix<double> >::transpose() const {
        if (ord2 == 0) {
            std::vector<matrix<double> > res(ord1,matrix<double>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i].transpose();
            return Expanded<matrix<double>>(res);
        } else {
            std::vector<std::vector<matrix<double> > > res(ord1, std::vector<matrix<double> >(ord2, matrix<double>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2].transpose();
            return Expanded<matrix<double>>(res);
        }
    }

template <> Expanded<matrix<complex> > Expanded<matrix<complex> >::transpose() const {
        if (ord2 == 0) {
            std::vector<matrix<complex> > res(ord1,matrix<complex>(data1[0].size_i(),data1[0].size_j(),0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i].transpose();
            return Expanded<matrix<complex>>(res);
        } else {
            std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data2[0][0].size_i(),data2[0][0].size_j(),0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2].transpose();
            return Expanded<matrix<complex>>(res);
        }
    }

template <> Expanded<double> Expanded<complex>::real() const {
    uint ord1 = getN1();
    uint ord2 = getN2();
    if (ord2 == 0) {
        std::vector<double> res(ord1, 0.);
        for (uint i = 0; i < ord1; i++)
            res[i] = getOrd(i).real();
        return Expanded<double>(res);
    } else {
        std::vector<std::vector<double> > res(ord1, std::vector<double>(ord2, 0.));
        for (uint i = 0; i < ord1; i++)
            for (uint j = 0; j < ord2; j++)
                res[i][j] = getOrd(i, j).real();
        return Expanded<double>(res);
    }
}

template <> Expanded<double> Expanded<complex>::abs2() const {
    return ((*this) * ((*this).conjugate())).real();
}

template <>
matrix<double> Expanded<matrix<double> >::Zero() const {
    uint id, jd;
    if(ord2 == 0) { 
        id=data1[0].size_i();
        jd=data1[0].size_j();
    } else {
        id=data2[0][0].size_i();
        jd=data2[0][0].size_j();        
    }
    return matrix<double>(id,jd,0.);
}

template <>
matrix<complex> Expanded<matrix<complex> >::Zero() const {
    uint id, jd;
    if(ord2 == 0) { 
        id=data1[0].size_i();
        jd=data1[0].size_j();
    } else {
        id=data2[0][0].size_i();
        jd=data2[0][0].size_j();        
    }
    return matrix<complex>(id,jd,0.);
}

template<>
vector<double> Expanded<vector<double> >::Zero() const {
    uint id;
    if(ord2 == 0) { 
        id=data1[0].size();
    } else {
        id=data2[0][0].size();
    }
    return vector<double>(id,0.);
}

template<>
vector<complex> Expanded<vector<complex> >::Zero() const {
    uint id;
    if(ord2 == 0) { 
        id=data1[0].size();
    } else {
        id=data2[0][0].size();
    }
    return vector<complex>(id,0.);
}

template<class T>
std::ostream& operator<<(std::ostream& output, const Expanded<T> x) {
    if(x.getN2()==0) {
        for (uint i=0;i<x.getN1();i++)
           output << "Order " << i << ": " << x.getOrd(i) << std::endl;
    } else {
        for (uint i=0;i<x.getN1();i++)
            for (uint j=0;j<x.getN2();j++)
                output << "Order (" << i << "," << j << "): " << x.getOrd(i,j) << std::endl;
    }
    return output;
}

/***************** Begin UnExpanded * Expanded *****************/

// scalar * X -> commute
    template <class T>
    Expanded<T> operator*(const double& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    template <class T>
    typename std::enable_if<
      isComp(T), Expanded<T> >::type
     operator*(const complex& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    Expanded<complex> operator*(const complex& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<vector<complex> > operator*(const complex& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue);
    }
    Expanded<matrix<complex> > operator*(const complex& ue, const Expanded<matrix<double> >& ex) {
        return (ex * ue);
    }
    
// X * scalar -> commute

    Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<vector<complex> > operator*(const vector<double>& ue, const Expanded<complex>& ex) {
        return (ex * ue);
    }
    Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<double>& ex) {
        return (ex * ue);
    }
    Expanded<matrix<complex> > operator*(const matrix<double>& ue, const Expanded<complex>& ex) {
        return (ex * ue);
    }

    template <class T>
    typename std::enable_if<
      isS(T), Expanded<vector<complex> > >::type
    operator*(const vector<complex>& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }
    template <class T>
    typename std::enable_if<
      isS(T), Expanded<matrix<complex> > >::type
    operator*(const matrix<complex>& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    // both vector -> commute
    
    Expanded<double> operator*(const vector<double>& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue);
    }
    template <class T, class Q>
    typename std::enable_if<
     isV(T) && isV(Q) && !(isVec(T,double) && isVec(Q,double)), Expanded<complex> >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex * ue);
    }

//  vector * matrix

    Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<matrix<double> >& ex) {
        return (ex.transpose() * ue);
    }

    template <class T, class Q>
    typename std::enable_if<
     isV(T) && isM(Q) && !(isVec(T,double) && isMat(Q,double)), Expanded<vector<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex.transpose() * ue);
    }

//  matrix * vector

    Expanded<vector<double> > operator*(const matrix<double>& ue, const Expanded<vector<double> >& ex) {
        return (ex * ue.transpose());
    }

    template <class T, class Q>
    typename std::enable_if<
     isM(T) && isV(Q) && !(isMat(T,double) && isVec(Q,double)), Expanded<vector<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return (ex * ue.transpose());
    }

//  matrix * matrix

    Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<matrix<double> >& ex) {
        return ((ex.transpose() * ue.transpose()).transpose());
    }

    template <class T, class Q>
    typename std::enable_if<
     isM(T) && isM(Q) && !(isMat(T,double) && isMat(Q,double)), Expanded<matrix<complex> > >::type
     operator*(const T& ue, const Expanded<Q>& ex) {
        return ((ex.transpose() * ue.transpose()).transpose());
    }

/***************** End UnExpanded * Expanded *****************/    

#endif /* EXPANDED_H */
