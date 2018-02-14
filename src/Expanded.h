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
            std::vector<T> res(up1, z.Zero());
            for (uint i = 0; i < up1; i++)
                for (uint j = 0; j <= i; j++) {
                    res[i] += data1[j] * z.getOrd(i - j);
                }
            return Expanded<T>(res);
        } else {
            std::vector<std::vector<T> > res(up1, std::vector<T>(up2, z.Zero()));
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

// T double, matrix<double> (expanded), Q matrix<double>, matrix<complex> 
    template <class Q>
    typename std::enable_if<
        (std::is_same<T, double>::value ||
         std::is_same<T, matrix<double> >::value) &&
        (std::is_same<Q,  matrix<complex> >::value ||
         std::is_same<Q,  matrix<double> >::value ), Expanded<Q> >::type
    operator*(const Q& z) const {
        uint id, jd;
        if (ord2 == 0) {
            id = z.size_i();
            jd = z.size_j();
            std::vector<Q> res(ord1, Q(id, jd, 0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<Q>(res);
        } else {
            id = z.size_i();
            jd = z.size_j();
            std::vector<std::vector<Q> >
                    res(ord2, std::vector<Q>(ord1, Q(id, jd, 0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<Q>(res);
        }
    }
    
// T double (expanded), Q double, complex 
    template <class Q>
    typename std::enable_if<
        std::is_same<T, double>::value &&
        (std::is_same<Q, double>::value ||
         std::is_same<Q, complex>::value), Expanded<Q> >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<Q> res(ord1, 0.);
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<Q>(res);
        } else {
            std::vector<std::vector<Q> >
                    res(ord2, std::vector<Q>(ord1, 0.));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<Q>(res);
        }
    }

// T complex (expanded), Q matrix<double>, matrix<complex> 
    template <class Q>
    typename std::enable_if<
        std::is_same<T, complex>::value &&
        (std::is_same<Q,  matrix<complex> >::value ||
         std::is_same<Q,  matrix<double> >::value), Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        uint id, jd;
        if (ord2 == 0) {
            id = z.size_i();
            jd = z.size_j();
            std::vector<matrix<complex> > res(ord1, matrix<complex>(id, jd, 0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<matrix<complex> >(res);
        } else {
            id = z.size_i();
            jd = z.size_j();
            std::vector<std::vector<matrix<complex> > >
                    res(ord2, std::vector<matrix<complex> >(ord1, matrix<complex>(id, jd, 0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<matrix<complex> >(res);
        }
    }
    
// T complex (expanded), Q double, complex 
    template <class Q>
    typename std::enable_if<
        std::is_same<T, complex>::value &&
        (std::is_same<Q, double>::value ||
         std::is_same<Q, complex>::value), Expanded<complex> >::type
    operator*(const Q& z) const {
        if (ord2 == 0) {
            std::vector<complex> res(ord1, 0.);
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<complex>(res);
        } else {
            std::vector<std::vector<complex> >
                    res(ord2, std::vector<complex>(ord1, 0.));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<complex>(res);
        }
    }

// T matrix<double> (expanded), Q complex 
    template <class Q>
    typename std::enable_if<
        std::is_same<T, matrix<double> >::value &&
        std::is_same<Q, complex>::value, Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        uint id, jd;
        if (ord2 == 0) {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<matrix<complex> > res(ord1, matrix<complex>(id, jd, 0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<matrix<complex> >(res);
        } else {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<std::vector<matrix<complex> > >
                    res(ord2, std::vector<matrix<complex> >(ord1, matrix<complex>(id, jd, 0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<matrix<complex> >(res);
        }
    }

// T matrix<double> (expanded), Q double 
    template <class Q>
    typename std::enable_if<
        std::is_same<T, matrix<double> >::value &&
        std::is_same<Q, double>::value, Expanded<T> >::type
    operator*(const Q& z) const {
        uint id, jd;
        if (ord2 == 0) {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<T> res(ord1, T(id, jd, 0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<T>(res);
        } else {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<std::vector<T> >
                    res(ord2, std::vector<T>(ord1, T(id, jd, 0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<T>(res);
        }
    }

// T matrix<complex> (expanded) Q double, complex, matrix<double>, matrix<complex>
    template <class Q>
    typename std::enable_if<
        std::is_same<T, matrix<complex> >::value  &&
       (std::is_same<Q,  double>::value ||
        std::is_same<Q,  matrix<double> >::value ||
        std::is_same<Q,  complex >::value ||
        std::is_same<Q,  matrix<complex> >::value
       ), Expanded<T> >::type
    operator*(const Q& z) const {
        uint id, jd;
        if (ord2 == 0) {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<T> res(ord1, T(id, jd, 0.));
            for (uint i = 0; i < ord1; i++)
                    res[i] = data1[i] * z;
            return Expanded<T>(res);
        } else {
            id = data1[0].size_i();
            jd = data1[0].size_j();
            std::vector<std::vector<T> >
                    res(ord2, std::vector<T>(ord1, T(id, jd, 0.)));
            for (uint i1 = 0; i1 < ord1; i1++)
                for (uint i2 = 0; i2 < ord2; i2++)
                            res[i1][i2] = data2[i1][i2] * z;
            return Expanded<T>(res);
        }
    }
    
    Expanded<double> real() const;

    Expanded<double> abs2() const;

    Expanded<T> conjugate() const;

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
           output << "Ordine " << i << ": " << x.getOrd(i) << std::endl;
    } else {
        for (uint i=0;i<x.getN1();i++)
            for (uint j=0;j<x.getN2();j++)
                output << "Ordine (" << i << "," << j << "): " << x.getOrd(i,j) << std::endl;
    }
    return output;
}

// Operator * between unexpanded and expanded

template <class T>
Expanded<T> operator*(const double& x1, const Expanded<T>& z2) {
    return (z2 * x1);
}
// Q complex T double, complex (expanded) 
template <class T>
typename std::enable_if<
     std::is_same<T, double>::value ||
     std::is_same<T, complex>::value, Expanded<complex> >::type
operator*(const complex& x1, const Expanded<T>& z2) {
         return (z2 * x1);
}
// Q complex T matrix<double>, matrix<complex> (expanded) 
template <class T>
typename std::enable_if<
     std::is_same<T, matrix<double> >::value ||
     std::is_same<T, matrix<complex> >::value, Expanded<matrix<complex>> >::type
operator*(const complex& x1, const Expanded<T>& z2) {
    return (z2 * x1);
}
// Q matrix<double> T double (expanded)
Expanded<matrix<double> > operator*(const matrix<double>& x1, const Expanded<double>& z2) {
    return (z2 * x1);
}

// Q matrix<double> T complex (expanded)
Expanded<matrix<complex> > operator*(const matrix<double>& x1, const Expanded<complex>& z2) {
    return (z2 * x1);
}

// Q matrix<double> T matrix<double> (expanded) 
Expanded<matrix<double> > operator*(const matrix<double>& x1, const Expanded<matrix<double> >& z2) {
        uint id, jd;
        if (z2.getN2() == 0) {
            id = x1.size_i();
            jd = x1.size_j();
            std::vector<matrix<double> > res(z2.getN1(), matrix<double>(id, jd, 0.));
            for (uint i = 0; i < z2.getN1(); i++)
                    res[i] = x1 * z2.getOrd(i);
            return Expanded<matrix<double> >(res);
        } else {
            id = x1.size_i();
            jd = x1.size_j();
            std::vector<std::vector<matrix<double> > >
                    res(z2.getN2(), std::vector<matrix<double> >(z2.getN1(), matrix<double>(id, jd, 0.)));
            for (uint i1 = 0; i1 < z2.getN1(); i1++)
                for (uint i2 = 0; i2 < z2.getN2(); i2++)
                            res[i1][i2] = x1 * z2.getOrd(i1,i2);
            return Expanded<matrix<double> >(res);
        }
    }
// Q matrix<double> T matrix<complex> (expanded)
// Q matrix<complex> T matrix<double> (expanded) 
// Q matrix<complex> T matrix<complex> (expanded) 
template<class T,class Q>
typename std::enable_if<
     (std::is_same<Q, matrix<double> >::value && std::is_same<T, matrix<complex> >::value )||
     (std::is_same<Q, matrix<complex> >::value && std::is_same<T, matrix<double> >::value )||
     (std::is_same<Q, matrix<complex> >::value && std::is_same<T, matrix<complex> >::value ), Expanded<matrix<complex> > >::type
operator*(const Q& x1, const Expanded<T>& z2) {
        uint id, jd;
        if (z2.getN2() == 0) {
            id = x1.size_i();
            jd = x1.size_j();
            std::vector<matrix<complex> > res(z2.getN1(), matrix<complex>(id, jd, 0.));
            for (uint i = 0; i < z2.getN1(); i++)
                    res[i] = x1 * z2.getOrd(i);
            return Expanded<matrix<complex> >(res);
        } else {
            id = x1.size_i();
            jd = x1.size_j();
            std::vector<std::vector<matrix<complex> > >
                    res(z2.getN2(), std::vector<matrix<complex> >(z2.getN1(), matrix<complex>(id, jd, 0.)));
            for (uint i1 = 0; i1 < z2.getN1(); i1++)
                for (uint i2 = 0; i2 < z2.getN2(); i2++)
                            res[i1][i2] = x1 * z2.getOrd(i1,i2);
            return Expanded<matrix<complex> >(res);
        }
    }
// Q matrix<complex> T double complex (expanded)
template <class T>
typename std::enable_if<
     std::is_same<T, double>::value ||
     std::is_same<T, complex>::value, Expanded<matrix<complex> > >::type
operator*(const matrix<complex>& x1, const Expanded<T>& z2) {
    return z2 * x1;
}

#endif /* EXPANDED_H */
