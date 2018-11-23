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

    Expanded<T>() {
        ord1 = 0;
        ord2 = 0;
    }

    Expanded<T>(std::vector<T> & dinp) {
        ord1 = 1;
        ord2 = dinp.size();
        data.push_back(dinp);
    }

    Expanded<T>(std::vector<std::vector<T> > & dinp) {
        ord1 = dinp.size();
        ord2 = dinp[0].size();
        data = dinp;
    }

    T Zero() const {
        return T();
    }

    template <class Q>
    typename std::enable_if<
    isSc(T, double) ||
    (isSc(T, complex) && isComp(Q)) ||
    (isMat(T, double) && !isS(Q)) ||
    (isMat(T, complex) && (isVec(Q, complex) || isMat(Q, complex))), Expanded<Q> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<Q>(res);
    }

    template <class Q>
    typename std::enable_if<
    (!isSc(T, double) && isSc(Q, double)) ||
    (!isS(T) && !isMat(T, double) && isMat(Q, double)) ||
    ((isMat(T, complex) || isVec(T, complex)) && isSc(Q, complex)) ||
    (isVec(T, complex) && isMat(Q, complex)), Expanded<T> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());
        std::vector<std::vector<T> > res(up1, std::vector<T>(up2, Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<T>(res);
    }

    template <class Q>
    typename std::enable_if<
    (isVec(T, double) && isSc(Q, complex)) ||
    (isVec(T, double) && isMat(Q, complex)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<vector<complex> > > res(up1, std::vector<vector<complex> >(up2, Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<vector<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    (isSc(T, complex) && isVec(Q, double)) ||
    (isMat(T, complex) && isVec(Q, double)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<vector<complex> > > res(up1, std::vector<vector<complex> >(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<vector<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    (isVec(T, double) && isVec(Q, complex)) ||
    (isVec(T, complex) && isVec(Q, double)) ||
    (isVec(T, complex) && isVec(Q, complex)), Expanded<complex> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<complex> > res(up1, std::vector<complex>(up2, 0.));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<complex>(res);
    }

    template <class Q>
    typename std::enable_if<
    (isMat(T, double) && isSc(Q, complex)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<matrix<complex> > > res(up1, std::vector<matrix<complex> >(up2, Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<matrix<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    (isSc(T, complex) && isMat(Q, double)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<matrix<complex> > > res(up1, std::vector<matrix<complex> >(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<matrix<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && isVec(Q, double), Expanded<double> >::type
    operator*(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<double> > res(up1, std::vector<double>(up2, 0.));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                for (uint j1 = 0; j1 <= i1; j1++)
                    for (uint j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i1 - j1, i2 - j2);
        return Expanded<double>(res);
    }

    // Operator +

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && isComp(T), Expanded<T> >::type
    operator+(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<T> > res(up1, std::vector<T>(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                res[i1][i2] = data[i1][i2] + z.getOrd(i1, i2);
        return Expanded<T>(res);
    }

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && !isComp(T), Expanded<Q> >::type
    operator+(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                res[i1][i2] = data[i1][i2] + z.getOrd(i1, i2);
        return Expanded<Q>(res);
    }

    // Operator -

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && isComp(T), Expanded<T> >::type
    operator-(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());
        std::vector<std::vector<T> > res(up1, std::vector<T>(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                res[i1][i2] = data[i1][i2] - z.getOrd(i1, i2);
        return Expanded<T>(res);
    }

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && !isComp(T), Expanded<Q> >::type
    operator-(const Expanded<Q>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        std::vector<std::vector<Q> > res(up1, std::vector<Q>(up2, z.Zero()));
        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                res[i1][i2] = data[i1][i2] - z.getOrd(i1, i2);
        return Expanded<Q>(res);
    }

    Expanded<T> operator-() const {
        std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = -data[i1][i2];
        return Expanded<T>(res);
    }

    bool operator==(const Expanded<T>& z) const {
        uint up1 = std::min(ord1, z.getN1());
        uint up2 = std::min(ord2, z.getN2());

        if(up1 != getN1() || up2 != getN2())
            throw std::runtime_error("Expanded::operator==: cannot compare objects expanded to diffent orders.");

        for (uint i1 = 0; i1 < up1; i1++)
            for (uint i2 = 0; i2 < up2; i2++)
                if(z.getOrd(i1, i2) != getOrd(i1, i2)) return(false);

        return true;
    }

    bool operator!=(const Expanded<T>& z) const {
        return !(z == *this);
    }

    /***************** End Expanded-Expanded *****************/

    /***************** Begin Expanded * UnExpanded *****************/

    // Return T.

    template <class Q>
    typename std::enable_if<
    (std::is_same<T, Q>::value && !isV(T)) ||
    isSc(Q, double) ||
    (!isS(T) && isMat(Q, double)) ||
    ((isMat(T, complex) || isVec(T, complex)) && isS(Q)) ||
    (isVec(T, complex) && isMat(Q, complex)) ||
    (isV(T) && isSc(Q, double)) ||
    (isMat(T, double) && isSc(Q, double)), Expanded<T> >::type
    operator*(const Q& z) const {
        std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<T>(res);
    }

    // Return double

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && isVec(Q, double), Expanded<double> >::type
    operator*(const Q& z) const {
        std::vector<std::vector<double> > res(ord1, std::vector<double>(ord2, 0.));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<double>(res);
    }

    // Return complex

    template <class Q>
    typename std::enable_if<
    (isSc(T, double) && isSc(Q, complex)) ||
    (isV(T) && isV(Q) && !(isVec(T, double) && isVec(Q, double))), Expanded<complex> >::type
    operator*(const Q& z) const {
        std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2, 0.));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<complex>(res);
    }
    // Return vector<double>

    template <class Q>
    typename std::enable_if<
    isVec(Q, double) && (isSc(T, double) || isMat(T, double)), Expanded<vector<double> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<double> > > res(ord1, std::vector<vector<double> >(ord2, vector<double>(z.size(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<vector<double> >(res);
    }


    // Return vector<complex>

    template <class Q>
    typename std::enable_if<
    (isVec(Q, complex) && (isS(T) || isM(T))) ||
    (isVec(Q, double) && (isSc(T, complex) || isMat(T, complex))), Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(z.size(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<vector<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && (isSc(Q, complex) || isMat(Q, complex)), Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(data[0][0].size(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<vector<complex> >(res);
    }

    // Return matrix<double>

    template <class Q>
    typename std::enable_if<
    (isSc(T, double) && isMat(Q, double)), Expanded<matrix<double> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<double> > > res(ord1,
                std::vector<matrix<double> >(ord2, matrix<double>(z.size_i(), z.size_j(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<matrix<double> >(res);
    }

    // Return matrix<complex>

    template <class Q>
    typename std::enable_if<
    (isS(T) && isMat(Q, complex)) ||
    (isSc(T, complex) && isMat(Q, double)) ||
    (isMat(T, double) && isMat(Q, complex)), Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res(ord1,
                std::vector<matrix<complex> >(ord2, matrix<complex>(z.size_i(), z.size_j(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<matrix<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<
    isMat(T, double) && isSc(Q, complex), Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] * z;
        return Expanded<matrix<complex> >(res);
    }

    /***************** End Expanded * UnExpanded *****************/

    /***************** Begin Expanded / UnExpanded-Scalar *****************/


    Expanded<T> operator/(const double& z) const {
        std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] / z;
        return Expanded<T>(res);
    }

    template <class Q>
    typename std::enable_if<isComp(T) && isSc(Q, complex), Expanded<T> >::type
    operator/(const Q& z) const {
        std::vector<std::vector<T> > res(ord1, std::vector<T>(ord2, Zero()));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] / z;
        return Expanded<T>(res);
    }

    template <class Q>
    typename std::enable_if<isSc(T, double) && isSc(Q, complex), Expanded<complex > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<complex> > res(ord1, std::vector<complex>(ord2, 0.));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] / z;
        return Expanded<complex>(res);
    }

    template <class Q>
    typename std::enable_if<isVec(T, double) && isSc(Q, complex), Expanded<vector<complex> > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res(ord1, std::vector<vector<complex> >(ord2, vector<complex>(data[0][0].size(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] / z;
        return Expanded<vector<complex> >(res);
    }

    template <class Q>
    typename std::enable_if<isMat(T, double) && isSc(Q, complex), Expanded<matrix<complex> > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res(ord1, std::vector<matrix<complex> >(ord2, matrix<complex>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
        for (uint i1 = 0; i1 < ord1; i1++)
            for (uint i2 = 0; i2 < ord2; i2++)
                res[i1][i2] = data[i1][i2] / z;
        return Expanded<matrix<complex> >(res);
    }

    /***************** End Expanded / UnExpanded-Scalar  *****************/
    Expanded<double> real() const;

    Expanded<double> abs2() const;

    Expanded<T> conjugate() const;

    Expanded<T> transpose() const;

    T& getAllOrd() {
        T res();
        for (uint i = 0; i < ord1; i++)
            for (uint j = 0; j < ord2; j++)
                res += data[i][j];
        return res;
    }

    T getOrd(uint j) const {
        return data[0][j];
    }
    
    const T& getOrd(uint i, uint j) const {
        return data[i][j];
    }

    void setOrd(uint j, T value) {
        data[0][j] = value;
    }

    void setOrd(uint i, uint j, T value) {
        data[i][j] = value;
    }
    
    template <class Q> typename std::enable_if<isMat(T, Q), void>::type 
    setSingleElem(uint i,uint j,uint h, uint k, Q x ) {
        data[i][j].assign(h,k,x);
    }

    template <class Q> typename std::enable_if<isVec(T, Q), void>::type 
    setSingleElem(uint i,uint j,uint h, Q x ) {
        data[i][j].assign(h,x);
    }

    const uint getN1() const {
        return ord1;
    }

    const uint getN2() const {
        return ord2;
    }

    friend std::ostream& operator<<(std::ostream& output, const Expanded<T> z) {
        for (uint i = 0; i < z.getN1(); i++)
            for (uint j = 0; j < z.getN2(); j++) {
                output << "Order (" << i << "," << j << "): " << z.getOrd(i, j);
                if (i != z.getN1() - 1 || j != z.getN2() - 1)
                    output << std::endl;
            }
        return output;
    }
    friend Expanded<T> operator*(const double& ue, const Expanded<T>& ex) {
        return(ex * ue);
    }

    private:
    uint ord1, ord2;
    std::vector<std::vector<T> > data;
};

// Template specialization

//  scalar  * expanded

template <class T>
typename std::enable_if<
isComp(T), Expanded<T> >::type
operator*(const complex& ue, const Expanded<T>& ex) {
    return (ex * ue);
}


//  vector  * expanded scalar

template <class T>
typename std::enable_if<
isS(T), Expanded<vector<complex> > >::type
operator*(const vector<complex>& ue, const Expanded<T>& ex) {
    return (ex * ue);
}

//  vector * vector

template <class T, class Q>
typename std::enable_if<
isV(T) && isV(Q) && !(isVec(T, double) && isVec(Q, double)), Expanded<complex> >::type
operator*(const T& ue, const Expanded<Q>& ex) {
    return (ex * ue);
}

//  vector * matrix

template <class T, class Q>
typename std::enable_if<
isV(T) && isM(Q) && !(isVec(T, double) && isMat(Q, double)), Expanded<vector<complex> > >::type
operator*(const T& ue, const Expanded<Q>& ex) {
    return (ex.transpose() * ue);
}

// matrix * expanded vector

template <class T, class Q>
typename std::enable_if<
isM(T) && isV(Q) && !(isMat(T, double) && isVec(Q, double)), Expanded<vector<complex> > >::type
operator*(const T& ue, const Expanded<Q>& ex) {
    return (ex * ue.transpose());
}

// matrix * expanded matrix

template <class T, class Q>
typename std::enable_if<
isM(T) && isM(Q) && !(isMat(T, double) && isMat(Q, double)), Expanded<matrix<complex> > >::type
operator*(const T& ue, const Expanded<Q>& ex) {
    return ((ex.transpose() * ue.transpose()).transpose());
}
// matrix * expanded scalar

template <class T>
typename std::enable_if<
isS(T), Expanded<matrix<complex> > >::type
operator*(const matrix<complex>& ue, const Expanded<T>& ex) {
    return (ex * ue);
}

template <> Expanded<complex> Expanded<complex>::conjugate() const;

template <> Expanded<matrix<double> > Expanded<matrix<double> >::transpose() const;

template <> Expanded<matrix<complex> > Expanded<matrix<complex> >::transpose() const;

template <> Expanded<double> Expanded<complex>::real() const;

template <> Expanded<double> Expanded<complex>::abs2() const;

template <> matrix<double> Expanded<matrix<double> >::Zero() const;

template <> matrix<complex> Expanded<matrix<complex> >::Zero() const;

template <> vector<double> Expanded<vector<double> >::Zero() const;

template<> vector<complex> Expanded<vector<complex> >::Zero() const;

/***************** Begin UnExpanded * Expanded *****************/

Expanded<complex> operator*(const complex& ue, const Expanded<double>& ex);

Expanded<vector<complex> > operator*(const complex& ue, const Expanded<vector<double> >& ex);

Expanded<matrix<complex> > operator*(const complex& ue, const Expanded<matrix<double> >& ex);

// X * scalar -> commute

Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<double>& ex);

Expanded<vector<complex> > operator*(const vector<double>& ue, const Expanded<complex>& ex);

Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<double>& ex);

Expanded<matrix<complex> > operator*(const matrix<double>& ue, const Expanded<complex>& ex);

// both vector -> commute

Expanded<double> operator*(const vector<double>& ue, const Expanded<vector<double> >& ex);

//  vector * matrix

Expanded<vector<double> > operator*(const vector<double>& ue, const Expanded<matrix<double> >& ex);

//  matrix * vector

Expanded<vector<double> > operator*(const matrix<double>& ue, const Expanded<vector<double> >& ex);

//  matrix * matrix

Expanded<matrix<double> > operator*(const matrix<double>& ue, const Expanded<matrix<double> >& ex);

/***************** End UnExpanded * Expanded *****************/

#endif /* EXPANDED_H */
