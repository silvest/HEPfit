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

#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "gslpp_complex.h"
#include "gslpp_vector_double.h"
#include "gslpp_vector_complex.h"
#include "gslpp_matrix_complex.h"
#include "gslpp_matrix_double.h"

using namespace gslpp;

#define isMat(X,M) (std::is_same<X,matrix<M> >::value)
#define isVec(X,V) (std::is_same<X,vector<V> >::value)
#define isSc(X,S)  (std::is_same<X,S>::value)
#define isS(X)   (isSc(X,double) || isSc(X,complex))
#define isV(X)   (isVec(X,double) || isVec(X,complex))
#define isM(X)   (isMat(X,double) || isMat(X,complex))
#define isComp(X) (isSc(X,complex) || isVec(X,complex) || isMat(X,complex))
#define sameType(X,Y) ( (isS(X) && isS(Y)) || (isV(X) && isV(Y)) || (isM(X) && isM(Y)))

/**
 * @class Expanded
 * @ingroup gslpp
 * @brief A template class for Taylor double expansion of several objects
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a template class that implement an expansion of the
 * object T in two generic expansion parameters. T can be double, gslpp::complex,
 * gslpp::vector<double>, gslpp::vector<gslpp::complex>, gslpp::matrix<double>,
 * and gslpp::matrix<gslpp::comlex>. The elements of the expansion are stored in
 * a private std::vector<std::vector<T> >. The variables int n1 and
 * std::vector<int> n2 contain the size of the outer and inner vectors
 * respectively. The variables minord1 and minord2 contain the starting powers
 * of the expansion in the two parameters. Operators +, -, *, /, ==, !=, << are
 * defined and the proper number of terms is retained in the result. 
 */
template <class T> class Expanded {
public:
    /**
     * @brief Empty constructor. All data initialized to zero.
     */
    Expanded<T>();

    /**
     * @brief Constructor of a single expansion.
     * @param[in] dinp a std::vector<T> with the elements of the expansion.
     * @param[in] minord2_i an int containing the starting power of the
     * expansion, default to 0. 
     */
    Expanded<T>(std::vector<T> & dinp, int minord2_i = 0);

    /**
     * @brief Constructor of a double expansion.
     * @param[in] dinp a std::vector<std::vector<T> > with the elements of the expansion.
     * @param[in] minord2_i an int containing the starting power of the inner
     * expansion, default to 0. 
     * @param[in] minord1_i an int containing the starting power of the outer
     * expansion, default to 0. 
     */
    Expanded<T>(std::vector<std::vector<T> > & dinp, int minord2_i = 0, int minord1_i = 0);

    template <class Q>
    typename std::enable_if<
    isSc(T, double) ||
    (isSc(T, complex) && isComp(Q)) ||
    (isMat(T, double) && !isS(Q)) ||
    (isMat(T, complex) && (isVec(Q, complex) || isMat(Q, complex))), Expanded<Q> >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<Q> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<Q>(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<Q>(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (!isSc(T, double) && isSc(Q, double)) ||
    (!isS(T) && !isMat(T, double) && isMat(Q, double)) ||
    ((isMat(T, complex) || isVec(T, complex)) && isSc(Q, complex)) ||
    (isVec(T, complex) && isMat(Q, complex)), Expanded<T> >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<T>(up2, Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<T>(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (isVec(T, double) && isSc(Q, complex)) ||
    (isVec(T, double) && isMat(Q, complex)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<vector<complex> > > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<vector<complex> >(up2, Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<vector<complex> >(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (isSc(T, complex) && isVec(Q, double)) ||
    (isMat(T, complex) && isVec(Q, double)), Expanded<vector<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<vector<complex> > > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<vector<complex> >(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<vector<complex> >(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (isVec(T, double) && isVec(Q, complex)) ||
    (isVec(T, complex) && isVec(Q, double)) ||
    (isVec(T, complex) && isVec(Q, complex)), Expanded<complex> >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<complex> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<complex>(up2, 0.));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<complex>(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (isMat(T, double) && isSc(Q, complex)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<matrix<complex> > > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<matrix<complex> >(up2, Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<matrix<complex> >(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    (isSc(T, complex) && isMat(Q, double)), Expanded<matrix<complex> > >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<matrix<complex> > > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<matrix<complex> >(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<matrix<complex> >(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && isVec(Q, double), Expanded<double> >::type
    operator*(const Expanded<Q>& z) const {
        int min1 = minord1 + z.getMin1();
        int min2 = minord2 + z.getMin2();
        int up1 = std::min(n1, z.getN1());
        std::vector<std::vector<double> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1), z.getN2().at(i1));
            res.push_back(std::vector<double>(up2, 0.));
            for (int i2 = 0; i2 < up2; i2++)
                for (int j1 = 0; j1 <= i1; j1++)
                    for (int j2 = 0; j2 <= i2; j2++)
                        res[i1][i2] += data[j1][j2] * z.getOrd(i2 - j2 + z.getMin2(), i1 - j1 + z.getMin1());
        }
        return Expanded<double>(res, min2, min1);
    }

    // Operator +

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && isComp(T), Expanded<T> >::type
    operator+(const Expanded<Q>& z) const {
        int min1 = std::min(minord1, z.getMin1());
        int min2 = std::min(minord2, z.getMin2());
        int up1 = std::min(n1 + minord1, z.getN1() + z.getMin1()) - min1;
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1) + minord2, z.getN2().at(i1) + z.getMin2()) - min2;
            res.push_back(std::vector<T>(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++) {
                if (i1 + min1 >= minord1 && i2 + min2 >= minord2) res[i1][i2] += data[i1 + min1 - minord1][i2 + min2 - minord2];
                if (i1 + min1 >= z.getMin1() && i2 + min2 >= z.getMin2()) res[i1][i2] += z.getOrd(i2 + min2, i1 + min1);
            }
        }
        return Expanded<T>(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && !isComp(T), Expanded<Q> >::type
    operator+(const Expanded<Q>& z) const {
        int min1 = std::min(minord1, z.getMin1());
        int min2 = std::min(minord2, z.getMin2());
        int up1 = std::min(n1 + minord1, z.getN1() + z.getMin1()) - min1;
        std::vector<std::vector<Q> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1) + minord2, z.getN2().at(i1) + z.getMin2()) - min2;
            res.push_back(std::vector<Q>(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++) {
                if (i1 + min1 >= minord1 && i2 + min2 >= minord2) res[i1][i2] += data[i1 + min1 - minord1][i2 + min2 - minord2];
                if (i1 + min1 >= z.getMin1() && i2 + min2 >= z.getMin2()) res[i1][i2] += z.getOrd(i2 + min2, i1 + min1);
            }
        }
        return Expanded<Q>(res, min2, min1);
    }

    // Operator -

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && isComp(T), Expanded<T> >::type
    operator-(const Expanded<Q>& z) const {
        int min1 = std::min(minord1, z.getMin1());
        int min2 = std::min(minord2, z.getMin2());
        int up1 = std::min(n1 + minord1, z.getN1() + z.getMin1()) - min1;
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1) + minord2, z.getN2().at(i1) + z.getMin2()) - min2;
            res.push_back(std::vector<T>(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++) {
                if (i1 + min1 >= minord1 && i2 + min2 >= minord2) res[i1][i2] += data[i1 + min1 - minord1][i2 + min2 - minord2];
                if (i1 + min1 >= z.getMin1() && i2 + min2 >= z.getMin2()) res[i1][i2] -= z.getOrd(i2 + min2, i1 + min1);
            }
        }
        return Expanded<T>(res, min2, min1);
    }

    template <class Q>
    typename std::enable_if<
    sameType(T, Q) && !isComp(T), Expanded<Q> >::type
    operator-(const Expanded<Q>& z) const {
        int min1 = std::min(minord1, z.getMin1());
        int min2 = std::min(minord2, z.getMin2());
        int up1 = std::min(n1 + minord1, z.getN1() + z.getMin1()) - min1;
        std::vector<std::vector<Q> > res;
        for (int i1 = 0; i1 < up1; i1++) {
            int up2 = std::min(n2.at(i1) + minord2, z.getN2().at(i1) + z.getMin2()) - min2;
            res.push_back(std::vector<Q>(up2, z.Zero()));
            for (int i2 = 0; i2 < up2; i2++) {
                if (i1 + min1 >= minord1 && i2 + min2 >= minord2) res[i1][i2] += data[i1 + min1 - minord1][i2 + min2 - minord2];
                if (i1 + min1 >= z.getMin1() && i2 + min2 >= z.getMin2()) res[i1][i2] -= z.getOrd(i2 + min2, i1 + min1);
            }
        }
        return Expanded<Q>(res, min2, min1);
    }

    Expanded<T> inverse() const {
        throw std::runtime_error("Expanded::inverse: method not implemented");
    }

    Expanded<T> operator-() const {
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<T>(n2.at(i1), Zero()));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = -data[i1][i2];
        }
        return Expanded<T>(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if< isSc(Q,double) , Expanded<T> >::type
    operator/(const Expanded<Q>& z) const {
        return (*this)*z.inverse();
    }

    template <class Q>
    typename std::enable_if< isSc(Q,complex) && isS(T) , Expanded<complex> >::type
    operator/(const Expanded<Q>& z) const {
        return (*this)*z.inverse();
    }

    template <class Q>
    typename std::enable_if< isSc(Q,complex) && isV(T) , Expanded<vector<complex> > >::type
    operator/(const Expanded<Q>& z) const {
        return (*this)*z.inverse();
    }

    template <class Q>
    typename std::enable_if< isSc(Q,complex) && isM(T) , Expanded<matrix<complex> > >::type
    operator/(const Expanded<Q>& z) const {
        return (*this)*z.inverse();
    }

    bool operator==(const Expanded<T>& z) const {

        if (minord1 != z.getMin1() || minord2 != z.getMin2() || n1 != z.getN1()) return false;

        for (int i1 = 0; i1 < n1; i1++)
        {
            if(n2.at(i1) != z.getN2().at(i1)) return false;
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                if (z.getOrd(i2 + z.getMin2(), i1 + z.getMin1()) != data[i1][i2]) return false;
        }

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
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<T>(n2.at(i1), Zero()));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<T>(res, minord2, minord1);
    }

    // Return double

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && isVec(Q, double), Expanded<double> >::type
    operator*(const Q& z) const {
        std::vector<std::vector<double> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<double>(n2.at(i1), 0.));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<double>(res, minord2, minord1);
    }

    // Return complex

    template <class Q>
    typename std::enable_if<
    (isSc(T, double) && isSc(Q, complex)) ||
    (isV(T) && isV(Q) && !(isVec(T, double) && isVec(Q, double))), Expanded<complex> >::type
    operator*(const Q& z) const {
        std::vector<std::vector<complex> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<complex>(n2.at(i1), 0.));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<complex>(res, minord2, minord1);
    }
    // Return vector<double>

    template <class Q>
    typename std::enable_if<
    isVec(Q, double) && (isSc(T, double) || isMat(T, double)), Expanded<vector<double> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<double> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<vector<double> >(n2.at(i1), vector<double>(z.size(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<vector<double> >(res, minord2, minord1);
    }


    // Return vector<complex>

    template <class Q>
    typename std::enable_if<
    (isVec(Q, complex) && (isS(T) || isM(T))) ||
    (isVec(Q, double) && (isSc(T, complex) || isMat(T, complex))), Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<vector<complex> >(n2.at(i1), vector<complex>(z.size(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<vector<complex> >(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<
    isVec(T, double) && (isSc(Q, complex) || isMat(Q, complex)), Expanded<vector<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<vector<complex> >(n2.at(i1), vector<complex>(data[0][0].size(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<vector<complex> >(res, minord2, minord1);
    }

    // Return matrix<double>

    template <class Q>
    typename std::enable_if<
    (isSc(T, double) && isMat(Q, double)), Expanded<matrix<double> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<double> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<matrix<double> >(n2.at(i1), matrix<double>(z.size_i(), z.size_j(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<matrix<double> >(res, minord2, minord1);
    }

    // Return matrix<complex>

    template <class Q>
    typename std::enable_if<
    (isS(T) && isMat(Q, complex)) ||
    (isSc(T, complex) && isMat(Q, double)) ||
    (isMat(T, double) && isMat(Q, complex)), Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<matrix<complex> >(n2.at(i1), matrix<complex>(z.size_i(), z.size_j(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<matrix<complex> >(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<
    isMat(T, double) && isSc(Q, complex), Expanded<matrix<complex> > >::type
    operator*(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<matrix<complex> >(n2.at(i1), matrix<complex>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] * z;
        }
        return Expanded<matrix<complex> >(res, minord2, minord1);
    }

    /***************** End Expanded * UnExpanded *****************/

    /***************** Begin Expanded / UnExpanded-Scalar *****************/


    Expanded<T> operator/(const double& z) const {
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<T>(n2.at(i1), Zero()));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] / z;
        }
        return Expanded<T>(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<isComp(T) && isSc(Q, complex), Expanded<T> >::type
    operator/(const Q& z) const {
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<T>(n2.at(i1), Zero()));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] / z;
        }
        return Expanded<T>(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<isSc(T, double) && isSc(Q, complex), Expanded<complex > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<complex> > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<complex>(n2.at(i1), 0.));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] / z;
        }
        return Expanded<complex>(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<isVec(T, double) && isSc(Q, complex), Expanded<vector<complex> > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<vector<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<vector<complex> >(n2.at(i1), vector<complex>(data[0][0].size(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] / z;
        }
        return Expanded<vector<complex> >(res, minord2, minord1);
    }

    template <class Q>
    typename std::enable_if<isMat(T, double) && isSc(Q, complex), Expanded<matrix<complex> > >::type
    operator/(const Q& z) const {
        std::vector<std::vector<matrix<complex> > > res;
        for (int i1 = 0; i1 < n1; i1++) {
            res.push_back(std::vector<matrix<complex> >(n2.at(i1), matrix<complex>(data[0][0].size_i(), data[0][0].size_j(), 0.)));
            for (int i2 = 0; i2 < n2.at(i1); i2++)
                res[i1][i2] = data[i1][i2] / z;
        }
        return Expanded<matrix<complex> >(res, minord2, minord1);
    }

   /***************** End Expanded / UnExpanded-Scalar  *****************/

    friend Expanded<T> operator*(const double& ue, const Expanded<T>& ex) {
        return (ex * ue);
    }

    /**
     * @brief Method to return the real part of a complex Expanded.
     * @return the real part of this Expanded
     */
    Expanded<double> real() const;

    /**
     * @brief Method to return the squared absolute value of a Expanded<complex>.
     * @return the squared absolute value of this Expanded<complex>
     */
    Expanded<double> abs2() const;

    /**
     * @brief Method to return the conjugate value of a Expanded<complex>.
     * @return the conjugate of this Expanded<complex>
     */
    Expanded<T> conjugate() const;

    /**
     * @brief Method to transpose a Expanded<matrix<double> > or a
     * Expanded<matrix<complex> >.
     * @return the transpose of this Expanded<matrix<*> >.
     */
    Expanded<T> transpose() const;

    /**
     * @brief Method to truncate a double expansion.
     * @param[in] n2max a std::vector<int> with the maximum order of the inner
     * expansion. n1max + 1 entries have to be provided.
     * @param[in] n1max a int with the maximum order of the outer expansion.
     * @return a new Expanded<T> with the truncated expansion.
     */
    Expanded<T> truncate(std::vector<int> n2max, int n1max) {
        if (n1max < minord1 || n1max >= minord1 + n1)
            throw std::runtime_error("Expanded::truncate(): order of the outer expansion not present in Expanded");
        std::vector<std::vector<T> > res;
        for (int i1 = 0; i1 <= n1max - minord1; i1++) {
            std::vector<T> tmp;
            if (n2max.at(i1) >= minord2 + n2.at(i1))
                throw std::runtime_error("Expanded::truncate(): order of the inner expansion not present in Expanded");
            for (int j1 = 0; j1 <= n2max.at(i1) - minord2; j1++)
                tmp.push_back(data[i1][j1]);
            res.push_back(tmp);
        }
        return Expanded<T>(res, minord2, minord1);
    }
    /**
     * @brief Method to truncate a single expansion.
     * @param[in] n2max a int with the maximum order of the expansion.
     * @return a new Expanded<T> with the truncated expansion.
     */
    Expanded<T> truncate(int n2max) {
        if(n1 >= 1 || minord1 > 0)
             throw std::runtime_error("Expanded::truncate(): wrong method, please use truncate(std::vector<int>, int)");
        return truncate(std::vector<int>(1, n2max), 0);
    }
    
    /**
     * @brief Method to sum a double expansion up to a given order.
     * @param[in] ord2 a std::vector<int> with the maximum order of the inner expansion in the sum.
     * @param[in] ord1 a int with the maximum order of the outer expansion in the sum.
     * @param[in] als2 a double with the value of the inner expansion parameter. Default to 1. 
     * @param[in] als1 a double with the value of the outer expansion parameter. Default to 1.
     * @return a instance of a class T containing the expansion summed to the requested orders.
     */
    T Series(std::vector<int> ord2, int ord1, double als2 = 1., double als1 = 1.) {
        T res;
        if (ord1 >= minord1)
        {
            for (int i = 0; i <= std::min(n1, ord1 - minord1); i++)
                {
                if (ord2.at(i) >= minord2)
                    for (int j = 0; j <= std::min(ord2.at(i) - minord2, n2.at(i)); j++)
                        res += data[i][j] * pow(als1, (double) (minord1 + i)) * pow(als2, (double) (minord2 + j));
                }
        }
        return res;
    }

    /**
     * @brief Method to sum an expansion up to a given order.
     * @param[in] ord2 a int with the maximum order of the expansion in the sum.
     * @param[in] als2 a double with the value of the expansion parameter. Default to 1. 
     * @return a instance of a class T containing the expansion summed to the requested orders.
     */
    T Series(int ord2, double als2 = 1.) {
        return Series(std::vector<int> (1, ord2), 0, als2);
    }

    /**
     * @brief Return an empty instance of class T.
     * @return An empty instance of class T.
     */
    T Zero() const {
        return T();
    }

    /**
     * @brief Get an element of a single expansion.
     * @param[in] j order of the requested element.
     * @return the j-th element of the expansion.
     */
    const T& getOrd(int j) const {
        return getOrd(j, minord1);
    }

    /**
     * @brief Get an element of a double expansion.
     * @param[in] j order of the requested element in the inner expansion.
     * @param[in] i order of the requested element in the outer expansion.
     * @return the (i-th, j-th) element of the double expansion.
     */
    const T& getOrd(int j, int i) const {
        checkOrd(j, i, "Expanded::getOrd(): order non present in Expanded");
        return data[i - minord1][j - minord2];
    }

    /**
     * @brief Set an element of a single expansion.
     * @param[in] j order of the element to be set.
     */
    void setOrd(int j, T value) {
        setOrd(j, minord1, value);
    }

    /**
     * @brief Set an element of a double expansion.
     * @param[in] j order of the element to be set in the inner expansion.
     * @param[in] i order of the element to be set in the outer expansion.
     */
    void setOrd(int j, int i, T value) {
        checkOrd(j, i, "Expanded::setOrd(): order non present in Expanded");
        data[i - minord1][j - minord2] = value;
    }

    void checkOrd(int j, int i, std::string s) const {
        if (i < minord1 || i >= minord1 + n1 || j < minord2 || j >= minord2 + n2.at(i - minord1))
            throw std::runtime_error(s);
    }
 
    /**
     * @brief Set an element of a matrix in a double Expanded<matrix<*> >
     * @param[in] j order in the inner expansion of the matrix whose element has to be set.
     * @param[in] i order in the outer expansion of the matrix whose element has to be set.
     * @param[in] h row index of the element to be set in the matrix.
     * @param[in] k column index of the element to be set in the matrix.
     * @param[in] x value to be assigned.
     */
    template <class Q> typename std::enable_if<isMat(T, Q), void>::type
    setMatrixElement(int j, int i, int h, int k, Q x) {
        checkOrd(j, i, "Expanded::setMatrixElement(): order non present in Expanded");
        data[i - minord1][j - minord2].assign(h, k, x);
    }
    /**
     * @brief Set an element of a matrix in a single Expanded<matrix<*> >
     * @param[in] j order in the inner expansion of the matrix whose element has to be set.
     * @param[in] h row index of the element to be set in the matrix.
     * @param[in] k column index of the element to be set in the matrix.
     * @param[in] x value to be assigned.
     */
    template <class Q> typename std::enable_if<isMat(T, Q), void>::type
    setMatrixElement(int j, int h, int k, Q x) {
        setMatrixElement(j, minord1, h, k, x);
    }

    /**
     * @brief Set an element of a matrix in a double Expanded<vector<*> >
     * @param[in] j order in the inner expansion of the vector whose element has to be set.
     * @param[in] i order in the outer expansion of the vector whose element has to be set.
     * @param[in] h index of the element to be set in the vector.
     * @param[in] x value to be assigned.
     */
    template <class Q> typename std::enable_if<isVec(T, Q), void>::type
    setVectorElement(int j, int i, int h, Q x) {
        checkOrd(j, i, "Expanded::setVectorElement(): order non present in Expanded");
        //THERE IS NO ASSIGN METHOD FOR VECTORS
        //data[i - minord1][j - minord2].assign(h, x);
        data[i - minord1][j - minord2](h) = x;
    }

    /**
     * @brief Set an element of a matrix in a single Expanded<vector<*> >
     * @param[in] j order in the inner expansion of the vector whose element has to be set.
     * @param[in] h index of the element to be set in the vector.
     * @param[in] x value to be assigned.
     */
    template <class Q> typename std::enable_if<isVec(T, Q), void>::type
    setVectorElement(int j, int h, Q x) {
        setVectorElement(j, minord1, h, x);
    }

    /**
     * @brief Get the minimum order of the outer expansion.
     * @return value of the minimum order of the outer expansion.
     */
    int getMin1() const {
        return minord1;
    }

    /**
     * @brief Get the minimum order of the inner expansion.
     * @return value of the minimum order of the inner expansion.
     */
    int getMin2() const {
        return minord2;
    }

    /**
     * @brief Get the number of terms in the outer expansion
     * @return value of the minimum order of the outer expansion.
     */
    int getN1() const {
        return n1;
    }

    /**
     * @brief Get the number of terms in the inner expansion
     * @return values of the minimum order of the inner expansion, contains n1
     * terms.
     */
    const std::vector<int>& getN2() const {
        return n2;
    }

    friend std::ostream& operator<<(std::ostream& output, const Expanded<T> z) {
        for (int i = 0; i < z.getN1(); i++) 
            for (int j = 0; j < z.getN2().at(i); j++) {
                output << "Order (" << z.getMin2() + j << "," << z.getMin1() + i << "): " << z.getOrd(j + z.getMin2(), i + z.getMin1());
                if (i != (z.getN1() - 1) || j != (z.getN2().at(i) - 1))
                    output << std::endl;
                }
        return output;
    }


private:
    /**
     * @brief Return the (i-th, j-th) coefficient of the expansion of 1/Expanded.
     * @return the (i-th, j-th) coefficient of the expansion of 1/Expanded.
     */
    T invCoeff(int i, int j) const {
        throw std::runtime_error("Expanded::invCoeff: method not implemented");
    }

    int minord1, minord2, n1;
    std::vector<int> n2;
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
