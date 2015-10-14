/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef XMATRICES_H
#define	XMATRICES_H

#include <gslpp.h>
#include <CKM.h>

class Xmatrices {
public:
    Xmatrices();

    void Update(const CKM &);
    
    gslpp::matrix<gslpp::complex> GetX1() const {
        return X1;
    }

    gslpp::matrix<gslpp::complex> GetX13() const {
        return X13;
    }

    gslpp::matrix<gslpp::complex> GetX2() const {
        return X2;
    }

    gslpp::matrix<gslpp::complex> GetX3() const {
        return X3;
    }

    gslpp::matrix<gslpp::complex> GetX4() const {
        return X4;
    }

    gslpp::matrix<gslpp::complex> GetX5() const {
        return X5;
    }

    gslpp::matrix<gslpp::complex> GetX6() const {
        return X6;
    }

    gslpp::matrix<gslpp::complex> GetX9() const {
        return X9;
    }
    
private:
    gslpp::matrix<gslpp::complex> X1, X2, X3, X4, X5, X6, X9, X13;
    double myrho, myeta, mylambda, myA;
};

#endif	/* XMATRICES_H */

