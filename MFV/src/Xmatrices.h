/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef XMATRICES_H
#define	XMATRICES_H

#include <gslpp_matrix_complex.h>
#include <CKM.h>

using namespace gslpp;

class Xmatrices {
public:
    Xmatrices();

    void Update(const CKM &);
    
    matrix<complex> GetX1() const {
        return X1;
    }

    matrix<complex> GetX13() const {
        return X13;
    }

    matrix<complex> GetX2() const {
        return X2;
    }

    matrix<complex> GetX3() const {
        return X3;
    }

    matrix<complex> GetX4() const {
        return X4;
    }

    matrix<complex> GetX5() const {
        return X5;
    }

    matrix<complex> GetX6() const {
        return X6;
    }

    matrix<complex> GetX9() const {
        return X9;
    }
    
private:
    matrix<complex> X1, X2, X3, X4, X5, X6, X9, X13;
    double myrho, myeta, mylambda, myA;
};

#endif	/* XMATRICES_H */

