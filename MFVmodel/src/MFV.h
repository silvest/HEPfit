/* 
 * File:   MFV.h
 * Author: silvest
 *
 * Created on September 24, 2010, 10:53 AM
 */

#ifndef MFV_H
#define	MFV_H
#include <stdlib.h>
#include <gslpp_complex.h>
#include <gslpp_vector_double.h>
#include <gslpp_vector_complex.h>
#include <gslpp_matrix_double.h>
#include <gslpp_matrix_complex.h>

class MFV {
public:
    MFV(const double mQtilde, const double mUtilde, const double mDtilde,
        const double Au, const double Ad, const double mLtilde,
        const double mEtilde, const double mNtilde, const double Ae,
        const double mhu, const double mhd, const double B, const double m1,
        const double m2, const double m3, const double muH,
        const gslpp::matrix<gslpp::complex>& Yu,
        const gslpp::matrix<gslpp::complex>& Yd,
        const gslpp::matrix<gslpp::complex>& Yl,
        const gslpp::matrix<gslpp::complex>& Ye);
    MFV(const MFV& orig);
    virtual ~MFV();
private:

};

#endif	/* MFV_H */

