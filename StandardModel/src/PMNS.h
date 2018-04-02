/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PMNS_H
#define PMNS_H

#include <math.h>
#include "gslpp.h"

/**
 * @class PMNS
 * @ingroup StandardModel
 * @brief A class for the %PMNS matrix elements.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining the %PMNS matrix and its elements with
 * the Wolfenstein and Gilman parameterizations. 
 */
class PMNS {
public:
    PMNS();

    /**
     * @brief A set method to calculate the PMNS matrix from PMNS parameters
     * @param[in] s12_v the PMNS parameter @f$ \sin\theta_{12} @f$
     * @param[in] s13_v the PMNS parameter @f$ \sin\theta_{13} @f$
     * @param[in] s23_v the PMNS parameter @f$ \sin\theta_{23} @f$
     * @param[in] delta_v the PMNS parameter @f$ \delta @f$
     * @param[in] alpha21_v the PMNS parameter @f$ \alpha_{21} @f$
     * @param[in] alpha21_v the PMNS parameter @f$ \alpha_{31} @f$
     */
    void computePMNS(double s12_v, double s13_v, double s23_v, double delta_v, double alpha21_v, double alpha31_v);

    /**
     * @brief A member for returning the PMNS matrix.
     * @return the CKM matrix
     */
    gslpp::matrix<gslpp::complex> getPMNS() const
    {
        return U;
    }
    // Gilman parameterization

    /**
     * @brief A member for returning the value of the sine of the PMNS parameter @f$ \theta_{12} @f$
     * @return the value of @f$ \sin\theta_{12} @f$
     */
    double gets12() const
    {
        return s12;
    }

    /**
     * @brief A member for returning the value of the sine of the PMNS parameter @f$ \theta_{13} @f$
     * @return the value of @f$ \sin\theta_{13} @f$
     */
    double gets13() const
    {
        return s13;
    }

    /**
     * @brief A member for returning the value of the sine of the PMNS parameter @f$ \theta_{23} @f$
     * @return the value of @f$ \sin\theta_{23} @f$
     */
    double gets23() const
    {
        return s23;
    }

    /**
     * @brief A member for returning the value of the cosine of the PMNS parameter @f$ \theta_{12} @f$
     * @return the value of @f$ \cos\theta_{12} @f$
     */
    double getc12() const
    {
        return c12;
    }

    /**
     * @brief A member for returning the value of the cosine of the PMNS parameter @f$ \theta_{23} @f$
     * @return the value of @f$ \cos\theta_{23} @f$
     */
    double getc23() const
    {
        return c23;
    }

    /**
     * @brief A member for returning the value of the cosine of the PMNS parameter @f$ \theta_{13} @f$
     * @return the value of @f$ \cos\theta_{13} @f$
     */
    double getc13() const
    {
        return c13;
    }

    /**
     * @brief A member for returning the value of the PMNS parameter @f$ \delta @f$
     * @return the value of @f$ \delta @f$
     */
    double getdelta() const
    {
        return delta;
    }

     /**
     * @brief A member for returning the value of the phase @f$ \alpha_{21} @f$
     * @return the value of @f$ \alpha_{21} @f$
     */
    double getalpha21()
    {
        return alpha21;
    }

    /**
     * @brief A member for returning the value of the phase @f$ \alpha_{31} @f$
     * @return the value of @f$ \alpha_{31} @f$
     */
    double getalpha31()
    {
        return alpha31;
    }

    //Complex values of PMNS elements

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{e1} @f$
     * @return the value of @f$ U_{e1} @f$
     */
    gslpp::complex getU_e1()
    {
        return U(0, 0);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{e2} @f$
     * @return the value of @f$ U_{e2} @f$
     */
    gslpp::complex getU_e2()
    {
        return U(0, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{e3} @f$
     * @return the value of @f$ U_{e3} @f$
     */
    gslpp::complex getU_e3()
    {
        return U(0, 2);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\mu 1} @f$
     * @return the value of @f$ U_{\mu 1} @f$
     */
    gslpp::complex getU_mu1()
    {
        return U(1, 0);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\mu 2} @f$
     * @return the value of @f$ U_{\mu 2} @f$
     */
    gslpp::complex getU_mu2()
    {
        return U(1, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\mu 3} @f$
     * @return the value of @f$ U_{\mu 3} @f$
     */
    gslpp::complex getU_mu3()
    {
        return U(1, 2);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\tau 1} @f$
     * @return the value of @f$ U_{\tau 1} @f$
     */
    gslpp::complex getU_tau1()
    {
        return U(2, 0);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\tau 2} @f$
     * @return the value of @f$ U_{\tau 2} @f$
     */
    gslpp::complex getU_tau2()
    {
        return U(2, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ U_{\tau 3} @f$
     * @return the value of @f$ U_{\tau 3} @f$
     */
    gslpp::complex getU_tau3()
    {
        return U(2, 2);
    }

    
private:
    double s12, s13, s23, delta, alpha21, alpha31;
    double c12, c23, c13;

    gslpp::matrix<gslpp::complex> U;

};

#endif	/* PMNS_H */

