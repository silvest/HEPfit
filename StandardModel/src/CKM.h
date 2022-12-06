/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CKM_H
#define CKM_H

#include <math.h>
#include "gslpp.h"

/**
 * @class CKM
 * @ingroup StandardModel
 * @brief A class for the %CKM matrix elements.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining the %CKM matrix and its elements with
 * the Wolfenstein and Gilman parameterizations. 
 */
class CKM {
public:
    CKM();

    /**
     * @brief A set method to calculate the CKM matrix from Wolfenstein parameters.
     * @param[in] Lambda_v the Wolfenstein parameter @f$ \lambda @f$
     * @param[in] A_v the Wolfenstein parameter @f$ A @f$
     * @param[in] Rho_v the Wolfenstein parameter @f$ \bar{\rho} @f$
     * @param[in] Rho_v the Wolfenstein parameter @f$ \bar{\eta} @f$
     */
    void computeCKMwithWolfenstein(double Lambda_v, double A_v, double Rho_v, double Eta_v);
    
    /**
     * @brief A set method to calculate the CKM matrix from CKM elements and @f$ \gamma @f$
     * @param[in] Vus_v the CKM element @f$ V_{us} @f$ (@f$ V_{ud} @f$ if the useVud flag is true)
     * @param[in] Vcb_v the CKM element @f$ V_{cb} @f$
     * @param[in] Vub_v the CKM element @f$ V_{ub} @f$
     * @param[in] gamma_v the CKM element @f$ \gamma @f$
     * @param[in] useVud (optional; if set to true, Vus_v is interpreted as the @f$ V_{ud} @f$ input value. Default: false
     */
    void computeCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v, bool useVud = false);

    /**
     * @brief A member for returning the CKM matrix.
     * @return the CKM matrix
     */
    gslpp::matrix<gslpp::complex> getCKM() const 
    {
        return V;
    }

    // Wolfenstein parameters
    
    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ \bar{\rho} @f$
     * @return the value of @f$ \bar{\rho} @f$
     */
    double getRhoBar() const 
    {
        return Rho;
    }
    
    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ \bar{\eta} @f$
     * @return the value of @f$ \bar{\eta} @f$
     */
    double getEtaBar() const
    {
        return Eta;
    }
    
    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ \lambda @f$
     * @return the value of @f$ \lambda @f$
     */
    double getLambda() const
    {
        return Lambda;
    }
    
    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ A @f$
     * @return the value of @f$ A @f$
     */
    double getA() const
    {
        return A;
    }

    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ \rho @f$
     * @return the value of @f$ \rho @f$
     */
    double getRho() const
    {
        return (s13 * cos(delta) / s12 / s23);
    }

    /**
     * @brief A member for returning the value of the Wolfenstein parameter @f$ \eta @f$
     * @return the value of @f$ \eta @f$
     */
    double getEta() const
    {
        return (s13 * sin(delta) / s12 / s23);
    }

    // Gilman parameterization

    /**
     * @brief A member for returning the value of the sine of the CKM parameter @f$ \theta_{12} @f$
     * @return the value of @f$ \sin\theta_{12} @f$
     */
    double gets12() const
    {
        return s12;
    }

    /**
     * @brief A member for returning the value of the sine of the CKM parameter @f$ \theta_{13} @f$
     * @return the value of @f$ \sin\theta_{13} @f$
     */
    double gets13() const
    {
        return s13;
    }

    /**
     * @brief A member for returning the value of the sine of the CKM parameter @f$ \theta_{23} @f$
     * @return the value of @f$ \sin\theta_{23} @f$
     */
    double gets23() const
    {
        return s23;
    }

    /**
     * @brief A member for returning the value of the cosine of the CKM parameter @f$ \theta_{12} @f$
     * @return the value of @f$ \cos\theta_{12} @f$
     */
    double getc12() const
    {
        return c12;
    }

    /**
     * @brief A member for returning the value of the cosine of the CKM parameter @f$ \theta_{23} @f$
     * @return the value of @f$ \cos\theta_{23} @f$
     */
    double getc23() const
    {
        return c23;
    }

    /**
     * @brief A member for returning the value of the cosine of the CKM parameter @f$ \theta_{13} @f$
     * @return the value of @f$ \cos\theta_{13} @f$
     */
    double getc13() const
    {
        return c13;
    }

    /**
     * @brief A member for returning the value of the CKM parameter @f$ \delta @f$
     * @return the value of @f$ \delta @f$
     */
    double getdelta() const
    {
        return delta;
    }

    // J_CP

    /**
     * @brief A member for returning the value of the Jarlskog determinant
     * @return the value of @f$ J @f$
     */
    double getJcp() const
    {
        return Eta * pow(A * pow(Lambda, 3), 2);
    }

    //Complex values of CKM elements

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{ud} @f$
     * @return the value of @f$ V_{ud} @f$
     */
    gslpp::complex getV_ud() const
    {
        return V(0, 0);
    }
    
    /**
     * @brief A member for returning the value of the CKM element @f$ V_{us} @f$
     * @return the value of @f$ V_{us} @f$
     */
    gslpp::complex getV_us() const
    {
        return V(0, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{ub} @f$
     * @return the value of @f$ V_{ub} @f$
     */
    gslpp::complex getV_ub() const
    {
        return V(0, 2);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{cd} @f$
     * @return the value of @f$ V_{cd} @f$
     */
    gslpp::complex getV_cd() const
    {
        return V(1, 0);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{cs} @f$
     * @return the value of @f$ V_{cs} @f$
     */
    gslpp::complex getV_cs() const
    {
        return V(1, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{cb} @f$
     * @return the value of @f$ V_{cb} @f$
     */
    gslpp::complex getV_cb() const
    {
        return V(1, 2);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{td} @f$
     * @return the value of @f$ V_{td} @f$
     */
    gslpp::complex getV_td() const
    {
        return V(2, 0);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{ts} @f$
     * @return the value of @f$ V_{ts} @f$
     */
    gslpp::complex getV_ts() const
    {
        return V(2, 1);
    }

    /**
     * @brief A member for returning the value of the CKM element @f$ V_{tb} @f$
     * @return the value of @f$ V_{tb} @f$
     */
    gslpp::complex getV_tb() const
    {
        return V(2, 2);
    }
    
     /**
     * @brief The %CKM angle @f$\beta@f$.
     * @return @f$\beta@f$ in radians
     */
    double computeBeta() const;

    /**
     * @brief The %CKM angle @f$\gamma@f$.
     * @return @f$\gamma@f$ in radians
     */
    double computeGamma() const;

    /**
     * @brief The %CKM angle @f$\alpha@f$.
     * @return @f$\alpha@f$ in radians
     */
    double computeAlpha() const;

    /**
     * @brief The %CKM angle @f$\beta_s@f$.
     * @return @f$\beta_s@f$ in radians
     */
    double computeBetas() const;

    /**
     * @brief The product of the %CKM elements @f$V_{td} V_{ts}^*@f$.
     * @return @f$V_{td} V_{ts}^*@f$
     */
    gslpp::complex computelamt() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cs}^*@f$.
     * @return @f$V_{cd} V_{cs}^*@f$
     */
    gslpp::complex computelamc() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{us}^*@f$.
     * @return @f$V_{ud} V_{us}^*@f$
     */
    gslpp::complex computelamu() const;

    /**
     * @brief The product of the %CKM elements @f$V_{td} V_{tb}^*@f$.
     * @return @f$V_{td} V_{tb}^*@f$
     */
    gslpp::complex computelamt_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cb}^*@f$.
     * @return @f$V_{cd} V_{cb}^*@f$
     */
    gslpp::complex computelamc_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{ub}^*@f$.
     * @return @f$V_{ud} V_{ub}^*@f$
     */
    gslpp::complex computelamu_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ts} V_{tb}^*@f$.
     * @return @f$V_{ts} V_{tb}^*@f$
     */
    gslpp::complex computelamt_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cs} V_{cb}^*@f$.
     * @return @f$V_{cs} V_{cb}^*@f$
     */
    gslpp::complex computelamc_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{us} V_{ub}^*@f$.
     * @return @f$V_{us} V_{ub}^*@f$
     */
    gslpp::complex computelamu_s() const;

    /**
     * @brief @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
     * @return @f$R_t@f$
     */
    double computeRt() const;

    /**
     * @brief @f$R_{ts}=|(V_{ts}V_{tb}^*)/(V_{cs}V_{cb}^*)|@f$.
     * @return @f$R_{ts}@f$
     */
    double computeRts() const;

    /**
     * @brief @f$R_b=|(V_{ud}V_{ub}^*)/(V_{ud}V_{ub}^*)|@f$.
     * @return @f$R_b@f$
     */
    double computeRb() const;
    
private:
    void computeCKMfromAngles();

    double Rho, Eta, Lambda, A; ///< The Wolfenstein parameters.
    double s12, s13, s23; ///< The sine of the three mixing angles
    double c12, c23, c13; ///< The cosine of the three mixing angles
    double delta; ///< The CP violating phase in the CKM matrix.

    gslpp::matrix<gslpp::complex> V; ///< The CKM matrix.

};

#endif	/* CKM_H */

