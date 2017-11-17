/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CKM_H
#define	CKM_H

#include <math.h>
#include <gslpp.h>

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
    ~CKM();

    void setWolfenstein(double Lambda_v, double A_v, double Rho_v, double Eta_v);
    void setCKM(double Vus_v, double Vcb_v, double Vub_v, double gamma_v);

    void setCKM(gslpp::matrix<gslpp::complex> V) {
        this->V = V;
    }

    gslpp::matrix<gslpp::complex> getCKM() const {
        return V;
    }

    // Wolfenstein parameters
    double getRho() const;
    double getEta() const;
    double getLambda() const;
    double getA() const;
    double getRhoNB() const;
    double getEtaNB() const;

    // Gilman parameterization
    double gets12() const;
    double gets13() const;
    double gets23() const;
    double getc12() const;
    double getc13() const;
    double getc23() const;
    double getdelta() const;

    // J_CP
    double getJcp() const;

    //Absolute values of CKM elements
    double getVud() const;
    double getVus() const;
    double getVub() const;
    double getVcd() const;
    double getVcs() const;
    double getVcb() const;
    double getVtd() const;
    double getVts() const;
    double getVtb() const;

    //Phases of CKM elements
    double getArgVud() const;
    double getArgVus() const;
    double getArgVub() const;
    double getArgVcd() const;
    double getArgVcs() const;
    double getArgVcb() const;
    double getArgVtd() const;
    double getArgVts() const;
    double getArgVtb() const;

    //Complex values of CKM elements
    gslpp::complex V_ud() const;
    gslpp::complex V_us() const;
    gslpp::complex V_ub() const;
    gslpp::complex V_cd() const;
    gslpp::complex V_cs() const;
    gslpp::complex V_cb() const;
    gslpp::complex V_td() const;
    gslpp::complex V_ts() const;
    gslpp::complex V_tb() const;

   
    // Angles
    double computeBeta() const;
    double computeGamma() const;
    double computeAlpha() const;
    double computeBetas() const;
    
    // Lambda_q
    gslpp::complex computelamt() const;
    gslpp::complex computelamc() const;
    gslpp::complex computelamu() const;
    
    gslpp::complex computelamt_d() const;
    gslpp::complex computelamc_d() const;
    gslpp::complex computelamu_d() const;
    
    gslpp::complex computelamt_s() const;
    gslpp::complex computelamc_s() const;
    gslpp::complex computelamu_s() const;
    
    // Sides
    double getRt() const;
    double getRts() const;
    double getRb() const;
    
private:
    void setCKMfromAngles();

    double Rho, Eta, Lambda, A;
    double s12, s13, s23, delta;
    double c12, c23, c13;

    gslpp::matrix<gslpp::complex> V;

};

#endif	/* CKM_H */

