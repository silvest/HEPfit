/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
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
    CKM(const CKM&);
    ~CKM();

    void setWolfenstein(double, double, double, double);
    void setCKM(double, double, double, double);

    void setCKM(gslpp::matrix<gslpp::complex> &);
    void getCKM(gslpp::matrix<gslpp::complex> &) const;

    // Wolfenstein parameters
    double getRho() const;
    double getEta() const;
    double getLambda() const;
    double getA() const;
    double getRhoNB();
    double getEtaNB();

    // Gilman parameterization
    double gets12();
    double gets13();
    double gets23();
    double getc12();
    double getc13();
    double getc23();
    double getdelta();

    // J_CP
    double getJcp();

    //Absolute values of CKM elements
    double getVud();
    double getVus();
    double getVub();
    double getVcd();
    double getVcs();
    double getVcb();
    double getVtd();
    double getVts();
    double getVtb();

    //Phases of CKM elements
    double getArgVud();
    double getArgVus();
    double getArgVub();
    double getArgVcd();
    double getArgVcs();
    double getArgVcb();
    double getArgVtd();
    double getArgVts();
    double getArgVtb();

    //Complex values of CKM elements
    gslpp::complex V_ud();
    gslpp::complex V_us();
    gslpp::complex V_ub();
    gslpp::complex V_cd();
    gslpp::complex V_cs();
    gslpp::complex V_cb();
    gslpp::complex V_td();
    gslpp::complex V_ts();
    gslpp::complex V_tb();

   
    // Angles
    double computeBeta();
    double computeGamma();
    double computeAlpha();
    double computeBetas();
    
    // Lambda_q
    gslpp::complex computelamt();
    gslpp::complex computelamc();
    gslpp::complex computelamu();
    
    gslpp::complex computelamt_d();
    gslpp::complex computelamc_d();
    gslpp::complex computelamu_d();
    
    gslpp::complex computelamt_s();
    gslpp::complex computelamc_s();
    gslpp::complex computelamu_s();
    
    // Sides
    double getRt();
    double getRts();
    double getRb();
    
private:
    double Rho, Eta, Lambda, A;
    double s12, s13, s23, delta;
    double c12, c23, c13;

    gslpp::complex Vud, Vcd, Vtd;
    gslpp::complex Vus, Vcs, Vts;
    gslpp::complex Vub, Vcb, Vtb;

};

#endif	/* CKM_H */

