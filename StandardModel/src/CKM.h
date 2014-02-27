/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CKM_H
#define	CKM_H

#include <math.h>
#include <gslpp.h>

using namespace gslpp;

/**
 * @class CKM
 * @ingroup StandardModel
 * @brief A class for the %CKM matrix elements.
 * @author SusyFit Collaboration
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
  
    void getCKM(matrix<complex> &) const;
  
    // Wolfenstein parameters
    double getRho() const;
    double getEta() const ;
    double getLambda() const;
    double getA() const ;
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
    complex V_ud();
    complex V_us();
    complex V_ub();
    complex V_cd();
    complex V_cs();
    complex V_cb();
    complex V_td();
    complex V_ts();
    complex V_tb();

    /*
    // Angles
    double computeBeta();
    double computeGamma();
    double computeAlpha();
    double computeBetas();
    
    // Lambda_q
    complex computelamt();
    complex computelamc();
    complex computelamu();
    
    complex computelamt_d();
    complex computelamc_d();
    complex computelamu_d();
    
    complex computelamt_s();
    complex computelamc_s();
    complex computelamu_s();
    
    // Sides
    double getRt();
    double getRts();
    double getRb();
    */
    
private:
    double Rho, Eta, Lambda, A;
    double s12, s13, s23, delta;
    double c12, c23, c13;
    
    complex Vud, Vcd, Vtd;
    complex Vus, Vcs, Vts;
    complex Vub, Vcb, Vtb;

};

#endif	/* CKM_H */

