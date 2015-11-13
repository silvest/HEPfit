/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SUSYMatching.h"
#include "SUSY.h"
#include <math.h>
#include <stdexcept>

SUSYMatching::SUSYMatching(const SUSY & SUSY_i) :

    StandardModelMatching(SUSY_i),
    mySUSY(SUSY_i),
    PV(true),

    mcdbd2(5, NDR, NLO),
    mcdbd2Hp(5, NDR, NLO),
    mcdbd2gg(5, NDR, NLO),
    mcdbd2ChiChi(5, NDR, NLO),
    mcdbd2Chi0Chi0(5, NDR, NLO),
    mcdbd2Chi0g(5, NDR, NLO),
    mcdbd2HpT(5, NDR, NLO),
    mcdbd2ggT(5, NDR, NLO),
    mcdbd2ChiChiT(5, NDR, NLO),
    mcdbd2Chi0Chi0T(5, NDR, NLO),
    mcdbd2Chi0gT(5, NDR, NLO),
    mcdbs2(5, NDR, NLO),
    mcdbs2Hp(5, NDR, NLO),
    mcdbs2gg(5, NDR, NLO),
    mcdbs2ChiChi(5, NDR, NLO),
    mcdbs2Chi0Chi0(5, NDR, NLO),
    mcdbs2Chi0g(5, NDR, NLO),
    mcdbs2HpT(5, NDR, NLO),
    mcdbs2ggT(5, NDR, NLO),
    mcdbs2ChiChiT(5, NDR, NLO),
    mcdbs2Chi0Chi0T(5, NDR, NLO),
    mcdbs2Chi0gT(5, NDR, NLO),
    mcdk2(5, NDR, NLO),
    mcdk2Hp(5, NDR, NLO),
    mcdk2gg(5, NDR, NLO),
    mcdk2ChiChi(5, NDR, NLO),
    mcdk2Chi0Chi0(5, NDR, NLO),
    mcdk2Chi0g(5, NDR, NLO),
    mcdk2HpT(5, NDR, NLO),
    mcdk2ggT(5, NDR, NLO),
    mcdk2ChiChiT(5, NDR, NLO),
    mcdk2Chi0Chi0T(5, NDR, NLO),
    mcdk2Chi0gT(5, NDR, NLO),
    mcdd2(5, NDR, NLO),
    mcdd2Hp(5, NDR, NLO),
    mcdd2gg(5, NDR, NLO),
    mcdd2ChiChi(5, NDR, NLO),
    mcdd2Chi0Chi0(5, NDR, NLO),
    mcdd2Chi0g(5, NDR, NLO),
    mcdd2HpT(5, NDR, NLO),
    mcdd2ggT(5, NDR, NLO),
    mcdd2ChiChiT(5, NDR, NLO),
    mcdd2Chi0Chi0T(5, NDR, NLO),
    mcdd2Chi0gT(5, NDR, NLO),
    mcDLij(2, NDR, LO),
    mcDLi3j(20, NDR, LO),
    mcmueconv(8, NDR, LO),
    mcgminus2mu(2, NDR, LO),
    mcbsg(8, NDR, NLO),
    mcbnlep(10, NDR, NLO, NLO_ew),
    mcbnlepCC(10, NDR, NLO),
    mcd1(10, NDR, NLO),
    mcd1Buras(10, NDR, NLO),
        
    myCKM(3, 3, 0.),
    myRu(6, 6, 0.),
    myRd(6, 6, 0.),
    myRl(6, 6, 0.),
    myRn(6, 6, 0.),
    mym_sn_sq(6, 0.),
    mym_se_sq(6, 0.),
    MChi0(4, 0.),
    myN(4, 4, 0.),
    MNeig(4, 0.),
    ON(4, 4, 0.),
    MChi(2, 0.),
    myV(2, 2, 0.),
    myU(2, 2, 0.),

    Eps_JCache(3,0.),
    Lambda0EpsYCache(3,3,0.),
    DeltaDL_Cache(3,3,0.),
    PHLRCache(3,3,0.),
    PHRLCache(3,3,0.),
    myCKM_cache(3, 3, 0.),
    VUDHH_cache(6, 6, 0.),
    DeltaMd_cache(3, 3, 0.),
    mySUSYMQ(6, 0.),

    Lepty(4, 6, 0.),
    Leptz(2, 3, 0.),
    Leptfa1(4, 6, 0.),
    Leptfa2(2, 3, 0.),
    Leptf1(4, 6, 0.),
    Leptf2(4, 6, 0.),
    Leptf3(2, 3, 0.),
    Leptf4(2, 3, 0.),
    CRlE(2, 3, 0.),
    CRlMU(2, 3, 0.),
    CRlTAU(2, 3, 0.),
    CRqUP(2, 3, 0.),
    CRqDOWN(2, 3, 0.),
    CLlE(2, 3, 0.),
    CLlMU(2, 3, 0.),
    CLlTAU(2, 3, 0.),
    CLqUP(2, 3, 0.),
    CLqDOWN(2, 3, 0.),
    NRlE(4, 6, 0.),
    NRlMU(4, 6, 0.),
    NRlTAU(4, 6, 0.),
    NRqUP(4, 6, 0.),
    NRqDOWN(4, 6, 0.),
    NLlE(4, 6, 0.),
    NLlMU(4, 6, 0.),
    NLlTAU(4, 6, 0.),
    NLqUP(4, 6, 0.),
    NLqDOWN(4, 6, 0.),
    AmpA1LN(4, 6, 0.),
    AmpA1RN(4, 6, 0.),
    AmpA1LC(2, 3, 0.),
    AmpA1RC(2, 3, 0.),
    AmpTauA1LN(4, 6, 0.),
    AmpTauA1RN(4, 6, 0.),
    AmpTauA1LC(2, 3, 0.),
    AmpTauA1RC(2, 3, 0.),
    AmpTEA1LN(4, 6, 0.),
    AmpTEA1RN(4, 6, 0.),
    AmpTEA1LC(2, 3, 0.),
    AmpTEA1RC(2, 3, 0.),
    AmpALN(4, 6, 0.),
    AmpARN(4, 6, 0.),
    AmpALC(2, 3, 0.),
    AmpARC(2, 3, 0.),
    AmpTauALN(4, 6, 0.),
    AmpTauARN(4, 6, 0.),
    AmpTauALC(2, 3, 0.),
    AmpTauARC(2, 3, 0.),
    AmpTEALN(4, 6, 0.),
    AmpTEARN(4, 6, 0.),
    AmpTEALC(2, 3, 0.),
    AmpTEARC(2, 3, 0.),

        TEhat(3, 3, 0.)
{
}


void SUSYMatching::updateSUSYParameters()
{
    myCKM = mySUSY.getVCKM();
    myRu = mySUSY.getRu();
    myRd = mySUSY.getRd();
    myRl = mySUSY.getRl();
    myRn = mySUSY.getRn();
    mym_sn_sq = mySUSY.getMsn2();
    mym_se_sq = mySUSY.getMse2();
    Q_S = mySUSY.getQ_SUSY();
    mu2R = Q_S * Q_S;
    tanb = mySUSY.getTanb();
    sinb = mySUSY.getSinb();
    cosb = mySUSY.getCosb();
    Als = mySUSY.Als(Q_S);
    Mg = mySUSY.getMGl();
    MChi0 = mySUSY.getMneu();
    MChi = mySUSY.getMch();
    v = mySUSY.v();
    v1 = mySUSY.v1();
    v2 = mySUSY.v2();
    gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    myN = mySUSY.getN();
    myV = mySUSY.getV();
    myU = mySUSY.getU();

}

/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
////////////// D0(x,y,z,t) and D2(x,y,z,t) limits

double SUSYMatching::D0N(double x, double y, double z, double t) {

    if ((fabs(z) < SUSYLEPS3) || (fabs(t) < SUSYLEPS3)) {
        return (0.);
    }
    else
        return (sqrt(z * t) * Dk(x, y, z, t, 0));

}

double SUSYMatching::D2LL0(double a, double b) {

    return ( (log(a) - log(b)) / (4. * (-a + b)) );

}


double SUSYMatching::DL0(double a, double b, double c, int k) {

    if (k == 2) {
        return ( (a * (-b + c) * log(a) + b * (a - c) * log(b)
                + (-a + b) * c * log(c)) / (4. * (a - b)*(-a + c)*(-b + c)));
    }
    else if (k == 0) {

        return (((-b + c) * log(a) + (a - c) * log(b) + (-a + b) * log(c)) /
                ((a - b)*(-a + c)*(-b + c)));
    }
    else {

        throw std::runtime_error("Error in DL0(a,b,c,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::DL(double a, double b, double c, int k) {
    if (k == 0) {
        if (fabs(a) < SUSYLEPS) {
            
            throw std::runtime_error("Error in DL(a,b,c,0) in SUSYMatching.cpp because the limit DL(0,b,c, k = 0) is singular");
        }
        else if (fabs(b) < SUSYLEPS) {

            return ( (-a + c + a * log(a) - a * log(c)) / (a * (a - c) * (a - c)));
        }
        else if (fabs(c) < SUSYLEPS) {

            return ( (-a + b + a * log(a) - a * log(b)) / (a * (a - b) * (a - b)));
        }
        else
            return( (a*a-b*c)/(a-b)/(a-b)/(a-c)/(a-c)*log(a) +
                     b/(a-b)/(a-b)/(c-b)*log(b) - c/(a-c)/(a-c)/(c-b)*log(c) 
                     - 1./(a-b)/(a-c) ) ;
    }
    else if (k == 2) {

        if (fabs(a) < SUSYLEPS) {

            return ( (-log(b) + log(c)) / (4. * (b - c)));
        }
        else if (fabs(b) < SUSYLEPS) {

            return ( (-a + c + c * log(a) - c * log(c)) / (4. * (a - c) * (a - c)));
        }
        else if (fabs(c) < SUSYLEPS) {

            return ( (-a + b + b * log(a) - b * log(b)) / (4. * (a - b) * (a - b)));
        }

            return ( a*(-2.*b*c + a*(b + c))/4./(a - b)/(a - b)/(a - c)/(a - c)*log(a)
                    + b*b/4./(a - b)/(a - b)/(c - b)*log(b) - c*c/4./(a - c)/(a - c)/
                    (c - b)*log(c) + a/4./(a - b)/(c - a) );
    }
    else {

        throw std::runtime_error("Error in DL(a,b,c,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::DLL(double a, double b, int k) {
    if (k == 0) {
        
        if (fabs(a) < SUSYLEPS) {

            throw std::runtime_error("Error in DLL(a,b,k = 0) in SUSYMatching.cpp because the limit DLL(0,b, k = 0) is singular");
        }

        if (fabs(b) < SUSYLEPS) {
            return (1. / (2. * a * a));
        }
        else
            return ( -(-a * a + b * b + 2 * a * b * log(a) - 2 * a * b * log(b)) /      // sign mistake found
                (2. * a * (a - b) * (a - b) * (a - b)));
    }
    else if (k == 2) {

        if (fabs(a) < SUSYLEPS) {

            throw std::runtime_error("Error in DLL(a,b,k = 2) in SUSYMatching.cpp because the limit DLL(0,b, k = 2) is singular");
        }

        if (fabs(b) < SUSYLEPS) {

            return ( -1. / (8. * a));
        }
        else
            return (a * a - 4 * a * b + 3 * b * b + 2 * b * b * log(a)
                - 2 * b * b * log(b)) / (8. * (-a + b) * (-a + b) * (-a + b));
    }
    else {

        throw std::runtime_error("Error in DLL(a,b,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::DLLp(double a, double b, int k) {
    if (k == 0) {
        
        if ( (fabs(a) < SUSYLEPS) || (fabs(b) < SUSYLEPS) ) {
            
            std::cout << "MChargini = " << mySUSY.getMch() << std::endl;
            std::cout << "MNeutralini = " << mySUSY.getMneu() << std::endl;
            std::cout << "MD2squarks = " << mySUSY.getMsd2() << std::endl;
            std::cout << "MU2squarks = " << mySUSY.getMsu2() << std::endl;
            std::cout << "Mgluino = " << mySUSY.getM3() << std::endl;
            std::cout << "MHp = " << mySUSY.getMHp() << std::endl;
            std::cout << "Mw_tree = " << mySUSY.Mw_tree() << std::endl;
            std::cout << "Mup(Q_S) = " << mySUSYMQ(0) << std::endl;
            std::cout << "Mdown(Q_S) = " << mySUSYMQ(1) << std::endl;
            std::cout << "Mc(Q_S) = " << mySUSYMQ(2) << std::endl;
            std::cout << "Ms(Q_S) = " << mySUSYMQ(3) << std::endl;
            std::cout << "Mtop(Q_S) = " << mySUSYMQ(4) << std::endl;
            std::cout << "Mb(Q_S) = " << mySUSYMQ(5) << std::endl;
            
            
            throw std::runtime_error("Error in DLLp function, because the limits D0(0,0,b,b) and D0(a,a,0,0) are singular "); 
        }
        return (-2 * a + 2 * b + (a + b) * log(a) - (a + b) * log(b))
                / ((a - b)*(a - b)*(a - b));
    }
    else if (k == 2) {
            
            if ( fabs(a) < SUSYLEPS ){
                
                return (-1. / (4. * b ));
            }
            else if ( fabs(b) < SUSYLEPS ){
                    
                return ( -1. / (4. * a) );
            }
        
        return (-a * a + b * b + 2 * a * b * log(a) - 2 * a * b * log(b))
                / (4. * (a - b)*(a - b)*(a - b));
    }
    else {

        throw std::runtime_error("Error in DLLp(a,b,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
    
}

double SUSYMatching::DLLL(double a, int k) {
    if (k == 0) {
        return (1. / (6. * a * a));
    }
    else if (k == 2) {
        return (-1. / (12. * a));
    }
    else {

        throw std::runtime_error("Error in DLLL(a,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////  Dk (k = 0,2) controls

double SUSYMatching::Dk(double x, double y, double z, double t, int k) {
   
    /// four variables equals
    if ((fabs(1. - y / x) < SUSYLEPS) && (fabs(1. - z / x) < SUSYLEPS) &&
            (fabs(1. - t / x) < SUSYLEPS)) {
        return DLLL(x, k);
    }

    /// three variables equals  
    if ((fabs(1. - y / x) < SUSYLEPS) && (fabs(1. - z / x) < SUSYLEPS)) {
        return DLL(x, t, k);
    }
    if ((fabs(1. - y / x) < SUSYLEPS) && (fabs(1. - t / x) < SUSYLEPS)) {
        return DLL(x, z, k);
    }
    if ((fabs(1. - z / x) < SUSYLEPS) && (fabs(1. - t / x) < SUSYLEPS)) {
        return DLL(x, y, k);
    }
    if ((fabs(1. - z / y) < SUSYLEPS) && (fabs(1. - t / y) < SUSYLEPS)) {
        return DLL(y, x, k);
    }

    /// pairs variables equals 
    if (fabs(1. - y / x) < SUSYLEPS) {
        if ((fabs(1. - t / z) < SUSYLEPS)) {
            return DLLp(x, z, k);
        }
        else
            return DL(x, z, t, k);
    }
    if (fabs(1. - z / x) < SUSYLEPS) {
        if (fabs(1. - t / y) < SUSYLEPS) {
            return DLLp(x, y, k);
        }
        else
            return DL(x, y, t, k);
    }
    if (fabs(1. - t / x) < SUSYLEPS) {
        if (fabs(1. - z / y) < SUSYLEPS) {
            return DLLp(x, y, k);
        }
        else
            return DL(x, z, y, k);
    }
    if ((fabs(1. - z / y) < SUSYLEPS)) {
        return DL(y, x, t, k);
    }
    if ((fabs(1. - t / y) < SUSYLEPS)) {
        return DL(y, z, x, k);
    }
    if ((fabs(1. - z / t) < SUSYLEPS)) {
        
        return (DL(x, y, z, k));
    }

    /// one variable equal to zero
    
    if (fabs(x) < SUSYLEPS) {

        return (DL0(y, z, t, k));
    }
    if (fabs(y) < SUSYLEPS) {

        return (DL0(x, z, t, k));
    }
    if (fabs(z) < SUSYLEPS) {

        return (DL0(x, y, t, k));
    }
    if (fabs(t) < SUSYLEPS) {

        return (DL0(x, y, z, k));
    }
    
    
    /// different variables 
    if (k == 0) {
        return x * log(x) / ((t - x)*(z - x)*(y - x)) + y * log(y)
                / ((t - y)*(z - y)*(x - y)) + z * log(z)
                / ((t - z)*(x - z)*(y - z)) +
                t * log(t) / ((x - t)*(z - t)*(y - t));
    }
    else if (k == 2) {
        return (((t * t * log(t)) / ((-t + x)*(-t + y)*(-t + z))
                + (x * x * log(x)) / ((t - x)*(-x + y)*(-x + z))
                + (y * y * log(y)) / ((t - y)*(x - y)*(-y + z))
                + (z * z * log(z)) / ((t - z)*(x - z)*(y - z)))/4.);
    }
    else {

        throw std::runtime_error("Error in Dk(x,y,z,t,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
    
}

////////////////////////////////////////////////////////////////////////////////
///// Ck limits

double SUSYMatching::CL(double a, double b, int k) {
    
    if (k == 0) {
        return ((a - b - b * log(a) + b * log(b)) / ((a - b)*(a - b)));

    }
    else if (k == 2) {

        return (-log(mu2R) + (a * (a - b) + a * (a - 2 * b) * log(a) + b * b * log(b)) /
                ((a - b)*(a - b)));
    }
    else {

        throw std::runtime_error("Error in CL(a,b,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::CLL(double a, int k) {
    
    if (k == 0) {
        return (1. / (2. * a));

    }
    else if (k == 2) {

        return (-log(mu2R) + 1.5 + log(a));
    }
    else {

        throw std::runtime_error("Error in CLL(a,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}


////////////////////////////////////////////////////////////////////////////////
///// Ck function

double SUSYMatching::Ck(double x, double y, double z,int k) {

    if ((fabs(1. - y / x) < SUSYLEPS)&(fabs(1. - z / x) < SUSYLEPS)) {

        return (CLL(x,k));

    }
    else if ((fabs(1. - y / x) < SUSYLEPS)) {

        return (CL(x,z,k));

    }
    else if ((fabs(1. - z / x) < SUSYLEPS)) {

        return (CL(x,y,k));

    }
    else if ((fabs(1. - z / y) < SUSYLEPS)) {

        return (CL(y,x,k));

    }

    if (k == 0) {
        return ((y * log(y / x)) / ((x - y)*(-y + z)) + (z * log(z / x)) /
                ((x - z)*(y - z)));
    }
    else if (k == 2) {
        return (-log(mu2R) + log(x) + y * y * log(y / x) / ((x - y)*(-y + z)) 
                + z * z * log(z / x) / ((x - z)*(y - z)));
    }
    else {

        throw std::runtime_error("Error in Ck(x,y,z,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
    
}

////////////////////////////////////////////////////////////////////////////////
////////////// Bk(0,x,y) functions and limits

double SUSYMatching::BL(double a, int k) {
    
    if (k == 0) {
        return (1. + log(a));
    }
    else if (k == 2) {
        return ((-log(mu2R) + 2. + log(a)) / 2.);
    }
    else {

        throw std::runtime_error("Error in BL(a,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::Bk(double x, double y, int k) {

    /** Limit of B0 function for zero x or zero y **/
    if((fabs(x) < SUSYLEPS2)&&(k==0)) return (log(y));
    if((fabs(y) < SUSYLEPS2)&&(k==0)) return (log(x));
    
    
    if (fabs((1. - y / x)) < SUSYLEPS2) {
        return (BL(x, k));

    }
    if (k == 0) {
        return (log(x) + (y * log(x / y)) / (x - y));
    }
    else if (k == 2) {

        return (1. / 4. + 0.5 * Ck(x, y, y, 2));
    }
    else {

        throw std::runtime_error("Error in Bk(x,y,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);

}
////////////////////////////////////////////////////////////////////////////////
//// Running quark masses to SUSY scale Q

void SUSYMatching::Comp_mySUSYMQ() {
    
        mySUSYMQ(0) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.UP).getMass_scale(),mySUSY.getQuarks(mySUSY.UP).getMass());
        mySUSYMQ(1) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.DOWN).getMass_scale(),mySUSY.getQuarks(mySUSY.DOWN).getMass());
        mySUSYMQ(2) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.CHARM).getMass());
        mySUSYMQ(3) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.STRANGE).getMass_scale(),mySUSY.getQuarks(mySUSY.STRANGE).getMass());
        mySUSYMQ(4) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.TOP).getMass());
        mySUSYMQ(5) = mySUSY.Mrun(Q_S,mySUSY.getQuarks(mySUSY.BOTTOM).getMass());
}

////////////////////////////////////////////////////////////////////////////////
////  Quark self-energies from Buras arXiv:hep-ph/0210145v2

void SUSYMatching::Comp_DeltaMd() {

    gslpp::vector<double> myMU2Squarks(6, 0.);
    gslpp::vector<double> myMD2Squarks(6, 0.);
    myMU2Squarks = mySUSY.getMsu2();
    myMD2Squarks = mySUSY.getMsd2();

    int k, l, I, J;//, i, j;
    
    for (J = 0; J < 3; J++) {
        for (I = 0; I < 3; I++) {

            gslpp::complex temp(0., 0., false);

            for (k = 0; k < 6; k++) {

                temp += -2. / 3. * Als / M_PI * Mg * myRd(k, J + 3).conjugate()
                        * myRd(k, I) * Bk(Mg * Mg, myMD2Squarks(k), 0);

            }
          
            for (k = 0; k < 6; k++) {
                for (l = 0; l < 4; l++) {

                    temp +=  1. / (16. * M_PI * M_PI) * VdDNR(J, k, l, 0).conjugate()
                            * VdDNL(I, k, l, 0)
                            * MChi0(l) * Bk(MChi0(l) * MChi0(l), myMD2Squarks(k), 0);
                }

            }
            
            for (k = 0; k < 6; k++) {
                for (l = 0; l < 2; l++) {

                    temp += 1. / (16. * M_PI * M_PI) * VdUCR(J, k, l, 0).conjugate() * VdUCL(I, k, l) *
                            MChi(l) * Bk(MChi(l) * MChi(l), myMU2Squarks(k), 0);
                }
            }
            
            DeltaMd_cache.assign(J, I, temp);
        }
    }
} 

gslpp::complex SUSYMatching::DeltaMd(int J, int I) {

    return (DeltaMd_cache(J,I));
}

////////////////////////////////////////////////////////////////////////////////
///// Epsilon_J and EpilonY_JI

void SUSYMatching::Comp_Eps_J() {
    
    for (int J = 0; J < 3; J++) {

        Eps_JCache.assign(J, DeltaMd(J, J) / (tanb * mySUSYMQ(2 * J + 1)));
    }
}

gslpp::complex SUSYMatching::Eps_J(int J) {

    return (Eps_JCache(J));

}

void SUSYMatching::Comp_Lambda0EpsY(){
    
    int I,J;
    double mtop = mySUSYMQ(4);
    
    for( I = 0;I < 3;I++){
        for( J = 0;J < 3;J++){
            
            Lambda0EpsYCache.assign(J,I,DeltaMd(J,I)/(tanb
            * mySUSYMQ(2 * J + 1)) /( 2 /(v2 * v2)  
            * mtop * mtop ));
        }
    }
}



gslpp::complex SUSYMatching::Lambda0EpsY(int J, int I){
    
    
    
    return (Lambda0EpsYCache(J,I));
    
}


////////////////////////////////////////////////////////////////////////////////
///// mySUSY_CKM
//// CKM matrix with large tan beta corrections

void SUSYMatching::Comp_DeltaDL() {

    int I, J;
    gslpp::complex C(0., 0., false);

    for (I = 0; I < 3; I++) {
        for (J = 0; J < 3; J++) {
            if (J != I) {
                double mdJ = mySUSYMQ(2 * J + 1);
                double mdI = mySUSYMQ(2 * I + 1);
                DeltaDL_Cache.assign(J, I, -(mdJ * DeltaMd(J, I) + DeltaMd(I, J).conjugate() * mdI) /
                        (mdJ * mdJ - mdI * mdI));
            } else
                if (J == I) {
                DeltaDL_Cache.assign(J, I, C);
            }

        }

    }
    
}

gslpp::complex SUSYMatching::DeltaDL(int J, int I) {

    return(DeltaDL_Cache(J,I));
}

gslpp::complex SUSYMatching::DeltaDR(int J, int I){
    
    gslpp::complex C(0., 0., false);
    C = DeltaDL_Cache(J,I).conjugate();
    
   return (C);  
}

void SUSYMatching::Comp_mySUSY_CKM() {

  
    
    gslpp::complex Delta_CKM_IJ(0.,0.,false);
    
    int l, I, J;

    for (I = 0; I < 3; I++) {
        for (J = 0; J < 3; J++) {
            
            for (l = 0; l < 3; l++) {

                Delta_CKM_IJ += myCKM(I, l) * DeltaDL(l, J);
                
            }
            
            myCKM_cache.assign(I, J, myCKM(I, J) - Delta_CKM_IJ);
            Delta_CKM_IJ.assign(0.,0.,0);
           
        }
    }    
}

gslpp::matrix<gslpp::complex> SUSYMatching::mySUSY_CKM() {

    return (myCKM_cache);

}
////////////////////////////////////////////////////////////////////////////////
///// VChiUdL element j,k,b <-> VdUCL element b,k,j

void SUSYMatching::Comp_VdUCL() {


    gslpp::complex VdUCL_bkj(0., 0., false);
    int l, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {


                for (l = 0; l < 3; l++) {
                    VdUCL_bkj += -gW * myRu(k, l) * myCKM(l, b) * myV(j, 0).conjugate()
                            + sqrt(2.) / v2 * mySUSYMQ(2 * l)
                            * myRu(k, l + 3) * myV(j, 1).conjugate() * myCKM(l, b);
                }

                VdUCL_cache[b][k][j] = VdUCL_bkj;
                VdUCL_bkj.assign(0., 0., 0);

            }
        }
    }
}

gslpp::complex SUSYMatching::VdUCL(int b, int k, int j) {

    return (VdUCL_cache[b][k][j]);

}


////////////////////////////////////////////////////////////////////////////////
/////  VdUCR element b,k,j <-> VChiUdR element j,k,b 
void SUSYMatching::Comp_VdUCR(int flag) {

    gslpp::complex VdUCR_bkj(0., 0., false);
    gslpp::complex Ydb(0., 0., false);
    gslpp::complex Ydp(0., 0., false);
    int l, p, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {
                for (l = 0; l < 3; l++) {
                    if(flag == 1) {
                        Ydb = mySUSYMQ(2 * b + 1) * sqrt(2.) / v1 /(1. + Eps_J(b)*tanb) ;
                        VdUCR_bkj += Ydb * myCKM(l, b) * myRu(k, l) * myU(j, 1);
                       
                        for (p = 0; p < 3; p++) {
                            Ydp = mySUSYMQ(2 * p + 1)* sqrt(2.) / v1 /(1. + Eps_J(p)* tanb) ;
                            VdUCR_bkj += (myCKM(l, p)*(Ydp * DeltaDR(p, b) 
                                    - DeltaDL(p, b) * Ydb) * myRu(k, l) * myU(j, 1));  
                        }
                              
                    }
                    else if(flag == 0){
                            Ydb = mySUSYMQ(2 * b + 1) * sqrt(2.) / v1;
                            VdUCR_bkj += Ydb * myRu(k, l) * myU(j, 1)
                                    * myCKM(l, b);
                    }
                    
                    else{throw std::runtime_error("Wrong flag assigned to vertex VdUCR_bkj å");}
                    
                }

                VdUCR_cache[b][k][j][flag] = VdUCR_bkj;
                VdUCR_bkj.assign(0., 0., 0);

            }
        }
    }
}


gslpp::complex SUSYMatching::VdUCR(int b, int k, int j, int flag) {
   
    return (VdUCR_cache[b][k][j][flag]);       
}


////////////////////////////////////////////////////////////////////////////////
///// VChiDdL element j,k,b <-> VdDNL element b,k,j

void SUSYMatching::Comp_VdDNL(int flag){
    
    gslpp::complex VdDNL_bkj(0., 0., false);
    /* tree-level cW2 */
    double cW2 = mySUSY.Mw_tree()*mySUSY.Mw_tree()/mySUSY.getMz()/mySUSY.getMz();
    /* SM value for cW2 in the on-shell scheme */
    //double cW2 = mySUSY.StandardModel::cW2();
    /* MSSM value for cW2 in the on-shell scheme */
    //double cW2 = mySUSY.cW2();
    double sW2 = 1.0 - cW2;
    double CosThetaW = sqrt(cW2);
    double SinThetaW = sqrt(sW2);
    int l, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {
                
                    if(flag==1) {

                            VdDNL_bkj += (-gW / sqrt(2.) * myRd(k, b) * (1. / 3.
                                    * SinThetaW / CosThetaW * myN(j, 0).conjugate() -
                                    myN(j, 1).conjugate()) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1) /
                                    (1. + Eps_J(b) * tanb)
                                    * myRd(k, b + 3) * myN(j, 2).conjugate());
                            
                            for (l = 0; l < 3; l++) {

                                // correzione vertici neutralini calcolate seguendo le indicazioni di Buras    

                                VdDNL_bkj += (-gW / sqrt(2.) * myRd(k, l) * (1. / 3.
                                        * SinThetaW / CosThetaW * myN(j, 0).conjugate() -
                                        myN(j, 1).conjugate()) - sqrt(2.) / v1
                                        * mySUSYMQ(2 * l + 1) /
                                        (1. + Eps_J(l)* tanb)
                                        * myRd(k, l + 3) * myN(j, 2).conjugate())
                                        *DeltaDL(l, b);
                            }
                        }

                    if(flag==0){
                            
                            VdDNL_bkj += -gW / sqrt(2.) * myRd(k, b) * (1. / 3.
                                    * SinThetaW / CosThetaW * myN(j, 0).conjugate() -
                                    myN(j, 1).conjugate()) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1)
                                    * myRd(k, b + 3) * myN(j, 2).conjugate();

                    }

                    VdDNL_cache[b][k][j][flag] = VdDNL_bkj;
                    VdDNL_bkj.assign(0., 0., 0);

            }

        }
    }
    
    if(flag != 0 && flag != 1) throw std::runtime_error("Error in Comp_VdDNL(flag) in SUSYMatching.cpp  ");
}

gslpp::complex SUSYMatching::VdDNL(int b, int k, int j, int flag) {


    return (VdDNL_cache[b][k][j][flag]);

}

////////////////////////////////////////////////////////////////////////////////
///// VdDNR element b,k,j <-> VChiDdR element j,k,b

void SUSYMatching::Comp_VdDNR(int flag) {

    gslpp::complex VdDNR_bkj(0., 0., false);
    /* tree-level cW2 */
    double cW2 = mySUSY.Mw_tree()*mySUSY.Mw_tree()/mySUSY.getMz()/mySUSY.getMz();
    /* SM value for cW2 in the on-shell scheme */
    //double cW2 = mySUSY.StandardModel::cW2();
    /* MSSM value for cW2 in the on-shell scheme */
    //double cW2 = mySUSY.cW2();
    double sW2 = 1.0 - cW2;
    
    double CosThetaW = sqrt(cW2);
    double SinThetaW = sqrt(sW2);
    int l, b, k, j;
    
    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {

                        if(flag==1){

                            VdDNR_bkj += -sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                    myRd(k, b + 3) * myN(j, 0) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1) /
                                    (1. + Eps_J(b) * tanb) * myRd(k, b) *
                                    myN(j, 2);

                            for (l = 0; l < 3; l++) {

                                // correzione vertici neutralini calcolate seguendo le indicazioni di Buras   

                                VdDNR_bkj += (-sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                        myRd(k, l + 3) * myN(j, 0) - sqrt(2.) / v1
                                        * mySUSYMQ(2 * l + 1) /
                                        (1. + Eps_J(l) * tanb) * myRd(k, l) *
                                        myN(j, 2)) * DeltaDR(l, b);
                            }
                        }

                        if(flag==0){
                            VdDNR_bkj += -sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                        myRd(k, b + 3) * myN(j, 0) - sqrt(2.) / v1
                                        * mySUSYMQ(2 * b + 1)* myRd(k, b) * myN(j, 2);
                        }
                   
                    VdDNR_cache[b][k][j][flag] = VdDNR_bkj;
                    VdDNR_bkj.assign(0., 0., 0);

            }
        }
    }
    
    if(flag != 0 && flag != 1) throw std::runtime_error("Error in Comp_VdDNR(flag) in SUSYMatching.cpp  "); 
    
}


gslpp::complex SUSYMatching::VdDNR(int b, int k, int j, int flag) {

    return (VdDNR_cache[b][k][j][flag]);
    
}


////////////////////////////////////////////////////////////////////////////////
/////  Feynmann rules of uDC vertex (D mixing)

void SUSYMatching::Comp_VuDCL() { // to be corrected


    gslpp::matrix<gslpp::complex> myCKM(3, 3, 0.);
    gslpp::complex VuDCL_bkj(0., 0., false);
    gslpp::complex YdI(0., 0., false);
    int I, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {

                for (I = 0; I < 3; I++) {
                    YdI = sqrt(2.) / v1 * mySUSYMQ(2 * I + 1); 
                    VuDCL_bkj += -(gW * myRd(k, I) * myU(j, 0).conjugate()
                            - YdI * myRd(k, I + 3) * myU(j, 1).conjugate()) *
                            myCKM(b, I).conjugate();
                }

                VuDCL_cache[b][k][j] = VuDCL_bkj;
                VuDCL_bkj.assign(0., 0., 0);

            }
        }
    }
}

gslpp::complex SUSYMatching::VuDCL(int b, int k, int j) {

    return (VuDCL_cache[b][k][j]);

}

void SUSYMatching::Comp_VuDCR() {


    gslpp::matrix<gslpp::complex> mySUSYCKM(3, 3, 0.);
    mySUSYCKM = mySUSY_CKM();
    gslpp::matrix<gslpp::complex> myV(2, 2, 0.);
    myV = mySUSY.getV();
    double Yub;
    gslpp::complex VuDCR_bkj(0., 0., false);
   
    int I, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {

                Yub = sqrt(2.) / v2 * mySUSYMQ(2 * b); // b is the up quark type index

                for (I = 0; I < 3; I++) {

                    VuDCR_bkj +=  mySUSYCKM(b, I).conjugate() * Yub * myRd(k, I) * myV(j, 1) ; 

                }

                VuDCR_cache[b][k][j] = VuDCR_bkj;
                VuDCR_bkj.assign(0., 0., 0);
            
                
            }
        }
    }
}

gslpp::complex SUSYMatching::VuDCR(int b, int k, int j) {

    return (VuDCR_cache[b][k][j]);
}

gslpp::complex SUSYMatching::VdUCL(int b, int k, int j, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdUCL(b, k, j));
    }
    else if (Dmixingflag == 1) {

        return (VuDCL(b, k, j));
    }
    else {
            
             throw std::runtime_error("Error in VdUCL(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  "); 
    }
}

gslpp::complex SUSYMatching::VdUCR(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdUCR(b, k, j, flag));
    }
    else if (Dmixingflag == 1) {

        return (VuDCR(b, k, j));
    }
    else {
            
             throw std::runtime_error("Error in VdUCR(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  "); 
    }
}

/* Vertices uUN from Buras arXiv:hep-ph/0210145v2 in SLHA convention  
     usefull in D - Dbar mixing */

void SUSYMatching::Comp_VuUN(){
 
    double TanThetaW = sqrt(mySUSY.sW2() / mySUSY.cW2());
    gslpp::complex VuUNL_bkj(0.,0.,false);
    gslpp::complex VuUNR_bkj(0.,0.,false);
    double Yub;
    
    int b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {

                Yub = sqrt(2.) / v2 * mySUSYMQ(2 * b);

                VuUNL_bkj += -1. / sqrt(2.) * gW * myRu(k, b) * (1. / (3. * TanThetaW) *
                        myN(j, 0).conjugate() + myN(j, 1).conjugate())
                        - Yub * myRu(k, b + 3) * myN(j, 3).conjugate();

                VuUNR_bkj += 2. * sqrt(2.) / 3. * gW * TanThetaW * myRu(k, b + 3) *
                        myN(j, 0) - Yub * myRu(k, b) * myN(j, 3);

                VuUNL_cache[b][k][j] = VuUNL_bkj;
                VuUNL_bkj.assign(0., 0., 0);
                VuUNR_cache[b][k][j] = VuUNR_bkj;
                VuUNR_bkj.assign(0.,0.,0);

            }
        }
    }    
}



gslpp::complex SUSYMatching::VuUN(int b, int k, int j, const std::string chirality) {


    if (chirality.compare("L") == 0) {

        return (VuUNL_cache[b][k][j]);

    }
    else if (chirality.compare("R") == 0) {

        return (VuUNR_cache[b][k][j]);

    }
    else {

        std::cout << " VuUN error in SUSYMatching.cpp" << std::endl;
        
    }

    return (EXIT_FAILURE);
    
}

/******************************************************************************/


gslpp::complex SUSYMatching::VdDNL(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdDNL(b, k, j, flag));

    }
    else if (Dmixingflag == 1) {

        return (VuUN(b, k, j, "L"));
    }
    else {

        throw std::runtime_error("Error in VdDNL(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::VdDNR(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdDNR(b, k, j, flag));

    }
    else if (Dmixingflag == 1) {

        return (VuUN(b, k, j, "R"));
    }
    else {

        throw std::runtime_error("Error in VdDNR(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
////
////////////////////////////////////////////////////////////////////////////////
//// PGLR, PGRL, PHLR, PHRL element j,i

gslpp::complex SUSYMatching::PGLR(int j, int i) {

    return (-sqrt(2.) / v * myCKM(j, i) *
            mySUSYMQ(2 * i + 1));
}

gslpp::complex SUSYMatching::PGRL(int j, int i) {

    return (sqrt(2.) / v * myCKM(j, i) *
            mySUSYMQ(2 * j));
}

void SUSYMatching::Comp_PHLR(){

    gslpp::complex PHLR(0., 0., false);
    int k;
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {



            PHLR += myCKM(j, i) *
                    mySUSYMQ(2 * i + 1) / (1. + Eps_J(i) * tanb);

            for (k = 0; k < 3; k++) {

                PHLR += myCKM(j, k) * mySUSYMQ(2 * k + 1) /
                        (1. + Eps_J(k) * tanb) * DeltaDR(k, i) - myCKM(j, k) *
                        DeltaDL(k, i) * mySUSYMQ(2 * i + 1) /
                        (1. + Eps_J(i) * tanb);
                
            }

            PHLR *= sqrt(2.) / mySUSY.v() * tanb;

            PHLRCache.assign(j,i,PHLR);
            PHLR.assign(0.,0.,0);
        }
    }     
}

gslpp::complex SUSYMatching::PHLR(int j, int i) {

    return(PHLRCache(j,i));
}


void SUSYMatching::Comp_VUDHH(){
    
    gslpp::complex VUDHijH(0., 0., false);
    gslpp::matrix<gslpp::complex> myTU(3, 3, 0.);
    gslpp::matrix<gslpp::complex> myTD(3, 3, 0.);
    gslpp::complex YuJ(0., 0., false);
    gslpp::complex YdI(0., 0., false);
    gslpp::matrix<gslpp::complex> mySUSYCKM(3, 3, 0.);
    mySUSYCKM = mySUSY_CKM();
    myTU = mySUSY.getTUhat();
    myTD = mySUSY.getTDhat();
    gslpp::matrix<gslpp::complex> ZH(2, 2, 0.);
    ZH.assign(0, 0, mySUSY.getSinb());
    ZH.assign(0, 1, -mySUSY.getCosb());
    ZH.assign(1, 0, mySUSY.getCosb());
    ZH.assign(1, 1, mySUSY.getSinb());
    int I, J, i, j;// ,l;

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {


            for (I = 0; I < 3; I++) {
                for (J = 0; J < 3; J++) {
                    

                    VUDHijH += v / sqrt(2.) * sqrt(2.) / v2 * mySUSYMQ(2 * J) * 
                            sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) /
                            (1 + Eps_J(I) * tanb) * mySUSYCKM(J, I) * myRd(j, I + 3).conjugate() *
                            myRu(i, J + 3) 
                           
                             + 1. / sqrt(2.) * (v1 * sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) /
                            (1 + Eps_J(I) * tanb) * sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) /
                            (1 + Eps_J(I) * tanb) * ZH(0, 0) +
                            v2 * sqrt(2.) / v2 * mySUSYMQ(2 * J) * sqrt(2.) / v2 * mySUSYMQ(2 * J)
                            * ZH(1, 0)) * mySUSYCKM(J, I) * myRd(j, I).conjugate() * myRu(i, J) 
                             
                                     
                            + (ZH(0, 0) * (mySUSY.getMuH()).conjugate() *
                            sqrt(2.) / v2 * mySUSYMQ(2 * J) * mySUSYCKM(J, I) +
                            ZH(1,0)*(myTU(J,0)*mySUSYCKM(0,I) + myTU(J,1)*mySUSYCKM(1,I) 
                            + myTU(J,2)*mySUSYCKM(2,I)))
                            * myRu(i, J + 3) * myRd(j, I).conjugate()
                            
                            + (ZH(0,0) * (myTD(I,0).conjugate() * mySUSYCKM(J,0) 
                            + myTD(I,1).conjugate() * mySUSYCKM(J,1) 
                            + myTD(I,2).conjugate() * mySUSYCKM(J,2))
                            + ZH(1,0) * mySUSY.getMuH() *  sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) /
                            (1 + Eps_J(I) * tanb) * mySUSYCKM(J,I)) * myRu(i,J) * 
                            myRd(j, I + 3).conjugate(); 
                                                            

                }
            }



            VUDHH_cache.assign(i, j, VUDHijH);
            VUDHijH.assign(0., 0., 0);
        }
    }
}

gslpp::complex SUSYMatching::VUDHH(int i, int j) {

    return (VUDHH_cache(i, j));
  
}


gslpp::complex SUSYMatching::DeltaFHL(int j, int i) {


    gslpp::complex DFHL_ji(0., 0., false);
    gslpp::vector<double> myMU2Squarks(6, 0.);
    gslpp::vector<double> myMD2Squarks(6, 0.);
    myMU2Squarks = mySUSY.getMsu2();
    myMD2Squarks = mySUSY.getMsd2();
    gslpp::vector<double> MChi0(4, 0.);
    gslpp::complex Yuj(0., 0., false);
    gslpp::complex Ydi(0., 0., false);
    int m, l;

    Yuj = sqrt(2.) / v2 * mySUSYMQ(2 * j);
    Ydi = sqrt(2.) / v1 * mySUSYMQ(2 * i + 1) /
            (1 + Eps_J(i) * tanb);

    for (m = 0; m < 6; m++) {
        for (l = 0; l < 6; l++) {

            DFHL_ji += VUDHH(m,l) *(-2 * Als / (3. * M_PI) * Mg 
                    * myRu(m, j + 3).conjugate() * myRd(l, i) *
                    Ck(Mg * Mg, myMU2Squarks(m), myMD2Squarks(l), 0)

                    - 1. / (16. * M_PI * M_PI) * Yuj * Ydi *  
                    myRu(m, j).conjugate() * myRd(l, i + 3) *
                    mySUSY.getMuH().conjugate() *
                    Ck(mySUSY.getMuH().abs2(), myMU2Squarks(m), myMD2Squarks(l), 0));

        }

    }
   return (DFHL_ji);
}


void SUSYMatching::Comp_PHRL(){
    
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {

            PHRLCache.assign(j, i, sqrt(2.) / (v * tanb) *
            mySUSYMQ(2 * j) * myCKM(j, i) + DeltaFHL(j, i));
        }
    }
}

gslpp::complex SUSYMatching::PHRL(int j, int i){
    

    return(PHRLCache(j,i));
}

gslpp::complex SUSYMatching::PLRk(int j, int i, int k){
    
    if(k == 0){
        return (PHLR(j,i));
    }
    else if(k == 1){
        return (PGLR(j,i));
    }
    else {
        throw std::runtime_error("Error in PLRk(j,i,k) in SUSYMatching.cpp  ");
    }
    
    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::PRLk(int j, int i, int k){
    
    if(k == 0){
        return (PHRL(j,i));
    }
    else if(k == 1){
        return (PGRL(j,i));
    }
    else {
        throw std::runtime_error("Error in PRLk(j,i,k) in SUSYMatching.cpp  ");
    }
}


/** Charged Higgs vertices genralized both for B and D mixing **/

gslpp::complex SUSYMatching::PRLk(int j, int i, int k, int Dmixingflag) {

    if (Dmixingflag == 0) {
        return (PRLk(j, i, k));

    }
    else if (Dmixingflag == 1) {
        return (PLRk(i, j, k).conjugate());
    }
    else {
        throw std::runtime_error("Error in PRLk(j,i,k,Dmixingflag) in SUSYMatching.cpp  ");
    }

}

gslpp::complex SUSYMatching::PLRk(int j, int i, int k, int Dmixingflag) {


    if (Dmixingflag == 0) {
        return (PLRk(j, i, k));

    }
    else if (Dmixingflag == 1) {
        return (PRLk(i, j, k).conjugate());
    }
    else {
        throw std::runtime_error("Error in PLRk(j,i,k,Dmixingflag) in SUSYMatching.cpp  ");
    }
    
    return (EXIT_FAILURE);
}
////////////////////////////////////////////////////////////////////////////////
////// Double Penguin Functions

gslpp::complex SUSYMatching::xdS(int S){
    
    double M2A = mySUSY.getMHa() * mySUSY.getMHa();
    double M2Z = mySUSY.getMz() * mySUSY.getMz();
    double tan2alpha = 2. * tanb / (1 - tanb * tanb) * (M2A + M2Z) / (M2A - M2Z); 
    double tana = (-1. + sqrt(1. + tan2alpha * tan2alpha)) / tan2alpha;
    double sina = tana / sqrt(1. + tana * tana);
    double cosa = 1. / sqrt(1 + tana * tana);
    
    gslpp::complex Cosa(cosa,0.,false);
    gslpp::complex Sina(sina,0.,false);
    gslpp::complex i(0.,1.,false);


    switch ( S ) {
        case 0:
            return (Cosa);

        case 1:
            return (-Sina);
            
        case 2:
            return ( i * mySUSY.getSinb());
            
        default:
            throw std::runtime_error("Error in xdS(S) in SUSYMatching.cpp  ");
            break;
    }
    
    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::xuS(int S){
    
    double M2A = mySUSY.getMHa() * mySUSY.getMHa();
    double M2Z = mySUSY.getMz() * mySUSY.getMz();
   
    double tan2alpha = 2. * tanb / (1 - tanb * tanb) * (M2A + M2Z) / (M2A - M2Z); 
    double tana = (-1. + sqrt(1. + tan2alpha * tan2alpha)) / tan2alpha;
    double sina = tana / sqrt(1. + tana * tana);
    double cosa = 1. / sqrt(1 + tana * tana);
    
    gslpp::complex Cosa(cosa,0.,false);
    gslpp::complex Sina(sina,0.,false);
    gslpp::complex i(0.,1.,false);
    
    switch ( S ) {
        case 0:
            return (Sina);

        case 1:
            return (Cosa);
            
        case 2:
            return (-i * mySUSY.getCosb());
            
        default:
            throw std::runtime_error("Error in xuS(S) in SUSYMatching.cpp  ");
            break;
    }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::XRLS(int J, int I, int S){
    
    if (J > I) {            
        return (mySUSYMQ(2 * J + 1) / (v1 * (1 + Eps_J(J) * tanb) *
                (1 + Eps_J(J) * tanb))  * Lambda0EpsY(J, I) * 2. / (v2 * v2) 
                * mySUSYMQ(4) * mySUSYMQ(4) *
                (xuS(S) - xdS(S) * tanb));
    }
    else if (J < I) {
        return (XLRS(I, J, S).conjugate());
    }
    else {
        throw std::runtime_error("Error in XRLS(J,I,S) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::XLRS(int J, int I, int S){


    double Y2ut = sqrt(2.) / v2 * mySUSYMQ(4);
    Y2ut *= Y2ut;
    gslpp::complex temp(0.,0.,false);
    temp = 1 + Eps_J(J) * tanb;
    gslpp::complex rJI(0.,0.,false);
        
    rJI = ((1. +  (Eps_J(J) + (Eps_J(I).conjugate() - Eps_J(J).conjugate()) *
            Lambda0EpsY(J, I) / Lambda0EpsY(I, J).conjugate()) * tanb) /
            (1 + Eps_J(I).conjugate() * tanb));   
    
    if (J > I) {    
        return (mySUSYMQ(2 * I + 1) / (v1 * temp.abs2()) *
                Lambda0EpsY(I, J).conjugate() * Y2ut * rJI *
                (xuS(S).conjugate() - xdS(S).conjugate() * tanb));
    }
    else if (J < I) {
        return (XRLS(I, J, S).conjugate());

    }
    else {
        throw std::runtime_error("Error in XLRS(J,I,S) in SUSYMatching.cpp  ");
    }
    
    return (EXIT_FAILURE);
    
}

////////////////////////////////////////////////////////////////////////////////
//// Charged Higgs contribution to Wilson coefficients of Delta F = 2 
//// in the Buras operator basis
//// Q_1 = (q \bar _L gamma_mu b_L) (q \bar _L gamma^mu b_L)
//// in the D - D \bar mixing q -> u and b -> c

gslpp::vector<gslpp::complex> SUSYMatching::CdF2dHp(int b, int q, int Dmixingflag) {

    gslpp::vector<double> M2S(3,0.);
    gslpp::vector<double> MQuarks(6,0.);
    int i;
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        for (i = 0; i < 6; i++) {
            MQuarks(i) += mySUSYMQ(i);
        }
    }
    else if (Dmixingflag == 1) {
        for (i = 0; i < 3; i++) {   
            
            MQuarks(2 * i) += mySUSYMQ(2 * i + 1);
            MQuarks(2 * i + 1) += mySUSYMQ(2 * i);
        }
        myCKM.transpose();
    }
    else {
        throw std::runtime_error("Error in Dmixingflag in SUSYMatching.cpp. Flag can be either 0 or 1  ");
    }

    M2S(0) = mySUSY.getMHl() * mySUSY.getMHl(); 
    M2S(1) = mySUSY.getMHh() * mySUSY.getMHh();
    M2S(2) = mySUSY.getMHa() * mySUSY.getMHa();

    gslpp::vector<gslpp::complex> VCLO(8, 0.);
    gslpp::complex CLO(0.,0.,false);
    
    
    int I, J, k, l, S, O;
    double M2W = mySUSY.Mw_tree();
    double M2H = mySUSY.getMHp();
    double M2I;
    double M2J;
    M2W *= M2W;
    M2H *= M2H;

    gslpp::vector<double> M2Hk(2, 0);
    M2Hk(0) = M2H;
    M2Hk(1) = M2W;
    
    int D = Dmixingflag;

    for (O = 1; O < 9; O++) {

        CLO.assign(0., 0., 0);

        if (O == 1) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {
                    
                        M2J = MQuarks(2 * J);
                        M2J *= M2J;
                        CLO +=  gW * gW / (32. * M_PI * M_PI) *
                                myCKM(I, q).conjugate() * myCKM(J, b) * PRLk(J, q, 0, D).conjugate() *
                                PRLk(I, b, 0, D)  * D0N(M2W, M2W, M2I, M2J)

                                -  1. / (32. * M_PI * M_PI) * PRLk(I, q, 0, D).conjugate() *
                                PRLk(J, q, 0, D).conjugate() * PRLk(I, b, 0, D) * PRLk(J, b, 0, D) *
                                Dk(M2W, M2W, M2I, M2J, 2)

                                -  1. / (16. * M_PI * M_PI) * PRLk(I, q, 1, D).conjugate() *
                                PRLk(J, q, 0, D).conjugate() * PRLk(I, b, 0, D) * PRLk(J, b, 1, D) *
                                Dk(M2W, M2W, M2I, M2J, 2);
                }
            }
        }
        else if (O == 2) {
            for (I = 1; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 1; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {
                   
                                    CLO +=  -1. / (32. * M_PI * M_PI) * PLRk(I, q, l, D).conjugate() *
                                    PLRk(J, q, k, D).conjugate() * PRLk(I, b, k, D) * PRLk(J, b, l, D) *
                                    D0N(M2Hk(k), M2Hk(l), M2I, M2J);
                                    
        
                        }
                    }
                }
            }
        }
        else if (O == 4) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {


                        CLO +=  gW * gW / (8. * M_PI * M_PI) *
                                myCKM(I, q).conjugate() * myCKM(J, b) *
                                PLRk(J, q, k, D).conjugate() * PLRk(I, b, k, D) *
                                Dk(M2W, M2Hk(k), M2I, M2J, 2);

                       
                        
                        for (l = 0; l < 2; l++) {

                            CLO +=  -1. / (16. * M_PI * M_PI) * PLRk(I, q, l, D).conjugate() *
                                    PRLk(J, q, k, D).conjugate() * PRLk(I, b, k, D) *
                                    PLRk(J, b, l, D) * D0N(M2Hk(k), M2Hk(l), M2I, M2J);
                        }
                    }
                }
            }
         
           /** The double Penguin contributions are calulated only for the B and K mixing **/

             if (D == 0) {
                for (S = 0; S < 3; S++) {

                    CLO += -XRLS(q, b, S) * XLRS(q, b, S) / M2S(S);

                }
            }
         
            /** end double Penguin contribution **/  

        }
        else if (O == 5) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);  
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {

                            CLO += 1. / 8. * PRLk(I, q, l, D).conjugate() *
                                    PLRk(J, q, k, D).conjugate() * PRLk(I, b, k, D) *
                                    PLRk(J, b, l, D) * Dk(M2Hk(k), M2Hk(l), M2I, M2J, 2);                        
                        }
                    }
                }
            }
        }
        else if (O == 6) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {

                            CLO += -1. / (32. * M_PI * M_PI) * PLRk(I, q, k, D).conjugate() *
                                    PLRk(J, q, l, D).conjugate() * PLRk(I, b, k, D) *
                                    PLRk(J, b, l, D) * Dk(M2Hk(k), M2Hk(l), M2I, M2J, 2);
                   
                        }
                    }
                }
            }
        }
        else if (O == 7) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {

                            CLO += -1. / (32. * M_PI * M_PI) * PRLk(I, q, l, D).conjugate() *
                                    PRLk(J, q, k, D).conjugate() * PLRk(I, b, k, D) * PLRk(J, b, l, D) 
                                    * D0N(M2Hk(k), M2Hk(l), M2I, M2J);      
                        
                        }
                    }
                }
            }
        }

        VCLO.assign(O - 1,CLO);
        
    }

    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Gluino contribution to Wilson coefficients of Delta F = 2 

gslpp::vector<gslpp::complex> SUSYMatching::CdF2dgg(int b, int q, int Dmixingflag) {

    gslpp::matrix<gslpp::complex> myR(6, 6, 0.);
    gslpp::vector<double> myM2Squarks(6, 0.);
    
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
        myR = mySUSY.getRd();
    }
    else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsu2();  
        
        /* in the D mixing Rd -> Ru^* , b -> c , q -> u */
        
        myR = mySUSY.getRu().hconjugate().transpose();
        
    }
    
    gslpp::complex CLO(0., 0., false);
    gslpp::vector<gslpp::complex> VCLO(8, 0.);


    double M2g = Mg*Mg;
    int h, k, O;

    
    for (O = 1; O < 9; O++) {
        
        CLO.assign(0., 0.,0);
        
        if (O == 1) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO +=  -Als * Als * 
                            myR(h, b) * myR(k, b)
                            * myR(h, q).conjugate() * myR(k, q).conjugate() * 
                            (1. / 9. * M2g *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) +
                            11. / 9. * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2));
                }
            }    
        }
        else if (O == 2) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 17. / 18. * M2g * myR(h, b)
                            * myR(k, b) * myR(h, q + 3).conjugate()
                            * myR(k, q + 3).conjugate() *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        }
        else if (O == 3) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += Als * Als * 1. / 6. * M2g * myR(h, b)
                            * myR(k, b) * myR(h, q + 3).conjugate()
                            * myR(k, q + 3).conjugate() *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        }
        else if (O == 4) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 7. / 3. * M2g * myR(h, b)
                            * myR(k, b + 3) * myR(h, q).conjugate()
                            * myR(k, q + 3).conjugate() *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) +
                            Als * Als * 2. / 9. *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2) 
                            * myR(h, b) * myR(k, b + 3) *
                            (6. * myR(h, q).conjugate() * myR(k, q + 3).conjugate()
                            + 11. * myR(k, q).conjugate() * myR(h, q + 3).conjugate())
                    ;
                }
            }
        }
        else if (O == 5) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 1. / 9. * M2g * myR(h, b)
                            * myR(k, b + 3) * myR(h, q).conjugate()
                            * myR(k, q + 3).conjugate() *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) +
                            Als * Als * 10. / 9. *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2) 
                            * myR(h, b) * myR(k, b + 3) *
                            (3. * myR(k, q).conjugate() * myR(h, q + 3).conjugate()
                            - 2. * myR(h, q).conjugate() * myR(k, q + 3).conjugate())
                            ;
                }
            }
        }
        else if (O == 6) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * myR(h, b + 3) * myR(k, b + 3)
                           * myR(h, q + 3).conjugate() * myR(k, q + 3).conjugate() *
                           (1. / 9. * M2g * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) -
                           11. / 9. * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2));
                }
            }
        }  
        else if (O == 7) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 17. / 18. * M2g * myR(h, b + 3)
                           * myR(k, b + 3) * myR(h, q).conjugate()
                           * myR(k, q).conjugate() *
                           Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        }
        else if (O == 8) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += Als * Als * 1. / 6. * M2g * myR(h, b + 3)
                           * myR(k, b + 3) * myR(h, q).conjugate()
                           * myR(k, q).conjugate() *
                           Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        }


          VCLO.assign(O - 1, CLO);
        
    }
        
    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Chargino contribution to Wilson coefficients of Delta F = 2

gslpp::vector<gslpp::complex> SUSYMatching::CdF2dChiChi(int b, int q, int Dmixingflag) {

    gslpp::vector<double> myM2Squarks(6, 0.);

    gslpp::complex CLO(0., 0., false);
    gslpp::vector<gslpp::complex> VCLO(8, 0.);
    int i, j, h, k, O;

    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsu2();
    }
    else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsd2();
    }

    int D = Dmixingflag;
    
    for (O = 1; O < 9; O++) {
        
        CLO.assign(0., 0.,0);
        
        if (O == 1) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i), MChi(j)
                                    * MChi(j), 2) * VdUCL(b, k, j, D) * VdUCL(b, h, i, D)
                                    * VdUCL(q, h, j, D).conjugate()
                                    * VdUCL(q, k, i, D).conjugate();
                        }
                    }
                }
            }
        }
        else if (O == 3) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) *
                                    VdUCR(q, k, j, 1, D).conjugate() *
                                    VdUCR(q, h, i, 1, D).conjugate() *
                                    VdUCL(b, k, i, D) * VdUCL(b, h, j, D) *
                                    D0N(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i),
                                    MChi(j) * MChi(j));
                        }
                    }
                }
            }
        }
        else if (O == 4) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += 1. / (8. * M_PI * M_PI) * VdUCL(q, h, i, D).conjugate() *
                                    VdUCR(q, k, j, 1, D).conjugate() * VdUCL(b, k, i, D) *
                                    VdUCR(b, h, j, 1, D) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i), MChi(j)
                                    * MChi(j), 2);
                        }
                    }
                }
            }
        }
        else if (O == 5) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (16. * M_PI * M_PI) * VdUCL(q, k, j, D).conjugate() *
                                    VdUCR(q, h, i, 1, D).conjugate() * VdUCL(b, k, i, D) *
                                    VdUCR(b, h, j, 1, D) * D0N(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i), MChi(j)
                                    * MChi(j));
                        }
                    }
                }
            }
        }
        else if (O == 6) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i), MChi(j)
                                    * MChi(j), 2) * VdUCR(q, k, i, 1, D).conjugate() *
                                    VdUCR(q, h, j, 1, D).conjugate() * VdUCR(b, h, i, 1, D) *
                                    VdUCR(b, k, j, 1, D);
                        }
                    }
                }
            }
        }
        else if (O == 8) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) * VdUCL(q, k, j, D).conjugate() *
                                    VdUCL(q, h, i, D).conjugate() * VdUCR(b, k, i, 1, D) *
                                    VdUCR(b, h, j, 1, D) *
                                    D0N(myM2Squarks(k), myM2Squarks(h),
                                    MChi(i) * MChi(i), MChi(j) * MChi(j));
                        }
                    }
                }
            }
        }
        
        VCLO.assign(O - 1, CLO);
    }

    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Neutralino contribution to Wilson coefficients of Delta F = 2
/**** In the D mixing b -> c , q -> u ****/

gslpp::vector<gslpp::complex> SUSYMatching::CdF2dChi0Chi0(int b, int q, int Dmixingflag) {  

    gslpp::vector<double> myM2Squarks(6, 0.);
    gslpp::vector<gslpp::complex> VCLO(8, 0.);
    gslpp::complex CLO(0., 0., false);
    int i, j, h, k, O;
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0

    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
    }
    else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsu2();
    }
    
    
    int D = Dmixingflag;
    
    for (O = 1; O < 9; O++) {

        CLO.assign(0., 0., 0);

        if (O == 1) {
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += - VdDNL(q, k, j, 1, D).conjugate()*VdDNL(b, k, i, 1, D)* 
                                    (1./ 32. / M_PI / M_PI * VdDNL(b, h, j, 1, D)
                                    * VdDNL(q, h, i, 1, D).conjugate()
                                    * Dk(myM2Squarks(k),myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) + MChi0(i) * MChi0(j) 
                                    / 64. / M_PI / M_PI * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i), MChi0(j)
                                    * MChi0(j), 0) * VdDNL(b, h, i, 1, D)*
                                    VdDNL(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        }
        else if (O == 2) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += MChi0(i) * MChi0(j) / 32. / M_PI / M_PI *
                                    Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNL(b, k, i, 1, D) * VdDNR(q, k, j, 1, D).conjugate()*
                                    VdDNL(b, h, i, 1, D) * VdDNR(q, h, j, 1, D).conjugate() ;
                        }
                    }
                }
            }
        }
        else if (O == 3) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -MChi0(i) * MChi0(j) / 32. / M_PI / M_PI *
                                    Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 0) * VdDNL(b, k, i, 1, D) *
                                    VdDNR(q, k, j, 1, D).conjugate() *
                                    (VdDNL(b, h, j, 1, D) * VdDNR(q, h, i, 1, D).conjugate() -
                                    VdDNL(b, h, i, 1, D) * VdDNR(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        }
        else if (O == 4) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += 1. / 8. / M_PI / M_PI * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) * VdDNR(b, k, i, 1, D) *
                                    VdDNL(q, k, j, 1, D).conjugate() * (
                                    VdDNL(b, h, j, 1, D) * (VdDNR(q, h, i, 1, D).conjugate()) +
                                    VdDNL(b, h, i, 1, D) * (VdDNR(q, h, j, 1, D).conjugate()));
                        }
                    }
                }
            }
        }
        else if (O == 5) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / 8. / M_PI / M_PI * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) * VdDNL(q, k, j, 1, D).conjugate() *
                                    VdDNR(q, h, j, 1, D).conjugate() * VdDNL(b, k, i, 1, D) *
                                    VdDNR(b, h, i, 1, D) - MChi0(i) * MChi0(j) / 16. / M_PI /
                                    M_PI * Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNL(q, h, i, 1, D).conjugate() * VdDNR(q, k, j, 1, D).conjugate() *
                                    VdDNL(b, h, j, 1, D) * VdDNR(b, k, i, 1, D);
                        }
                    }
                }
            }
        }
        else if (O == 6) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += - VdDNR(q, k, j, 1, D).conjugate()*VdDNR(b, k, i, 1, D)* 
                                    (1./ 32. / M_PI / M_PI * VdDNR(b, h, j, 1, D)
                                    * VdDNR(q, h, i, 1, D).conjugate() 
                                    * Dk(myM2Squarks(k),myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) + MChi0(i) * MChi0(j) 
                                    / 64. / M_PI / M_PI * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i), MChi0(j)
                                    * MChi0(j), 0) * VdDNR(b, h, i, 1, D)*
                                    VdDNR(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        }
        else if (O == 7) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += MChi0(i) * MChi0(j) / 32. / M_PI / M_PI *
                                    Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNR(b, k, i, 1, D) * VdDNL(q, k, j, 1, D).conjugate()*
                                    VdDNR(b, h, i, 1, D) * VdDNL(q, h, j, 1, D).conjugate() ;
                        }
                    }
                }
            }

        }
        else if (O == 8) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                           CLO += -MChi0(i) * MChi0(j) / 32. / M_PI / M_PI *
                                    Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 0) * VdDNR(b, k, i, 1, D) *
                                    VdDNL(q, k, j, 1, D).conjugate() *
                                    (VdDNR(b, h, j, 1, D) * VdDNL(q, h, i, 1, D).conjugate() -
                                    VdDNR(b, h, i, 1, D) * VdDNL(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        } 

        VCLO.assign(O - 1, CLO);

    }
    return (VCLO);

}

////////////////////////////////////////////////////////////////////////////////
//// Neutralino - Gluino contribution to Wilson coefficients of Delta F = 2

gslpp::vector<gslpp::complex> SUSYMatching::CdF2dChi0g(int b, int q, int Dmixingflag) {   


    gslpp::matrix<gslpp::complex> myR(6, 6, 0.);
    gslpp::vector<double> myM2Squarks(6, 0.);
    gslpp::vector<gslpp::complex> VCLO(8, 0.);
    gslpp::complex CLO(0., 0., false);
    int i, h, k, O;
    double M2g = Mg*Mg;


    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0

    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
        myR = mySUSY.getRd();
    }
    else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsu2();
        
        /* in the D mixing Rd -> Ru^* , b -> c , q -> u */
        
        myR = mySUSY.getRu().hconjugate().transpose();
    }
 
    int D = Dmixingflag;
    
    
    for (O = 1; O < 9; O++) {

        CLO.assign(0., 0., 0);

        if (O == 1) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                CLO += -Als * 2. / 3. / 4. / M_PI * Dk(myM2Squarks(k), myM2Squarks(h),
                        MChi0(i) * MChi0(i), M2g, 2) * myR(k, q).conjugate()
                        * myR(h, b) * VdDNL(q, h, i, 1, D).conjugate() *
                        VdDNL(b, k, i, 1, D) - Als / 4. / M_PI * MChi0(i) * Mg / 6. *
                        Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i)
                        , M2g, 0) * (myR(h, q).conjugate() * myR(k, q).conjugate()
                        * VdDNL(b, h, i, 1, D) *  VdDNL(b, k, i, 1, D) 
                        + myR(h, b) * myR(k, b) * VdDNL(q, h, i, 1, D).conjugate() *
                        VdDNL(q, k, i, 1, D).conjugate());

                    }
                }
            }
        }
        else if (O == 2) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO +=  Als /4. / M_PI * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (3. * myR(h, b) *
                                myR(k, q + 3).conjugate() * VdDNL(b, k, i, 1, D) * 
                                VdDNR(q, h, i, 1, D).conjugate()
                                + myR(k, b) * myR(h, b) *
                                VdDNR(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate() +
                                myR(k, q + 3).conjugate() * myR(h, q + 3).conjugate() * 
                                VdDNL(b, k, i, 1, D) * VdDNL(b, h, i, 1, D));

                    }
                }
            }
        }
        else if (O == 3) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO +=  -Als / 4. / M_PI * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b) * myR(k, q + 3).conjugate() * VdDNL(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() - myR(k, b) *
                                myR(h, b) * VdDNR(q, k, i, 1, D).conjugate() *
                                VdDNR(q, h, i, 1, D).conjugate() - myR(k, q + 3).conjugate() *
                                myR(h, q + 3).conjugate() * VdDNL(b, k, i, 1, D) *
                                VdDNL(b, h, i, 1, D));

                    }
                }
            }
        }
        else if (O == 4) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO +=  -Als / 4. / M_PI * 2. / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 2) * (
                                myR(h, b) * myR(k, q).conjugate() * VdDNR(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() +
                                myR(h, b + 3) *
                                myR(k, q + 3).conjugate() * VdDNL(b, k, i, 1, D) * 
                                VdDNL(q, h, i, 1, D).conjugate()
                                - myR(h, b) * myR(k, b + 3) *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate()
                                - myR(h, q).conjugate() * myR(k, q + 3).conjugate() * 
                                VdDNL(b, k, i, 1, D) * VdDNR(b, h, i, 1, D) 
                                - 3. * myR(k, b + 3) * myR(h, b) * 
                                VdDNL(q, h, i, 1, D).conjugate() *VdDNR(q, k, i, 1, D).conjugate() 
                                - 3. * myR(k, q + 3).conjugate() *
                                myR(h, q).conjugate() * VdDNL(b, h, i, 1, D) * VdDNR(b, k, i, 1, D)) +
                                Als / 4. / M_PI * MChi0(i) * Mg *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)* MChi0(i), M2g, 0) * 
                                (myR(h, b) * myR(k, q + 3).conjugate() *
                                VdDNR(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(h, b + 3) * myR(k, q).conjugate() *
                                VdDNL(b, k, i, 1, D) * VdDNR(q, h, i, 1, D).conjugate());

                    }
                }
            }
        }
        else if (O == 5) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO +=  Als / 4. / M_PI * 2. / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 2)*(
                                3. * myR(h, b) * myR(k, q).conjugate() * VdDNR(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() + 3. * myR(h, b + 3) *
                                myR(k, q + 3).conjugate() * VdDNL(b, k, i, 1, D) * 
                                VdDNL(q, h, i, 1, D).conjugate()
                                - 3. * myR(h, b) * myR(k, b + 3) *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate()
                                - 3. * myR(h, q).conjugate() * myR(k, q + 3).conjugate() * 
                                VdDNL(b, k, i, 1, D) * VdDNR(b, h, i, 1, D) - myR(k, b + 3) *
                                myR(h, b) * VdDNL(q, h, i, 1, D).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate() - myR(k, q + 3).conjugate() *
                                myR(h, q).conjugate() * VdDNL(b, h, i, 1, D) * VdDNR(b, k, i, 1, D)) -
                                Als / 4. / M_PI * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0)* (myR(h, b) * myR(k, q + 3).conjugate() *
                                VdDNR(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(h, b + 3) * myR(k, q).conjugate() *
                                VdDNL(b, k, i, 1, D) * VdDNR(q, h, i, 1, D).conjugate());

                    }
                }
            }
        }
        else if (O == 6) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als * 2. / 3. / 4. / M_PI * Dk(myM2Squarks(k), myM2Squarks(h),
                                MChi0(i) * MChi0(i), M2g, 2) * myR(k, q + 3).conjugate()
                                * myR(h, b + 3) * VdDNR(q, h, i, 1, D).conjugate() *
                                VdDNR(b, k, i, 1, D) - Als / 4. / M_PI * MChi0(i) * Mg / 6. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i)
                                , M2g, 0) * (myR(h, q + 3).conjugate() * myR(k, q + 3).conjugate()
                                * VdDNR(b, h, i, 1, D) *  VdDNR(b, k, i, 1, D) 
                                + myR(h, b + 3) * myR(k, b + 3) * VdDNR(q, h, i, 1, D).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate());

                    }
                }
            }
        }
        else if (O == 7) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO +=  Als /4. / M_PI * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (3. * myR(h, b + 3) *
                                myR(k, q).conjugate() * VdDNR(b, k, i, 1, D) * 
                                VdDNL(q, h, i, 1, D).conjugate()
                                + myR(k, b + 3) * myR(h, b + 3) *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(k, q).conjugate() * myR(h, q).conjugate() * 
                                VdDNR(b, k, i, 1, D) * VdDNR(b, h, i, 1, D));

                    }
                }
            }
        }
        else if (O == 8) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als / 4. / M_PI * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b + 3) * myR(k, q).conjugate() * VdDNR(b, k, i, 1, D) *
                                VdDNL(q, h, i, 1, D).conjugate() - myR(k, b + 3) *
                                myR(h, b + 3) * VdDNL(q, k, i, 1, D).conjugate() *
                                VdDNL(q, h, i, 1, D).conjugate() - myR(k, q).conjugate() *
                                myR(h, q).conjugate() * VdDNR(b, k, i, 1, D) *
                                VdDNR(b, h, i, 1, D));

                    }
                }
            }
        }
        
        VCLO.assign(O - 1, CLO);

    }
    return (VCLO);
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////// 


/******************************************************************************/
std::vector<WilsonCoefficient>& SUSYMatching::CMdbd2() {

    int i;
    gslpp::vector<gslpp::complex> CdF2dHpT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dggT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChigT(8, 0.);

    if (mySUSY.IsFlag_h()) CdF2dHpT = CdF2dHp(0, 2, 0);
    if (mySUSY.IsFlag_ch()) CdF2dChiChiT = CdF2dChiChi(0, 2, 0);
    if (mySUSY.IsFlag_g()) CdF2dggT = CdF2dgg(0, 2, 0);
    if (mySUSY.IsFlag_ne()) CdF2dChi0Chi0T = CdF2dChi0Chi0(0, 2, 0);
    if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) CdF2dChigT = CdF2dChi0g(0, 2, 0);
    
    vmdbd2 = StandardModelMatching::CMdbd2();
    
    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdbd2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdbd2Hp.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbd2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdbd2.push_back(mcdbd2Hp);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {

                mcdbd2gg.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbd2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdbd2.push_back(mcdbd2gg);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdbd2ChiChi.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbd2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdbd2.push_back(mcdbd2ChiChi);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdbd2Chi0Chi0.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbd2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0Chi0);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdbd2Chi0g.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                mcdbd2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0g);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdbd2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdbd2HpT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbd2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2HpT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {

                mcdbd2ggT.setMu(Q_S);
                for(i = 0; i < 3 ; i++){
                mcdbd2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2ggT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdbd2ChiChiT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbd2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2ChiChiT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {

                mcdbd2Chi0Chi0T.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbd2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0Chi0T);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");
    }
    switch (mcdbd2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdbd2Chi0gT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbd2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0gT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbd2()");

    }
    
    return (vmdbd2);
}

/******************************************************************************/

std::vector<WilsonCoefficient>& SUSYMatching::CMdbs2() {
    int i;
    gslpp::vector<gslpp::complex> CdF2dHpT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dggT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFlag_h()) CdF2dHpT = CdF2dHp(1, 2, 0);
    if (mySUSY.IsFlag_ch()) CdF2dChiChiT = CdF2dChiChi(1, 2, 0);
    if (mySUSY.IsFlag_g()) CdF2dggT = CdF2dgg(1, 2, 0);  
    if (mySUSY.IsFlag_ne()) CdF2dChi0Chi0T = CdF2dChi0Chi0(1, 2, 0);
    if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) CdF2dChigT = CdF2dChi0g(1, 2, 0);
    
    vmdbs2 = StandardModelMatching::CMdbs2();

    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdbs2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdbs2Hp.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbs2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdbs2.push_back(mcdbs2Hp);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdbs2gg.setMu(Q_S);
                for(i = 0; i < 5 ; i++){
                mcdbs2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdbs2.push_back(mcdbs2gg);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdbs2ChiChi.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbs2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdbs2.push_back(mcdbs2ChiChi);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdbs2Chi0Chi0.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbs2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0Chi0);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdbs2Chi0g.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdbs2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0g);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
  
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdbs2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdbs2HpT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbs2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2HpT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdbs2ggT.setMu(Q_S);
                for(i = 0; i < 3 ; i++){
                mcdbs2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2ggT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdbs2ChiChiT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbs2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2ChiChiT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdbs2Chi0Chi0T.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbs2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0Chi0T);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    switch (mcdbs2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdbs2Chi0gT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdbs2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0gT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdbs2()");
    }
    
    
    return (vmdbs2);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<WilsonCoefficient>& SUSYMatching::CMdk2() {

    int i;
    gslpp::vector<gslpp::complex> CdF2dHpT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dggT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFlag_h()) CdF2dHpT = CdF2dHp(0, 1, 0);
    if (mySUSY.IsFlag_ch()) CdF2dChiChiT = CdF2dChiChi(0, 1, 0);
    if (mySUSY.IsFlag_g()) CdF2dggT = CdF2dgg(0, 1, 0);
    if (mySUSY.IsFlag_ne()) CdF2dChi0Chi0T = CdF2dChi0Chi0(0, 1, 0);
    if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) CdF2dChigT = CdF2dChi0g(0, 1, 0);
    
    vmdk2 = StandardModelMatching::CMdk2();

    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdk2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdk2Hp.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdk2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdk2.push_back(mcdk2Hp);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdk2gg.setMu(Q_S);
                for(i = 0; i < 5 ; i++){
                mcdk2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdk2.push_back(mcdk2gg);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdk2ChiChi.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdk2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdk2.push_back(mcdk2ChiChi);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdk2Chi0Chi0.setMu(Q_S);    
                for (i = 0; i < 5; i++) {
                    mcdk2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdk2.push_back(mcdk2Chi0Chi0);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdk2Chi0g.setMu(Q_S);    
                for (i = 0; i < 5; i++) {
                    mcdk2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdk2.push_back(mcdk2Chi0g);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdk2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdk2HpT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdk2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2HpT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdk2ggT.setMu(Q_S);
                for(i = 0; i < 3 ; i++){
                mcdk2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2ggT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdk2ChiChiT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdk2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2ChiChiT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdk2Chi0Chi0T.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdk2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdk2.push_back(mcdk2Chi0Chi0T);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    switch (mcdk2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdk2Chi0gT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdk2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2Chi0gT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdk2()");
    }
    
    
    return (vmdk2);
}


////////////////////////////////////////////////////////////////////////////////

 std::vector<WilsonCoefficient>& SUSYMatching::CMdd2(){

    int i;
    gslpp::vector<gslpp::complex> CdF2dHpT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dggT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<gslpp::complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFlag_h()) CdF2dHpT = CdF2dHp(1, 0, 1);
    if (mySUSY.IsFlag_g()) CdF2dggT = CdF2dgg(1, 0, 1);
    if (mySUSY.IsFlag_ch()) CdF2dChiChiT = CdF2dChiChi(1, 0, 1); 
    if (mySUSY.IsFlag_ne()) CdF2dChi0Chi0T = CdF2dChi0Chi0(1, 0, 1);
    if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) CdF2dChigT = CdF2dChi0g(1, 0, 1);
    
    vmdd2 = StandardModelMatching::CMdd2();
    
    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
   
     switch (mcdd2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdd2Hp.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdd2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdd2.push_back(mcdd2Hp);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
     }
    switch (mcdd2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdd2gg.setMu(Q_S);
                for(i = 0; i < 5 ; i++){
                mcdd2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdd2.push_back(mcdd2gg);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdd2ChiChi.setMu(Q_S);
                for (i = 0; i < 5; i++) {
                    mcdd2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdd2.push_back(mcdd2ChiChi);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdd2Chi0Chi0.setMu(Q_S);    
                for (i = 0; i < 5; i++) {
                    mcdd2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdd2.push_back(mcdd2Chi0Chi0);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdd2Chi0g.setMu(Q_S);    
                for (i = 0; i < 5; i++) {
                    mcdd2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdd2.push_back(mcdd2Chi0g);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdd2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_h()) {
                mcdd2HpT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdd2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2HpT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_g()) {
                mcdd2ggT.setMu(Q_S);
                for(i = 0; i < 3 ; i++){
                mcdd2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2ggT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ch()) {
                mcdd2ChiChiT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdd2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2ChiChiT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFlag_ne()) {
                mcdd2Chi0Chi0T.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdd2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdd2.push_back(mcdd2Chi0Chi0T);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    switch (mcdd2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFlag_g()) || (mySUSY.IsFlag_ne())) {
                mcdd2Chi0gT.setMu(Q_S);
                for (i = 0; i < 3; i++) {
                    mcdd2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2Chi0gT);
            }
            break;
        case FULLNLO:
        case NNLO:
        case FULLNNLO:
        default:
            throw std::runtime_error("Error in SUSYMatching::CMdd2()");
    }
    
    return (vmdd2);
    
}

/******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
/*******************************************************************************
 * Wilson coefficients misiak base for b -> s gamma                            * 
 * operator basis: - magnetic moment operator Q7                               *          
 *                                                                             * 
 * ****************************************************************************/

gslpp::complex SUSYMatching::EpsPrime(int J, int I){
   
    gslpp::matrix<gslpp::complex> myCKM(3, 3, 0.);
    myCKM = mySUSY_CKM();
    gslpp::complex YuJ(0., 0., false);

    YuJ = sqrt(2.) / v2 * mySUSYMQ(2 * J);
    
    
    return (-DeltaFHL(J,I) / (tanb * myCKM(J,I) * YuJ));
}
double SUSYMatching::F7k(double x, int k) {

    if (k == 1) {

        if (fabs(x - 1.) < SUSYLEPS) return (-5. / 48.);

        return ((x * (7 - 5 * x - 8 * x * x)) / (24. * (-1 + x)*(-1 + x)*(-1 + x)) +
                x * x * (-2 + 3 * x) * log(x) / (4. * (-1 + x)*(-1 + x)*(-1 + x)*(-1 + x)));

    }
    else if (k == 2) {

        if (fabs(x - 1.) < SUSYLEPS) return (-7. / 36.);

        return (((3 - 5 * x) * x) / (12. * (-1 + x)*(-1 + x)) +
                (x * (-2 + 3 * x) * log(x)) / (6. * (-1 + x)*(-1 + x)*(-1 + x)));

    }
    else throw std::runtime_error("Error in F7k "); 

}

gslpp::vector <gslpp::complex> SUSYMatching::CalcC7(int b, int q) {

    gslpp::vector<gslpp::complex> VCLO(5, 0.);
    gslpp::complex CLO(0., 0., false);
    double m2top = mySUSYMQ(4) * mySUSYMQ(4);
    double M2Hp = mySUSY.getMHp() * mySUSY.getMHp();
    double M2W = mySUSY.Mw_tree() * mySUSY.Mw_tree();

    gslpp::vector<double> myMU2Squarks(6, 0.);
    myMU2Squarks = mySUSY.getMsu2();
    int j, k;


    VCLO.assign(0, F7k(m2top / M2W, 0) + (Eps_J(2) - EpsPrime(2, 2)) / (1 + Eps_J(2)
            * tanb) * F7k(m2top / M2W, 2));

    VCLO.assign(1, 1. / (3. * tanb * tanb) * F7k(m2top / M2Hp, 1) +
            F7k(m2top / M2Hp, 2) - (EpsPrime(2, 1) + Eps_J(2)) / (1. + Eps_J(2) * tanb) *
            tanb * F7k(m2top / M2Hp, 2));

    for (j = 0; j < 2; j++) {
        for (k = 0; k < 6; k++) {
            CLO += 1. / (3. * gW * gW * myCKM(2, 1).conjugate() * myCKM(2, 2)) *
                    M2W / (MChi(j) * MChi(j)) * (-1. / 2. * VdUCL(q, k, j).conjugate()
                    * VdUCL(b, k, j).conjugate()
                    * F7k(myMU2Squarks(k) / (MChi(j) * MChi(j)), 1) +
                    VdUCL(q, k, j).conjugate() * VdUCR(b, k, j,1) *
                    MChi(j) / mySUSYMQ(5) * 6. * 
                    F7k(myMU2Squarks(k) / (MChi(j) * MChi(j)), 2));
        }
    }

    VCLO.assign(2, CLO);

    return (VCLO);
    
}


 std::vector<WilsonCoefficient>& SUSYMatching::CMbsg(){

    vmcbsg = StandardModelMatching::CMbsg();

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("StandardModel::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(mySUSY.getMuw());

    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
           
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}

/* LEPTON FLAVOUR */

void SUSYMatching::NeutralinoRemixing() {

    gslpp::matrix<gslpp::complex> MN_tmp(4, 4, 0.);
    gslpp::matrix<gslpp::complex> Mdiag(4, 4, 0.);
    for (int a=0;a<4;a++) {
        Mdiag.assign(a,a,MChi0(a));
    }
//    std::cout<<"Old eigenvalues "<<MChi0<<std::endl;
//    std::cout<<"Old mixing matrix "<<myN<<std::endl;
    MN_tmp = myN.transpose() * Mdiag * myN;
    gslpp::matrix<gslpp::complex> cmplxONT(4, 4, 0.);
    MN_tmp.eigensystem(cmplxONT, MNeig);
    for (int a=0;a<4;a++) {
        for (int b=0;b<4;b++) {
            ON.assign(a,b,cmplxONT(b,a).real());
        }
    }
//    std::cout<<"New eigenvalues"<<MNeig<<std::endl;
//    std::cout<<"New mixing matrix "<<ON<<std::endl;
}

    gslpp::vector<gslpp::complex> SUSYMatching::AFunctions(int li_to_lj) {
        //gamma penguin contributions

    NeutralinoRemixing();

    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double piconst = 1.0/(32.0 * pi * pi);
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    sinb = mySUSY.getSinb();
    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);

    gslpp::vector<gslpp::complex> AFunctions(4, 0.);

    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlE.assign(a, x, - (g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 0) + (mE/cdenn)*ON(a, 2)*myRl(x, 3)));
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            NRlTAU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 2) + (mTAU/cdenn)*ON(a, 2)*myRl(x, 5)));
            //     RL + RR TYPE MI
            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2)*myRl(x, 0) + 2.0*ON(a, 0)*ttw*myRl(x, 3)));
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2)*myRl(x, 2) + 2.0*ON(a, 0)*ttw*myRl(x, 5)));
//    Commented expressions might be useful for complex neutralino mixing matrices
//            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2).conjugate()*myRl(x, 0) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 3)));
//            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2).conjugate()*myRl(x, 1) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 4)));
//            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2).conjugate()*myRl(x, 2) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 5)));
        }
    }

    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlE.assign(a, x, - (g2*myV(a, 0)*myRn(x, 0)));
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            CRlTAU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 2)));
            //     LR-TYPE
            CLlE.assign(a, x, g2*mE/cdenc*myU(a, 1).conjugate()*myRn(x, 0));
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
            CLlTAU.assign(a, x, g2*mTAU/cdenc*myU(a, 1).conjugate()*myRn(x, 2));
        }
    }

    //     Definition of y and z - remember they are dimensionless quantities.
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            Lepty.assign(a, x, MNeig(a) * MNeig(a) / mym_se_sq(x) );
        }
    }

    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            Leptz.assign(a, x, MChi(a) * MChi(a) / mym_sn_sq(x) );
        }
    }

    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            if (fabs(1.0 - Lepty(a, x)) > 0.01) {
                Leptfa1.assign(a, x, (1.0/mym_se_sq(x))*(1.0/pow((1.0 - Lepty(a, x)),4.0))*
                                (2.0 - 9.0*Lepty(a, x) + 18.0*pow(Lepty(a, x),2.0) - 11.0*pow(Lepty(a, x),3.0)
                                + 6.0*pow(Lepty(a, x),3.0)*log(Lepty(a, x))) );
                Leptf1.assign(a, x, ((1.0 - 6.0*Lepty(a, x) + 3.0 * pow(Lepty(a, x),2.0) + 
                                      2.0*pow(Lepty(a, x),3.0) - 6.0*pow(Lepty(a,x),2.0)*log(Lepty(a, x))))/
                                    (6.0 * pow((1.0 - Lepty(a,x)),4.0)) );
                Leptf2.assign(a, x, (1.0 - pow(Lepty(a, x),2.0) +  2.0 * Lepty(a, x) * log(Lepty(a, x)))/
                                    (pow((1.0-Lepty(a, x)),3.0)));
            }
            else {
                Leptfa1.assign(a, x, (3.0/2.0)*(1.0/mym_se_sq(x)));
                Leptf1.assign(a, x, 1.0/12.0 - (Lepty(a, x) - 1.0)/30.0);
                Leptf2.assign(a, x, 1.0/3.0 - (Lepty(a, x) - 1.0)/6.0);
            }
        }
    }

    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            if(fabs(1.0-Leptz(a, x)) > 0.01) {
                Leptfa2.assign(a, x, (1.0/mym_sn_sq(x))*(1.0/pow(1.0 - Leptz(a,x),4.0))*
                                (16.0 - 45.0*Leptz(a,x) + 36.0*pow(Leptz(a,x),2.0) - 7.0*pow(Leptz(a,x),3.0)
                                + 6.0*(2.0 - 3.0*Leptz(a,x))*log(Leptz(a,x))) );
                Leptf3.assign(a, x, ((2.0 + 3.0*Leptz(a, x) - 6.0*pow(Leptz(a, x),2.0) 
                                      + pow(Leptz(a, x),3.0) + 6.0*Leptz(a, x)*log(Leptz(a, x)))/
		                     (6.0*pow((1.0 - Leptz(a, x)),4.0))) );
                Leptf4.assign(a, x, ((-3.0 + 4.0*Leptz(a, x) - pow(Leptz(a, x),2.0)
                                      - 2.0*log(Leptz(a, x)))/
                                     pow((1.0 - Leptz(a, x)),3.0)) );
            }
            else {
                Leptfa2.assign(a, x, (-9.0/2.0)*(1.0/mym_sn_sq(x)) );
                Leptf3.assign(a, x, 1.0/12.0 - (Leptz(a, x) - 1.0)/20.0 );
                Leptf4.assign(a, x, 2.0/3.0 - (Leptz(a, x) - 1.0)/2.0 );
            }
        }
    }

    if (li_to_lj == 1) // mu -> (3)e
    {
        //      Neutralino contributions
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                AmpA1RN.assign(a, x, (piconst/18.0)*NLlE(a,x)*NLlMU(a,x).conjugate()*Leptfa1(a,x) );
                AmpA1LN.assign(a, x, (piconst/18.0)*NRlE(a,x)*NRlMU(a,x).conjugate()*Leptfa1(a,x) );
                AmpARN.assign(a, x, piconst*(NRlE(a, x) * NRlMU(a, x).conjugate() * Leptf1(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                             + NLlE(a, x) * NLlMU(a, x).conjugate() * (mE/mMU) * Leptf1(a, x)
//...until here
                                             + NRlE(a, x) * NLlMU(a, x).conjugate() * (MNeig(a)/mMU) * Leptf2(a, x)) /mym_se_sq(x) );
                AmpALN.assign(a, x, piconst*(NLlE(a, x) * NLlMU(a, x).conjugate() * Leptf1(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                             + NRlE(a, x) * NRlMU(a, x).conjugate() * (mE/mMU) * Leptf1(a, x)
//...until here
                                             + NLlE(a, x) * NRlMU(a, x).conjugate() * (MNeig(a)/mMU) * Leptf2(a, x)) /mym_se_sq(x) );
            }
        }
        gslpp::complex A1RN = 0.0;
        gslpp::complex A1LN = 0.0;
        gslpp::complex ARN = 0.0;
        gslpp::complex ALN = 0.0;
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                A1RN = A1RN + AmpA1RN(a,x);
                A1LN = A1LN + AmpA1LN(a,x);
                ARN = ARN + AmpARN(a,x);
                ALN = ALN + AmpALN(a,x);
            }
        }
        //      Chargino contributions
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                AmpA1RC.assign(a, x, -(piconst/18.0)*CLlE(a,x)*CLlMU(a,x).conjugate()*Leptfa2(a,x) );
                AmpA1LC.assign(a, x, -(piconst/18.0)*CRlE(a,x)*CRlMU(a,x).conjugate()*Leptfa2(a,x) );
                AmpARC.assign(a, x, -(piconst/mym_sn_sq(x)) * (CRlE(a, x)*CRlMU(a, x).conjugate() * Leptf3(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                               + CLlE(a, x) * CLlMU(a, x).conjugate() * (mE/mMU) * Leptf3(a, x)
//...until here
                                                               + CRlE(a, x) * CLlMU(a, x).conjugate() * (MChi(a)/mMU) * Leptf4(a, x)) );
                AmpALC.assign(a, x, -(piconst/mym_sn_sq(x)) * (CLlE(a, x) * CLlMU(a, x).conjugate() * Leptf3(a, x) 
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                               + CRlE(a, x) * CRlMU(a, x).conjugate() * (mE/mMU) * Leptf3(a, x)
//...until here
                                                               + CLlE(a, x) * CRlMU(a, x).conjugate() * (MChi(a)/mMU) * Leptf4(a, x)) );
            }
        }
        gslpp::complex A1RC = 0.0;
        gslpp::complex A1LC = 0.0;
        gslpp::complex ARC = 0.0;
        gslpp::complex ALC = 0.0;
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                A1RC = A1RC + AmpA1RC(a,x);
                A1LC = A1LC + AmpA1LC(a,x);
                ARC = ARC + AmpARC(a,x);
                ALC = ALC + AmpALC(a,x);
            }
        }
        //     write AR and AL into a vector for mu->(3)e
        AFunctions.assign(0, A1RN + A1RC );
        AFunctions.assign(1, A1LN + A1LC );
        AFunctions.assign(2, ARN + ARC );
        AFunctions.assign(3, ALN + ALC );

    }

    if (li_to_lj == 2) // tau -> (3)mu
    {
        //      Neutralino contributions
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                AmpTauA1RN.assign(a, x, (piconst/18.0)*NLlMU(a,x)*NLlTAU(a,x).conjugate()*Leptfa1(a,x) );
                AmpTauA1LN.assign(a, x, (piconst/18.0)*NRlMU(a,x)*NRlTAU(a,x).conjugate()*Leptfa1(a,x) );
                AmpTauARN.assign(a, x, piconst * (NRlMU(a, x) * NRlTAU(a, x).conjugate() * Leptf1(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                  + NLlMU(a, x) * NLlTAU(a, x).conjugate() * (mMU/mTAU) * Leptf1(a, x)
//...until here
                                                  + NRlMU(a, x) * NLlTAU(a, x).conjugate() * (MNeig(a)/mTAU) * Leptf2(a, x)) /mym_se_sq(x) );
                AmpTauALN.assign(a, x, piconst * (NLlMU(a, x) * NLlTAU(a, x).conjugate() * Leptf1(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                  + NRlMU(a, x) * NRlTAU(a, x).conjugate() * (mMU/mTAU) * Leptf1(a, x)
//...until here
                                                  + NLlMU(a, x) * NRlTAU(a, x).conjugate() * (MNeig(a)/mTAU) * Leptf2(a, x)) /mym_se_sq(x) );
            }
        }
        gslpp::complex TauA1RN = 0.0;
        gslpp::complex TauA1LN = 0.0;
        gslpp::complex TauARN = 0.0;
        gslpp::complex TauALN = 0.0;
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                TauA1RN = TauA1RN + AmpTauA1RN(a,x);
                TauA1LN = TauA1LN + AmpTauA1LN(a,x);
                TauARN = TauARN + AmpTauARN(a,x);
                TauALN = TauALN + AmpTauALN(a,x);
            }
        }
        //      Chargino contributions
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                AmpTauA1RC.assign(a, x, -(piconst/18.0)*CLlMU(a,x)*CLlTAU(a,x).conjugate()*Leptfa2(a,x) );
                AmpTauA1LC.assign(a, x, -(piconst/18.0)*CRlMU(a,x)*CRlTAU(a,x).conjugate()*Leptfa2(a,x) );
                AmpTauARC.assign(a, x, -piconst / mym_sn_sq(x) * (CRlMU(a, x) * CRlTAU(a, x).conjugate() * Leptf3(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                                  + CLlMU(a, x) * CLlTAU(a, x).conjugate() * (mMU/mTAU) * Leptf3(a, x)
//...until here
                                                                  + CRlMU(a, x) * CLlTAU(a, x).conjugate() * (MChi(a)/mTAU) * Leptf4(a, x)) );
                AmpTauALC.assign(a, x, -piconst / mym_sn_sq(x) * (CLlMU(a, x) * CLlTAU(a, x).conjugate() * Leptf3(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                                  + CRlMU(a, x) * CRlTAU(a, x).conjugate() * (mMU/mTAU) * Leptf3(a, x)
//...until here
                                                                  + CLlMU(a, x) * CRlTAU(a, x).conjugate() * (MChi(a)/mTAU) * Leptf4(a, x)) );
            }
        }
        gslpp::complex TauA1RC = 0.0;
        gslpp::complex TauA1LC = 0.0;
        gslpp::complex TauARC = 0.0;
        gslpp::complex TauALC = 0.0;
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                TauA1RC = TauA1RC + AmpTauA1RC(a,x);
                TauA1LC = TauA1LC + AmpTauA1LC(a,x);
                TauARC = TauARC + AmpTauARC(a,x);
                TauALC = TauALC + AmpTauALC(a,x);
            }
        }
        //     write AR and AL into a vector for tau->(3)mu
        AFunctions.assign(0, TauA1RC + TauA1RN );
        AFunctions.assign(1, TauA1LC + TauA1LN );
        AFunctions.assign(2, TauARC + TauARN );
        AFunctions.assign(3, TauALC + TauALN );
    }

    if (li_to_lj == 3) // tau -> (3)e
    {
        //      Neutralino contributions
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                AmpTEA1RN.assign(a, x, (piconst/18.0)*NLlE(a,x)*NLlTAU(a,x).conjugate()*Leptfa1(a,x) );
                AmpTEA1LN.assign(a, x, (piconst/18.0)*NRlE(a,x)*NRlTAU(a,x).conjugate()*Leptfa1(a,x) );
                AmpTEARN.assign(a, x, piconst * (NRlE(a, x) * NRlTAU(a, x).conjugate() * Leptf1(a,x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                 + NLlE(a, x) * NLlTAU(a, x).conjugate() * (mE/mTAU) * Leptf1(a,x)
//...until here
                                                 + NRlE(a, x) * NLlTAU(a, x).conjugate() * (MNeig(a)/mTAU) * Leptf2(a, x)) / mym_se_sq(x) );
                AmpTEALN.assign(a, x, piconst * (NLlE(a, x) * NLlTAU(a, x).conjugate() * Leptf1(a, x) 
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                 + NRlE(a, x) * NRlTAU(a, x).conjugate() * (mE/mTAU) * Leptf1(a, x)
//...until here
                                                 + NLlE(a, x) * NRlTAU(a, x).conjugate() * (MNeig(a)/mTAU) * Leptf2(a, x)) / mym_se_sq(x) );
            }
        }
        gslpp::complex TEA1RN = 0.0;
        gslpp::complex TEA1LN = 0.0;
        gslpp::complex TEARN = 0.0;
        gslpp::complex TEALN = 0.0;
        for (int a=0;a<4;a++) {
            for (int x=0;x<6;x++) {
                TEA1RN = TEA1RN + AmpTEA1RN(a,x);
                TEA1LN = TEA1LN + AmpTEA1LN(a,x);
                TEARN = TEARN + AmpTEARN(a,x);
                TEALN = TEALN + AmpTEALN(a,x);
            }
        }
        //      Chargino contributions
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                AmpTEA1RC.assign(a, x, -(piconst/18.0)*CLlE(a,x)*CLlTAU(a,x).conjugate()*Leptfa2(a,x) );
                AmpTEA1LC.assign(a, x, -(piconst/18.0)*CRlE(a,x)*CRlTAU(a,x).conjugate()*Leptfa2(a,x) );
                AmpTEARC.assign(a, x, -piconst / mym_sn_sq(x) * (CRlE(a, x) * CRlTAU(a, x).conjugate() * Leptf3(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                                 + CLlE(a, x) * CLlTAU(a, x).conjugate() * (mE/mTAU) * Leptf3(a, x)
//...until here
                                                                 + CRlE(a, x) * CLlTAU(a, x).conjugate() *(MChi(a)/mTAU) * Leptf4(a, x)) );
                AmpTEALC.assign(a, x, -piconst / mym_sn_sq(x) * (CLlE(a, x) * CLlTAU(a, x).conjugate() * Leptf3(a, x)
// The following contribution is absent in PRD 53.2442 (Hisano et al.), but appears in PRD 73.055003 (Arganda & Herrero)...
                                                                 + CRlE(a, x) * CRlTAU(a, x).conjugate() * (mE/mTAU) * Leptf3(a, x)
//...until here
                                                                 + CLlE(a, x) * CRlTAU(a, x).conjugate() * (MChi(a)/mTAU) * Leptf4(a, x)) );
            }
        }
        gslpp::complex TEA1RC = 0.0;
        gslpp::complex TEA1LC = 0.0;
        gslpp::complex TEARC = 0.0;
        gslpp::complex TEALC = 0.0;
        for (int a=0;a<2;a++) {
            for (int x=0;x<3;x++) {
                TEA1RC = TEA1RC + AmpTEA1RC(a, x);
                TEA1LC = TEA1LC + AmpTEA1LC(a,x);
                TEARC = TEARC + AmpTEARC(a,x);
                TEALC = TEALC + AmpTEALC(a,x);
            }
        }
        //     write AR and AL into a vector for tau->(3)e
        AFunctions.assign(0, TEA1RC + TEA1RN );
        AFunctions.assign(1, TEA1LC + TEA1LN );
        AFunctions.assign(2, TEARC + TEARN );
        AFunctions.assign(3, TEALC + TEALN );

    }

    return(AFunctions);
}

gslpp::vector<gslpp::complex> SUSYMatching::BFunctions(int li_to_lj) {
        //box diagram contributions

    NeutralinoRemixing();

    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();

    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);
    double alph = mySUSY.getAle();

    gslpp::vector<gslpp::complex> BFunctions(8, 0.);

    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlE.assign(a, x, - (g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 0) + (mE/cdenn)*ON(a, 2)*myRl(x, 3)));
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            NRlTAU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 2) + (mTAU/cdenn)*ON(a, 2)*myRl(x, 5)));
            //     RL + RR TYPE MI
            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2)*myRl(x, 0) + 2.0*ON(a, 0)*ttw*myRl(x, 3)));
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2)*myRl(x, 2) + 2.0*ON(a, 0)*ttw*myRl(x, 5)));
//    Commented expressions might be useful for complex neutralino mixing matrices
//            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2).conjugate()*myRl(x, 0) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 3)));
//            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2).conjugate()*myRl(x, 1) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 4)));
//            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2).conjugate()*myRl(x, 2) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 5)));
        }
    }

    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlE.assign(a, x, - (g2*myV(a, 0)*myRn(x, 0)));
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            CRlTAU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 2)));
            //     LR-TYPE
            CLlE.assign(a, x, g2*mE/cdenc*myU(a, 1).conjugate()*myRn(x, 0));
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
            CLlTAU.assign(a, x, g2*mTAU/cdenc*myU(a, 1).conjugate()*myRn(x, 2));
        }
    }

    if (li_to_lj == 1) // mu -> 3e
    {
        //      Neutralino contributions
//    J4n(a,b,x,t)=/*J4n(a,b,x,t)*/
//                                                                    *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
//    I4n(a,b,x,t)=/*I4n(a,b,x,t)*/
//                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi)
        gslpp::complex B1nRMu3E = 0.0;
        gslpp::complex B2nRMu3E = 0.0;
        gslpp::complex B3nRMu3E = 0.0;
        gslpp::complex B4nRMu3E = 0.0;
        gslpp::complex B1nLMu3E = 0.0;
        gslpp::complex B2nLMu3E = 0.0;
        gslpp::complex B3nLMu3E = 0.0;
        gslpp::complex B4nLMu3E = 0.0;
        for (int a=0;a<4;a++) {
            for (int b=0;b<4;b++) {
                for (int x=0;x<6;x++) {
                    for (int t=0;t<6;t++) {
                        B1nRMu3E = B1nRMu3E + (1.0/(4.0*pi*alph))*(0.5*NLlMU(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)/*J4n(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                   +MNeig(a)*MNeig(b)*NLlMU(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nRMu3E = B2nRMu3E + (1.0/(4.0*pi*alph))*(0.25*(NLlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x)
                                                                         +NLlMU(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x)
                                                                         -NLlMU(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x))/*J4n(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                   -0.5*MNeig(a)*MNeig(b)*NLlMU(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nRMu3E = B3nRMu3E + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NLlMU(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x)
                                                                                     +0.5*NLlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x))/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nRMu3E = B4nRMu3E + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NLlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B1nLMu3E = B1nLMu3E + (1.0/(4.0*pi*alph))*(0.5*NRlMU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)/*J4n(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                   +MNeig(a)*MNeig(b)*NRlMU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nLMu3E = B2nLMu3E + (1.0/(4.0*pi*alph))*(0.25*(NRlMU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x)
                                                                         +NRlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x)
                                                                         -NRlMU(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x))/*J4n(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                   -0.5*MNeig(a)*MNeig(b)*NRlMU(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nLMu3E = B3nLMu3E + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NRlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x)
                                                                                     +0.5*NRlMU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x))/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nLMu3E = B4nLMu3E + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NRlMU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
        //      Chargino contributions
//    J4c(a,b,x,t)=/*J4c(a,b,x,t)*/
//                                                                    *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
//    I4c(a,b,x,t)=/*I4c(a,b,x,t)*/
//                                                                    *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi)
        gslpp::complex B1cRMu3E = 0.0;
        gslpp::complex B2cRMu3E = 0.0;
        gslpp::complex B3cRMu3E = 0.0;
        gslpp::complex B1cLMu3E = 0.0;
        gslpp::complex B2cLMu3E = 0.0;
        gslpp::complex B3cLMu3E = 0.0;
        for (int a=0;a<2;a++) {
            for (int b=0;b<2;b++) {
                for (int x=0;x<3;x++) {
                    for (int t=0;t<3;t++) {
                        B1cRMu3E = B1cRMu3E + (1.0/(8.0*pi*alph))*CLlMU(a,x)*CLlE(a,t)*CLlE(b,t)*CLlE(b,x)/*J4c(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi);
                        B2cRMu3E = B2cRMu3E + (1.0/(4.0*pi*alph))*(0.25*CLlMU(a,x)*CLlE(a,t)*CRlE(b,t)*CRlE(b,x)/*J4c(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                   -0.5*MChi(a)*MChi(b)*CLlMU(a,x)*CRlE(a,t)*CRlE(b,t)*CLlE(b,x)/*I4c(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cRMu3E = B3cRMu3E + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CLlMU(a,x)*CRlE(a,t)*CLlE(b,t)*CRlE(b,x)/*I4c(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                        B1cLMu3E = B1cLMu3E + (1.0/(8.0*pi*alph))*(CRlMU(a,x)*CRlE(a,t)*CRlE(b,t)*CRlE(b,x)/*J4c(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi));
                        B2cLMu3E = B2cLMu3E + (1.0/(4.0*pi*alph))*(0.25*CRlMU(a,x)*CRlE(a,t)*CLlE(b,t)*CLlE(b,x)/*J4c(a,b,x,t)*/
                                                                    *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                   -0.5*MChi(a)*MChi(b)*CRlMU(a,x)*CLlE(a,t)*CLlE(b,t)*CRlE(b,x)/*I4c(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cLMu3E = B3cLMu3E + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CRlMU(a,x)*CLlE(a,t)*CRlE(b,t)*CLlE(b,x)/*I4c(a,b,x,t)*/
                                                                    *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
                
//     write BR and BL into a vector for mu->3e
        BFunctions.assign(0, B1nRMu3E + B1cRMu3E );
        BFunctions.assign(1, B1nLMu3E + B1cLMu3E );
        BFunctions.assign(2, B2nRMu3E + B2cRMu3E );
        BFunctions.assign(3, B2nLMu3E + B2cLMu3E );
        BFunctions.assign(4, B3nRMu3E + B3cRMu3E );
        BFunctions.assign(5, B3nLMu3E + B3cLMu3E );
        BFunctions.assign(6, B4nRMu3E );
        BFunctions.assign(7, B4nLMu3E );
    }
    if (li_to_lj == 2) // tau -> 3mu
    {
        //      Neutralino contributions
        gslpp::complex B1nRTau3Mu = 0.0;
        gslpp::complex B2nRTau3Mu = 0.0;
        gslpp::complex B3nRTau3Mu = 0.0;
        gslpp::complex B4nRTau3Mu = 0.0;
        gslpp::complex B1nLTau3Mu = 0.0;
        gslpp::complex B2nLTau3Mu = 0.0;
        gslpp::complex B3nLTau3Mu = 0.0;
        gslpp::complex B4nLTau3Mu = 0.0;
        for (int a=0;a<4;a++) {
            for (int b=0;b<4;b++) {
                for (int x=0;x<6;x++) {
                    for (int t=0;t<6;t++) {
                        B1nRTau3Mu = B1nRTau3Mu + (1.0/(4.0*pi*alph))*(0.5*NLlTAU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)/*J4n(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                       +MNeig(a)*MNeig(b)*NLlTAU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nRTau3Mu = B2nRTau3Mu + (1.0/(4.0*pi*alph))*(0.25*(NLlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)
                                                                             +NLlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NRlMU(b,x)
                                                                             -NLlTAU(a,x)*NRlMU(a,t)*NRlMU(b,t)*NLlMU(b,x))/*J4n(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                       -0.5*MNeig(a)*MNeig(b)*NLlTAU(a,x)*NRlMU(a,t)*NRlMU(b,t)*NLlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nRTau3Mu = B3nRTau3Mu + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NLlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NRlMU(b,x)
                                                                                         +0.5*NLlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NRlMU(b,x))/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nRTau3Mu = B4nRTau3Mu + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NLlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B1nLTau3Mu = B1nLTau3Mu + (1.0/(4.0*pi*alph))*(0.5*NRlTAU(a,x)*NRlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)/*J4n(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                       +MNeig(a)*MNeig(b)*NRlTAU(a,x)*NRlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nLTau3Mu = B2nLTau3Mu + (1.0/(4.0*pi*alph))*(0.25*(NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)
                                                                           +NRlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NLlMU(b,x)
                                                                           -NRlTAU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NRlMU(b,x))/*J4n(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                     -0.5*MNeig(a)*MNeig(b)*NRlTAU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NRlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nLTau3Mu = B3nLTau3Mu + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NRlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NLlMU(b,x)
                                                                                       +0.5*NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x))/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nLTau3Mu = B4nLTau3Mu + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)/*I4n(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
        //      Chargino contributions
        gslpp::complex B1cRTau3Mu = 0.0;
        gslpp::complex B2cRTau3Mu = 0.0;
        gslpp::complex B3cRTau3Mu = 0.0;
        gslpp::complex B1cLTau3Mu = 0.0;
        gslpp::complex B2cLTau3Mu = 0.0;
        gslpp::complex B3cLTau3Mu = 0.0;
        for (int a=0;a<2;a++) {
            for (int b=0;b<2;b++) {
                for (int x=0;x<3;x++) {
                    for (int t=0;t<3;t++) {
                        B1cRTau3Mu = B1cRTau3Mu + (1.0/(8.0*pi*alph))*CLlTAU(a,x)*CLlMU(a,t)*CLlMU(b,t)*CLlMU(b,x)/*J4c(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi);
                        B2cRTau3Mu = B2cRTau3Mu + (1.0/(4.0*pi*alph))*(0.25*CLlTAU(a,x)*CLlMU(a,t)*CRlMU(b,t)*CRlMU(b,x)/*J4c(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                       -0.5*MChi(a)*MChi(b)*CLlTAU(a,x)*CRlMU(a,t)*CRlMU(b,t)*CLlMU(b,x)/*I4c(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cRTau3Mu = B3cRTau3Mu + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CLlTAU(a,x)*CRlMU(a,t)*CLlMU(b,t)*CRlMU(b,x)/*I4c(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                        B1cLTau3Mu = B1cLTau3Mu + (1.0/(8.0*pi*alph))*CRlTAU(a,x)*CRlMU(a,t)*CRlMU(b,t)*CRlMU(b,x)/*J4c(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi);
                        B2cLTau3Mu = B2cLTau3Mu + (1.0/(4.0*pi*alph))*(0.25*CRlTAU(a,x)*CRlMU(a,t)*CLlMU(b,t)*CLlMU(b,x)/*J4c(a,b,x,t)*/
                                                                        *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                       -0.5*MChi(a)*MChi(b)*CRlTAU(a,x)*CLlMU(a,t)*CLlMU(b,t)*CRlMU(b,x)/*I4c(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cLTau3Mu = B3cLTau3Mu + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CRlTAU(a,x)*CLlMU(a,t)*CRlMU(b,t)*CLlMU(b,x)/*I4c(a,b,x,t)*/
                                                                        *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
        //     write BR and BL into a vector for tau->3mu
        BFunctions.assign(0, B1nRTau3Mu + B1cRTau3Mu );
        BFunctions.assign(1, B1nLTau3Mu + B1cLTau3Mu );
        BFunctions.assign(2, B2nRTau3Mu + B2cRTau3Mu );
        BFunctions.assign(3, B2nLTau3Mu + B2cLTau3Mu );
        BFunctions.assign(4, B3nRTau3Mu + B3cRTau3Mu );
        BFunctions.assign(5, B3nLTau3Mu + B3cLTau3Mu );
        BFunctions.assign(6, B4nRTau3Mu );
        BFunctions.assign(7, B4nLTau3Mu );
    }
    if (li_to_lj == 3) // tau -> 3e
    {
        //      Neutralino contributions
        gslpp::complex B1nRTau3E = 0.0;
        gslpp::complex B2nRTau3E = 0.0;
        gslpp::complex B3nRTau3E = 0.0;
        gslpp::complex B4nRTau3E = 0.0;
        gslpp::complex B1nLTau3E = 0.0;
        gslpp::complex B2nLTau3E = 0.0;
        gslpp::complex B3nLTau3E = 0.0;
        gslpp::complex B4nLTau3E = 0.0;
        for (int a=0;a<4;a++) {
            for (int b=0;b<4;b++) {
                for (int x=0;x<6;x++) {
                    for (int t=0;t<6;t++) {
                        B1nRTau3E = B1nRTau3E + (1.0/(4.0*pi*alph))*(0.5*NLlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)/*J4n(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                     +MNeig(a)*MNeig(b)*NLlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nRTau3E = B2nRTau3E + (1.0/(4.0*pi*alph))*(0.25*(NLlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x)
                                                                           +NLlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x)
                                                                           -NLlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x))/*J4n(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                     -0.5*MNeig(a)*MNeig(b)*NLlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nRTau3E = B3nRTau3E + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NLlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x)
                                                                                       +0.5*NLlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x))/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nRTau3E = B4nRTau3E + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NLlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B1nLTau3E = B1nLTau3E + (1.0/(4.0*pi*alph))*(0.5*NRlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)/*J4n(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                     +MNeig(a)*MNeig(b)*NRlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B2nLTau3E = B2nLTau3E + (1.0/(4.0*pi*alph))*(0.25*(NRlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x)
                                                                           +NRlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x)
                                                                           -NRlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x))/*J4n(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(4.0*pi*pi)
                                                                     -0.5*MNeig(a)*MNeig(b)*NRlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi));
                        B3nLTau3E = B3nLTau3E + (1.0/(4.0*pi*alph))*MNeig(a)*MNeig(b)*(NRlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x)
                                                                                       +0.5*NRlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x))/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                        B4nLTau3E = B4nLTau3E + (1.0/(32.0*pi*alph))*MNeig(a)*MNeig(b)*NRlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x)/*I4n(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
        //      Chargino contributions
        gslpp::complex B1cRTau3E = 0.0;
        gslpp::complex B2cRTau3E = 0.0;
        gslpp::complex B3cRTau3E = 0.0;
        gslpp::complex B1cLTau3E = 0.0;
        gslpp::complex B2cLTau3E = 0.0;
        gslpp::complex B3cLTau3E = 0.0;
        for (int a=0;a<2;a++) {
            for (int b=0;b<2;b++) {
                for (int x=0;x<3;x++) {
                    for (int t=0;t<3;t++) {
                        B1cRTau3E = B1cRTau3E + (1.0/(8.0*pi*alph))*CLlTAU(a,x)*CLlE(a,t)*CLlE(b,t)*CLlE(b,x)/*J4c(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi);
                        B2cRTau3E = B2cRTau3E + (1.0/(4.0*pi*alph))*(0.25*CLlTAU(a,x)*CLlE(a,t)*CRlE(b,t)*CRlE(b,x)/*J4c(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                     -0.5*MChi(a)*MChi(b)*CLlTAU(a,x)*CRlE(a,t)*CRlE(b,t)*CLlE(b,x)/*I4c(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cRTau3E = B3cRTau3E + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CLlTAU(a,x)*CRlE(a,t)*CLlE(b,t)*CRlE(b,x)/*I4c(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                        B1cLTau3E = B1cLTau3E + (1.0/(8.0*pi*alph))*CRlTAU(a,x)*CRlE(a,t)*CRlE(b,t)*CRlE(b,x)/*J4c(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi);
                        B2cLTau3E = B2cLTau3E + (1.0/(4.0*pi*alph))*(0.25*CRlTAU(a,x)*CRlE(a,t)*CLlE(b,t)*CLlE(b,x)/*J4c(a,b,x,t)*/
                                                                      *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(4.0*pi*pi)
                                                                     -0.5*MChi(a)*MChi(b)*CRlTAU(a,x)*CLlE(a,t)*CLlE(b,t)*CRlE(b,x)/*I4c(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi));
                        B3cLTau3E = B3cLTau3E + (1.0/(4.0*pi*alph))*MChi(a)*MChi(b)*CRlTAU(a,x)*CLlE(a,t)*CRlE(b,t)*CLlE(b,x)/*I4c(a,b,x,t)*/
                                                                      *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(t))/(16.0*pi*pi);
                    }
                }
            }
        }
        //     write BR and BL into a vector for tau->3e
        BFunctions.assign(0, B1nRTau3E + B1cRTau3E );
        BFunctions.assign(1, B1nLTau3E + B1cLTau3E );
        BFunctions.assign(2, B2nRTau3E + B2cRTau3E );
        BFunctions.assign(3, B2nLTau3E + B2cLTau3E );
        BFunctions.assign(4, B3nRTau3E + B3cRTau3E );
        BFunctions.assign(5, B3nLTau3E + B3cLTau3E );
        BFunctions.assign(6, B4nRTau3E );
        BFunctions.assign(7, B4nLTau3E );
    }
//    if (li_to_lj == 4) // tau -> muee
//    {
//        //     write BR and BL into a vector for tau->muee
//        BFunctions.assign(0, 0. );
//        BFunctions.assign(1, 0. );
//        BFunctions.assign(2, 0. );
//        BFunctions.assign(3, 0. );
//        BFunctions.assign(4, 0. );
//        BFunctions.assign(5, 0. );
//        BFunctions.assign(6, 0. );
//        BFunctions.assign(7, 0. );
//    }

    return(BFunctions);
}

int SUSYMatching::delta_ab(int a, int b) {
    if(a == b) return 1;
    else return 0;
}

////gslpp::complex SUSYMatching::PVC11(double m12, double m22, double m32) {
////    return 0.5*(-(m12-m22)*(m12-m32)*(m22-m32)*m12
////                +m12*m12*m22*(2.0*m12-m22)*log(m12/m22) +m12*m12*m32*(-2.0*m12+m32)*log(m12/m32)
////                +m22*m32*(-2.0*m12+m22)*(-2.0*m12+m32)*log(m22/m32))/((m12-m22)*(m12-m22)*(m12-m32)*(m12-m32)*(m22-m32));
////}
////
////gslpp::complex SUSYMatching::PVC12(double m12, double m22, double m32) {
////    return 0.5*((m12-m22)*(m12-m32)*(m22-m32)*m32 +m22*m22*m32*(2.0*m12-m32)*log(m22/m32)
////                +m12*m12*(m22*m22*log(m12/m22)+m32*(-2.0*m22+m32)*log(m12/m32)))/((m12-m22)*(m12-m32)*(m12-m32)*(m22-m32)*(m22-m32));
////}

gslpp::vector<gslpp::complex> SUSYMatching::BHFunctions(int li_to_lj) {
        //Higgs penguin contributions from PhysRevD.73.055003


    // To do:
    // Are the trilinear couplings from (A24) correct?
    // Should we update the missing parameters in the update function (sinalpha, A_l, ...)?

    NeutralinoRemixing();

    double MZ = mySUSY.getMz();
    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double piconst = 1.0/(32.0 * pi * pi);
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    double cos2b = 2.0*cosb*cosb-1.0;
    gslpp::complex sina = mySUSY.getSaeff();
    gslpp::complex cosa = sqrt(1.0-sina*sina);
    gslpp::complex sinapb = sina*cosb+cosa*sinb;
    gslpp::complex cosapb = cosa*cosb-sina*sinb;
    gslpp::complex cosbma = cosb*cosa+sinb*sina;
    gslpp::complex sinbma = sinb*cosa-cosb*sina;
    double mh = mySUSY.getMHl();
    double mH = mySUSY.getMHh();
    double mA = mySUSY.getMHa();
    gslpp::complex M1 = mySUSY.getM1();
    gslpp::complex M2 = mySUSY.getM2();
    gslpp::complex muH = mySUSY.getMuH();
    TEhat = mySUSY.getTEhat();

    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);
    double alph = mySUSY.getAle();

    gslpp::vector<gslpp::complex> BHFunctions(4, 0.);

    std::cout<<"MZ"<<MZ<<std::endl;
    std::cout<<"MW"<<MW<<std::endl;
    std::cout<<"sw2"<<sw2<<std::endl;
    std::cout<<"mE"<<mE<<std::endl;
    std::cout<<"mMU"<<mMU<<std::endl;
    std::cout<<"mTAU"<<mTAU<<std::endl;
    std::cout<<"sina"<<sina<<std::endl;
    std::cout<<"tanb"<<tanb<<std::endl;
    std::cout<<"mh"<<mh<<std::endl;
    std::cout<<"mH"<<mH<<std::endl;
    std::cout<<"mA"<<mA<<std::endl;
    std::cout<<"M1"<<M1<<std::endl;
    std::cout<<"M2"<<M2<<std::endl;
    std::cout<<"muH"<<muH<<std::endl;
    std::cout<<"TEhat"<<TEhat<<std::endl;
    std::cout<<"g2"<<g2<<std::endl;
    std::cout<<"alph"<<alph<<std::endl;
    std::cout<<"myRl"<<myRl<<std::endl;
    std::cout<<"ON"<<ON<<std::endl;
    std::cout<<"myRn"<<myRn<<std::endl;
    std::cout<<"myV"<<myV<<std::endl;
    std::cout<<"myU"<<myU<<std::endl;
    std::cout<<"mym_se_sq"<<mym_se_sq<<std::endl;
    std::cout<<"mym_sn_sq"<<mym_sn_sq<<std::endl;
    std::cout<<"MNeig"<<MNeig<<std::endl;
    std::cout<<"MChi"<<MChi<<std::endl;
gslpp::complex HpengMuEEENL0temp=0.;


    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlE.assign(a, x, - (g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 0) + (mE/cdenn)*ON(a, 2)*myRl(x, 3)));
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            NRlTAU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 2) + (mTAU/cdenn)*ON(a, 2)*myRl(x, 5)));
            //     RL + RR TYPE MI
            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2)*myRl(x, 0) + 2.0*ON(a, 0)*ttw*myRl(x, 3)));
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2)*myRl(x, 2) + 2.0*ON(a, 0)*ttw*myRl(x, 5)));
//    Commented expressions might be useful for complex neutralino mixing matrices
//            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2).conjugate()*myRl(x, 0) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 3)));
//            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2).conjugate()*myRl(x, 1) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 4)));
//            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2).conjugate()*myRl(x, 2) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 5)));
        }
    }

    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlE.assign(a, x, - (g2*myV(a, 0)*myRn(x, 0)));
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            CRlTAU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 2)));
            //     LR-TYPE
            CLlE.assign(a, x, g2*mE/cdenc*myU(a, 1).conjugate()*myRn(x, 0));
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
            CLlTAU.assign(a, x, g2*mTAU/cdenc*myU(a, 1).conjugate()*myRn(x, 2));
        }
    }

    gslpp::vector<gslpp::complex> sigma1(3, 0.);
    gslpp::vector<gslpp::complex> sigma2(3, 0.);
    gslpp::vector<gslpp::complex> sigma3(3, 0.);
    gslpp::vector<gslpp::complex> sigma4(3, 0.);
    gslpp::vector<gslpp::complex> sigma5(3, 0.);
    sigma1.assign(0, sina);
    sigma1.assign(1, -cosa);
    sigma1.assign(2, sinb*complex::i());
    sigma2.assign(0, cosa);
    sigma2.assign(1, sina);
    sigma2.assign(2, -cosb*complex::i());
    sigma3.assign(0, sinapb);
    sigma3.assign(1, -cosapb);
    sigma3.assign(2, 0.);
    sigma4.assign(0, -sina);
    sigma4.assign(1, cosa);
    sigma4.assign(2, 0.);
    sigma5.assign(0, -cosbma);
    sigma5.assign(1, sinbma);
    sigma5.assign(2, cos2b*complex::i());

    gslpp::matrix<gslpp::complex> Qpp(4, 4, 0.), Rpp(4, 4, 0.);
    gslpp::matrix<gslpp::complex> DL0(4, 4, 0.), DR0(4, 4, 0.), DL1(4, 4, 0.), DR1(4, 4, 0.), DL2(4, 4, 0.), DR2(4, 4, 0.);
    for (int a=0;a<4;a++) {
        for (int b=0;b<4;b++) {
            Qpp.assign(a, b, 0.5*(ON(a,2)*(ON(b,1)-ttw*ON(b,0))+ON(b,2)*(ON(a,1)-ttw*ON(a,0))) );
            Rpp.assign(a, b, (M2.conjugate()*ON(a,1)*ON(b,1) +M1.conjugate()*ON(a,0)*ON(b,0) -muH.conjugate()*(ON(a,2)*ON(b,3)+ON(a,3)*ON(b,2)))/(2.0*MW) );
            DL0.assign(b, a, -g2/sinb * (Qpp(a,b).conjugate()*sigma5(0) -Rpp(a,b).conjugate()*sigma2(0) +MNeig(a)/(2.0*MW)*sigma2(0)*delta_ab(a,b)) );
            DR0.assign(b, a, DL0(b,a).conjugate() );
            DL1.assign(b, a, -g2/sinb * (Qpp(a,b).conjugate()*sigma5(1) -Rpp(a,b).conjugate()*sigma2(1) +MNeig(a)/(2.0*MW)*sigma2(1)*delta_ab(a,b)) );
            DR1.assign(b, a, DL1(b,a).conjugate() );
            DL2.assign(b, a, -g2/sinb * (Qpp(a,b).conjugate()*sigma5(2) -Rpp(a,b).conjugate()*sigma2(2) +MNeig(a)/(2.0*MW)*sigma2(2)*delta_ab(a,b)) );
            DR2.assign(b, a, DL2(b,a).conjugate() );
        }
    }

    gslpp::matrix<gslpp::complex> Qch(2, 2, 0.), Rch(2, 2, 0.);
    gslpp::matrix<gslpp::complex> WL0(2, 2, 0.), WR0(2, 2, 0.), WL1(2, 2, 0.), WR1(2, 2, 0.), WL2(2, 2, 0.), WR2(2, 2, 0.);
    for (int a=0;a<2;a++) {
        for (int b=0;b<2;b++) {
            Qch.assign(a, b, myU(a,1)*myV(b,0)/sqrt(2.0) );
            Rch.assign(a, b, (M2.conjugate()*myU(a,0)*myV(b,0) +muH.conjugate()*myU(a,1)*myV(b,1))/(2.0*MW) );
            WL0.assign(a, b, -g2/sinb * (Qch(b,a)*sigma5(0).conjugate() -Rch(b,a)*sigma2(0).conjugate() +MChi(a)/(2.0*MW)*sigma2(0).conjugate()*delta_ab(a,b)) );
            WR0.assign(a, b, WL0(b,a).conjugate() );
            WL1.assign(a, b, -g2/sinb * (Qch(b,a)*sigma5(1).conjugate() -Rch(b,a)*sigma2(1).conjugate() +MChi(a)/(2.0*MW)*sigma2(1).conjugate()*delta_ab(a,b)) );
            WR1.assign(a, b, WL1(b,a).conjugate() );
            WL2.assign(a, b, -g2/sinb * (Qch(b,a)*sigma5(2).conjugate() -Rch(b,a)*sigma2(2).conjugate() +MChi(a)/(2.0*MW)*sigma2(2).conjugate()*delta_ab(a,b)) );
            WR2.assign(a, b, WL2(b,a).conjugate() );
        }
    }

    gslpp::vector<gslpp::complex> gLLE(3, 0.), gRRE(3, 0.), gLRE(3, 0.), gRLE(3, 0.);
    gslpp::vector<gslpp::complex> gLLMU(3, 0.), gRRMU(3, 0.), gLRMU(3, 0.), gRLMU(3, 0.);
    gslpp::vector<gslpp::complex> gLLTAU(3, 0.), gRRTAU(3, 0.), gLRTAU(3, 0.), gRLTAU(3, 0.);
    gslpp::vector<gslpp::complex> gLLNU(3, 0.);
    for (int p=0;p<3;p++) {
        gLLE.assign(p, MZ/ctw*sigma3(p)*(0.5-sw2) + mE*mE/(MW*cosb)*sigma4(p));
        gLLMU.assign(p, MZ/ctw*sigma3(p)*(0.5-sw2) + mMU*mMU/(MW*cosb)*sigma4(p));
        gLLTAU.assign(p, MZ/ctw*sigma3(p)*(0.5-sw2) + mTAU*mTAU/(MW*cosb)*sigma4(p));
        gRRE.assign(p, MZ/ctw*sigma3(p)*sw2 + mE*mE/(MW*cosb)*sigma4(p));
        gRRMU.assign(p, MZ/ctw*sigma3(p)*sw2 + mMU*mMU/(MW*cosb)*sigma4(p));
        gRRTAU.assign(p, MZ/ctw*sigma3(p)*sw2 + mTAU*mTAU/(MW*cosb)*sigma4(p));
        gLRE.assign(p, (-sigma1(p)*TEhat(0,0)/mE*v1/sqrt(2.0)-sigma2(p).conjugate()*muH)*mE/(2.0*MW*cosb));
        gLRMU.assign(p, (-sigma1(p)*TEhat(1,1)/mMU*v1/sqrt(2.0)-sigma2(p).conjugate()*muH)*mMU/(2.0*MW*cosb));
        gLRTAU.assign(p, (-sigma1(p)*TEhat(2,2)/mTAU*v1/sqrt(2.0)-sigma2(p).conjugate()*muH)*mTAU/(2.0*MW*cosb));
        gRLE.assign(p, gLRE(p).conjugate());
        gRLMU.assign(p, gLRMU(p).conjugate());
        gRLTAU.assign(p, gLRTAU(p).conjugate());
        gLLNU.assign(p, -0.5*MZ/ctw*sigma3(p));
    }

    //Note that the dependence on the Rl elements here is different than the one in PhysRevD.73.055003
    //due to the different choice of flavour basis.
    //The replacement rules are:
    //AH basis (PhysRevD.73.055003) -> SLHA basis
    //  (from 1 to 6)               -> (from 0 to 5)
    //  Rl(x,1)                     -> Rl(0,x)
    //  Rl(x,2)                     -> Rl(3,x)
    //  Rl(x,3)                     -> Rl(1,x)
    //  Rl(x,4)                     -> Rl(4,x)
    //  Rl(x,5)                     -> Rl(2,x)
    //  Rl(x,6)                     -> Rl(5,x)
    gslpp::matrix<gslpp::complex> Gl0(6, 6, 0.), Gl1(6, 6, 0.), Gl2(6, 6, 0.);
    for (int x=0;x<6;x++) {
        for (int y=0;y<6;y++) {
            Gl0.assign(x, y, -g2*( gLLE(0)*myRl(x,0).conjugate()*myRl(y,0) +gRRE(0)*myRl(x,3).conjugate()*myRl(y,3) +gLRE(0)*myRl(x,0).conjugate()*myRl(y,3) +gRLE(0)*myRl(x,3).conjugate()*myRl(y,0)
                                  +gLLMU(0)*myRl(x,1).conjugate()*myRl(y,1) +gRRMU(0)*myRl(x,4).conjugate()*myRl(y,4) +gLRMU(0)*myRl(x,1).conjugate()*myRl(y,4) +gRLMU(0)*myRl(x,4).conjugate()*myRl(y,1)
                                  +gLLTAU(0)*myRl(x,2).conjugate()*myRl(y,2) +gRRTAU(0)*myRl(x,5).conjugate()*myRl(y,5) +gLRTAU(0)*myRl(x,2).conjugate()*myRl(y,5) +gRLTAU(0)*myRl(x,5).conjugate()*myRl(y,2)));
            Gl1.assign(x, y, -g2*( gLLE(1)*myRl(x,0).conjugate()*myRl(y,0) +gRRE(1)*myRl(x,3).conjugate()*myRl(y,3) +gLRE(1)*myRl(x,0).conjugate()*myRl(y,3) +gRLE(1)*myRl(x,3).conjugate()*myRl(y,0)
                                  +gLLMU(1)*myRl(x,1).conjugate()*myRl(y,1) +gRRMU(1)*myRl(x,4).conjugate()*myRl(y,4) +gLRMU(1)*myRl(x,1).conjugate()*myRl(y,4) +gRLMU(1)*myRl(x,4).conjugate()*myRl(y,1)
                                  +gLLTAU(1)*myRl(x,2).conjugate()*myRl(y,2) +gRRTAU(1)*myRl(x,5).conjugate()*myRl(y,5) +gLRTAU(1)*myRl(x,2).conjugate()*myRl(y,5) +gRLTAU(1)*myRl(x,5).conjugate()*myRl(y,2)));
            Gl2.assign(x, y, -g2*( gLLE(2)*myRl(x,0).conjugate()*myRl(y,0) +gRRE(2)*myRl(x,3).conjugate()*myRl(y,3) +gLRE(2)*myRl(x,0).conjugate()*myRl(y,3) +gRLE(2)*myRl(x,3).conjugate()*myRl(y,0)
                                  +gLLMU(2)*myRl(x,1).conjugate()*myRl(y,1) +gRRMU(2)*myRl(x,4).conjugate()*myRl(y,4) +gLRMU(2)*myRl(x,1).conjugate()*myRl(y,4) +gRLMU(2)*myRl(x,4).conjugate()*myRl(y,1)
                                  +gLLTAU(2)*myRl(x,2).conjugate()*myRl(y,2) +gRRTAU(2)*myRl(x,5).conjugate()*myRl(y,5) +gLRTAU(2)*myRl(x,2).conjugate()*myRl(y,5) +gRLTAU(2)*myRl(x,5).conjugate()*myRl(y,2)));
        }
    }

    gslpp::matrix<gslpp::complex> Gnu0(3, 3, 0.), Gnu1(3, 3, 0.), Gnu2(3, 3, 0.);
    for (int x=0;x<3;x++) {
        for (int y=0;y<3;y++) {
            Gnu0.assign(x, y, -g2*( gLLNU(0)*myRn(0,x).conjugate()*myRn(0,y) +gLLNU(0)*myRn(1,x).conjugate()*myRn(1,y) +gLLNU(0)*myRn(2,x).conjugate()*myRn(2,y)) );
            Gnu1.assign(x, y, -g2*( gLLNU(1)*myRn(0,x).conjugate()*myRn(0,y) +gLLNU(1)*myRn(1,x).conjugate()*myRn(1,y) +gLLNU(1)*myRn(2,x).conjugate()*myRn(2,y)) );
            Gnu2.assign(x, y, -g2*( gLLNU(2)*myRn(0,x).conjugate()*myRn(0,y) +gLLNU(2)*myRn(1,x).conjugate()*myRn(1,y) +gLLNU(2)*myRn(2,x).conjugate()*myRn(2,y)) );
        }
    }

    gslpp::vector<gslpp::complex> SRE(3, 0.), SLE(3, 0.), SRMU(3, 0.), SLMU(3, 0.), SRTAU(3, 0.), SLTAU(3, 0.);
    for (int p=0;p<3;p++) {
        SRE.assign(p, g2*mE/(2.0*MW*cosb) * sigma1(p));
        SLE.assign(p, g2*mE/(2.0*MW*cosb) * sigma1(p).conjugate());
        SRMU.assign(p, g2*mMU/(2.0*MW*cosb) * sigma1(p));
        SLMU.assign(p, g2*mMU/(2.0*MW*cosb) * sigma1(p).conjugate());
        SRTAU.assign(p, g2*mTAU/(2.0*MW*cosb) * sigma1(p));
        SLTAU.assign(p, g2*mTAU/(2.0*MW*cosb) * sigma1(p).conjugate());
    }

    if (li_to_lj == 1) // mu -> 3e
    {
        //     Neutralino contributions
        gslpp::complex HpengMuEEENR0 = 0.0;
        gslpp::complex HpengMuEEENL0 = 0.0;
        gslpp::complex HpengMuEEENR1 = 0.0;
        gslpp::complex HpengMuEEENL1 = 0.0;
        gslpp::complex HpengMuEEENR2 = 0.0;
        gslpp::complex HpengMuEEENL2 = 0.0;
        for (int x=0;x<6;x++) {
            for (int a=0;a<4;a++) {
                for (int b=0;b<4;b++) {
                    //     h R contribution
                    HpengMuEEENR0 = HpengMuEEENR0 - 2.0*piconst*(NRlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NLlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR0(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL0(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NLlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR0(a,b)*NLlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                    //     h L contribution
    if(x==0&&a==2) std::cout<<"HpengMuEEENR0(3,1)_1b = "<< - 2.0*piconst*(NRlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NLlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR0(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL0(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NLlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR0(a,b)*NLlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))<<std::endl;
                    HpengMuEEENL0 = HpengMuEEENL0 - 2.0*piconst*(NLlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NRlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL0(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR0(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NRlE(a,x)*DR0(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL0(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL0(a,b)*NRlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                    //     H R contribution
                    HpengMuEEENR1 = HpengMuEEENR1 - 2.0*piconst*(NRlE(a,x)*DL1(a,b)*NLlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NLlE(a,x)*DR1(a,b)*NRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR1(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL1(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NLlE(a,x)*DL1(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR1(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                          -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR1(a,b)*NLlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                    //     H L contribution
                    HpengMuEEENL1 = HpengMuEEENL1 - 2.0*piconst*(NLlE(a,x)*DR1(a,b)*NRlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NRlE(a,x)*DL1(a,b)*NLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL1(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR1(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NRlE(a,x)*DR1(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL1(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                          -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL1(a,b)*NRlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                    //     A R contribution
                    HpengMuEEENR2 = HpengMuEEENR2 - 2.0*piconst*(NRlE(a,x)*DL2(a,b)*NLlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NLlE(a,x)*DR2(a,b)*NRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR2(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL2(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NLlE(a,x)*DL2(a,b)*NLlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR2(a,b)*NRlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                          -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DR2(a,b)*NLlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                    //     A L contribution
                    HpengMuEEENL2 = HpengMuEEENL2 - 2.0*piconst*(NLlE(a,x)*DR2(a,b)*NRlMU(b,x).conjugate()*(PV.B0(1.,0.,MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            -mym_se_sq(x)*PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mMU*mMU*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                            +mE*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))))
                                                                 +NRlE(a,x)*DL2(a,b)*NLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                    -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NRlE(a,x)*DL2(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(b)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DR2(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(b)*PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                 +NRlE(a,x)*DR2(a,b)*NRlMU(b,x).conjugate()*mE*MNeig(a)*(PV.C11(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                         -PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL2(a,b)*NLlMU(b,x).conjugate()*mMU*MNeig(a)*(PV.C12(mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))
                                                                                                                          -PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b)))
                                                                 +NLlE(a,x)*DL2(a,b)*NRlMU(b,x).conjugate()*MNeig(a)*MNeig(b)*(-PV.C0(0.,mym_se_sq(x),MNeig(a)*MNeig(a),MNeig(b)*MNeig(b))));
                }
                for (int y=0;y<6;y++) {
                    //     h R contribution
                    HpengMuEEENR0 = HpengMuEEENR0 - 2.0*piconst*Gl0(x,y)*(-NLlE(a,x)*NLlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NRlE(a,x)*NRlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NRlE(a,x)*NLlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
                    //     h L contribution
                    HpengMuEEENL0 = HpengMuEEENL0 - 2.0*piconst*Gl0(x,y)*(-NRlE(a,x)*NRlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NLlE(a,x)*NLlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NLlE(a,x)*NRlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
    if(x==0&&a==2) std::cout<<"HpengMuEEENR0(3,1)_2y = "<< - 2.0*piconst*Gl0(x,y)*(-NLlE(a,x)*NLlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NRlE(a,x)*NRlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NRlE(a,x)*NLlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))<<std::endl;
                    //     H R contribution
                    HpengMuEEENR1 = HpengMuEEENR1 - 2.0*piconst*Gl1(x,y)*(-NLlE(a,x)*NLlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NRlE(a,x)*NRlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NRlE(a,x)*NLlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
                    //     H L contribution
                    HpengMuEEENL1 = HpengMuEEENL1 - 2.0*piconst*Gl1(x,y)*(-NRlE(a,x)*NRlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NLlE(a,x)*NLlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NLlE(a,x)*NRlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
                    //     A R contribution
                    HpengMuEEENR2 = HpengMuEEENR2 - 2.0*piconst*Gl2(x,y)*(-NLlE(a,x)*NLlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NRlE(a,x)*NRlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NRlE(a,x)*NLlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
                    //     A L contribution
                    HpengMuEEENL2 = HpengMuEEENL2 - 2.0*piconst*Gl2(x,y)*(-NRlE(a,x)*NRlMU(a,y).conjugate()*mE*(PV.C11(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                                                                -PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)))
                                                                          -NLlE(a,x)*NLlMU(a,y).conjugate()*mMU*PV.C12(MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y))
                                                                          -NLlE(a,x)*NRlMU(a,y).conjugate()*MNeig(a)*PV.C0(0.,MNeig(a)*MNeig(a),mym_se_sq(x),mym_se_sq(y)));
                }
                //     h R contribution
                HpengMuEEENR0 = HpengMuEEENR0 - 2.0*piconst*(SRMU(0)/(mE*mE-mMU*mMU)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SRE(0)/(mMU*mMU-mE*mE)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));
                //     h L contribution
                HpengMuEEENL0 = HpengMuEEENL0 - 2.0*piconst*(SLMU(0)/(mE*mE-mMU*mMU)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SLE(0)/(mMU*mMU-mE*mE)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));
    if(x==0&&a==2) std::cout<<"HpengMuEEENR0(3,1)_3 = "<< - 2.0*piconst*(SRMU(0)/(mE*mE-mMU*mMU)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SRE(0)/(mMU*mMU-mE*mE)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))))<<std::endl;
    if(x==0&&a==2) std::cout<<"SRMU(0) = "<< SRMU(0)<<std::endl;
    if(x==0&&a==2) std::cout<<"NRlE(3,1) = "<< NRlE(a,x)<<std::endl;
    if(x==0&&a==2) std::cout<<"NRlMU(3,1) = "<< NRlMU(a,x)<<std::endl;
    if(x==0&&a==2) std::cout<<"NLlE(3,1) = "<< NLlE(a,x)<<std::endl;
    if(x==0&&a==2) std::cout<<"NLlMU(3,1) = "<< NLlMU(a,x)<<std::endl;
    if(x==0&&a==2) std::cout<<"B0(3,1) = "<< PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))<<std::endl;
    if(x==0&&a==2) std::cout<<"B1(3,1) = "<< PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))<<std::endl;
                
    if(x==0&&a==2) std::cout<<"HpengMuEEENR0(3,1)_3_1 = "<< - 2.0*piconst*(SRMU(0)/(mE*mE-mMU*mMU)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))))<<std::endl;
    if(x==0&&a==2) std::cout<<"HpengMuEEENR0(3,1)_3_2 = "<< - 2.0*piconst*(SRE(0)/(mMU*mMU-mE*mE)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))))<<std::endl;
                //     H R contribution
                HpengMuEEENR1 = HpengMuEEENR1 - 2.0*piconst*(SRMU(1)/(mE*mE-mMU*mMU)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SRE(1)/(mMU*mMU-mE*mE)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));
                //     H L contribution
                HpengMuEEENL1 = HpengMuEEENL1 - 2.0*piconst*(SLMU(1)/(mE*mE-mMU*mMU)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SLE(1)/(mMU*mMU-mE*mE)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));
                //     A R contribution
                HpengMuEEENR2 = HpengMuEEENR2 - 2.0*piconst*(SRMU(2)/(mE*mE-mMU*mMU)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SRE(2)/(mMU*mMU-mE*mE)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));
                //     A L contribution
                HpengMuEEENL2 = HpengMuEEENL2 - 2.0*piconst*(SLMU(2)/(mE*mE-mMU*mMU)*(-NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NRlE(a,x)*NRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x)))
                                                             +SLE(2)/(mMU*mMU-mE*mE)*(-NRlE(a,x)*NRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NRlE(a,x)*NLlMU(a,x).conjugate()*mMU*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      -NLlE(a,x)*NLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))
                                                                                      +NLlE(a,x)*NRlMU(a,x).conjugate()*mE*MNeig(a)*PV.B0(1.,0.,MNeig(a)*MNeig(a),mym_se_sq(x))));

                std::cout<<"HpengMuEEENR0("<<a<<","<<x<<") = "<< HpengMuEEENR0-HpengMuEEENL0temp<<std::endl;
                HpengMuEEENL0temp=HpengMuEEENR0;

            }
        }
    std::cout<<"HpengMuEEENR1 = "<< HpengMuEEENR1<<std::endl;
    std::cout<<"HpengMuEEENL1 = "<< HpengMuEEENL1<<std::endl;
    std::cout<<"HpengMuEEENR2 = "<< HpengMuEEENR2<<std::endl;
    std::cout<<"HpengMuEEENL2 = "<< HpengMuEEENL2<<std::endl;

        //     summing up the h (0), H (1) and A (2) parts
        gslpp::complex B2HiggsnR = (-0.5*HpengMuEEENR0*SLE(0)/(mh*mh)-0.5*HpengMuEEENR1*SLE(1)/(mH*mH)-0.5*HpengMuEEENR2*SLE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B2HiggsnL = (-0.5*HpengMuEEENL0*SRE(0)/(mh*mh)-0.5*HpengMuEEENL1*SRE(1)/(mH*mH)-0.5*HpengMuEEENL2*SRE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B3HiggsnR = (HpengMuEEENR0*SRE(0)/(mh*mh)+HpengMuEEENR1*SRE(1)/(mH*mH)+HpengMuEEENR2*SRE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B3HiggsnL = (HpengMuEEENL0*SLE(0)/(mh*mh)+HpengMuEEENL1*SLE(1)/(mH*mH)+HpengMuEEENL2*SLE(2)/(mA*mA))/(4.0*pi*alph);

        //     Chargino contributions
        gslpp::complex HpengMuEEECR0 = 0.0;
        gslpp::complex HpengMuEEECL0 = 0.0;
        gslpp::complex HpengMuEEECR1 = 0.0;
        gslpp::complex HpengMuEEECL1 = 0.0;
        gslpp::complex HpengMuEEECR2 = 0.0;
        gslpp::complex HpengMuEEECL2 = 0.0;
        for (int x=0;x<3;x++) {
            for (int a=0;a<2;a++) {
                for (int b=0;b<2;b++) {
                    //     h R contribution
                    HpengMuEEECR0 = HpengMuEEECR0 - 2.0*piconst*(CRlE(a,x)*WL0(a,b)*CLlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CLlE(a,x)*WR0(a,b)*CRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR0(a,b)*CLlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL0(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CLlE(a,x)*WL0(a,b)*CLlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR0(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR0(a,b)*CLlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                    //     h L contribution
                    HpengMuEEECL0 = HpengMuEEECL0 - 2.0*piconst*(CLlE(a,x)*WR0(a,b)*CRlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CRlE(a,x)*WL0(a,b)*CLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL0(a,b)*CRlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR0(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CRlE(a,x)*WR0(a,b)*CRlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL0(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL0(a,b)*CRlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                    //     H R contribution
                    HpengMuEEECR1 = HpengMuEEECR1 - 2.0*piconst*(CRlE(a,x)*WL1(a,b)*CLlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CLlE(a,x)*WR1(a,b)*CRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR1(a,b)*CLlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL1(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CLlE(a,x)*WL1(a,b)*CLlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR1(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR1(a,b)*CLlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                    //     H L contribution
                    HpengMuEEECL1 = HpengMuEEECL1 - 2.0*piconst*(CLlE(a,x)*WR1(a,b)*CRlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CRlE(a,x)*WL1(a,b)*CLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL1(a,b)*CRlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR1(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CRlE(a,x)*WR1(a,b)*CRlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL1(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL1(a,b)*CRlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                    //     A R contribution
                    HpengMuEEECR2 = HpengMuEEECR2 - 2.0*piconst*(CRlE(a,x)*WL2(a,b)*CLlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CLlE(a,x)*WR2(a,b)*CRlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR2(a,b)*CLlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL2(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CLlE(a,x)*WL2(a,b)*CLlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR2(a,b)*CRlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WR2(a,b)*CLlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                    //     A L contribution
                    HpengMuEEECL2 = HpengMuEEECL2 - 2.0*piconst*(CLlE(a,x)*WR2(a,b)*CRlMU(b,x).conjugate()*(PV.B0(1.,0.,MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            -mym_sn_sq(x)*PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mMU*mMU*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                            +mE*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))))
                                                                 +CRlE(a,x)*WL2(a,b)*CLlMU(b,x).conjugate()*mMU*mE*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                    -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CRlE(a,x)*WL2(a,b)*CRlMU(b,x).conjugate()*mE*MChi(b)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WR2(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(b)*PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                 +CRlE(a,x)*WR2(a,b)*CRlMU(b,x).conjugate()*mE*MChi(a)*(PV.C11(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                        -PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL2(a,b)*CLlMU(b,x).conjugate()*mMU*MChi(a)*(PV.C12(mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))
                                                                                                                         -PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b)))
                                                                 +CLlE(a,x)*WL2(a,b)*CRlMU(b,x).conjugate()*MChi(a)*MChi(b)*(-PV.C0(0.,mym_sn_sq(x),MChi(a)*MChi(a),MChi(b)*MChi(b))));
                }
                for (int y=0;y<3;y++) {
                    //     h R contribution
                    HpengMuEEECR0 = HpengMuEEECR0 - 2.0*piconst*Gnu0(x,y)*(-CLlE(a,x)*CLlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CRlE(a,x)*CRlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CRlE(a,x)*CLlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                    //     h L contribution
                    HpengMuEEECL0 = HpengMuEEECL0 - 2.0*piconst*Gnu0(x,y)*(-CRlE(a,x)*CRlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CLlE(a,x)*CLlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CLlE(a,x)*CRlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                    //     H R contribution
                    HpengMuEEECR1 = HpengMuEEECR1 - 2.0*piconst*Gnu1(x,y)*(-CLlE(a,x)*CLlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CRlE(a,x)*CRlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CRlE(a,x)*CLlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                    //     H L contribution
                    HpengMuEEECL1 = HpengMuEEECL1 - 2.0*piconst*Gnu1(x,y)*(-CRlE(a,x)*CRlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CLlE(a,x)*CLlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CLlE(a,x)*CRlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                    //     A R contribution
                    HpengMuEEECR2 = HpengMuEEECR2 - 2.0*piconst*Gnu2(x,y)*(-CLlE(a,x)*CLlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CRlE(a,x)*CRlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CRlE(a,x)*CLlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                    //     A L contribution
                    HpengMuEEECL2 = HpengMuEEECL2 - 2.0*piconst*Gnu2(x,y)*(-CRlE(a,x)*CRlMU(a,y).conjugate()*mE*(PV.C11(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                                                                 -PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)))
                                                                           -CLlE(a,x)*CLlMU(a,y).conjugate()*mMU*PV.C12(MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y))
                                                                           -CLlE(a,x)*CRlMU(a,y).conjugate()*MChi(a)*PV.C0(0.,MChi(a)*MChi(a),mym_sn_sq(x),mym_sn_sq(y)));
                }
                //     h R contribution
                HpengMuEEECR0 = HpengMuEEECR0 - 2.0*piconst*(SRMU(0)/(mE*mE-mMU*mMU)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SRE(0)/(mMU*mMU-mE*mE)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
                //     h L contribution
                HpengMuEEECL0 = HpengMuEEECL0 - 2.0*piconst*(SLMU(0)/(mE*mE-mMU*mMU)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SLE(0)/(mMU*mMU-mE*mE)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
                //     H R contribution
                HpengMuEEECR1 = HpengMuEEECR1 - 2.0*piconst*(SRMU(1)/(mE*mE-mMU*mMU)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SRE(1)/(mMU*mMU-mE*mE)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
                //     H L contribution
                HpengMuEEECL1 = HpengMuEEECL1 - 2.0*piconst*(SLMU(1)/(mE*mE-mMU*mMU)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SLE(1)/(mMU*mMU-mE*mE)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
                //     A R contribution
                HpengMuEEECR2 = HpengMuEEECR2 - 2.0*piconst*(SRMU(2)/(mE*mE-mMU*mMU)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SRE(2)/(mMU*mMU-mE*mE)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
                //     A L contribution
                HpengMuEEECL2 = HpengMuEEECL2 - 2.0*piconst*(SLMU(2)/(mE*mE-mMU*mMU)*(-CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mE*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CRlE(a,x)*CRlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x)))
                                                             +SLE(2)/(mMU*mMU-mE*mE)*(-CRlE(a,x)*CRlMU(a,x).conjugate()*mMU*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CRlE(a,x)*CLlMU(a,x).conjugate()*mMU*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      -CLlE(a,x)*CLlMU(a,x).conjugate()*mE*mMU*PV.B1(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))
                                                                                      +CLlE(a,x)*CRlMU(a,x).conjugate()*mE*MChi(a)*PV.B0(1.,0.,MChi(a)*MChi(a),mym_sn_sq(x))));
           }
        }
    std::cout<<"HpengMuEEECR0 = "<< HpengMuEEECR0<<std::endl;
    std::cout<<"HpengMuEEECL0 = "<< HpengMuEEECL0<<std::endl;
    std::cout<<"HpengMuEEECR1 = "<< HpengMuEEECR1<<std::endl;
    std::cout<<"HpengMuEEECL1 = "<< HpengMuEEECL1<<std::endl;
    std::cout<<"HpengMuEEECR2 = "<< HpengMuEEECR2<<std::endl;
    std::cout<<"HpengMuEEECL2 = "<< HpengMuEEECL2<<std::endl;

        //     summing up the h (0), H (1) and A (2) parts
        gslpp::complex B2HiggscR = (-0.5*HpengMuEEECR0*SLE(0)/(mh*mh)-0.5*HpengMuEEECR1*SLE(1)/(mH*mH)-0.5*HpengMuEEECR2*SLE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B2HiggscL = (-0.5*HpengMuEEECL0*SRE(0)/(mh*mh)-0.5*HpengMuEEECL1*SRE(1)/(mH*mH)-0.5*HpengMuEEECL2*SRE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B3HiggscR = (HpengMuEEECR0*SRE(0)/(mh*mh)+HpengMuEEECR1*SRE(1)/(mH*mH)+HpengMuEEECR2*SRE(2)/(mA*mA))/(4.0*pi*alph);
        gslpp::complex B3HiggscL = (HpengMuEEECL0*SLE(0)/(mh*mh)+HpengMuEEECL1*SLE(1)/(mH*mH)+HpengMuEEECL2*SLE(2)/(mA*mA))/(4.0*pi*alph);

        //     write B2H and B3H into a vector for mu->3e
        BHFunctions.assign(0, B2HiggsnR+B2HiggscR );
        BHFunctions.assign(1, B2HiggsnL+B2HiggscL );
        BHFunctions.assign(2, B3HiggsnR+B3HiggscR );
        BHFunctions.assign(3, B3HiggsnL+B3HiggscL );
    }
//    if (li_to_lj == 2) // tau -> 3mu
//    {
//    }
//    if (li_to_lj == 3) // tau -> 3e
//    {
//    }
//    if (li_to_lj == 4) // tau -> muee
//    {
//    }

    return(BHFunctions);
}

gslpp::vector<gslpp::complex> SUSYMatching::DFunctions() {
        //box contributions to mu->e conversion

    NeutralinoRemixing();

    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mUP = mySUSY.getQuarks(QCD::UP).getMass();
    double mDOWN = mySUSY.getQuarks(QCD::DOWN).getMass();
    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);
    double alph = mySUSY.getAle();

    gslpp::vector<gslpp::complex> DFunctions(4, 0.);
    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlE.assign(a, x, - (g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 0) + (mE/cdenn)*ON(a, 2)*myRl(x, 3)));
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            NRqUP.assign(a, x, -(g2t)*((ON(a, 1)+ON(a, 0)*ttw/3.0)*myRu(x, 0)+mUP/(MW*sinb)*ON(a, 3)*myRu(x, 3)));
            NRqDOWN.assign(a, x, -(g2t)*((-ON(a, 1)+ON(a, 0)*ttw/3.0)*myRd(x, 0)+mDOWN/(MW*cosb)*ON(a, 2)*myRd(x, 3)));
            //     RL + RR TYPE MI
            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2)*myRl(x, 0) + 2.0*ON(a, 0)*ttw*myRl(x, 3)));
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
            NLqUP.assign(a, x, -(g2t)*(mUP/(MW*sinb)*ON(a, 3)*myRu(x, 0)-4.0*ON(a, 0)*myRu(x, 3)*ttw/3.0));
            NLqDOWN.assign(a, x, -(g2t)*(mDOWN/(MW*cosb)*ON(a, 2)*myRd(x, 0)+2.0*ON(a, 0)*myRd(x, 3)*ttw/3.0));
        }
    }
    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlE.assign(a, x, - (g2*myV(a, 0)*myRn(x, 0)));
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            CRqUP.assign(a, x, g2*(-myU(a, 0)*myRd(x, 0) + mDOWN/cdenc*myU(a, 1)*myRd(x, 3)));
            CRqDOWN.assign(a, x, g2*(-myV(a, 0)*myRu(x, 0) + mUP/(sqrt(2.0)*MW*sinb)*myV(a, 1)*myRu(x, 3)));
            //     LR-TYPE
            CLlE.assign(a, x, g2*mE/cdenc*myU(a, 1).conjugate()*myRn(x, 0));
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
            CLqUP.assign(a, x, g2*mUP/(sqrt(2.0)*MW*sinb)*myV(a, 1)*myRd(x, 0));
            CLqDOWN.assign(a, x, g2*mDOWN/cdenc*myU(a, 1)*myRu(x, 0));
        }
    }

    //     Neutralino contributions
    gslpp::complex DunR = 0.0;
    gslpp::complex DunL = 0.0;
    gslpp::complex DdnR = 0.0;
    gslpp::complex DdnL = 0.0;
    for (int a=0;a<4;a++) {
        for (int b=0;b<4;b++) {
            for (int x=0;x<6;x++) {
                for (int y=0;y<6;y++) {
                    DunR = DunR + (0.125*(NLlMU(a,x).conjugate()*NLlE(b,x)*NLqUP(a,y)*NLqUP(b,y).conjugate()
                                          -NLlMU(a,x).conjugate()*NLlE(b,x)*NRqUP(a,y).conjugate()*NRqUP(b,y))
                                         *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(4.0*pi*pi)
                                   -0.25*MNeig(a)*MNeig(b)*(NLlMU(a,x).conjugate()*NLlE(b,x)*NRqUP(a,y)*NRqUP(b,y).conjugate()
                                                            -NLlMU(a,x).conjugate()*NLlE(b,x)*NLqUP(a,y).conjugate()*NLqUP(b,y))
                                                           *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
                    DunL = DunL + (0.125*(NRlMU(a,x).conjugate()*NRlE(b,x)*NRqUP(a,y)*NRqUP(b,y).conjugate()
                                          -NRlMU(a,x).conjugate()*NRlE(b,x)*NLqUP(a,y).conjugate()*NLqUP(b,y))
                                         *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(4.0*pi*pi)
                                   -0.25*MNeig(a)*MNeig(b)*(NRlMU(a,x).conjugate()*NRlE(b,x)*NLqUP(a,y)*NLqUP(b,y).conjugate()
                                                            -NRlMU(a,x).conjugate()*NRlE(b,x)*NRqUP(a,y).conjugate()*NRqUP(b,y))
                                                           *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
                    DdnR = DdnR + (0.125*(NLlMU(a,x).conjugate()*NLlE(b,x)*NLqDOWN(a,y)*NLqDOWN(b,y).conjugate()
                                          -NLlMU(a,x).conjugate()*NLlE(b,x)*NRqDOWN(a,y).conjugate()*NRqDOWN(b,y))
                                         *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(4.0*pi*pi)
                                   -0.25*MNeig(a)*MNeig(b)*(NLlMU(a,x).conjugate()*NLlE(b,x)*NRqDOWN(a,y)*NRqDOWN(b,y).conjugate()
                                                            -NLlMU(a,x).conjugate()*NLlE(b,x)*NLqDOWN(a,y).conjugate()*NLqDOWN(b,y))
                                                           *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
                    DdnL = DdnL + (0.125*(NRlMU(a,x).conjugate()*NRlE(b,x)*NRqDOWN(a,y)*NRqDOWN(b,y).conjugate()
                                          -NRlMU(a,x).conjugate()*NRlE(b,x)*NLqDOWN(a,y).conjugate()*NLqDOWN(b,y))
                                         *PV.D00(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(4.0*pi*pi)
                                   -0.25*MNeig(a)*MNeig(b)*(NRlMU(a,x).conjugate()*NRlE(b,x)*NLqDOWN(a,y)*NLqDOWN(b,y).conjugate()
                                                            -NRlMU(a,x).conjugate()*NRlE(b,x)*NRqDOWN(a,y).conjugate()*NRqDOWN(b,y))
                                                           *PV.D0(0., 0., MNeig(a)*MNeig(a), MNeig(b)*MNeig(b), mym_se_sq(x), mym_se_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
                }
            }
        }
    }

    //     Chargino contributions
    gslpp::complex DucR = 0.0;
    gslpp::complex DucL = 0.0;
    gslpp::complex DdcR = 0.0;
    gslpp::complex DdcL = 0.0;
    for (int a=0;a<2;a++) {
        for (int b=0;b<2;b++) {
            for (int x=0;x<3;x++) {
                for (int y=0;y<3;y++) {
                    DucL = DucL + (-0.125*CRlMU(a,x).conjugate()*CRlE(b,x)*CLqUP(a,y).conjugate()*CLqUP(b,y)
                                   *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(y))/(4.0*pi*pi)
                                   +0.25*MChi(a)*MChi(b)*CRlMU(a,x).conjugate()*CRlE(b,x)*CRqUP(a,y).conjugate()*CRqUP(b,y)
                                   *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
                    DdcL = DdcL + (0.125*CRlMU(a,x).conjugate()*CRlE(b,x)*CRqDOWN(a,y)*CRqDOWN(b,y).conjugate()
                                   *PV.D00(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(y))/(4.0*pi*pi)
                                   -0.25*MChi(a)*MChi(b)*CRlMU(a,x).conjugate()*CRlE(b,x)*CLqDOWN(a,y)*CLqDOWN(b,y).conjugate()
                                   *PV.D0(0., 0., MChi(a)*MChi(a), MChi(b)*MChi(b), mym_sn_sq(x), mym_sn_sq(y))/(16.0*pi*pi))/(4.0*pi*alph);
// And the R contributions???
                }
            }
        }
    }

    //     write DuR, DuL, DdR and DdL into a vector for mu->e conversion
    DFunctions.assign(0, 0.0 );  //DuR
    DFunctions.assign(1, 0.0 );  //DuL
    DFunctions.assign(2, 0.0 );  //DdR
    DFunctions.assign(3, 0.0 );  //DdL

    return(DFunctions);
}

gslpp::vector<gslpp::complex> SUSYMatching::FFunctions(int li_to_lj) {
        //Z penguin contributions

    NeutralinoRemixing();

    double MZ = mySUSY.getMz();
    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double piconst = 1.0/(32.0 * pi * pi);
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();

    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);

    gslpp::vector<gslpp::complex> FFunctions(4, 0.);
//      Note: Here contributions proportional to the Yukawa couplings are set to 
//      zero. In particular this leads to Right Chargino Contributions to be zero. 

    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlE.assign(a, x, - (g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 0) + (mE/cdenn)*ON(a, 2)*myRl(x, 3)));
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            NRlTAU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 2) + (mTAU/cdenn)*ON(a, 2)*myRl(x, 5)));
            //     RL + RR TYPE MI
            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2)*myRl(x, 0) + 2.0*ON(a, 0)*ttw*myRl(x, 3)));
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2)*myRl(x, 2) + 2.0*ON(a, 0)*ttw*myRl(x, 5)));
//    Commented expressions might be useful for complex neutralino mixing matrices
//            NLlE.assign(a, x, -(g2t)*((mE/cdenn)*ON(a, 2).conjugate()*myRl(x, 0) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 3)));
//            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2).conjugate()*myRl(x, 1) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 4)));
//            NLlTAU.assign(a, x, -(g2t)*((mTAU/cdenn)*ON(a, 2).conjugate()*myRl(x, 2) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 5)));
        }
    }

    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlE.assign(a, x, - (g2*myV(a, 0)*myRn(x, 0)));
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            CRlTAU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 2)));
            //     LR-TYPE
            CLlE.assign(a, x, g2*mE/cdenc*myU(a, 1).conjugate()*myRn(x, 0));
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
            CLlTAU.assign(a, x, g2*mTAU/cdenc*myU(a, 1).conjugate()*myRn(x, 2));
        }
    }

//      Neutralino contributions

    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            Lepty.assign(a, x, ( MNeig(a) * MNeig(a) ) / mym_se_sq(x) );
        }
    }

    for (int x=0;x<6;x++) {
        for (int a=0;a<4;a++) {
            for (int b=0;b<4;b++) {
                if (a != b && std::abs(Lepty(a,x)-Lepty(b,x)) > 0.01 && std::abs(1.0-Lepty(a,x)) > 0.01 && std::abs(1.0-Lepty(b,x)) > 0.01) {
                    Leptfzn[x][a][b] = log(Lepty(a,x)) + 1.0/(Lepty(a,x)-Lepty(b,x))*
                                     ( (Lepty(a,x)*Lepty(a,x)*log(Lepty(a,x)))/(1.0-Lepty(a,x)) - (Lepty(b,x)*Lepty(b,x)*log(Lepty(b,x)))/(1.0-Lepty(b,x)) );
                    Leptgzn[x][a][b] = (MNeig(a)*MNeig(b)/mym_se_sq(x))/(Lepty(a,x)-Lepty(b,x))*
                                     ( (Lepty(a,x)*log(Lepty(a,x)))/(1.0-Lepty(a,x)) - (Lepty(b,x)*log(Lepty(b,x)))/(1.0-Lepty(b,x)) );
                }
                else if (a != b && std::abs(Lepty(a,x)-Lepty(b,x)) > 0.01 && std::abs(1.0-Lepty(a,x)) > 0.01 && std::abs(1.0-Lepty(b,x)) <= 0.01) {
                    Leptfzn[x][a][b] = (-2.0*log(Lepty(a,x))*Lepty(a,x)+log(Lepty(a,x))+Lepty(a,x)-1.0)/(Lepty(a,x)-1.0)/(Lepty(a,x)-1.0);
                    Leptgzn[x][a][b] = -(MNeig(a)*MNeig(b)/mym_se_sq(x))*(log(Lepty(a,x))*Lepty(a,x)-Lepty(a,x)+1.0)/(Lepty(a,x)-1.0)/(Lepty(a,x)-1.0);
                }
                else if (a != b && std::abs(Lepty(a,x)-Lepty(b,x)) > 0.01 && std::abs(1.0-Lepty(b,x)) > 0.01 && std::abs(1.0-Lepty(a,x)) <= 0.01) {
                    Leptfzn[x][a][b] = (-log(Lepty(b,x))*Lepty(b,x)*Lepty(b,x)+Lepty(b,x)-1.0)/(Lepty(b,x)-1.0)/(Lepty(b,x)-1.0); 
                    Leptgzn[x][a][b] = -(MNeig(a)*MNeig(b)/mym_se_sq(x))*(log(Lepty(b,x))*Lepty(b,x)-Lepty(b,x)+1.0)/(Lepty(b,x)-1.0)/(Lepty(b,x)-1.0);
                }
                else if ((a == b || std::abs(Lepty(a,x)-Lepty(b,x)) <= 0.01) && std::abs(1.0-Lepty(a,x)) > 0.01 && std::abs(1.0-Lepty(b,x)) > 0.01) {
                    Leptfzn[x][a][b] = (-Lepty(a,x)*Lepty(a,x)+Lepty(a,x)+log(Lepty(a,x)))/(Lepty(a,x)-1.0)/(Lepty(a,x)-1.0);
                    Leptgzn[x][a][b] = (MNeig(b)*MNeig(b)/mym_se_sq(x))*(-Lepty(a,x)+log(Lepty(a,x))+1.0)/(Lepty(a,x)-1.0)/(Lepty(a,x)-1.0);
                }
                else {
                    Leptfzn[x][a][b] = -3.0/2.0-(Lepty(a,x)-1.0)/3.0;
                    Leptgzn[x][a][b] = -MNeig(a)*MNeig(b)/(2.0*mym_se_sq(x));
                }
            }
        }
    }

//      Chargino contributions
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            Leptz.assign(a, x, MChi(a)*MChi(a)/mym_sn_sq(x) );
        }
    }

    for (int x=0;x<3;x++) {
        for (int a=0;a<2;a++) {
            for (int b=0;b<2;b++) {
                if (a != b && std::abs(Leptz(a,x)-Leptz(b,x)) > 0.01 && std::abs(1.0-Leptz(a,x)) > 0.01 && std::abs(1.0-Leptz(b,x)) > 0.01) {
                    Leptfzc[x][a][b] = log(Leptz(a,x)) + 1.0/(Leptz(a,x)-Leptz(b,x))*
                                     ( (Leptz(a,x)*Leptz(a,x)*log(Leptz(a,x)))/(1.0-Leptz(a,x)) - (Leptz(b,x)*Leptz(b,x)*log(Leptz(b,x)))/(1.0-Leptz(b,x)) );
                    Leptgzc[x][a][b] = (MChi(a)*MChi(b)/mym_sn_sq(x))/(Leptz(a,x)-Leptz(b,x))*
                                     ( (Leptz(a,x)*log(Leptz(a,x)))/(1.0-Leptz(a,x)) - (Leptz(b,x)*log(Leptz(b,x)))/(1.0-Leptz(b,x)) );
                }
                else if (a != b && std::abs(Leptz(a,x)-Leptz(b,x)) > 0.01 && std::abs(1.0-Leptz(a,x)) > 0.01 && std::abs(1.0-Leptz(b,x)) <= 0.01) {
                    Leptfzc[x][a][b] = (-2.0*log(Leptz(a,x))*Leptz(a,x)+Leptz(a,x)+log(Leptz(a,x))-1.0)/(Leptz(a,x)-1.0)/(Leptz(a,x)-1.0);
                    Leptgzc[x][a][b] = -(MChi(a)*MChi(b)/mym_sn_sq(x))*(log(Leptz(a,x))*Leptz(a,x)-Leptz(a,x)+1.0)/(Leptz(a,x)-1.0)/(Leptz(a,x)-1.0);
                }
                else if (a != b && std::abs(Leptz(a,x)-Leptz(b,x)) > 0.01 && std::abs(1.0-Leptz(b,x)) > 0.01 && std::abs(1.0-Leptz(a,x)) <= 0.01) {
                    Leptfzc[x][a][b] = (-log(Leptz(b,x))*Leptz(b,x)*Leptz(b,x)+Leptz(b,x)-1.0)/(Leptz(b,x)-1.0)/(Leptz(b,x)-1.0);
                    Leptgzc[x][a][b] = -(MChi(a)*MChi(b)/mym_sn_sq(x))*(log(Leptz(b,x))*Leptz(b,x)-Leptz(b,x)+1.0)/(Leptz(b,x)-1.0)/(Leptz(b,x)-1.0);
                }
                else if ((a == b || std::abs(Leptz(a,x)-Leptz(b,x)) <= 0.01) && std::abs(1.0-Leptz(a,x)) > 0.01 && std::abs(1.0-Leptz(b,x)) > 0.01) {
                    Leptfzc[x][a][b] = (-Leptz(a,x)*Leptz(a,x)+Leptz(a,x)+log(Leptz(a,x)))/(Leptz(a,x)-1.0)/(Leptz(a,x)-1.0);
                    Leptgzc[x][a][b] = (MChi(b)*MChi(b)/mym_sn_sq(x))*((-Leptz(a,x)+log(Leptz(a,x))+1.0))/(Leptz(a,x)-1.0)/(Leptz(a,x)-1.0);
                }
                else {
                    Leptfzc[x][a][b] = -3.0/2.0-(Leptz(a,x)-1.0)/3.0;
                    Leptgzc[x][a][b] = -MChi(a)*MChi(b)/(2.0*mym_sn_sq(x));
                }
            }
        }
    }

    if (li_to_lj == 1) // mu -> 3e
    {
        //     Neutralino contributions
        gslpp::complex ZpengMuEEENR = 0.0;
        gslpp::complex ZpengMuEEENL = 0.0;
        for (int x=0;x<6;x++) {
            for (int a=0;a<4;a++) {
                for (int b=0;b<4;b++) {
                    ZpengMuEEENR = ZpengMuEEENR - NLlE(a,x)*NLlMU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                    ZpengMuEEENL = ZpengMuEEENL + NRlE(a,x)*NRlMU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                }
            }
        }
        //     Chargino contributions
        gslpp::complex ZpengMuEEEC = 0.0;
        for (int x=0;x<3;x++) {
            for (int a=0;a<2;a++) {
                for (int b=0;b<2;b++) {
                    ZpengMuEEEC = ZpengMuEEEC + CRlE(a,x)*CRlMU(b,x)*piconst*
                                                    (myU(a,1)*myU(b,1)*gslpp::complex(std::real(Leptgzc[x][a][b]),std::imag(Leptgzc[x][a][b]))
                                                     -myV(a,1)*myV(b,1)*gslpp::complex(std::real(Leptfzc[x][a][b]),std::imag(Leptfzc[x][a][b]))/2.0);
                }
            }
        }
        //     write FRR, FRL, FLR and FLL into a vector for mu->3e
        FFunctions.assign(0, ZpengMuEEENR/(MZ*MZ*ctw*ctw) );                                 //FRR
        FFunctions.assign(1, ZpengMuEEENR*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );                  //FRL
        FFunctions.assign(2, (ZpengMuEEENL + ZpengMuEEEC)/(MZ*MZ*ctw*ctw) );                 //FLR
        FFunctions.assign(3, (ZpengMuEEENL + ZpengMuEEEC)*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );  //FLL
    }
    if (li_to_lj == 2) // tau -> 3mu
    {
        //     Neutralino contributions
        gslpp::complex ZpengTauMuMuMuNR = 0.0;
        gslpp::complex ZpengTauMuMuMuNL = 0.0;
        for (int x=0;x<6;x++) {
            for (int a=0;a<4;a++) {
                for (int b=0;b<4;b++) {
                    ZpengTauMuMuMuNR = ZpengTauMuMuMuNR - NLlMU(a,x)*NLlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                    ZpengTauMuMuMuNL = ZpengTauMuMuMuNL + NRlMU(a,x)*NRlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                }
            }
        }
        //     Chargino contributions
        gslpp::complex ZpengTauMuMuMuC = 0.0;
        for (int x=0;x<3;x++) {
            for (int a=0;a<2;a++) {
                for (int b=0;b<2;b++) {
                    ZpengTauMuMuMuC = ZpengTauMuMuMuC + CRlMU(a,x)*CRlTAU(b,x)*piconst*
                                                    (myU(a,1)*myU(b,1)*gslpp::complex(std::real(Leptgzc[x][a][b]),std::imag(Leptgzc[x][a][b]))
                                                     -myV(a,1)*myV(b,1)*gslpp::complex(std::real(Leptfzc[x][a][b]),std::imag(Leptfzc[x][a][b]))/2.0);
                }
            }
        }
        //     write FRR, FRL, FLR and FLL into a vector for tau->3mu
        FFunctions.assign(0, ZpengTauMuMuMuNR/(MZ*MZ*ctw*ctw) );                                     //FRR
        FFunctions.assign(1, ZpengTauMuMuMuNR*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );                      //FRL
        FFunctions.assign(2, (ZpengTauMuMuMuNL + ZpengTauMuMuMuC)/(MZ*MZ*ctw*ctw) );                 //FLR
        FFunctions.assign(3, (ZpengTauMuMuMuNL + ZpengTauMuMuMuC)*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );  //FLL
    }
    if (li_to_lj == 3) // tau -> 3e
    {
        //      Neutralino contributions
        gslpp::complex ZpengTauEEENR = 0.0;
        gslpp::complex ZpengTauEEENL = 0.0;
        for (int x=0;x<6;x++) {
            for (int a=0;a<4;a++) {
                for (int b=0;b<4;b++) {
                    ZpengTauEEENR = ZpengTauEEENR - NLlE(a,x)*NLlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                    ZpengTauEEENL = ZpengTauEEENL + NRlE(a,x)*NRlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                }
            }
        }
        //      Chargino contributions
        gslpp::complex ZpengTauEEEC = 0.0;
        for (int x=0;x<3;x++) {
            for (int a=0;a<2;a++) {
                for (int b=0;b<2;b++) {
                    ZpengTauEEEC = ZpengTauEEEC + CRlE(a,x)*CRlTAU(b,x)*piconst*
                                                    (myU(a,1)*myU(b,1)*gslpp::complex(std::real(Leptgzc[x][a][b]),std::imag(Leptgzc[x][a][b]))
                                                     -myV(a,1)*myV(b,1)*gslpp::complex(std::real(Leptfzc[x][a][b]),std::imag(Leptfzc[x][a][b]))/2.0);
                }
            }
        }
        //     write FRR, FRL, FLR and FLL into a vector for tau->3e
        FFunctions.assign(0, ZpengTauEEENR/(MZ*MZ*ctw*ctw) );                                  //FRR
        FFunctions.assign(1, ZpengTauEEENR*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );                   //FRL
        FFunctions.assign(2, (ZpengTauEEENL + ZpengTauEEEC)/(MZ*MZ*ctw*ctw) );                 //FLR
        FFunctions.assign(3, (ZpengTauEEENL + ZpengTauEEEC)*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );  //FLL
    }
    if (li_to_lj == 4) // tau -> muee
    {
        //     Neutralino contributions
        gslpp::complex ZpengTauMuEENR = 0.0;
        gslpp::complex ZpengTauMuEENL = 0.0;
        for (int x=0;x<6;x++) {
            for (int a=0;a<4;a++) {
                for (int b=0;b<4;b++) {
                    ZpengTauMuEENR = ZpengTauMuEENR - NLlMU(a,x)*NLlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                    ZpengTauMuEENL = ZpengTauMuEENL + NRlMU(a,x)*NRlTAU(b,x)*piconst*
                                                    (ON(a,2)*ON(b,2)-ON(a,3)*ON(b,3))*gslpp::complex(std::real(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]),std::imag(Leptfzn[x][a][b]+2.0*Leptgzn[x][a][b]));
                }
            }
        }
        //     Chargino contributions
        gslpp::complex ZpengTauMuEEC = 0.0;
        for (int x=0;x<3;x++) {
            for (int a=0;a<2;a++) {
                for (int b=0;b<2;b++) {
                    ZpengTauMuEEC = ZpengTauMuEEC + CRlMU(a,x)*CRlTAU(b,x)*piconst*
                                                    (myU(a,1)*myU(b,1)*gslpp::complex(std::real(Leptgzc[x][a][b]),std::imag(Leptgzc[x][a][b]))
                                                     -myV(a,1)*myV(b,1)*gslpp::complex(std::real(Leptfzc[x][a][b]),std::imag(Leptfzc[x][a][b]))/2.0);
                }
            }
        }
        //     write FRR, FRL, FLR and FLL into a vector for tau->muee
        FFunctions.assign(0, ZpengTauMuEENR/(MZ*MZ*ctw*ctw) );                                   //FRR
        FFunctions.assign(1, ZpengTauMuEENR*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );                    //FRL
        FFunctions.assign(2, (ZpengTauMuEENL + ZpengTauMuEEC)/(MZ*MZ*ctw*ctw) );                 //FLR
        FFunctions.assign(3, (ZpengTauMuEENL + ZpengTauMuEEC)*(-0.5+sw2)/(MZ*MZ*sw2*ctw*ctw) );  //FLL
    }

    return(FFunctions);
}

gslpp::vector<gslpp::complex> SUSYMatching::gminus2mu() {

    NeutralinoRemixing();

    double MW = mySUSY.Mw_tree();
    double pi = M_PI;
    double piconst = 1.0/(32.0 * pi * pi);
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double stw = sqrt(sw2);
    double ctw = sqrt(1.0 - sw2);
    double ttw = stw/ctw;
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();

    double cdenc = sqrt(2.0)*MW*cosb;
    double cdenn = MW*cosb;
    double g2 = gW;
    double g2t = g2/sqrt(2.0);

    gslpp::vector<gslpp::complex> gminus2mu(2, 0.);

    //     Neutralino-Fermion-Sfermion couplings
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            //     LL + RL TYPE MI
            NRlMU.assign(a, x, -(g2t)*((-ON(a, 1) - ON(a, 0)*ttw)*myRl(x, 1) + (mMU/cdenn)*ON(a, 2)*myRl(x, 4)));
            //     RL + RR TYPE MI
            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2)*myRl(x, 1) + 2.0*ON(a, 0)*ttw*myRl(x, 4)));
//    Commented expressions might be useful for complex neutralino mixing matrices
//            NLlMU.assign(a, x, -(g2t)*((mMU/cdenn)*ON(a, 2).conjugate()*myRl(x, 1) + 2.0*ON(a, 0).conjugate()*ttw*myRl(x, 4)));
        }
    }

    //     Chargino-Fermion-Sfermion couplings
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            //     LL-TYPE
            CRlMU.assign(a, x, - (g2*myV(a, 0)*myRn(x, 1)));
            //     LR-TYPE
            CLlMU.assign(a, x, g2*mMU/cdenc*myU(a, 1).conjugate()*myRn(x, 1));
        }
    }

    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            Lepty.assign(a, x, MNeig(a) * MNeig(a) / mym_se_sq(x) );
        }
    }

    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            Leptz.assign(a, x, MChi(a) * MChi(a) / mym_sn_sq(x) );
        }
    }

    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            if (fabs(1.0 - Lepty(a, x)) > 0.01) {
                Leptf1.assign(a, x, ((1.0 - 6.0*Lepty(a, x) + 3.0 * pow(Lepty(a, x),2.0) + 
                                      2.0*pow(Lepty(a, x),3.0) - 6.0*pow(Lepty(a,x),2.0)*log(Lepty(a, x))))/
                                    (6.0 * pow((1.0 - Lepty(a,x)),4.0)) );
                Leptf2.assign(a, x, (1.0 - pow(Lepty(a, x),2.0) +  2.0 * Lepty(a, x) * log(Lepty(a, x)))/
                                    (pow((1.0-Lepty(a, x)),3.0)));
            }
            else {
                Leptf1.assign(a, x, 1.0/12.0 - (Lepty(a, x) - 1.0)/30.0);
                Leptf2.assign(a, x, 1.0/3.0 - (Lepty(a, x) - 1.0)/6.0);
            }
        }
    }

    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            if(fabs(1.0-Leptz(a, x)) > 0.01) {
                Leptf3.assign(a, x, ((2.0 + 3.0*Leptz(a, x) - 6.0*pow(Leptz(a, x),2.0) 
                                      + pow(Leptz(a, x),3.0) + 6.0*Leptz(a, x)*log(Leptz(a, x)))/
		                     (6.0*pow((1.0 - Leptz(a, x)),4.0))) );
                Leptf4.assign(a, x, ((-3.0 + 4.0*Leptz(a, x) - pow(Leptz(a, x),2.0)
                                      - 2.0*log(Leptz(a, x)))/
                                     pow((1.0 - Leptz(a, x)),3.0)) );
            }
            else {
                Leptf3.assign(a, x, 1.0/12.0 - (Leptz(a, x) - 1.0)/20.0 );
                Leptf4.assign(a, x, 2.0/3.0 - (Leptz(a, x) - 1.0)/2.0 );
            }
        }
    }

    //      Neutralino contributions
    gslpp::complex g2ARN = 0.0;
    gslpp::complex g2ALN = 0.0;
    for (int a=0;a<4;a++) {
        for (int x=0;x<6;x++) {
            g2ARN = g2ARN -mMU*mMU*piconst*(4.0*NRlMU(a,x)*NRlMU(a,x).conjugate()*Leptf1(a,x)
                                            +2.0*NRlMU(a,x)*NLlMU(a,x).conjugate()*(MNeig(a)/mMU)*Leptf2(a,x))/mym_se_sq(x);
            g2ALN = g2ALN -mMU*mMU*piconst*(4.0*NLlMU(a,x)*NLlMU(a,x).conjugate()*Leptf1(a,x)
                                            +2.0*NLlMU(a,x)*NRlMU(a,x).conjugate()*(MNeig(a)/mMU)*Leptf2(a,x))/mym_se_sq(x);
        }
    }

    //      Chargino contributions
    gslpp::complex g2ARC = 0.0;
    gslpp::complex g2ALC = 0.0;
    for (int a=0;a<2;a++) {
        for (int x=0;x<3;x++) {
            g2ARC = g2ARC +mMU*mMU*piconst*(4.0*CRlMU(a,x)*CRlMU(a,x).conjugate()*Leptf3(a,x)
                                            +2.0*CRlMU(a,x)*CLlMU(a,x).conjugate()*(MChi(a)/mMU)*Leptf4(a,x))/mym_sn_sq(x);
            g2ALC = g2ALC +mMU*mMU*piconst*(4.0*CLlMU(a,x)*CLlMU(a,x).conjugate()*Leptf3(a,x)
                                            +2.0*CLlMU(a,x)*CRlMU(a,x).conjugate()*(MChi(a)/mMU)*Leptf4(a,x))/mym_sn_sq(x);
        }
    }

    //     write R and L contributions to the muon g-2 into a vector
    gminus2mu.assign(0, g2ARN + g2ARC );    //g-2_muR
    gminus2mu.assign(1, g2ALN + g2ALC );    //g-2_muL

    return(gminus2mu);
}

gslpp::vector<gslpp::complex> SUSYMatching::C7_Lepton(int li_to_lj) {

    gslpp::vector<gslpp::complex> C7(2, 0.);

    //     write C7 and C7' into a vector
    C7.assign(0, -0.5*AFunctions(li_to_lj)(2) );
    C7.assign(1, -0.5*AFunctions(li_to_lj)(3) );
    return(C7);
}

gslpp::vector<gslpp::complex> SUSYMatching::C9_Lepton(int li_to_lj) {

    double MZ = mySUSY.getMz();
    double MW = mySUSY.Mw_tree();
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double cw2 = 1.0 - sw2;

    gslpp::vector<gslpp::complex> C9(2, 0.);

        //     write C9 and C9' into a vector
        C9.assign(0, AFunctions(li_to_lj)(1)
                     + (1.0-1.0/(4.0*sw2))/(cw2*MZ*MZ)*FFunctions(li_to_lj)(1)
                     + (BFunctions(li_to_lj)(1)+BFunctions(li_to_lj)(3))/2.0 );
        C9.assign(1, AFunctions(li_to_lj)(0)
                     + (1.0-1.0/(4.0*sw2))/(cw2*MZ*MZ)*FFunctions(li_to_lj)(0)
                     + (BFunctions(li_to_lj)(0)+BFunctions(li_to_lj)(2))/2.0 );

    return(C9);
}

gslpp::vector<gslpp::complex> SUSYMatching::C10_Lepton(int li_to_lj) {

    double MZ = mySUSY.getMz();
    double MW = mySUSY.Mw_tree();
    double sw2 = mySUSY.StandardModel::sW2(MW);
    double cw2 = 1.0 - sw2;

    gslpp::vector<gslpp::complex> C10(2, 0.);

        //     write C10 and C10' into a vector
        C10.assign(0, 1.0/(4.0*sw2*cw2*MZ*MZ)*FFunctions(li_to_lj)(1)
                     + (BFunctions(li_to_lj)(3)-BFunctions(li_to_lj)(1))/2.0 );
        C10.assign(1, 1.0/(4.0*sw2*cw2*MZ*MZ)*FFunctions(li_to_lj)(0)
                     + (BFunctions(li_to_lj)(0)-BFunctions(li_to_lj)(2))/2.0 );

    return(C10);
}

gslpp::vector<gslpp::complex> SUSYMatching::CS_Lepton(int li_to_lj) {

    gslpp::vector<gslpp::complex> CS(2, 0.);

        //     write CS and CS' into a vector
        CS.assign(0, BFunctions(li_to_lj)(4)/2.0 );
        CS.assign(1, BFunctions(li_to_lj)(5)/2.0 );

    return(CS);
}

gslpp::vector<gslpp::complex> SUSYMatching::CP_Lepton(int li_to_lj) {

    gslpp::vector<gslpp::complex> CP(2, 0.);

        //     write CP and CP' into a vector
        CP.assign(0, BFunctions(li_to_lj)(4)/2.0 );
        CP.assign(1, -BFunctions(li_to_lj)(5)/2.0 );

    return(CP);
}

gslpp::vector<gslpp::complex> SUSYMatching::CT_Lepton(int li_to_lj) {

    gslpp::vector<gslpp::complex> CT(2, 0.);

        //     write CT and CT5 into a vector
        CT.assign(0, (BFunctions(li_to_lj)(6)+BFunctions(li_to_lj)(7))/2.0 );
        CT.assign(1, (BFunctions(li_to_lj)(6)+BFunctions(li_to_lj)(7))/2.0 );

    return(CT);
}

std::vector<WilsonCoefficient>& SUSYMatching::CMDLij(int li_lj) {

    vmcDLij = StandardModelMatching::CMDLij(li_lj);

    gslpp::vector<gslpp::complex> C7vector=C7_Lepton(li_lj);
    switch (mcDLij.getOrder()) {
        case LO:
            mcDLij.setCoeff(0, C7vector(0), LO);
            mcDLij.setCoeff(1, C7vector(1), LO);
            break;
        case NNLO:
        case NLO:
        default:
            std::stringstream out;
            out << mcDLij.getOrder();
            throw std::runtime_error("SSUSYlMatching::CMDLij(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }

    vmcDLij.push_back(mcDLij);
    return(vmcDLij);

}

std::vector<WilsonCoefficient>& SUSYMatching::CMDLi3j(int li_lj) {

    vmcDLi3j = StandardModelMatching::CMDLi3j(li_lj);

    gslpp::vector<gslpp::complex> Avector=AFunctions(li_lj);
    gslpp::vector<gslpp::complex> Bvector=BFunctions(li_lj);
    gslpp::vector<gslpp::complex> BHvector=BHFunctions(li_lj);
    gslpp::vector<gslpp::complex> Fvector=FFunctions(li_lj);
    switch (mcDLi3j.getOrder()) {
        case LO:
//            mcDLi3j.setCoeff(0, C7_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(1, C7_Lepton(li_lj)(1), LO);
//            mcDLi3j.setCoeff(2, C9_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(3, C9_Lepton(li_lj)(1), LO);
//            mcDLi3j.setCoeff(4, C10_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(5, C10_Lepton(li_lj)(1), LO);
//            mcDLi3j.setCoeff(6, CS_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(7, CS_Lepton(li_lj)(1), LO);
//            mcDLi3j.setCoeff(8, CP_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(9, CP_Lepton(li_lj)(1), LO);
//            mcDLi3j.setCoeff(10, CT_Lepton(li_lj)(0), LO);
//            mcDLi3j.setCoeff(11, CT_Lepton(li_lj)(1), LO);
            mcDLi3j.setCoeff(0, Avector(0), LO);      //A1R
            mcDLi3j.setCoeff(1, Avector(1), LO);      //A1L
            mcDLi3j.setCoeff(2, Avector(2), LO);      //A2R
            mcDLi3j.setCoeff(3, Avector(3), LO);      //A2L
            mcDLi3j.setCoeff(4, Bvector(0), LO);      //B1R
            mcDLi3j.setCoeff(5, Bvector(1), LO);      //B1L
            mcDLi3j.setCoeff(6, Bvector(2), LO);      //B2R
            mcDLi3j.setCoeff(7, Bvector(3), LO);      //B2L
            mcDLi3j.setCoeff(8, Bvector(4), LO);      //B3R
            mcDLi3j.setCoeff(9, Bvector(5), LO);      //B3L
            mcDLi3j.setCoeff(10, Bvector(6), LO);     //B4R
            mcDLi3j.setCoeff(11, Bvector(7), LO);     //B4L
            mcDLi3j.setCoeff(12, BHvector(0), LO);    //BH2R
            mcDLi3j.setCoeff(13, BHvector(1), LO);    //BH2L
            mcDLi3j.setCoeff(14, BHvector(2), LO);    //BH3R
            mcDLi3j.setCoeff(15, BHvector(3), LO);    //BH3L
            mcDLi3j.setCoeff(16, Fvector(0), LO);     //FRR
            mcDLi3j.setCoeff(17, Fvector(1), LO);     //FRL
            mcDLi3j.setCoeff(18, Fvector(2), LO);     //FLR
            mcDLi3j.setCoeff(19, Fvector(3), LO);     //FLL
            break;
        case NNLO:
        case NLO:
        default:
            std::stringstream out;
            out << mcDLi3j.getOrder();
            throw std::runtime_error("SSUSYlMatching::CMDLi3j(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }

    vmcDLi3j.push_back(mcDLi3j);
    return(vmcDLi3j);

}

std::vector<WilsonCoefficient>& SUSYMatching::CMmueconv() {

    double MW = mySUSY.Mw_tree();
    double sw2 = mySUSY.StandardModel::sW2(MW);

    vmcmueconv = StandardModelMatching::CMmueconv();

    gslpp::vector<gslpp::complex> Avector=AFunctions(1);
    gslpp::vector<gslpp::complex> Dvector=DFunctions();
    gslpp::vector<gslpp::complex> Fvector=FFunctions(1);
    switch (mcmueconv.getOrder()) {
        case LO:
            mcmueconv.setCoeff(0, Avector(0), LO);                                        //A1R
            mcmueconv.setCoeff(1, Avector(1), LO);                                        //A1L
            mcmueconv.setCoeff(2, Avector(2), LO);                                        //A2R
            mcmueconv.setCoeff(3, Avector(3), LO);                                        //A2L
            mcmueconv.setCoeff(4, Dvector(0)+(1.0/(4.0*sw2)-2.0/3.0)*Fvector(0), LO);     //DbaruR
            mcmueconv.setCoeff(5, Dvector(1)+(1.0/(4.0*sw2)-2.0/3.0)*Fvector(2), LO);     //DbaruL
            mcmueconv.setCoeff(6, Dvector(2)+(-1.0/(4.0*sw2)+1.0/3.0)*Fvector(0), LO);    //DbardR
            mcmueconv.setCoeff(7, Dvector(3)+(-1.0/(4.0*sw2)+1.0/3.0)*Fvector(2), LO);    //DbardL
            break;
        case NNLO:
        case NLO:
        default:
            std::stringstream out;
            out << mcmueconv.getOrder();
            throw std::runtime_error("SSUSYlMatching::CMmueconv(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }

    vmcmueconv.push_back(mcmueconv);
    return(vmcmueconv);

}

std::vector<WilsonCoefficient>& SUSYMatching::CMgminus2mu() {

    vmcgminus2mu = StandardModelMatching::CMmueconv();

    gslpp::vector<gslpp::complex> gminus2muvector=gminus2mu();
    switch (mcgminus2mu.getOrder()) {
        case LO:
            mcgminus2mu.setCoeff(0, gminus2muvector(0), LO);  //g-2_muR
            mcgminus2mu.setCoeff(1, gminus2muvector(1), LO);  //g-2_muL
            break;
        case NNLO:
        case NLO:
        default:
            std::stringstream out;
            out << mcgminus2mu.getOrder();
            throw std::runtime_error("SSUSYlMatching::CMgminus2mu(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }

    vmcgminus2mu.push_back(mcgminus2mu);
    return(vmcgminus2mu);

}
