/* 
 * File:   SUSYMatching.cpp
 * Author: girardi_mac
 * 
 * Created on 14 maggio 2012, 11.53
 */

#include "SUSYMatching.h"
#include "SUSY.h"
#include <math.h>
#include <stdexcept>

SUSYMatching::SUSYMatching(const SUSY & SUSY_i) : StandardModelMatching(SUSY_i), mySUSY(SUSY_i), 
mcdbd2Hp(5, NDR, NLO), mcdbd2gg(5, NDR, NLO), mcdbd2ChiChi(5, NDR, NLO),
mcdbd2Chi0Chi0(5, NDR, NLO), mcdbd2Chi0g(5, NDR, NLO),
mcdbd2HpT(5, NDR, NLO), mcdbd2ggT(5, NDR, NLO), mcdbd2ChiChiT(5, NDR, NLO),
mcdbd2Chi0Chi0T(5, NDR, NLO), mcdbd2Chi0gT(5, NDR, NLO),
mcdbs2Hp(5, NDR, NLO), mcdbs2gg(5, NDR, NLO), mcdbs2ChiChi(5, NDR, NLO),
mcdbs2Chi0Chi0(5, NDR, NLO), mcdbs2Chi0g(5, NDR, NLO),
mcdbs2HpT(5, NDR, NLO), mcdbs2ggT(5, NDR, NLO), mcdbs2ChiChiT(5, NDR, NLO),
mcdbs2Chi0Chi0T(5, NDR, NLO), mcdbs2Chi0gT(5, NDR, NLO),
mcdk2Hp(5, NDR, NLO), mcdk2gg(5, NDR, NLO), mcdk2ChiChi(5, NDR, NLO),
mcdk2Chi0Chi0(5, NDR, NLO), mcdk2Chi0g(5, NDR, NLO),
mcdk2HpT(5, NDR, NLO), mcdk2ggT(5, NDR, NLO), mcdk2ChiChiT(5, NDR, NLO),
mcdk2Chi0Chi0T(5, NDR, NLO), mcdk2Chi0gT(5, NDR, NLO),        
mcdd2Hp(5, NDR, NLO), mcdd2gg(5, NDR, NLO), mcdd2ChiChi(5, NDR, NLO),
mcdd2Chi0Chi0(5, NDR, NLO), mcdd2Chi0g(5, NDR, NLO),
mcdd2HpT(5, NDR, NLO), mcdd2ggT(5, NDR, NLO), mcdd2ChiChiT(5, NDR, NLO),
mcdd2Chi0Chi0T(5, NDR, NLO), mcdd2Chi0gT(5, NDR, NLO),                
mcdbd2(5, NDR, NLO), mcdbs2(5, NDR, NLO),
mcdd2(5, NDR, NLO), mcdk2(5, NDR, NLO),
mcbsg(10, NDR, NLO), mcbnlep(10, NDR, NLO, NLO_ew),
mcbnlepCC(10, NDR, NLO), mcd1(10, NDR, NLO),
mcd1Buras(10, NDR, NLO), myCKM_cache(3, 3, 0.), DeltaMd_cache(3, 3, 0.),
mySUSYMQ(6, 0.), VUDHH_cache(6, 6, 0.) {

    swa = 0.;
    swb = 0.;
    swc = 0.;
    xcachea = 0.;
    xcacheb = 0.;
    xcachec = 0.;

    for (int j = 0; j < 10; j++) {
        CWbsgArrayLO[j] = 0.;
        CWbsgArrayNLO[j] = 0.;
        CWD1ArrayLO[j] = 0.;
        CWD1ArrayNLO[j] = 0.;
        CWbnlepArrayLOqcd[j] = 0.;
        CWbnlepArrayNLOqcd[j] = 0.;
        CWbnlepArrayLOew[j] = 0.;
        CWbnlepArrayNLOew[j] = 0.;
    };
}


/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
////////////// D0(x,y,z,t) and D2(x,y,z,t) limits

double SUSYMatching::D0N(double x, double y, double z, double t) {

    if ((fabs(z) < SUSYLEPS) || (fabs(t) < SUSYLEPS)) {

        return (0.);
    } else
        return (sqrt(z * t) * Dk(x, y, z, t, 0));

}

double SUSYMatching::D2LL0(double a, double b) {

    return ( (log(a) - log(b)) / (4. * (-a + b)) );

}


double SUSYMatching::DL0(double a, double b, double c, int k) {

    if (k == 2) {
        return ( (a * (-b + c) * log(a) + b * (a - c) * log(b)
                + (-a + b) * c * log(c)) / (4. * (a - b)*(-a + c)*(-b + c)));
    } else if (k == 1) {

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
        return ((-a * a * b + (a * a + b * b) * c - b * c * c) * log(a)
                + (-a + c)*((a - b)*(-b + c) + b * (-a + c) * log(b))
                - (a - b)*(a - b) * c * log(c)) / ((a - b)*(a - b)*(-a
                + c)*(-a + c)*(-b + c));
    } else
        if (k == 2) {

        if (a < SUSYLEPS) {

            return ( (-log(b) + log(c)) / (4. * (b - c)));
        } else
            if (b < SUSYLEPS) {

            return ( (-a + c + c * log(a) - c * log(c)) / (4. * (a - c) * (a - c)));
        } else
            if (c < SUSYLEPS) {

            return ( (-a + b + b * log(a) - b * log(b)) / (4. * (a - b) * (a - b)));
        }
        return (a * (-b + c)*(a * b + (a - 2. * b) * c) * log(a)
                + (-a + c)*(a * (a - b)*(-b + c)
                + b * b * (-a + c) * log(b))
                - (a - b)*(a - b) * c * c * log(c)) /
                (4. * (a - b)*(a - b)*(-a + c)*(-a + c)*(-b + c));
    }
    
    else {

        throw std::runtime_error("Error in DL(a,b,c,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::DLL(double a, double b, int k) {
    if (k == 0) {

        if (b < SUSYLEPS) {
            return (1. / (2. * a * a));
        } else
            return ( (-a * a + b * b + 2 * a * b * log(a) - 2 * a * b * log(b)) /
                (2. * a * (a - b)*(a - b)*(a - b)));
    } else
        if (k == 2) {

        if (b < SUSYLEPS) {

            return ( -1. / (8. * a));
        } else
            return (a * a - 4 * a * b + 3 * b * b + 2 * b * b * log(a)
                - 2 * b * b * log(b)) / (8. * (-a + b)*(-a + b)*(-a + b));
    }
    
    else {

        throw std::runtime_error("Error in DLL(a,b,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::DLLp(double a, double b, int k) {
    if (k == 0) {
        
        if ( (a < SUSYLEPS) || (b < SUSYLEPS) ) {
            
            throw std::runtime_error("Error in DLLp function, because the limit D0(0,0,b,b) and D0(a,a,0,0) are singular "); 
        }
        return (-2 * a + 2 * b + (a + b) * log(a) - (a + b) * log(b))
                / ((a - b)*(a - b)*(a - b));
    } else
        if (k == 2) {
            
            if ( a < SUSYLEPS ){
                
                return (-1. / (4. * b ));
            } else 
                if ( b < SUSYLEPS ){
                    
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
    } else
        if (k == 2) {
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
        } else
            return DL(x, z, t, k);
    }
    if (fabs(1. - z / x) < SUSYLEPS) {
        if (fabs(1. - t / y) < SUSYLEPS) {
            return DLLp(x, y, k);
        } else
            return DL(x, y, t, k);
    }
    if (fabs(1. - t / x) < SUSYLEPS) {
        if (fabs(1. - z / y) < SUSYLEPS) {
            return DLLp(x, y, k);
        } else
            return DL(x, z, y, k);
    }
    if ((fabs(1. - z / y) < SUSYLEPS)) {
        return DL(y, x, t, k);
    }
    if ((fabs(1. - t / y) < SUSYLEPS)) {
        return DL(y, z, x, k);
    }
    if (( fabs(1. - z / t) < SUSYLEPS)) {
        
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
    } else
        if (k == 2) {
        return ((t * t * log(t)) / ((-t + x)*(-t + y)*(-t + z))
                + (x * x * log(x)) / ((t - x)*(-x + y)*(-x + z))
                + (y * y * log(y)) / ((t - y)*(x - y)*(-y + z))
                + (z * z * log(z)) / ((t - z)*(x - z)*(y - z))) / 4.;
    }

    else {

        throw std::runtime_error("Error in Dk(x,y,z,t,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);

}

////////////////////////////////////////////////////////////////////////////////
///// Ck limits

double SUSYMatching::CL(double a, double b, int k) {

    double mu2R = mySUSY.GetQ() * mySUSY.GetQ();
    
    if (k == 0) {
        return ((a - b - b * log(a) + b * log(b)) / ((a - b)*(a - b)));

    } else
        if (k == 2) {

        return (-log(mu2R) + (a * (a - b) + a * (a - 2 * b) * log(a) + b * b * log(b)) /
                ((a - b)*(a - b)));
    }

    else {

        throw std::runtime_error("Error in CL(a,b,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::CLL(double a, int k) {

    double mu2R = mySUSY.GetQ() * mySUSY.GetQ();
    
    if (k == 0) {
        return (1. / (2. * a));

    } else
        if (k == 2) {

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

    double mu2R = mySUSY.GetQ() * mySUSY.GetQ();

    if ((fabs(1. - y / x) < SUSYLEPS)&(fabs(1. - z / x) < SUSYLEPS)) {

        return (CLL(x,k));

    } else
        if ((fabs(1. - y / x) < SUSYLEPS)) {

        return (CL(x,z,k));

    } else
        if ((fabs(1. - z / x) < SUSYLEPS)) {

        return (CL(x,y,k));

    } else
        if ((fabs(1. - z / y) < SUSYLEPS)) {

        return (CL(y,x,k));

    }

    if (k == 0) {
        return ((y * log(y / x)) / ((x - y)*(-y + z)) + (z * log(z / x)) /
                ((x - z)*(y - z)));
    } else
        if (k == 2) {
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

    double mu2R = mySUSY.GetQ() * mySUSY.GetQ();
    
    if (k == 0) {
        return (1. + log(a));
    } else
        if (k == 2) {
        return ((-log(mu2R) + 2. + log(a)) / 2.);
    }
    else {

        throw std::runtime_error("Error in BL(a,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

double SUSYMatching::Bk(double x, double y, int k) {

    if (fabs((1. - y / x)) < SUSYLEPS) {
        return (BL(x, k));

    }
    if (k == 0) {
        return (log(x) + (y * log(x / y)) / (x - y));
    } else
        if (k == 2) {

        return (1. / 4. + 0.5 * Ck(x, y, y, 2));
    } else {

        throw std::runtime_error("Error in Bk(x,y,k) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);

}
////////////////////////////////////////////////////////////////////////////////
//// Running quark masses to SUSY scale Q

void SUSYMatching::Comp_mySUSYMQ() {

    double Q = mySUSY.GetQ();
    
    
        mySUSYMQ(0) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.UP).getMass_scale(),mySUSY.getQuarks(mySUSY.UP).getMass());
        mySUSYMQ(1) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.DOWN).getMass_scale(),mySUSY.getQuarks(mySUSY.DOWN).getMass());
        mySUSYMQ(2) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.CHARM).getMass());
        mySUSYMQ(3) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.STRANGE).getMass_scale(),mySUSY.getQuarks(mySUSY.STRANGE).getMass());
        mySUSYMQ(4) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.TOP).getMass());
        mySUSYMQ(5) = mySUSY.Mrun(Q,mySUSY.getQuarks(mySUSY.BOTTOM).getMass());
}

////////////////////////////////////////////////////////////////////////////////
////  Quark self-energies from Buras arXiv:hep-ph/0210145v2

void SUSYMatching::Comp_DeltaMd() {

    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::vector<double> myMU2Squarks(6, 0.);
    gslpp::vector<double> myMD2Squarks(6, 0.);
    gslpp::vector<double> MChi(2, 0.);
    gslpp::vector<double> MChi0(4, 0.);
    double Q = mySUSY.GetQ();
    double Als = mySUSY.Als(Q);
    double Mg = mySUSY.getM3();

    MChi = mySUSY.getMch();
    MChi0 = mySUSY.getMneu();
    myMU2Squarks = mySUSY.getMsu2();
    myMD2Squarks = mySUSY.getMsu2();
    myRd = mySUSY.getRd();
    int k, l, I, J;

    for (J = 0; J < 3; J++) {
        for (I = 0; I < 3; I++) {

            complex temp(0., 0., false);

            for (k = 0; k < 6; k++) {

                temp += -2. / 3. * Als / M_PI * Mg * myRd(k, J + 3)
                        * myRd(k, I).conjugate() * Bk(Mg * Mg, myMD2Squarks(k), 0);

                //                + 1. / 3. * Als / M_PI * myRd(k, J + 3) * myRd(k, I + 3).conjugate()
                //                * Bk(Mg * Mg, myMD2Squarks(k), 1)
                //                * mySUSYMQ(2 * I + 1).getMass()
                //
                //                + 1. / 3. * Als / M_PI * myRd(k, J) * myRd(k, I).conjugate() *
                //                Bk(Mg * Mg, myMD2Squarks(k), 1) * mySUSYMQ(2 * J + 1).getMass();

            }
            
//            double prova = temp.real();
          
            for (k = 0; k < 6; k++) {
                for (l = 0; l < 4; l++) {

                    temp += 1. / (16. * M_PI * M_PI) * VdDNR(J, k, l, 0).conjugate()
                            * VdDNL(I, k, l, 0)
                            * MChi0(l) * Bk(MChi0(l) * MChi0(l), myMD2Squarks(k), 0);

                    //                    + 1. / (32. * M_PI * M_PI) * VdDNR(J, k, l, 0).conjugate()
                    //                    * VdDNR(I, k, l, 0) * Bk(MChi0(l) * MChi0(l), myMD2Squarks(k), 1)
                    //                    * mySUSYMQ(2 * I + 1).getMass()
                    //
                    //                    + 1. / (32. * M_PI * M_PI) * VdDNL(J, k, l, 0).conjugate() *
                    //                    VdDNL(I, k, l, 0) * Bk(MChi0(l) * MChi0(l), myMD2Squarks(k), 1) *
                    //                    mySUSYMQ(2 * J + 1).getMass();

                }

            }
            
//            double prova2 = temp.real();
            for (k = 0; k < 6; k++) {
                for (l = 0; l < 2; l++) {

                    temp += 1. / (16. * M_PI * M_PI) * VdUCR(J, k, l, 0).conjugate() * VdUCL(I, k, l) *
                            MChi(l) * Bk(MChi(l) * MChi(l), myMU2Squarks(k), 0);

                    //                    + 1. / (32. * M_PI * M_PI) * VdUCR(J, k, l, 0).conjugate() *
                    //                    VdUCR(I, k, l, 0) * Bk(MChi(l) * MChi(l), myMU2Squarks(k), 1)
                    //                    * mySUSYMQ(2 * I + 1).getMass()
                    //
                    //                    + 1. / (32. * M_PI * M_PI) * VdUCL(J, k, l).conjugate() *
                    //                    VdUCL(I, k, l) * Bk(MChi(l) * MChi(l), myMU2Squarks(k), 1);
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

gslpp::complex SUSYMatching::Eps_J(int J) {

  
    return (DeltaMd(J, J) / (mySUSY.getTanb()
            * mySUSYMQ(2 * J + 1)));

}

gslpp::complex SUSYMatching::Lambda0EpsY(int J, int I){
    
    double v2 = mySUSY.v() * mySUSY.getSinb();
    double mtop = mySUSYMQ(4);
    
    return (DeltaMd(J,I)/(mySUSY.getTanb()
            * mySUSYMQ(2 * J + 1) * 2 /(v2 * v2)  
            * mtop * mtop ));
    
}


////////////////////////////////////////////////////////////////////////////////
///// mySUSY_CKM
//// CKM matrix with large tan beta corrections

gslpp::complex SUSYMatching::DeltaDL(int J, int I) {

    complex C(0., 0., false);

    if (J != I) {

        double mdJ = mySUSYMQ(2 * J + 1);
        double mdI = mySUSYMQ(2 * I + 1);
              
        return (-(mdJ * DeltaMd(J, I) + DeltaMd(I, J).conjugate() * mdI) /
                (mdJ * mdJ - mdI * mdI));
    } else
        if (J == I) {

        return (C);
    }
    else {

        throw std::runtime_error("Error in DeltaDL(J,I) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::DeltaDR(int J, int I){
    
   return (DeltaDL(I,J).conjugate()); 
   
}

void SUSYMatching::Comp_mySUSY_CKM() {


    gslpp::matrix<complex> myCKM(3, 3, 0.);
    gslpp::matrix<complex> myTempCKM(3, 3, 0.);
    mySUSY.getCKM().getCKM(myCKM);
    mySUSY.getCKM().getCKM(myTempCKM);
 
    gslpp::matrix<complex> DeltaCKM(3, 3, 0.);
    
    complex Delta_CKM_IJ(0.,0.,false);
    
    int l, I, J;

    for (I = 0; I < 3; I++) {
        for (J = 0; J < 3; J++) {
            
            complex prova(0.,0.,false);
            


            for (l = 0; l < 3; l++) {

                Delta_CKM_IJ += -myTempCKM(I, l) * DeltaDL(l, J);
                
                prova += -myTempCKM(I, l) * DeltaDL(l, J);
                
            }
            
            
            
            myCKM_cache.assign(I, J, myCKM(I, J) + Delta_CKM_IJ);
            
            /*** test - lines ***/
            
            DeltaCKM.assign(I, J, Delta_CKM_IJ.abs() / myCKM(I, J).abs() * 100);
            
            /** end - test */
            
            Delta_CKM_IJ.assign(0.,0.,0);
            
           
        }
    }
    
    /** test - lines  **/
    
    //std::cout << "Delta_CKM % = " << DeltaCKM << std::endl;
    
    /** end - test */
    
}

gslpp::matrix<complex> SUSYMatching::mySUSY_CKM() {

    return (myCKM_cache);

}
////////////////////////////////////////////////////////////////////////////////
///// VChiUdL element j,k,b <-> VdUCL element b,k,j

void SUSYMatching::Comp_VdUCL() {


    complex VdUCL_bkj(0., 0., false);
    gslpp::matrix<complex> myRu(6, 6, 0.);
    gslpp::matrix<complex> myV(2, 2, 0.);
    gslpp::matrix<complex> myCKM(3, 3, 0.);

    mySUSY.getCKM().getCKM(myCKM);

    int l;
    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    double v2 = mySUSY.v() * mySUSY.getSinb();
    myRu = mySUSY.getRu();
    myV = mySUSY.getV();



    int b, k, j;

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

    complex VdUCR_bkj(0., 0., false);
    gslpp::matrix<complex> myCKM(3, 3, 0.);
    double tanb = mySUSY.getTanb();

    mySUSY.getCKM().getCKM(myCKM);

    gslpp::matrix<complex> myRu(6, 6, 0.);
    myRu = mySUSY.getRu();
    gslpp::matrix<complex> myU(2, 2, 0.);
    myU = mySUSY.getU();

    double v1 = mySUSY.v() * mySUSY.getCosb();
    complex Mdb(0., 0., false);
    complex Mdp(0., 0., false);
    int l, p, b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {


                for (l = 0; l < 3; l++) {

                    switch (flag) {
                        case 1:

                            Mdb += mySUSYMQ(2 * b + 1) /
                                    (1. + Eps_J(b) * tanb);

                            VdUCR_bkj += Mdb * myCKM(l, b);

                            for (p = 0; p < 3; p++) {

                                Mdp += mySUSYMQ(2 * p + 1) /
                                        (1. + Eps_J(p) * tanb);

                                VdUCR_bkj += myCKM(l, p) * (Mdp * DeltaDR(p, b) -
                                        DeltaDL(p, b) * Mdb);
                            }

                            VdUCR_bkj *= sqrt(2.) / v1 * myRu(k, l) * myU(j, 1);

                        case 0:
                            VdUCR_bkj += sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1) * myRu(k, l) * myU(j, 1)
                                    * myCKM(l, b);
                    }
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
    
    
    complex VdDNL_bkj(0., 0., false);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::matrix<complex> myN(4, 4, 0.);


    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    double v1 = mySUSY.v() * mySUSY.getCosb();
    double CosThetaW = sqrt(mySUSY.cW2());
    double SinThetaW = sqrt(mySUSY.sW2());
    int l;

    myRd = mySUSY.getRd();
    myN = mySUSY.getN();

    
    int b, k, j;
    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {
                

                    switch (flag) {

                        case 1:

                            VdDNL_bkj += -gW / sqrt(2.) * myRd(k, b).conjugate() * (1. / 3.
                                    * SinThetaW / CosThetaW * myN(j, 1).conjugate() -
                                    myN(j, 2).conjugate()) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1) /
                                    (1 + Eps_J(b) * mySUSY.getTanb())
                                    * myRd(k, b + 3).conjugate() * myN(j, 3).conjugate();
                            for (l = 0; l < 3; l++) {

                                // correzione vertici neutralini calcolate seguendo le indicazioni di Buras    

                                VdDNL_bkj += (-gW / sqrt(2.) * myRd(k, l).conjugate() * (1. / 3.
                                        * SinThetaW / CosThetaW * myN(j, 1).conjugate() -
                                        myN(j, 2).conjugate()) - sqrt(2.) / v1
                                        * mySUSYMQ(2 * l + 1) /
                                        (1 + Eps_J(l) * mySUSY.getTanb())
                                        * myRd(k, l + 3).conjugate() * myN(j, 3).conjugate()) *
                                        DeltaDL(l, b);
                            }

                        case 0:
                            VdDNL_bkj += -gW / sqrt(2.) * myRd(k, b).conjugate() * (1. / 3.
                                    * SinThetaW / CosThetaW * myN(j, 1).conjugate() -
                                    myN(j, 2).conjugate()) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1)
                                    * myRd(k, b + 3).conjugate() * myN(j, 3).conjugate();

                    }

                    VdDNL_cache[b][k][j][flag] = VdDNL_bkj;
                    VdDNL_bkj.assign(0., 0., 0);

            }

        }
    }
}

gslpp::complex SUSYMatching::VdDNL(int b, int k, int j, int flag) {


    return (VdDNL_cache[b][k][j][flag]);

}

////////////////////////////////////////////////////////////////////////////////
///// VdDNR element b,k,j <-> VChiDdR element j,k,b

void SUSYMatching::Comp_VdDNR(int flag) {

    complex VdDNR_bkj(0., 0., false);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::matrix<complex> myN(4, 4, 0.);
    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    double v1 = mySUSY.v() * mySUSY.getCosb();
    double CosThetaW = sqrt(mySUSY.cW2());
    double SinThetaW = sqrt(mySUSY.sW2());
    int l;
    myRd = mySUSY.getRd();
    myN = mySUSY.getN();


    int b, k, j;
    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {
                

                    switch (flag) {

                        case 1:

                            VdDNR_bkj += -sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                    myRd(k, b + 3).conjugate() * myN(j, 1) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1) /
                                    (1 + Eps_J(b) * mySUSY.getTanb()) * myRd(k, b).conjugate() *
                                    myN(j, 3);

                            for (l = 0; l < 3; l++) {

                                // correzione vertici neutralini calcolate seguendo le indicazioni di Buras   

                                VdDNR_bkj += (-sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                        myRd(k, l + 3).conjugate() * myN(j, 1) - sqrt(2.) / v1
                                        * mySUSYMQ(2 * l + 1) /
                                        (1 + Eps_J(l) * mySUSY.getTanb()) * myRd(k, l).conjugate() *
                                        myN(j, 3)) * DeltaDR(l, b);
                            }

                        case 0:
                            VdDNR_bkj += -sqrt(2.) / 3. * gW * SinThetaW / CosThetaW *
                                    myRd(k, b + 3).conjugate() * myN(j, 1) - sqrt(2.) / v1
                                    * mySUSYMQ(2 * b + 1)
                                    * myRd(k, b).conjugate() * myN(j, 3);

                    }



                    VdDNR_cache[b][k][j][flag] = VdDNR_bkj;
                    VdDNR_bkj.assign(0., 0., 0);

            }
        }
    }

}


gslpp::complex SUSYMatching::VdDNR(int b, int k, int j, int flag) {

    return (VdDNR_cache[b][k][j][flag]);
    
}

////////////////////////////////////////////////////////////////////////////////
/////  Feynmann rules of uDC vertex (D mixing)

void SUSYMatching::Comp_VuDCL() {


    gslpp::matrix<complex> myCKM(3, 3, 0.);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::matrix<complex> myU(2, 2, 0.);
    complex VuDCL_bkj(0., 0., false);
    complex YdI(0., 0., false);
    myCKM = mySUSY_CKM();

    int I;
    double tanb = mySUSY.getTanb();
    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    double v1 = mySUSY.v() * mySUSY.getCosb();
    myRd = mySUSY.getRd();
    myU = mySUSY.getU();


    int b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {



                for (I = 0; I < 3; I++) {

                    YdI = sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) / (1 + Eps_J(I) * tanb);
                    VuDCL_bkj += -(gW * myRd(k, I).conjugate() * myU(j, 0).conjugate()
                            - YdI * myRd(k, I + 3).conjugate() * myU(j, 1).conjugate()) *
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


    gslpp::matrix<complex> myCKM(3, 3, 0.);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::matrix<complex> myV(2, 2, 0.);
    myRd = mySUSY.getRd();
    myV = mySUSY.getV();
    complex Yub(0., 0., false);
    complex VuDCR_bkj(0., 0., false);
    myCKM = mySUSY_CKM();
    int I;
    double v2 = mySUSY.v() * mySUSY.getSinb();


    int b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 2; j++) {



                Yub = sqrt(2.) / v2 * mySUSYMQ(2 * b); // b is the up quark type index

                for (I = 0; I < 3; I++) {

                    VuDCR_bkj += Yub * myRd(k, I).conjugate() * myV(j, 1) * myCKM(b, I).conjugate();

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
    } else
        if (Dmixingflag == 1) {

        return (VuDCL(b, k, j));
    } else {
            
             throw std::runtime_error("Error in VdUCL(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  "); 
        }
}

gslpp::complex SUSYMatching::VdUCR(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdUCR(b, k, j, flag));
    } else
        if (Dmixingflag == 1) {

        return (VuDCR(b, k, j));
    } else {
            
             throw std::runtime_error("Error in VdUCR(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  "); 
        }
}

/* Vertices uUN from Buras arXiv:hep-ph/0210145v2 in SLHA convention  
     usefull in D - Dbar mixing */

void SUSYMatching::Comp_VuUN(){
 
    double TanThetaW = sqrt(mySUSY.sW2() / mySUSY.cW2());
    double v = mySUSY.v();
    double v2 = v * mySUSY.getSinb();
    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    gslpp::matrix<complex> myN(4, 4, 0.);
    gslpp::matrix<complex> myRu(6, 6, 0.);
    myN = mySUSY.getN();
    myRu = mySUSY.getRu();
    complex VuUNL_bkj(0.,0.,false);
    complex VuUNR_bkj(0.,0.,false);
    double Yub;
    
    int b, k, j;

    for (b = 0; b < 3; b++) {
        for (k = 0; k < 6; k++) {
            for (j = 0; j < 4; j++) {

                Yub = sqrt(2.) / v2 * mySUSYMQ(2 * b);


                VuUNL_bkj += (-1. / sqrt(2.) * gW * myRu(k, b) * (1. / (3. * TanThetaW) *
                        myN(j, 0).conjugate() + myN(j, 1).conjugate())
                        - Yub * myRu(k, b + 3) * myN(j, 3).conjugate());


                VuUNR_bkj += (2. * sqrt(2.) / 3. * gW * TanThetaW * myRu(k, b + 3) *
                        myN(j, 0) - Yub * myRu(k, b) * myN(j, 3));

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

    } else if (chirality.compare("R") == 0) {

        return (VuUNR_cache[b][k][j]);

    } else {

        std::cout << " VuUN error in SUSYMatching.cpp" << std::endl;
        
    }

    return (EXIT_FAILURE);
    
}

/******************************************************************************/


gslpp::complex SUSYMatching::VdDNL(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdDNL(b, k, j, flag));

    } else if (Dmixingflag == 1) {

        return (VuUN(b, k, j, "L"));
    } else {

        throw std::runtime_error("Error in VdDNL(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::VdDNR(int b, int k, int j, int flag, int Dmixingflag) {

    if (Dmixingflag == 0) {

        return (VdDNR(b, k, j, flag));

    } else if (Dmixingflag == 1) {

        return (VuUN(b, k, j, "R"));
    } else {

        throw std::runtime_error("Error in VdDNR(b,k,j,flag,Dmixingflag) in SUSYMatching.cpp  ");
    }

    return (EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
////
////////////////////////////////////////////////////////////////////////////////
//// PGLR, PGRL, PHLR, PHRL element j,i

gslpp::complex SUSYMatching::PGLR(int j, int i) {

    gslpp::matrix<complex> myCKM(3, 3, 0.);
    mySUSY.getCKM().getCKM(myCKM);

    return (-sqrt(2.) / mySUSY.v() * myCKM(j, i) *
            mySUSYMQ(2 * i + 1));
}

gslpp::complex SUSYMatching::PGRL(int j, int i) {

    gslpp::matrix<complex> myCKM(3, 3, 0.);
    mySUSY.getCKM().getCKM(myCKM);
    
    return (sqrt(2.) / mySUSY.v() * myCKM(j, i) *
            mySUSYMQ(2 * i));

}

gslpp::complex SUSYMatching::PHLR(int j, int i) {

    gslpp::matrix<complex> myCKM(3, 3, 0.);
    gslpp::complex PHLR(0., 0., false);
    mySUSY.getCKM().getCKM(myCKM);
    int k;
    
    PHLR +=  myCKM(j, i) *
            mySUSYMQ(2 * i + 1) / (1. + Eps_J(i) * mySUSY.getTanb());

    for( k = 0;k < 3; k++){
        
        PHLR += myCKM(j,k) * mySUSYMQ(2 * k + 1) / 
                (1. + Eps_J(k) * mySUSY.getTanb()) * DeltaDR(k,i) - myCKM(j,k) * 
                DeltaDL(k,i) * mySUSYMQ(2 * i + 1) / 
                (1. + Eps_J(i) * mySUSY.getTanb());
        
    }
    
    PHLR *= sqrt(2.) / mySUSY.v() * mySUSY.getTanb();
    
    return(PHLR);

}


void SUSYMatching::Comp_VUDHH(){
    
    complex VUDHijH(0., 0., false);
    gslpp::matrix<complex> myRu(6, 6, 0.);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::matrix<complex> myTU(3, 3, 0.);
    gslpp::matrix<complex> myTD(3, 3, 0.);
    double v = mySUSY.v();
    double tanb = mySUSY.getTanb();
    double v1 = v * mySUSY.getCosb();
    double v2 = v * mySUSY.getSinb();
    complex YuJ(0., 0., false);
    complex YdI(0., 0., false);
    gslpp::matrix<complex> myCKM(3, 3, 0.);
    myCKM = mySUSY_CKM();
    myRu = mySUSY.getRu();
    myRd = mySUSY.getRd();
    myTU = mySUSY.GetTU();
    myTD = mySUSY.GetTD();
    gslpp::matrix<complex> ZH(2, 2, 0.);
    ZH.assign(0, 0, mySUSY.getSinb());
    ZH.assign(0, 1, -mySUSY.getCosb());
    ZH.assign(1, 0, mySUSY.getCosb());
    ZH.assign(1, 1, mySUSY.getSinb());
    int I, J, l, i, j;

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {


            for (I = 0; I < 3; I++) {
                for (J = 0; J < 3; J++) {
                    YuJ = sqrt(2.) / v2 * mySUSYMQ(2 * J);
                    YdI = sqrt(2.) / v1 * mySUSYMQ(2 * I + 1) /
                            (1 + Eps_J(I) * tanb);

                    VUDHijH += v / sqrt(2.) * YuJ * YdI * myCKM(J, I) * myRd(j, I + 3) *
                            myRu(i, J + 3) + 1. / sqrt(2.) * (v1 * YdI * YdI * ZH(0, 0) +
                            v2 * YuJ * YuJ * ZH(1, 0)) * myCKM(J, I) * myRd(j, I)
                            * myRu(i, J) + ZH(0, 0) * mySUSY.getMuH().conjugate() *
                            YuJ * myCKM(J, I) * myRu(i, J + 3) * myRd(j, I) +
                            ZH(1, 0) * mySUSY.getMuH() * YdI * myCKM(J, I) * myRu(i, J) *
                            myRd(j, I + 3);
                    for (l = 0; l < 3; l++) {

                        VUDHijH += ZH(1, 0) * myTU(l, J) * myCKM(l, I) * myRu(i, J + 3) *
                                myRd(j, I) + ZH(0, 0) * myTD(I, l).conjugate() * myCKM(J, l) *
                                myRu(i, J) * myRd(j, I + 3);

                    }

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


    complex DFHL_ji(0., 0., false);
    gslpp::matrix<complex> myRu(6, 6, 0.);
    gslpp::matrix<complex> myRd(6, 6, 0.);
    gslpp::vector<double> myMU2Squarks(6, 0.);
    gslpp::vector<double> myMD2Squarks(6, 0.);
    myMU2Squarks = mySUSY.getMsu2();
    myMD2Squarks = mySUSY.getMsd2();
    gslpp::vector<double> MChi0(4, 0.);
    myRu = mySUSY.getRu();
    myRd = mySUSY.getRd();
    double Q = mySUSY.GetQ();
    double Als = mySUSY.Als(Q);
    double Mg = mySUSY.getM3();
    complex Yuj(0., 0., false);
    complex Ydi(0., 0., false);
    double v = mySUSY.v();
    double v1 = v * mySUSY.getCosb();
    double v2 = v * mySUSY.getSinb();
    double tanb = mySUSY.getTanb();
    int m, l;

    Yuj = sqrt(2.) / v2 * mySUSYMQ(2 * j);
    Ydi = sqrt(2.) / v1 * mySUSYMQ(2 * i + 1) /
            (1 + Eps_J(i) * tanb);

    for (m = 0; m < 6; m++) {
        for (l = 0; l < 6; l++) {

            DFHL_ji += VUDHH(m, l)*(-2 * Als / (3. * M_PI) * Mg
                    * myRu(m, j + 3).conjugate() * myRd(l, i).conjugate() *
                    Ck(Mg * Mg, myMU2Squarks(m), myMD2Squarks(l), 0) +

                    1. / (16. * M_PI * M_PI) * Yuj * Ydi *
                    myRu(m, j).conjugate() * myRd(l, i + 3).conjugate() *
                    mySUSY.getMuH().conjugate() *
                    Ck(mySUSY.getMuH().abs2(), myMU2Squarks(m), myMD2Squarks(l), 0));

        }

    }
   return (DFHL_ji);
}


gslpp::complex SUSYMatching::PHRL(int j, int i){
    
    return (sqrt(2.)/(mySUSY.v() * mySUSY.getTanb()) * mySUSYMQ(2 * j) *
            mySUSY_CKM()(j,i) + DeltaFHL(j,i));
    
}

gslpp::complex SUSYMatching::PLRk(int j, int i, int k){
    
    if(k == 0){
        return (PHLR(j,i));
    }else
        if(k == 1){
            return (PGLR(j,i));
        } else {
            
             throw std::runtime_error("Error in PLRk(j,i,k) in SUSYMatching.cpp  "); 
        }
    
    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::PRLk(int j, int i, int k){
    
    if(k == 0){
        return (PHRL(j,i));
    } else
        if(k == 1){
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

    } else
        if (Dmixingflag == 1) {

        return (PLRk(i, j, k).conjugate());
    } else {
            
             throw std::runtime_error("Error in PRLk(j,i,k,Dmixingflag) in SUSYMatching.cpp  "); 
        }

}

gslpp::complex SUSYMatching::PLRk(int j, int i, int k, int Dmixingflag) {


    if (Dmixingflag == 0) {

        return (PLRk(j, i, k));

    } else
        if (Dmixingflag == 1) {

        return (PRLk(i, j, k).conjugate());
        } else {
            
             throw std::runtime_error("Error in PLRk(j,i,k,Dmixingflag) in SUSYMatching.cpp  "); 
        }
    
    return (EXIT_FAILURE);
}
////////////////////////////////////////////////////////////////////////////////
////// Double Penguin Functions



gslpp::complex SUSYMatching::xdS(int S){
    
    double M2A = mySUSY.getMHa() * mySUSY.getMHa();
    double M2Z = mySUSY.getMz() * mySUSY.getMz();
    double tanb = mySUSY.getTanb();

    
    
    double tan2alpha = 2. * tanb / (1 - tanb * tanb) * (M2A + M2Z) / (M2A - M2Z); 
    double tana = (-1. + sqrt(1. + tan2alpha * tan2alpha)) / tan2alpha;
    double sina = tana / sqrt(1. + tana * tana);
    double cosa = 1. / sqrt(1 + tana * tana);
    
    complex Cosa(cosa,0.,false);
    complex Sina(sina,0.,false);
    complex i(0.,1.,false);


    if (S == 0) {
        return (Cosa);
    }
    if (S == 1) {
        return (-Sina);
    }
    if (S == 2) {
        return ( i * mySUSY.getSinb());
    } else {
            
             throw std::runtime_error("Error in xdS(S) in SUSYMatching.cpp  "); 
        }
    
    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::xuS(int S){
    
    double M2A = mySUSY.getMHa() * mySUSY.getMHa();
    double M2Z = mySUSY.getMz() * mySUSY.getMz();
    double tanb = mySUSY.getTanb();
   
    double tan2alpha = 2. * tanb / (1 - tanb * tanb) * (M2A + M2Z) / (M2A - M2Z); 
    double tana = (-1. + sqrt(1. + tan2alpha * tan2alpha)) / tan2alpha;
    double sina = tana / sqrt(1. + tana * tana);
    double cosa = 1. / sqrt(1 + tana * tana);
    
    complex Cosa(cosa,0.,false);
    complex Sina(sina,0.,false);
    complex i(0.,1.,false);
    

    if (S == 0) {
        return (Sina);
    }
    if (S == 1) {
        return (Cosa);
    }
    if (S == 2) {
        return (-i * mySUSY.getCosb());
    } else {
            
             throw std::runtime_error("Error in xuS(S) in SUSYMatching.cpp  "); 
        }
    
    return (EXIT_FAILURE);
    
}

gslpp::complex SUSYMatching::XRLS(int J, int I, int S){
    
    double v = mySUSY.v();
    double v1 = v * mySUSY.getCosb();
    double v2 = v * mySUSY.getSinb();
    double tanb = mySUSY.getTanb();
    double Y2ut = sqrt(2.) / v2 * mySUSYMQ(4);
    Y2ut *= Y2ut;
 
    if (J > I) {            
        return (mySUSYMQ(2 * J + 1) / (v1 * (1 + Eps_J(J) * tanb) *
                (1 + Eps_J(J) * tanb)) * Lambda0EpsY(J, I) * Y2ut *
                (xuS(S) - xdS(S) * tanb));
    } else
        if (J < I) {

        return (XLRS(I, J, S).conjugate());
    } else {
            
             throw std::runtime_error("Error in XRLS(J,I,S) in SUSYMatching.cpp  "); 
        }

    return (EXIT_FAILURE);
}

gslpp::complex SUSYMatching::XLRS(int J, int I, int S){
    
    double v = mySUSY.v();
    double v1 = v * mySUSY.getCosb();
    double v2 = v * mySUSY.getSinb();
    double tanb = mySUSY.getTanb();
    double Y2ut = sqrt(2.) / v2 * mySUSYMQ(4);
    Y2ut *= Y2ut;
    complex temp(0.,0.,false);
    temp = 1 + Eps_J(J) * tanb;
    complex rJI(0.,0.,false);
    

    rJI = ((1. + (Eps_J(J) + (Eps_J(I).conjugate() - Eps_J(J).conjugate()) *
            Lambda0EpsY(J, I) / Lambda0EpsY(I, J).conjugate()) * tanb) /
            (1 + Eps_J(I).conjugate() * tanb));   
    
    if (J > I) {    

        return (mySUSYMQ(2 * I + 1) / (v1 * temp.abs2()) * 
                Lambda0EpsY(I, J).conjugate() * Y2ut * rJI *
                (xuS(S).conjugate() - xdS(S).conjugate() * tanb));
    } else
        if (J < I) {

        return (XRLS(I, J, S).conjugate());

    } else {
           
             throw std::runtime_error("Error in XLRS(J,I,S) in SUSYMatching.cpp  "); 
        
        }
    
    return (EXIT_FAILURE);
    
}

////////////////////////////////////////////////////////////////////////////////
//// Charged Higgs contribution to Wilson coefficients of Delta F = 2 
//// in the Buras operator basis
//// Q_1 = (q \bar _L gamma_mu b_L) (q \bar _L gamma^mu b_L)
//// in the D - D \bar mixing q -> u and b -> c

gslpp::vector<complex> SUSYMatching::CdF2dHp(int b, int q, int Dmixingflag) {
    
    gslpp::matrix<complex> myCKM(3, 3, 0.);
    myCKM = mySUSY_CKM();
    gslpp::matrix<complex> CKMeff(3, 3, 0.);
    mySUSY.getCKM().getCKM(CKMeff);
    gslpp::vector<double> M2S(3,0.);
    gslpp::vector<double> MQuarks(6,0.);
    int i;
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        for (i = 0; i < 6; i++) {
            MQuarks(i) += mySUSYMQ(i);
        }
    } else if (Dmixingflag == 1) {
        for (i = 0; i < 3; i++) {   
            
            MQuarks(2 * i) += mySUSYMQ(2 * i + 1);
            MQuarks(2 * i + 1) += mySUSYMQ(2 * i);
        }
        myCKM.transpose();
    }
    
   
    //double M2Z = mySUSY.getMz() * mySUSY.getMz();

    M2S(0) = mySUSY.getMHl() * mySUSY.getMHl(); 
    M2S(1) = mySUSY.getMHh() * mySUSY.getMHh();
    M2S(2) = mySUSY.getMHa() * mySUSY.getMHa();

    gslpp::vector<complex> VCLO(8, 0.);
    complex CLO(0.,0.,false);
    
    
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
                        CLO += 1. / 4. * mySUSY.getGF() * M2W / (sqrt(2) * M_PI * M_PI) *
                                myCKM(I, q).conjugate() * myCKM(J, b) * PRLk(J, q, 0, D).conjugate() *
                                PRLk(I, b, 0, D) * D0N(M2W, M2H, M2I, M2J)

                                - 1. / (32. * M_PI * M_PI) * PRLk(I, q, 0, D).conjugate() *
                                PRLk(J, q, 0, D).conjugate() * PRLk(I, b, 0, D) * PRLk(J, b, 0, D) *
                                Dk(M2H, M2H, M2I, M2J, 2)

                                - 1. / (16. * M_PI * M_PI) * PRLk(I, q, 1, D).conjugate() *
                                PRLk(J, q, 0, D).conjugate() * PRLk(I, b, 0, D) * PRLk(J, b, 1, D) *
                                Dk(M2W, M2H, M2I, M2J, 2);
                }
            }
        } else if (O == 2) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {

                            CLO += -1. / (32. * M_PI * M_PI) * PLRk(I, q, l, D).conjugate() *
                                    PLRk(J, q, k, D).conjugate() * PRLk(I, b, k, D) * PRLk(J, b, l, D) *
                                    D0N(M2Hk(k), M2Hk(l), M2I, M2J);
                        }
                    }
                }
            }
        } else if (O == 4) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {


                        CLO += mySUSY.getGF() * M2W / (sqrt(2) * M_PI * M_PI) *
                                myCKM(I, q).conjugate() * myCKM(J, b) *
                                PLRk(J, q, k, D).conjugate() * PLRk(I, b, k, D) *
                                Dk(M2W, M2Hk(k), M2I, M2J, 2);

                        //std::cout << "CLO = " << CLO << std::endl;
                        
                        for (l = 0; l < 2; l++) {

                            CLO += -1. / (16. * M_PI * M_PI) * PLRk(I, q, l, D).conjugate() *
                                    PRLk(J, q, k, D).conjugate() * PRLk(I, b, k, D) *
                                    PLRk(J, b, l, D) * D0N(M2Hk(k), M2Hk(l), M2I, M2J);
                        }
                    }
                }
            }
         
            /** The double Penguin contributions are calulated only for the B and K mixing **/

            if (D != 0) {
                for (S = 0; S < 3; S++) {

                    CLO += -XRLS(q, b, S) * XLRS(q, b, S) / M2S(S);


                    //                complex temp(0.,0.,false);
                    //                temp +=  - XRLS(q, b, S) * XLRS(q, b, S) / M2S(S);

                }
            }
            
            /** end double Penguin contribution **/

        } else if (O == 5) {
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
        } else if (O == 6) {
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
        } else if (O == 7) {
            for (I = 0; I < 3; I++) {
                M2I = MQuarks(2 * I);
                M2I *= M2I;

                for (J = 0; J < 3; J++) {

                    M2J = MQuarks(2 * J);
                    M2J *= M2J;

                    for (k = 0; k < 2; k++) {
                        for (l = 0; l < 2; l++) {

                            CLO += -1. / (32. * M_PI * M_PI) * PRLk(I, q, l, D).conjugate() *
                                    PRLk(J, q, k, D).conjugate() * PLRk(I, b, k, D) * PLRk(J, b, l, D) *
                                    D0N(M2Hk(k), M2Hk(l), M2I, M2J);      
                        
                        }
                    }
                }
            }
        }
        
       
//        double CLOreal = CLO.real();
//        double CLOimm = CLO.imag();
        VCLO.assign(O - 1,CLO);
    }

    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Gluino contribution to Wilson coefficients of Delta F = 2 

gslpp::vector<complex> SUSYMatching::CdF2dgg(int b, int q, int Dmixingflag) {

    double Q = mySUSY.GetQ();

    gslpp::matrix<complex> myR(6, 6, 0.);
    gslpp::vector<double> myM2Squarks(6, 0.);
    
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
        myR = mySUSY.getRd();
    } else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsu2();  
        
        /* in the D mixing Rd -> Ru^* , b -> c , q -> u */
        
        myR = mySUSY.getRu().hconjugate().transpose();
        
    }
    
    complex CLO(0., 0., false);
    gslpp::vector<complex> VCLO(8, 0.);
    double Als = mySUSY.Als(Q);
    double Mg = mySUSY.getM3();
    double M2g = Mg*Mg;
    int h, k, O;


    

    for (O = 1; O < 9; O++) {
        
        CLO.assign(0., 0.,0);
        
        if (O == 1) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * myR(h, b).conjugate() * myR(k, b).conjugate()
                            * myR(h, q) * myR(k, q)*(1. / 9. * M2g *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) -
                            11. / 9. * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2));
                }
            }
        } else if (O == 2) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 17. / 18. * M2g * myR(h, b).conjugate()
                            * myR(k, b).conjugate() * myR(h, q + 3)
                            * myR(k, q + 3)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        } else if (O == 3) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += Als * Als * 1. / 6. * M2g * myR(h, b).conjugate()
                            * myR(k, b).conjugate() * myR(h, q + 3)
                            * myR(k, q + 3)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        } else if (O == 4) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 7. / 3. * M2g * myR(h, b).conjugate()
                            * myR(k, b + 3).conjugate() * myR(h, q)
                            * myR(k, q + 3)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) +
                            Als * Als * 2. / 9. *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2) *
                            myR(h, b).conjugate() * myR(k, b + 3).conjugate() *
                            (6. * myR(h, q) * myR(k, q + 3) + 11. * myR(k, q)
                            * myR(h, q + 3));
                }
            }
        } else if (O == 5) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 1. / 9. * M2g * myR(h, b).conjugate()
                            * myR(k, b + 3).conjugate() * myR(h, q)
                            * myR(k, q + 3)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) +
                            Als * Als * 10. / 9. *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2) *
                            myR(h, b).conjugate() * myR(k, b + 3).conjugate() *
                            (3. * myR(k, q) * myR(h, q + 3) - 2. * myR(h, q)
                            * myR(k, q + 3));
                }
            }
        } else if (O == 6) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * myR(h, b + 3).conjugate() * myR(k, b + 3).conjugate()
                            * myR(h, q + 3) * myR(k, q + 3)*(1. / 9. * M2g *
                            Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0) -
                            11. / 9. * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 2));
                }
            }
        } else if (O == 7) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += -Als * Als * 17. / 18. * M2g * myR(h, b).conjugate()
                            * myR(k, b + 3).conjugate() * myR(h, q)
                            * myR(k, q)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        } else if (O == 8) {
            for (h = 0; h < 6; h++) {
                for (k = 0; k < 6; k++) {

                    CLO += Als * Als * 1. / 6. * M2g * myR(h, b + 3).conjugate()
                            * myR(k, b + 3).conjugate() * myR(h, q)
                            * myR(k, q)
                            * Dk(myM2Squarks(h), myM2Squarks(k), M2g, M2g, 0);
                }
            }
        }

//        double CLOreal = CLO.real();
//        double CLOimm = CLO.imag();
        VCLO.assign(O - 1, CLO);
    }
    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Chargino contribution to Wilson coefficients of Delta F = 2

gslpp::vector<complex> SUSYMatching::CdF2dChiChi(int b, int q, int Dmixingflag) {

    gslpp::vector<double> myM2Squarks(6, 0.);
    gslpp::vector<double> MChi(2, 0.);
    complex CLO(0., 0., false);
    gslpp::vector<complex> VCLO(8, 0.);
    int i, j, h, k, O;

    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsu2();
    } else if (Dmixingflag == 1) {
        myM2Squarks = mySUSY.getMsd2();
    }
  
    MChi = mySUSY.getMch();
    
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
        } else if (O == 3) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) *
                                    VdUCR(q, k, j, 1, D).conjugate() *
                                    VdUCR(q, h, i, 1, D).conjugate() *
                                    VdUCL(b, k, i, D) * VdUCL(b, h, j, D) *
                                    MChi(i) * MChi(j) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i),
                                    MChi(j) * MChi(j), 0);
                        }
                    }
                }
            }
        } else if (O == 4) {
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
        } else if (O == 5) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (16. * M_PI * M_PI) * VdUCL(q, k, j, D).conjugate() *
                                    VdUCR(q, h, i, 1, D).conjugate() * VdUCL(b, h, j, D) *
                                    VdUCR(b, h, j, 1, D) * MChi(i) * MChi(j) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi(i) * MChi(i), MChi(j)
                                    * MChi(j), 0);
                        }
                    }
                }
            }
        } else if (O == 6) {
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
        } else if (O == 8) {
            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (32. * M_PI * M_PI) * VdUCL(q, k, j, D).conjugate() *
                                    VdUCL(q, h, i, D).conjugate() * VdUCR(b, k, i, 1, D) *
                                    VdUCR(b, h, j, 1, D) * MChi(i) * MChi(j) *
                                    Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi(i) * MChi(i), MChi(j) * MChi(j), 0);
                        }
                    }
                }
            }
        }

//        double CLOreal = CLO.real();
//        double CLOimm = CLO.imag();
        VCLO.assign(O - 1, CLO);

    }


    return (VCLO);
}

////////////////////////////////////////////////////////////////////////////////
//// Neutralino contribution to Wilson coefficients of Delta F = 2
/**** In the D mixing b -> c , q -> u ****/

gslpp::vector<complex> SUSYMatching::CdF2dChi0Chi0(int b, int q, int Dmixingflag) {  

    gslpp::vector<double> myM2Squarks(6, 0.);
    gslpp::vector<double> MChi0(4, 0.);
    gslpp::vector<complex> VCLO(8, 0.);
    complex CLO(0., 0., false);
    int i, j, h, k, O;

    MChi0 = mySUSY.getMneu();
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0
    
    

    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
    } else if (Dmixingflag == 1) {
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

                            CLO += VdDNL(b, k, i, 1, D) * VdDNL(q, k, j, 1, D).conjugate()
                                    * (-1. / (32. * M_PI * M_PI) * VdDNL(b, h, j, 1, D)
                                    * VdDNL(q, h, i, 1, D).conjugate() * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) - MChi0(i) * MChi0(j)
                                    / (64. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i), MChi0(j)
                                    * MChi0(j), 0) * VdDNL(b, h, i, 1, D)
                                    * VdDNL(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        } else if (O == 2) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += MChi0(i) * MChi0(j) / (32. * M_PI * M_PI) *
                                    Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNL(b, k, i, 1, D) * VdDNR(q, k, j, 1, D).conjugate() *
                                    VdDNL(b, h, i, 1, D) * VdDNR(q, h, j, 1, D).conjugate() ;
                        }
                    }
                }
            }
        } else if (O == 3) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -MChi0(i) * MChi0(j) / (32. * M_PI * M_PI) *
                                    Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 0) * VdDNL(b, k, i, 1, D) *
                                    VdDNR(q, k, j, 1, D).conjugate() *
                                    (VdDNL(b, h, j, 1, D) * VdDNR(q, h, i, 1, D).conjugate() -
                                    VdDNL(b, h, i, 1, D) * VdDNR(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        } else if (O == 4) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += 1. / (8. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) * VdDNR(b, k, i, 1, D) *
                                    VdDNL(q, k, j, 1, D).conjugate() * (
                                    VdDNL(b, h, j, 1, D) * VdDNR(q, h, i, 1, D).conjugate() +
                                    VdDNL(b, h, i, 1, D) * VdDNR(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        } else if (O == 5) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += -1. / (8. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) * VdDNL(q, k, j, 1, D).conjugate() *
                                    VdDNR(q, h, j, q, D).conjugate() * VdDNL(b, k, i, 1, D) *
                                    VdDNR(b, h, i, 1, D) - MChi0(i) * MChi0(j) / (16. * M_PI *
                                    M_PI) * Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNL(q, h, i, 1, D).conjugate() * VdDNR(q, k, j, 1, D).conjugate() *
                                    VdDNL(b, h, j, 1, D) * VdDNR(b, k, i, 1, D);
                        }
                    }
                }
            }
        } else if (O == 6) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += VdDNR(b, k, i, 1, D) * VdDNR(q, k, j, 1, D).conjugate()
                                    * (-1. / (32. * M_PI * M_PI) * VdDNR(b, h, j, 1, D)
                                    * VdDNR(q, h, i, 1, D).conjugate() * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i),
                                    MChi0(j) * MChi0(j), 2) - MChi0(i) * MChi0(j)
                                    / (64. * M_PI * M_PI) * Dk(myM2Squarks(k),
                                    myM2Squarks(h), MChi0(i) * MChi0(i), MChi0(j)
                                    * MChi0(j), 0) * VdDNR(b, h, i, 1, D)
                                    * VdDNR(q, h, j, 1, D).conjugate());
                        }
                    }
                }
            }
        } else if (O == 7) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                            CLO += MChi0(i) * MChi0(j) / (32. * M_PI * M_PI) *
                                    Dk(myM2Squarks(k), myM2Squarks(h),
                                    MChi0(i) * MChi0(i), MChi0(j) * MChi0(j), 0) *
                                    VdDNR(b, k, i, 1, D) * VdDNL(q, k, j, 1, D).conjugate() *
                                    VdDNR(b, h, i, 1, D) * VdDNL(q, h, j, 1, D).conjugate();
                        }
                    }
                }
            }

        } else if (O == 8) {

            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    for (h = 0; h < 6; h++) {
                        for (k = 0; k < 6; k++) {

                           CLO += -MChi0(i) * MChi0(j) / (32. * M_PI * M_PI) *
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

//        double CLOreal = CLO.real();
//        double CLOimm = CLO.imag();
        VCLO.assign(O - 1, CLO);

    }
    return (VCLO);

}

////////////////////////////////////////////////////////////////////////////////
//// Neutralino - Gluino contribution to Wilson coefficients of Delta F = 2

gslpp::vector<complex> SUSYMatching::CdF2dChi0g(int b, int q, int Dmixingflag) {   

    double Q = mySUSY.GetQ();

    gslpp::matrix<complex> myR(6, 6, 0.);

    gslpp::vector<double> myM2Squarks(6, 0.);
    gslpp::vector<double> MChi0(4, 0.);
    gslpp::vector<complex> VCLO(8, 0.);
    complex CLO(0., 0., false);
    int i, h, k, O;
    double Mg = mySUSY.getM3();
    double M2g = Mg*Mg;
    double Als = mySUSY.Als(Q);

    MChi0 = mySUSY.getMneu();
    
    // Set the D - Dbar mixing flag
    // in D - Dbar mixing the flag = 1 otherwise the flag = 0

    if (Dmixingflag == 0) {
        myM2Squarks = mySUSY.getMsd2();
        myR = mySUSY.getRd();
    } else if (Dmixingflag == 1) {
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

                CLO += -Als * 2. / (3. * 4. * M_PI ) * Dk(myM2Squarks(k), myM2Squarks(h),
                                MChi0(i) * MChi0(i), M2g, 2) * myR(k, q)
                                * myR(h, b).conjugate() * VdDNL(q, h, i, 1, D).conjugate() *
                                VdDNL(b, k, i, 1, D) - Als /( 4. * M_PI ) * MChi0(i) * Mg / 6. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i)
                                , M2g, 0) * (myR(h, q) * myR(k, q) * VdDNL(b, h, i, 1, D) *
                                VdDNL(b, k, i, 1, D) + myR(h, b).conjugate()
                                * myR(k, b).conjugate() * VdDNL(q, h, i, 1, D).conjugate() *
                                VdDNL(q, k, i, 1, D).conjugate());

                    }
                }
            }
        } else if (O == 2) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += Als / (4. * M_PI) * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (3 * myR(h, b).conjugate() *
                                myR(k, q + 3) * VdDNL(b, k, i, 1, D) * VdDNR(q, h, i, 1, D).conjugate()
                                + myR(k, b).conjugate() * myR(h, b).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate() +
                                myR(k, q + 3) * myR(h, q + 3) * VdDNL(b, k, i, 1, D) *
                                VdDNL(b, h, i, 1, D));

                    }
                }
            }
        } else if (O == 3) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als / (4. * M_PI) * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b).conjugate() * myR(k, q + 3) * VdDNL(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() - myR(k, b).conjugate() *
                                myR(h, b).conjugate() * VdDNR(q, k, i, 1, D).conjugate() *
                                VdDNR(q, h, i, 1, D).conjugate() - myR(k, q + 3) *
                                myR(h, q + 3) * VdDNL(b, k, i, 1, D) * VdDNL(b, h, i, 1, D));

                    }
                }
            }
        } else if (O == 4) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als / (4. * M_PI) * 2. / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 2) * (
                                myR(h, b).conjugate() * myR(k, q) * VdDNR(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() + myR(h, b + 3).conjugate() *
                                myR(k, q + 3) * VdDNL(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate()
                                - myR(h, b).conjugate() * myR(k, b + 3).conjugate() *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate()
                                - myR(h, q) * myR(k, q + 3) * VdDNL(b, k, i, 1, D)
                                * VdDNR(b, h, i, 1, D) - 3. * myR(k, b + 3).conjugate() *
                                myR(h, b).conjugate() * VdDNL(q, h, i, 1, D).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate() - 3. * myR(k, q + 3) *
                                myR(h, q) * VdDNL(b, h, i, 1, D) * VdDNR(b, k, i, 1, D)) +
                                Als / (4. * M_PI) * MChi0(i) * Mg *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b).conjugate() * myR(k, q + 3) *
                                VdDNR(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(h, b + 3).conjugate() * myR(h, q) *
                                VdDNL(b, k, i, 1, D) * VdDNR(q, h, i, 1, D).conjugate());

                    }
                }
            }
        } else if (O == 5) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += Als / (4. * M_PI) * 2. / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 2) * (
                                3. * myR(h, b).conjugate() * myR(k, q) * VdDNR(b, k, i, 1, D) *
                                VdDNR(q, h, i, 1, D).conjugate() + 3. * myR(h, b + 3).conjugate() *
                                myR(k, q + 3) * VdDNL(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate()
                                - 3. * myR(h, b).conjugate() * myR(k, b + 3).conjugate() *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNR(q, h, i, 1, D).conjugate()
                                - 3. * myR(h, q) * myR(k, q + 3) * VdDNL(b, k, i, 1, D)
                                * VdDNR(b, h, i, 1, D) - myR(k, b + 3).conjugate() *
                                myR(h, b).conjugate() * VdDNL(q, h, i, 1, D).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate() - myR(k, q + 3) *
                                myR(h, q) * VdDNL(b, h, i, 1, D) * VdDNR(b, k, i, 1, D)) -
                                Als / (4. * M_PI) * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b).conjugate() * myR(k, q + 3) *
                                VdDNR(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(h, b + 3).conjugate() * myR(h, q) *
                                VdDNL(b, k, i, 1, D) * VdDNR(q, h, i, 1, D).conjugate());

                    }
                }
            }
        } else if (O == 6) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als * 2. / (3. * 4. * M_PI) * Dk(myM2Squarks(k), myM2Squarks(h),
                                MChi0(i) * MChi0(i), M2g, 2) * myR(k, q + 3)
                                * myR(h, b + 3).conjugate() * VdDNR(q, h, i, 1, D).conjugate() *
                                VdDNR(b, k, i, 1, D) - Als / (4. * M_PI) * MChi0(i) * Mg / 6. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i) * MChi0(i)
                                , M2g, 0) * (myR(h, q + 3) * myR(k, q + 3) * VdDNR(b, h, i, 1, D) *
                                VdDNR(b, k, i, 1, D) + myR(h, b + 3).conjugate()
                                * myR(k, b + 3).conjugate() * VdDNR(q, h, i, 1, D).conjugate() *
                                VdDNR(q, k, i, 1, D).conjugate());

                    }
                }
            }
        } else if (O == 7) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += Als / (4. * M_PI) * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (3. * myR(h, b + 3).conjugate() *
                                myR(k, q) * VdDNR(b, k, i, 1, D) * VdDNL(q, h, i, 1, D).conjugate()
                                + myR(k, b + 3).conjugate() * myR(h, b + 3).conjugate() *
                                VdDNL(q, k, i, 1, D).conjugate() * VdDNL(q, h, i, 1, D).conjugate() +
                                myR(k, q) * myR(h, q) * VdDNR(b, k, i, 1, D) *
                                VdDNR(b, h, i, 1, D));

                    }
                }
            }
        } else if (O == 8) {
            for (i = 0; i < 4; i++) {
                for (h = 0; h < 6; h++) {
                    for (k = 0; k < 6; k++) {

                        CLO += -Als / (4. * M_PI) * MChi0(i) * Mg / 3. *
                                Dk(myM2Squarks(k), myM2Squarks(h), MChi0(i)
                                * MChi0(i), M2g, 0) * (
                                myR(h, b + 3).conjugate() * myR(k, q) * VdDNR(b, k, i, 1, D) *
                                VdDNL(q, h, i, 1, D).conjugate() - myR(k, b + 3).conjugate() *
                                myR(h, b + 3).conjugate() * VdDNL(q, k, i, 1, D).conjugate() *
                                VdDNL(q, h, i, 1, D).conjugate() - myR(k, q) *
                                myR(h, q) * VdDNR(b, k, i, 1, D) * VdDNR(b, h, i, 1, D));

                    }
                }
            }
        }
        
//        double CLOreal = CLO.real();
//        double CLOimm = CLO.imag();
        VCLO.assign(O - 1, CLO);
    }
    return (VCLO);
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////// 

/******************************************************************************/
const std::vector<WilsonCoefficient>& SUSYMatching::CMdbd2() {

    double Q = mySUSY.GetQ();
    int i;
    gslpp::vector<complex> CdF2dHpT(8, 0.);
    gslpp::vector<complex> CdF2dggT(8, 0.);
    gslpp::vector<complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<complex> CdF2dChigT(8, 0.);

    if (mySUSY.IsFh()) CdF2dHpT = CdF2dHp(0, 2, 0);
    if (mySUSY.IsFChi()) CdF2dChiChiT = CdF2dChiChi(0, 2, 0);
    if (mySUSY.IsFg()) CdF2dggT = CdF2dgg(0, 2, 0);
    if (mySUSY.IsFChi0()) CdF2dChi0Chi0T = CdF2dChi0Chi0(0, 2, 0);
    if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) CdF2dChigT = CdF2dChi0g(0, 2, 0);
    
    vmdbd2 = StandardModelMatching::CMdbd2();   
                                                
    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdbd2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdbd2Hp.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbd2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdbd2.push_back(mcdbd2Hp);
            }
            break;
    }
    switch (mcdbd2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {

                mcdbd2gg.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbd2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdbd2.push_back(mcdbd2gg);
            }
            break;
    }
    switch (mcdbd2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdbd2ChiChi.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbd2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdbd2.push_back(mcdbd2ChiChi);
            }
            break;
    }
    switch (mcdbd2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdbd2Chi0Chi0.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbd2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0Chi0);
            }
            break;
    }
    switch (mcdbd2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdbd2Chi0g.setMu(Q);
                for (i = 0; i < 5; i++) {
                mcdbd2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0g);
            }
            break;

    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdbd2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdbd2HpT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbd2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2HpT);
            }
            break;
    }
    switch (mcdbd2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {

                mcdbd2ggT.setMu(Q);
                for(i = 0; i < 3 ; i++){
                mcdbd2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2ggT);
            }
            break;
    }
    switch (mcdbd2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdbd2ChiChiT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbd2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2ChiChiT);
            }
            break;
    }
    switch (mcdbd2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {

                mcdbd2Chi0Chi0T.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbd2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0Chi0T);
            }
            break;
    }
    switch (mcdbd2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdbd2Chi0gT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbd2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdbd2.push_back(mcdbd2Chi0gT);
            }
            break;

    }
    
    return (vmdbd2);
}

/******************************************************************************/

const std::vector<WilsonCoefficient>& SUSYMatching::CMdbs2() {

    double Q = mySUSY.GetQ();
    int i;
    gslpp::vector<complex> CdF2dHpT(8, 0.);
    gslpp::vector<complex> CdF2dggT(8, 0.);
    gslpp::vector<complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFh()) CdF2dHpT = CdF2dHp(1, 2, 0);
    if (mySUSY.IsFChi()) CdF2dChiChiT = CdF2dChiChi(1, 2, 0);
    if (mySUSY.IsFg()) CdF2dggT = CdF2dgg(1, 2, 0);  
    if (mySUSY.IsFChi0()) CdF2dChi0Chi0T = CdF2dChi0Chi0(1, 2, 0);
    if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) CdF2dChigT = CdF2dChi0g(1, 2, 0);
    
    vmdbs2 = StandardModelMatching::CMdbs2();

    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdbs2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdbs2Hp.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbs2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdbs2.push_back(mcdbs2Hp);
            }
            break;
    }
    switch (mcdbs2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdbs2gg.setMu(Q);
                for(i = 0; i < 5 ; i++){
                mcdbs2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdbs2.push_back(mcdbs2gg);
            }
            break;
    }
    switch (mcdbs2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdbs2ChiChi.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbs2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdbs2.push_back(mcdbs2ChiChi);
            }
            break;
    }
    switch (mcdbs2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdbs2Chi0Chi0.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbs2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0Chi0);
            }
            break;
    }
    switch (mcdbs2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdbs2Chi0g.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdbs2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0g);
            }
            break;
    }
  
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdbs2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdbs2HpT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbs2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2HpT);
            }
            break;
    }
    switch (mcdbs2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdbs2ggT.setMu(Q);
                for(i = 0; i < 3 ; i++){
                mcdbs2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2ggT);
            }
            break;
    }
    switch (mcdbs2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdbs2ChiChiT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbs2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2ChiChiT);
            }
            break;
    }
    switch (mcdbs2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdbs2Chi0Chi0T.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbs2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0Chi0T);
            }
            break;
    }
    switch (mcdbs2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdbs2Chi0gT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdbs2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdbs2.push_back(mcdbs2Chi0gT);
            }
            break;
    }
    
    
    return (vmdbs2);
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<WilsonCoefficient>& SUSYMatching::CMdk2() {



    double Q = mySUSY.GetQ();
    int i;
    gslpp::vector<complex> CdF2dHpT(8, 0.);
    gslpp::vector<complex> CdF2dggT(8, 0.);
    gslpp::vector<complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFh()) CdF2dHpT = CdF2dHp(0, 1, 0);
    if (mySUSY.IsFChi()) CdF2dChiChiT = CdF2dChiChi(0, 1, 0);
    if (mySUSY.IsFg()) CdF2dggT = CdF2dgg(0, 1, 0);
    if (mySUSY.IsFChi0()) CdF2dChi0Chi0T = CdF2dChi0Chi0(0, 1, 0);
    if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) CdF2dChigT = CdF2dChi0g(0, 1, 0);
    
    vmdk2 = StandardModelMatching::CMdk2();

    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
    
    switch (mcdk2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdk2Hp.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdk2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdk2.push_back(mcdk2Hp);
            }
            break;
    }
    switch (mcdk2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdk2gg.setMu(Q);
                for(i = 0; i < 5 ; i++){
                mcdk2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdk2.push_back(mcdk2gg);
            }
            break;
    }
    switch (mcdk2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdk2ChiChi.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdk2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdk2.push_back(mcdk2ChiChi);
            }
            break;
    }
    switch (mcdk2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdk2Chi0Chi0.setMu(Q);    
                for (i = 0; i < 5; i++) {
                    mcdk2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdk2.push_back(mcdk2Chi0Chi0);
            }
            break;
    }
    switch (mcdk2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdk2Chi0g.setMu(Q);    
                for (i = 0; i < 5; i++) {
                    mcdk2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdk2.push_back(mcdk2Chi0g);
            }
            break;
    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdk2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdk2HpT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdk2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2HpT);
            }
            break;
    }
    switch (mcdk2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdk2ggT.setMu(Q);
                for(i = 0; i < 3 ; i++){
                mcdk2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2ggT);
            }
            break;
    }
    switch (mcdk2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdk2ChiChiT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdk2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2ChiChiT);
            }
            break;
    }
    switch (mcdk2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdk2Chi0Chi0T.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdk2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdk2.push_back(mcdk2Chi0Chi0T);
            }
            break;
    }
    switch (mcdk2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdk2Chi0gT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdk2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdk2.push_back(mcdk2Chi0gT);
            }
            break;
    }
    
    
    return (vmdk2);
}


////////////////////////////////////////////////////////////////////////////////

const std::vector<WilsonCoefficient>& SUSYMatching::CMdd2(){

    double Q = mySUSY.GetQ();
    int i;
    gslpp::vector<complex> CdF2dHpT(8, 0.);
    gslpp::vector<complex> CdF2dggT(8, 0.);
    gslpp::vector<complex> CdF2dChiChiT(8, 0.);
    gslpp::vector<complex> CdF2dChi0Chi0T(8, 0.);
    gslpp::vector<complex> CdF2dChigT(8, 0.);
    
    if (mySUSY.IsFh()) CdF2dHpT = CdF2dHp(1, 0, 1);
    if (mySUSY.IsFg()) CdF2dggT = CdF2dgg(1, 0, 1);
    if (mySUSY.IsFChi()) CdF2dChiChiT = CdF2dChiChi(1, 0, 1); 
    if (mySUSY.IsFChi0()) CdF2dChi0Chi0T = CdF2dChi0Chi0(1, 0, 1);
    if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) CdF2dChigT = CdF2dChi0g(1, 0, 1);
    
    vmdd2 = StandardModelMatching::CMdd2();
    
    /** Wilson coefficients of operator Q_1,2,3,4,5 **/
   
     switch (mcdd2Hp.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdd2Hp.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdd2Hp.setCoeff(i, CdF2dHpT(i), LO);
                }
                vmdd2.push_back(mcdd2Hp);
            }
            break;
    }
    switch (mcdd2gg.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdd2gg.setMu(Q);
                for(i = 0; i < 5 ; i++){
                mcdd2gg.setCoeff(i, CdF2dggT(i), LO);
                }
                vmdd2.push_back(mcdd2gg);
            }
            break;
    }
    switch (mcdd2ChiChi.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdd2ChiChi.setMu(Q);
                for (i = 0; i < 5; i++) {
                    mcdd2ChiChi.setCoeff(i, CdF2dChiChiT(i), LO);
                }
                vmdd2.push_back(mcdd2ChiChi);
            }
            break;
    }
    switch (mcdd2Chi0Chi0.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdd2Chi0Chi0.setMu(Q);    
                for (i = 0; i < 5; i++) {
                    mcdd2Chi0Chi0.setCoeff(i, CdF2dChi0Chi0T(i), LO);
                }
                vmdd2.push_back(mcdd2Chi0Chi0);
            }
            break;
    }
    switch (mcdd2Chi0g.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdd2Chi0g.setMu(Q);    
                for (i = 0; i < 5; i++) {
                    mcdd2Chi0g.setCoeff(i, CdF2dChigT(i), LO);
                }
                vmdd2.push_back(mcdd2Chi0g);
            }
            break;
    }
    
    /** Wilson coefficients of operator Q_1,2,3 tilde **/
    
    switch (mcdd2HpT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFh()) {
                mcdd2HpT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdd2HpT.setCoeff(i, CdF2dHpT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2HpT);
            }
            break;
    }
    switch (mcdd2ggT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFg()) {
                mcdd2ggT.setMu(Q);
                for(i = 0; i < 3 ; i++){
                mcdd2ggT.setCoeff(i, CdF2dggT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2ggT);
            }
            break;
    }
    switch (mcdd2ChiChiT.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi()) {
                mcdd2ChiChiT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdd2ChiChiT.setCoeff(i, CdF2dChiChiT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2ChiChiT);
            }
            break;
    }
    switch (mcdd2Chi0Chi0T.getOrder()) {
        case NLO:
        case LO:
            if (mySUSY.IsFChi0()) {
                mcdd2Chi0Chi0T.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdd2Chi0Chi0T.setCoeff(i, CdF2dChi0Chi0T(i + 5), LO);
                }
                vmdd2.push_back(mcdd2Chi0Chi0T);
            }
            break;
    }
    switch (mcdd2Chi0gT.getOrder()) {
        case NLO:
        case LO:
            if ((mySUSY.IsFg()) || (mySUSY.IsFChi0())) {
                mcdd2Chi0gT.setMu(Q);
                for (i = 0; i < 3; i++) {
                    mcdd2Chi0gT.setCoeff(i, CdF2dChigT(i + 5), LO);
                }
                vmdd2.push_back(mcdd2Chi0gT);
            }
            break;
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
   
    gslpp::matrix<complex> myCKM(3, 3, 0.);
    myCKM = mySUSY_CKM();
    complex YuJ(0., 0., false);
    double v = mySUSY.v();
    double v2 = v * mySUSY.getSinb();

    YuJ = sqrt(2.) / v2 * mySUSYMQ(2 * J);
    
    
    return (-DeltaFHL(J,I) / (mySUSY.getTanb() * myCKM(J,I) * YuJ));
}
double SUSYMatching::F7k(double x, int k) {

    if (k == 1) {

        if (fabs(x - 1.) < SUSYLEPS) return (-5. / 48.);

        return ((x * (7 - 5 * x - 8 * x * x)) / (24. * (-1 + x)*(-1 + x)*(-1 + x)) +
                x * x * (-2 + 3 * x) * log(x) / (4. * (-1 + x)*(-1 + x)*(-1 + x)*(-1 + x)));

    } else if (k == 2) {

        if (fabs(x - 1.) < SUSYLEPS) return (-7. / 36.);

        return (((3 - 5 * x) * x) / (12. * (-1 + x)*(-1 + x)) +
                (x * (-2 + 3 * x) * log(x)) / (6. * (-1 + x)*(-1 + x)*(-1 + x)));

    }
    else throw std::runtime_error("Error in F7k "); 

}

gslpp::vector <complex> SUSYMatching::CalcC7(int b, int q) {

    gslpp::vector<complex> VCLO(5, 0.);
    complex CLO(0., 0., false);
    double m2top = mySUSYMQ(4) * mySUSYMQ(4);
    double M2Hp = mySUSY.getMHp() * mySUSY.getMHp();
    double M2W = mySUSY.Mw_tree() * mySUSY.Mw_tree();
    gslpp::matrix<complex> myCKM(3, 3, 0.);
    mySUSY.getCKM().getCKM(myCKM);
    double gW = sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    gslpp::vector<double> MChi(2, 0.);
    MChi = mySUSY.getMch();
    gslpp::vector<double> myMU2Squarks(6, 0.);
    myMU2Squarks = mySUSY.getMsu2();
    int j, k;


    VCLO.assign(0, F7k(m2top / M2W, 0) + (Eps_J(2) - EpsPrime(2, 2)) / (1 + Eps_J(2)
            * mySUSY.getTanb()) * F7k(m2top / M2W, 2));

    VCLO.assign(1, 1. / (3. * mySUSY.getTanb() * mySUSY.getTanb()) * F7k(m2top / M2Hp, 1) +
            F7k(m2top / M2Hp, 2) - (EpsPrime(2, 1) + Eps_J(2)) / (1. + Eps_J(2) * mySUSY.getTanb()) *
            mySUSY.getTanb() * F7k(m2top / M2Hp, 2));

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


/*** FUNZIONE - TEST  ***/
    
void SUSYMatching::Test() {
    
    gslpp::matrix<complex> myRu = mySUSY.getRu();
    gslpp::matrix<double> myRur(6, 6, 0.);
    gslpp::matrix<double> myRui(6, 6, 0.);
    gslpp::matrix<complex> myRd = mySUSY.getRd();
    gslpp::matrix<double> myRdr(6, 6, 0.);
    gslpp::matrix<double> myRdi(6, 6, 0.);
    
     
      
    StandardModelMatching::CMdbd2();
    
    
    int i, j;
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {

            double tempRur = myRu(i, j).real();
            double tempRui = myRu(i, j).imag();
            double tempRdr = myRd(i, j).real();
            double tempRdi = myRd(i, j).imag();
            myRur(i, j) += tempRur ;
            myRui(i, j) += tempRui ;
            myRdr(i, j) += tempRdr ;
            myRdi(i, j) += tempRdi ;

            
        }
    }
    std::cout << " " << std::endl;
    std::cout << "CdF2Hp B_d = "  << CdF2dHp(0, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2ChiChi B_d = "  << CdF2dChiChi(0, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2gg B_d = "  << CdF2dgg(0, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0Chi0 B_d = "  << CdF2dChi0Chi0(0, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0g B_d = "  << CdF2dChi0g(0, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    StandardModelMatching::CMdbs2();
    std::cout << " " << std::endl;
    std::cout << "CdF2Hp B_s = "  << CdF2dHp(1, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2ChiChi B_s = "  << CdF2dChiChi(1, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2gg B_s = "  << CdF2dgg(1, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0Chi0 B_s = "  << CdF2dChi0Chi0(1, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0g B_s = "  << CdF2dChi0g(1, 2, 0) << std::endl;
    std::cout << " " << std::endl;
    
    std::cout << "EPSILON K C_i " << std::endl;
    std::cout << "CdF2Hp K = "  << CdF2dHp(0, 1, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2ChiChi K = "  << CdF2dChiChi(0, 1, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2gg K = "  << CdF2dgg(0, 1, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0Chi0 K = "  << CdF2dChi0Chi0(0, 1, 0) << std::endl;
    std::cout << " " << std::endl;
    std::cout << "CdF2Chi0g K = "  << CdF2dChi0g(0, 1, 0) << std::endl;
    std::cout << " " << std::endl;
    
    
    
    
    
    std::cout << "Rur = " << myRur << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Rui = " << myRui << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Rdr = " << myRdr << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Rdi = " << myRdi << std::endl;
    std::cout << " " << std::endl;
    gslpp::vector<double> myM2US(6,0.),myM2DS(6.,0);
    myM2US = mySUSY.getMsu2();
    myM2DS = mySUSY.getMsd2();
    for(i = 0;i < 6; i++) {
        myM2US(i) /= sqrt(myM2US(i));
        myM2DS(i) /= sqrt(myM2DS(i));
    }
    std::cout << "MUSquarks = " << myM2US << std::endl;
    std::cout << " " << std::endl;
    std::cout << "MDSquarks = " << myM2DS << std::endl;
    std::cout << " " << std::endl;
    std::cout << "M2Chi = " << mySUSY.getMch() << std::endl;
    std::cout << " " << std::endl;
    std::cout << "M2Chi0 = " << mySUSY.getMneu() << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Vckm_input = " << mySUSY.getVCKM() << std::endl;
    std::cout << "mtop[Q] = " << mySUSYMQ(4) << std::endl;
    std::cout << "mup[Q] = " << mySUSYMQ(0) << std::endl;
    std::cout << "mcharm[Q] = " << mySUSYMQ(2) << std::endl;

    gslpp::matrix<complex> Temp(6,6,0.);
    gslpp::vector<double> MTemp = mySUSY.getMsu2();
    
    for(i = 0; i < 6; i++) 
        Temp.assign(i,i,MTemp(i));
    
    std::cout << "M2utilde a mano = " <<
            mySUSY.getRu().hconjugate() * Temp * mySUSY.getRu() << std::endl;
   
    
    std::cout << " CKM_eff CKM_eff^{dag} = " << mySUSY.getVCKM() * mySUSY.getVCKM().hconjugate() << std::endl;
    std::cout << " CKM_eff = " << mySUSY_CKM() << std::endl;
    std::cout << " CKM = " << mySUSY.getVCKM() << std::endl;
    
    std::cout << "End Test " << std::endl;
    
    
}
    
/**************************************/









const std::vector<WilsonCoefficient>& SUSYMatching::CMbsg(){


  
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
