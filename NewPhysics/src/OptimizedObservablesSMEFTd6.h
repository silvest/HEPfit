/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OPTIMIZEDOBSERVABLESSMEFTD6_H
#define OPTIMIZEDOBSERVABLESSMEFTD6_H

class StandardModel;
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class eeWW : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    eeWW(const StandardModel& SM_i);
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    
    
protected:
    
    void updateParameters();
    
    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double Gamma_H, Gamma_W, Gamma_Z, Gamma_T;
    double ymtau, ymt, ymb;
    double aS, Gf, aEWM1;
    double MH, MZ, MTA, MT, MB;
    double conjg__CKM3x3, CKM3x3, conjg__CKM1x1;
    double MZ2, MZ4, MH2, aEW, MW, sqrt__aEW, ee, MW2;
    double sw2, cw, sqrt__sw2, sw, g1, gw, vev, vev2, lam;
    double yb, yt, ytau, muH, ee2, cw2;
    std::complex<double> complexi, I1x33, I2x33, I3x33, I4x33;
    // Model parameters dependent on aS
    double sqrt__aS, G, G2; 
    // Model couplings independent of aS
    std::complex<double> GC_3, GC_51, GC_53, GC_59, GC_100;
    
private:
    
    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual std::string name() const 
    {
        return "e+ e- > mu+ vm u~ d (sm)";
    }

    virtual int code() const 
    {
        return 1;
    }

    const std::vector<double> & getMasses() const 
    {
        return mME;
    }

    // Get and set momenta for matrix element evaluation
    std::vector < double * > getMomenta()
    {
        return p;
    }
    
    void setMomenta(std::vector < double * > & momenta)
    {
        p = momenta;
    }
    
    void setInitial(int inid1, int inid2)
    {
        id1 = inid1; 
        id2 = inid2;
    }

    // Get matrix element vector
    const double * getMatrixElements() const 
    {
        return matrix_element;
    }

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 6; 
    static const int nprocesses = 1; 
    static const int nwavefuncs = 11; 
    static const int namplitudes = 3;
    
    int perm[nexternal];

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    
    std::complex<double> w[nwavefuncs][18]; 
    std::complex<double> amp[namplitudes]; 
    
    double matrix_1_epem_wpwm_wp_mupvm_wm_uxd(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // vector with external particle masses
    std::vector<double> mME; 

    // vector with momenta (to be changed each event)
    std::vector <double * > p; 
    
    // Initial particle ids
    int id1, id2;
    
    double Sgn(double e, double f);

    void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double> fo[6]);

    void sxxxxx(double p[4], int nss, std::complex<double> sc[3]);

    void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double> fi[6]);

    void txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double> fi[18]);

    void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double> v[6]);

    void VVV1_0(std::complex<double> V1[], std::complex<double> V2[],
            std::complex<double> V3[], std::complex<double> COUP, std::complex<double>& vertex);

    void FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);

    void FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);
    void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
            W3, std::complex<double> V3[]);

    void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, double M3, double W3, std::complex<double> V3[]);

    void FFV2_1(std::complex<double> F2[], std::complex<double> V3[],
            std::complex<double> COUP, double M1, double W1, std::complex<double> F1[]);

    void FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
            & vertex);

};

#endif /* OPTIMIZEDOBSERVABLESSMEFTD6_H */

