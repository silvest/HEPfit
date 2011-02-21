/* 
 * File:   ZFitterTest-2.cpp
 * Author: mishima
 *
 * Created on Feb 21, 2011, 5:00:52 AM
 */

/*
 *  Test for Table 6.1 in hep-ph/0507146
 *
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "../src/ZFitter.h"


// GeV^{-2} --> nb
const double GeVminus2_to_nb = pow(10.0, -6.0)
                               / pow(10.0, -28.0)
                               / pow(299792458.0, -2.0)
                               / pow(6.58211899*pow(10.0,-22.0), -2.0)
                               * pow(10.0, 9.0);


int main(int argc, char** argv) {
 
    /* set parameters */
    double mZ = 91.1875;
    double mt = 175.0;
    double mh = 150.0;
    double alphas = 0.118;
    double dalpha5h = 0.02758;
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter */
    ZFitter ZF(mZ, mt, mh, alphas, dalpha5h, vtb,
               mu_constituent, md_constituent);

    int indexFermion;

    /* sqrt(s) = mZ */
    double sqrt_s = ZF.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* print input parameters */
    ZF.printInputs();

    /* set flags */
    // default: AMT4=4, ALEM=3
    int AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
        BOXD = 1, CONV = 1, FINR = 1, FOT2 = 3, GAMS = 1,
        DIAG = 1, INTF = 1, BARB = 2, PART = 0, POWR = 1,
        PRNT = 0, ALEM = 2, QCDC = 3, VPOL = 1, WEAK = 1,
        FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
        CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
        ISPP = 2, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
        IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
        ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
        DSWW = 0;
    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    ZF.setAllFlags(flags, 1);

    /* set cuts */
    int ICUT[12];
    double ACOL[12], EMIN[12], S_PR[12], ANG0[12], ANG1[12], SPP[12];
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        ICUT[indexFermion] = -1;
        ACOL[indexFermion] = 0.0;
        EMIN[indexFermion] = 0.0;
        S_PR[indexFermion] = 0.01*s;
        ANG0[indexFermion] = 0.0;
        ANG1[indexFermion] = 180.0;
        SPP[indexFermion]  = 0.25*s;
    }
    ICUT[11] = 0;
    ACOL[11] = 15.0;
    EMIN[11] = 10.0;
    S_PR[11] = 0.0;
    ANG0[11] = 35.0;
    ANG1[11] = 145.0;
    SPP[11]  = 0.0;
    ZF.setAllCuts(ICUT, ACOL, EMIN, S_PR, ANG0, ANG1, SPP, 1);
    std::cout << std::endl;

    /* calculate EW common blocks */
    ZF.calcCommonBlocks();
    std::cout << std::endl;

    std::cout << "###  see Table 6.1 in hep-ph/0507146  ###" << std::endl;
    double mW, Gamma_W, sw2;
    double Gamma_inv, Gamma_had, Gamma_total;
    double sigma0_e, sigma0_mu, sigma0_tau, sigma0_had;
    double R0_e, R0_mu, R0_tau, R0_b, R0_c, R0_s;
    double A_e, A_mu, A_tau, A_b, A_c, A_s;
    double AFB0_e, AFB0_mu, AFB0_tau, AFB0_b, AFB0_c, AFB0_s;
    double s2teff_e, s2teff_mu, s2teff_tau, s2teff_b, s2teff_c, s2teff_s;
    ZF.calcPO(&mW, &Gamma_W, &sw2,
              &Gamma_inv, &Gamma_had, &Gamma_total,
              &sigma0_e, &sigma0_mu, &sigma0_tau, &sigma0_had,
              &R0_e, &R0_mu, &R0_tau, &R0_b, &R0_c, &R0_s,
              &A_e, &A_mu, &A_tau, &A_b, &A_c, &A_s,
              &AFB0_e, &AFB0_mu, &AFB0_tau, &AFB0_b, &AFB0_c, &AFB0_s,
              &s2teff_e, &s2teff_mu, &s2teff_tau, &s2teff_b, &s2teff_c, &s2teff_s);
    std::cout << std::setw(15) << "m_W [GeV]" << std::setw(13)<< mW << std::endl
              << std::setw(15) << "Gamma_W [GeV]" << std::setw(13) << Gamma_W << std::endl
              << std::setw(15) << "sin^2(th_W)" << std::setw(13) << sw2 << std::endl
              << std::setw(15) << "Gamma_inv [GeV]" << std::setw(13) << Gamma_inv << std::endl
              << std::setw(15) << "Gamma_had [GeV]" << std::setw(13) << Gamma_had << std::endl
              << std::setw(15) << "Gamma_Z [GeV]" << std::setw(13) << Gamma_total << std::endl
              << std::setw(15) << "sigma0_e [nb]" << std::setw(13)
              << sigma0_e*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_mu [nb]" << std::setw(13)
              << sigma0_mu*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_tau [nb]" << std::setw(13)
              << sigma0_tau*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_had [nb]" << std::setw(13)
              << sigma0_had*GeVminus2_to_nb << std::endl
              << std::setw(15) << "R0_e" << std::setw(13) << R0_e << std::endl
              << std::setw(15) << "R0_mu" << std::setw(13) << R0_mu << std::endl
              << std::setw(15) << "R0_tau" << std::setw(13) << R0_tau << std::endl
              << std::setw(15) << "R0_b" << std::setw(13) << R0_b << std::endl
              << std::setw(15) << "R0_c" << std::setw(13) << R0_c << std::endl
              << std::setw(15) << "R0_s" << std::setw(13) << R0_s << std::endl
              << std::setw(15) << "A_e" << std::setw(13) << A_e << std::endl
              << std::setw(15) << "A_mu" << std::setw(13) << A_mu << std::endl
              << std::setw(15) << "A_tau" << std::setw(13) << A_tau << std::endl
              << std::setw(15) << "A_b" << std::setw(13) << A_b << std::endl
              << std::setw(15) << "A_c" << std::setw(13) << A_c << std::endl
              << std::setw(15) << "A_s" << std::setw(13) << A_s << std::endl
              << std::setw(15) << "AFB0_e" << std::setw(13) << AFB0_e << std::endl
              << std::setw(15) << "AFB0_mu" << std::setw(13) << AFB0_mu << std::endl
              << std::setw(15) << "AFB0_tau" << std::setw(13) << AFB0_tau << std::endl
              << std::setw(15) << "AFB0_b" << std::setw(13) << AFB0_b << std::endl
              << std::setw(15) << "AFB0_c" << std::setw(13) << AFB0_c << std::endl
              << std::setw(15) << "AFB0_s" << std::setw(13) << AFB0_s << std::endl
              << std::setw(15) << "sin^2(teff_e)" << std::setw(13) << s2teff_e << std::endl
              << std::setw(15) << "sin^2(teff_mu)" << std::setw(13) << s2teff_mu << std::endl
              << std::setw(15) << "sin^2(teff_tau)" << std::setw(13) << s2teff_tau << std::endl
              << std::setw(15) << "sin^2(teff_b)" << std::setw(13) << s2teff_b << std::endl
              << std::setw(15) << "sin^2(teff_c)" << std::setw(13) << s2teff_c << std::endl
              << std::setw(15) << "sin^2(teff_s)" << std::setw(13) << s2teff_s << std::endl
              << std::endl;
    
    return (EXIT_SUCCESS);
}

