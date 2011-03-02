/* 
 * File:   EWphysicsTest.cpp
 * Author: mishima
 *
 * Created on Mar 1, 2011, 4:29:34 AM
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ZFitter.h"
#include "EWphysics.h"


int main(int argc, char** argv) {

    /* set parameters */
    double mZ_i = 91.1875;
    double mHl_i = 150.0;
    double alsMz_i = 0.118;
    double dalpha5h_i = 0.02758;
    double GF_i = 1.16637e-5;
    double ale_i = 1.0/137.0359895;
    //double aleMz_i = ale_i/(1.0 - 0.031419 - dalpha5h_i + 0.00007);
    double aleMz_i = 0.00775497; /* used in ZFitterTest-2 */
    double mu_i = 0.03;
    double md_i = 0.03;
    double mc_i = 1.25;
    double ms_i = 0.110;
    double mt_i = 175.0;
    double mb_i = 4.2;
    double me_i = 0.000510999;
    double mmu_i = 0.105658;
    double mtau_i = 1.77705;
    //
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter with the above parameters (essential)
     * set flags and cuts to be their defalut values */
    ZFitter ZF(mZ_i, mt_i, mHl_i, alsMz_i, dalpha5h_i, 
               vtb, mu_constituent, md_constituent);

    /* sqrt(s) = mZ */
    double sqrt_s = ZF.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* set flags (not necessary if using the default flags) */
    // default: AMT4=4, ALEM=3
    ZF.flag("AMT4", 6);
    ZF.flag("ALEM", 2);
 
    /* calculate EW common blocks (necessary before computing observables) */
    ZF.calcCommonBlocks();
    std::cout << std::endl;

    std::cout << "#######################  ZFitter class  "
              << "#######################" << std::endl << std::endl;

    /* print intermediate results */
    ZF.printIntermediateResults();
    std::cout << std::endl;

    /* calcualte and print EW precision obserbables */
    ZF.printPO();

    gslpp::complex gZf_i[10];
    gslpp::complex rhoZf_i[10];
    double Delta_r_i;

    int indexFermion;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        gZf_i[indexFermion].real() = ZF.getCommonARVEFZ(indexFermion);
        gZf_i[indexFermion].imag() = ZF.getCommonAIVEFZ(indexFermion);
        rhoZf_i[indexFermion].real() = ZF.getCommonAROTFZ(indexFermion);
        rhoZf_i[indexFermion].imag() = ZF.getCommonAIROFZ(indexFermion);
    }
    Delta_r_i = ZF.Delta_r();

    EWphysics EWP(gZf_i, rhoZf_i, Delta_r_i,
                  mu_i, md_i, mc_i, ms_i, mt_i, mb_i, me_i, mmu_i, mtau_i,
                  mZ_i, mHl_i, alsMz_i, GF_i, ale_i, aleMz_i);


    std::cout << "##############  This test program (EWphysics class)  "
              << "###############" << std::endl << std::endl;

    std::cout << "--- Inputs from the ZFitter class to the EWphysics class ---"
              << std::endl;
    std::cout << std::setw(12) << "Re[g_Z^f]"
              << std::setw(13) << "Im[g_Z^f]"
              << std::setw(13) << "Re[rho_Z^f]"
              << std::setw(13) << "Im[rho_Z^f]"
              << std::endl;
    for (indexFermion = 0; indexFermion < 10; indexFermion++) {
        std::cout << std::setw(12) << EWP.getGZf(indexFermion).real()
                  << std::setw(13) << EWP.getGZf(indexFermion).imag()
                  << std::setw(12) << EWP.getRhoZf(indexFermion).real()
                  << std::setw(14) << EWP.getRhoZf(indexFermion).imag()
                  << std::endl;
    }
    std::cout << std::endl;
    std::cout << "  Delta r = " << EWP.getDelta_r() << std::endl
              << std::endl;

    std::cout << "--- Outputs from EWphysics class ---"
              << std::endl;

    // GeV^{-2} --> nb
    const double GeVminus2_to_nb = pow(10.0, -6.0)
                                   / pow(10.0, -28.0)
                                   / pow(299792458.0, -2.0)
                                   / pow(6.58211899 * pow(10.0, -22.0), -2.0)
                                   * pow(10.0, 9.0);

    std::cout << std::setw(15) << "m_W [GeV]" << std::setw(13)
              << EWP.mW() << std::endl
              << std::setw(15) << "Gamma_W [GeV]" << std::setw(13)
              << EWP.Gamma_W() << std::endl
              << std::setw(15) << "sin^2(th_W)" << std::setw(13)
              << EWP.sw2() << std::endl
              << std::setw(15) << "sin^2(teff_e)" << std::setw(13)
              << EWP.s2teff_f("e") << std::endl
              << std::setw(15) << "sin^2(teff_mu)" << std::setw(13)
              << EWP.s2teff_f("mu") << std::endl
              << std::setw(15) << "sin^2(teff_tau)" << std::setw(13)
              << EWP.s2teff_f("tau") << std::endl
              << std::setw(15) << "sin^2(teff_b)" << std::setw(13)
              << EWP.s2teff_f("b") << std::endl
              << std::setw(15) << "sin^2(teff_c)" << std::setw(13)
              << EWP.s2teff_f("c") << std::endl
              << std::setw(15) << "sin^2(teff_s)" << std::setw(13)
              << EWP.s2teff_f("s") << std::endl
              << std::setw(15) << "Gamma_inv [GeV]" << std::setw(13)
              << EWP.Gamma_inv() << std::endl
              << std::setw(15) << "Gamma_had [GeV]" << std::setw(13)
              << EWP.Gamma_had() << std::endl
              << std::setw(15) << "Gamma_Z [GeV]" << std::setw(13)
              << EWP.Gamma_Z() << std::endl
              << std::setw(15) << "sigma0_e [nb]" << std::setw(13)
              << EWP.sigma0_l("e")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_mu [nb]" << std::setw(13)
              << EWP.sigma0_l("mu")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_tau [nb]" << std::setw(13)
              << EWP.sigma0_l("tau")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_had [nb]" << std::setw(13)
              << EWP.sigma0_had()*GeVminus2_to_nb << std::endl
              << std::setw(15) << "R0_e" << std::setw(13)
              << EWP.R0_l("e") << std::endl
              << std::setw(15) << "R0_mu" << std::setw(13)
              << EWP.R0_l("mu") << std::endl
              << std::setw(15) << "R0_tau" << std::setw(13)
              << EWP.R0_l("tau") << std::endl
              << std::setw(15) << "R0_b" << std::setw(13)
              << EWP.R0_q("b") << std::endl
              << std::setw(15) << "R0_c" << std::setw(13)
              << EWP.R0_q("c") << std::endl
              << std::setw(15) << "R0_s" << std::setw(13)
              << EWP.R0_q("s") << std::endl
              << std::setw(15) << "A_e" << std::setw(13)
              << EWP.A_f("e") << std::endl
              << std::setw(15) << "A_mu" << std::setw(13)
              << EWP.A_f("mu") << std::endl
              << std::setw(15) << "A_tau" << std::setw(13)
              << EWP.A_f("tau") << std::endl
              << std::setw(15) << "A_b" << std::setw(13)
              << EWP.A_f("b") << std::endl
              << std::setw(15) << "A_c" << std::setw(13)
              << EWP.A_f("c") << std::endl
              << std::setw(15) << "A_s" << std::setw(13)
              << EWP.A_f("s") << std::endl
              << std::setw(15) << "AFB0_e" << std::setw(13)
              << EWP.AFB0_f("e") << std::endl
              << std::setw(15) << "AFB0_mu" << std::setw(13)
              << EWP.AFB0_f("mu") << std::endl
              << std::setw(15) << "AFB0_tau" << std::setw(13)
              << EWP.AFB0_f("tau") << std::endl
              << std::setw(15) << "AFB0_b" << std::setw(13)
              << EWP.AFB0_f("b") << std::endl
              << std::setw(15) << "AFB0_c" << std::setw(13)
              << EWP.AFB0_f("c") << std::endl
              << std::setw(15) << "AFB0_s" << std::setw(13)
              << EWP.AFB0_f("s") << std::endl
              << std::endl;

    std::cout << "  epsilon1 = " << EWP.obliqueEpsilon1() << std::endl
              << "  epsilon2 = " << EWP.obliqueEpsilon2() << std::endl
              << "  epsilon3 = " << EWP.obliqueEpsilon3() << std::endl
              << "         S = " << EWP.obliqueS() << std::endl
              << "         T = " << EWP.obliqueT() << std::endl
              << "         U = " << EWP.obliqueU() << std::endl
              << std::endl;


    return (EXIT_SUCCESS);
}

