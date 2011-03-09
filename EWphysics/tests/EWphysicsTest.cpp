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
#include <ZFitter.h>
#include "EWphysics.h"
#include <StandardModel.h>


int main(int argc, char** argv) {

    /*----- set parameters (same as in ZFitterTest-2) -----*/

    double mZ_i = 91.1875;
    double alsMz_i = 0.118;
    double GF_i = 1.16637e-5;
    double ale_i = 1.0/137.0359895;
    double aleMz_i = 0.007754972109; /* used in ZFitterTest-2 */
    double mcMz_i = 0.563817;
    double mbMz_i = 2.81944;
    double mt_i = 175.0;
    double me_i = 0.000510999;
    double mmu_i = 0.105658;
    double mtau_i = 1.77705;
    //
    /* additional inputs for ZFitter class */
    double mHl_i = 150.0;
    double dAle5Mz_i = 0.02758;
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    //
    /* additional inputs for StandardModel class */
    double mc_i = 1.5;
    double mb_i = 4.2;

    std::cout << std::setprecision(6);
    std::cout << "Inputs:" << std::endl;
    std::cout << "     m_Z      = " << mZ_i << std::endl;
    std::cout << "    m_Hl      = " << mHl_i << std::endl;
    std::cout << "     m_t      = " << mt_i << std::endl;
    std::cout << " alpha_s(m_Z) = " << alsMz_i << std::endl;
    std::cout << "      GF      = " << GF_i << std::endl;
    std::cout << "   alpha(0)   = " << ale_i << std::endl;
    std::cout << "   alpha(mZ)  = " << aleMz_i << std::endl;
    std::cout << "    m_c(m_Z)  = " << mcMz_i << std::endl;
    std::cout << "    m_b(m_Z)  = " << mbMz_i << std::endl;
    std::cout << "     m_e      = " << me_i << std::endl;
    std::cout << "     m_mu     = " << mmu_i << std::endl;
    std::cout << "     m_tau    = " << mtau_i << std::endl;
    std::cout << std::endl;
    std::cout << std::setprecision(6);

    /* Constuct a SM object from the above inputs */
    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
    std::vector<double> bpars_i;
    Parameters Par;
    //
    Par.Set("AlsMz", alsMz_i);
    Par.Set("Mz", mZ_i);
    Par.Set("mup", 0.003);
    Par.Set("mdown", 0.007);
    Par.Set("mcharm", mc_i);
    Par.Set("mstrange", 0.1);
    Par.Set("mtop", mt_i);
    Par.Set("mbottom", mb_i);
    Par.Set("mu1_qcd", mt_i);
    Par.Set("mu2_qcd", mb_i);
    Par.Set("mu3_qcd", mc_i);
    Par.Set("MBd", 0.0);
    Par.Set("MBs", 0.0);
    Par.Set("MBp", 0.0);
    Par.Set("MK0", 0.0);
    Par.Set("MKp", 0.0);
    Par.Set("FBd", 0.0);
    Par.Set("BBd", bpars_i);
    //
    Par.Set("GF", GF_i);
    Par.Set("mneutrino_1", 0.0);
    Par.Set("mneutrino_2", 0.0);
    Par.Set("mneutrino_3", 0.0);
    Par.Set("melectron", me_i);
    Par.Set("mmu", mmu_i);
    Par.Set("mtau", mtau_i);
    Par.Set("VCKM", VCKM_i);
    Par.Set("UPMNS", UPMNS_i);
    Par.Set("ale", ale_i);
    Par.Set("dAle5Mz", dAle5Mz_i);
    Par.Set("mHl", mHl_i);
    Par.Set("muw", 90.0);

    StandardModel SM(Par);

    /* TESTS */
    std::cout << "TESTS for the StandardModel class: " << std::endl;
    std::cout << "  me= " << SM.getMass(SM.ELECTRON) << std::endl;
    std::cout << "  mc(mc)= " << SM.getMass(SM.CHARM) << std::endl;
    std::cout << "  mb(mb)= " << SM.getMass(SM.BOTTOM) << std::endl;
    std::cout << "  1/getAle()= " << 1.0/SM.getAle() << std::endl;
    std::cout << "  getDAle5Mz()= " << SM.getDAle5Mz() << std::endl;
    std::cout << "  dAleLepMz()= " << SM.dAleLepMz() << std::endl;
    std::cout << "  dAleTopMz()= " << SM.dAleTopMz() << std::endl;
    std::cout << "  1/aleMz()= " << 1.0/SM.aleMz() << std::endl;
    std::cout << "  mc(mZ)= " << SM.mcMz() << std::endl;
    std::cout << "  mb(mZ)= " << SM.mbMz() << std::endl;
    std::cout << std::endl;
    

    /*-----------------------------------------------------*/
 
    std::cout << "#######################  ZFitter class  "
              << "#######################" << std::endl << std::endl;

    /* call a constructor of ZFitter with the above parameters (essential)
     * set flags and cuts to be their defalut values */
    ZFitter ZF(mZ_i, mt_i, mHl_i, alsMz_i, dAle5Mz_i,
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

    /* print intermediate results */
    ZF.printIntermediateResults();
    std::cout << std::endl;

    /* calcualte and print EW precision obserbables */
    ZF.printPO();


    /*-----------------------------------------------------*/

    std::cout << "#######################  EWphysics class  "
              << "#######################" << std::endl << std::endl;

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


#ifndef ZFITTER
    ////////////////////////////////////////////////////////////
    /* using the SM object */

    EWphysics EWP(&SM);

    
#else
    ////////////////////////////////////////////////////////////
    /* not using the SM object */

    EWphysics EWP(gZf_i, rhoZf_i, Delta_r_i,
                  mcMz_i, mbMz_i, mt_i, me_i, mmu_i, mtau_i,
                  mZ_i, alsMz_i, GF_i, ale_i, aleMz_i);

    std::cout << "--- Inputs from the ZFitter class to the EWphysics class ---"
              << std::endl << std::endl;
    std::cout << std::setw(10) << "Channel"
              << std::setw(12) << "Re[g_Z^f]"
              << std::setw(13) << "Im[g_Z^f]"
              << std::setw(13) << "Re[rho_Z^f]"
              << std::setw(13) << "Im[rho_Z^f]"
              << std::endl;
    for (indexFermion = 0; indexFermion < 10; indexFermion++) {
        std::cout << std::setw(10) << EWP.flavour_int_to_st(indexFermion)
                  << std::setw(12) << ZF.getCommonARVEFZ(indexFermion)
                  << std::setw(13) << ZF.getCommonAIVEFZ(indexFermion)
                  << std::setw(12) << ZF.getCommonAROTFZ(indexFermion)
                  << std::setw(14) << ZF.getCommonAIROFZ(indexFermion)
                  << std::endl;
    }
    std::cout << std::endl;
    std::cout << "  Delta r = " << ZF.Delta_r() << std::endl << std::endl;

    ////////////////////////////////////////////////////////////
#endif

    
    std::cout << "--- Outputs from EWphysics class ---"
              << std::endl << std::endl;

    std::cout << "Partial widths of the Z boson" << std::endl;
    for (indexFermion = 0; indexFermion < 10; indexFermion++) {
        std::cout << std::setw(10) << EWP.flavour_int_to_st(indexFermion)
                  << std::setw(12) << EWP.Gamma_f(indexFermion)
                  << std::endl;
    }
    std::cout << std::endl;

    // GeV^{-2} --> nb
    const double GeVminus2_to_nb = pow(10.0, -6.0)
                                   / pow(10.0, -28.0)
                                   / pow(299792458.0, -2.0)
                                   / pow(6.58211899 * pow(10.0, -22.0), -2.0)
                                   * pow(10.0, 9.0);

    std::cout << std::setw(15) << "m_W [GeV]" << std::setw(13)
              << EWP.mW() << std::endl
              //<< std::setw(15) << "Gamma_W [GeV]" << std::setw(13)
              //<< EWP.Gamma_W() << std::endl
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

    std::cout << " * Small differences in the patial widths Gamma_q() "
              << "originates from the C_{04} term, which is absent in ZFITTER"
              << std::endl;

    return (EXIT_SUCCESS);
}

