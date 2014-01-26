/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSM_Output.h"

EWSM_Output::EWSM_Output(const EWSM& EWSM_in)
: myEWSM(EWSM_in), SM(myEWSM.SM)
{
}


void EWSM_Output::outputEachDeltaR(const double Mw_i) const
{
    std::cout << "Mw_SM          = " << myEWSM.Mw_SM() << std::endl;
    std::cout << "DeltaR_SM()    = " << myEWSM.DeltaR_SM() << std::endl;
    std::cout << "DeltaRbar_SM() = " << myEWSM.DeltaRbar_SM() << std::endl;
    std::cout << "Mw(input)      = " << Mw_i << std::endl;

    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;

    double DeltaRho[EWSM::orders_EW_size];
    DeltaRho[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaRho(Mw_i);

    double DeltaR_rem[EWSM::orders_EW_size];
    DeltaR_rem[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaR_rem(Mw_i);
    DeltaR_rem[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaR_rem(Mw_i);
    DeltaR_rem[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaR_rem(Mw_i);
    DeltaR_rem[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaR_rem(Mw_i);
    DeltaR_rem[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaR_rem(Mw_i);
    DeltaR_rem[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaR_rem(Mw_i);

    double f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)*sW2_TMP*cW2_TMP/M_PI/SM.getAle();
    //f_AlphaToGF = 1.0; /* for test */
    double DeltaRho_sum = f_AlphaToGF*DeltaRho[EWSM::EW1]
                          + f_AlphaToGF*DeltaRho[EWSM::EW1QCD1]
                          + f_AlphaToGF*DeltaRho[EWSM::EW1QCD2]
                          + pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]
                          + pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2QCD1]
                          + pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3];
    double DeltaRho_G = f_AlphaToGF*DeltaRho[EWSM::EW1];

    if (myEWSM.schemeMw == myEWSM.NORESUM) {

        f_AlphaToGF = 1.0;
        DeltaRho[EWSM::EW1QCD2] *= f_AlphaToGF;
        DeltaRho[EWSM::EW2QCD1] *= pow(f_AlphaToGF, 2.0);
        DeltaRho[EWSM::EW3] *= pow(f_AlphaToGF, 3.0);

        // Full EW one-loop contribution (without the full DeltaAlphaL5q)
        double DeltaR_EW1 = - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW1] + DeltaR_rem[EWSM::EW1];

        // Full EW two-loop contribution with reducible corrections
        double DeltaR_EW2_rem = myEWSM.myApproximateFormulae->DeltaR_TwoLoopEW_rem(Mw_i);

        // EW two-loop irreducible contributions with large-mt expansion
        double DeltaR_EW2_old_red = myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q()
                                    - 2.0*cW2_TMP/sW2_TMP*myEWSM.DeltaAlphaL5q()*DeltaRho[EWSM::EW1]
                                    + pow(cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW1], 2.0);
        double DeltaR_EW2_old_irred = - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW2] + DeltaR_rem[EWSM::EW2];

        // Delta r, including the full EW two-loop contribution
        double deltaR = myEWSM.DeltaAlphaL5q();
        for (int j=0; j<EWSM::orders_EW_size; ++j) {
            deltaR += - cW2_TMP/sW2_TMP*DeltaRho[(EWSM::orders_EW)j];
            deltaR += DeltaR_rem[(EWSM::orders_EW)j];
        }
        deltaR -= -cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW2];
        deltaR -= DeltaR_rem[EWSM::EW2];
        deltaR += myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() + 2.0*myEWSM.DeltaAlphaL5q()*DeltaR_EW1 + DeltaR_EW2_rem;

        std::cout << "(1+dr) - 1        =  " << deltaR << std::endl;
        std::cout << "  EW1             =  " << myEWSM.DeltaAlphaL5q() + DeltaR_EW1 << std::endl;
        std::cout << "    DeltaAlphaL5q =  " << myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    dR            = " << DeltaR_EW1 << std::endl;
        std::cout << "  EW1QCD1         =  " << - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW1QCD1] + DeltaR_rem[EWSM::EW1QCD1] << std::endl;
        std::cout << "  EW2(full)       =  " << DeltaR_EW2_rem + myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() + 2.0*myEWSM.DeltaAlphaL5q()*DeltaR_EW1 << std::endl;
        std::cout << "    dAle*dAle     =  " << myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    2*dAle*dR     = " << 2.0*myEWSM.DeltaAlphaL5q()*DeltaR_EW1 << std::endl;
        std::cout << "    others        =  " << DeltaR_EW2_rem << std::endl;
        std::cout << "  EW1QCD2         =  " << - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW1QCD2] + DeltaR_rem[EWSM::EW1QCD2] << std::endl;
        std::cout << "  EW2QCD1         = " << - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW2QCD1] + DeltaR_rem[EWSM::EW2QCD1] << std::endl;
        std::cout << "  EW3             = " << - cW2_TMP/sW2_TMP*DeltaRho[EWSM::EW3] + DeltaR_rem[EWSM::EW3] << std::endl;
        std::cout << "  EW2(old,irreducible) = " << DeltaR_EW2_old_irred << std::endl;
        std::cout << "  EW2(old,red+irred)   = " << DeltaR_EW2_old_red + DeltaR_EW2_old_irred << std::endl;
        std::cout << "  EW2(old,red+irred-dAle*dAle-2*dAle*dR) = "
                  << DeltaR_EW2_old_red + DeltaR_EW2_old_irred
                     - myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q()
                     - 2.0*myEWSM.DeltaAlphaL5q()*DeltaR_EW1 << std::endl;

    } else if (myEWSM.schemeMw == myEWSM.OMSI) {

        // R = 1/(1 - Delta r)
        double R = 1.0/( 1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)
                   /(1.0 - myEWSM.DeltaAlphaL5q() - DeltaR_rem[EWSM::EW1] - DeltaR_rem[EWSM::EW1QCD1] - DeltaR_rem[EWSM::EW2] );

        std::cout << "1/(1-dr) - 1 (exact)                 = " << R - 1.0 << std::endl;
        std::cout << "     --> dr = " << 1.0 - 1.0/R << std::endl;

        // each contribution
        double DeltaR_EW1 = myEWSM.DeltaAlphaL5q() - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1] + DeltaR_rem[EWSM::EW1];
        double DeltaR_EW1QCD1 = - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD1] + DeltaR_rem[EWSM::EW1QCD1];
        double DeltaR_EW2 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]
                            + DeltaR_rem[EWSM::EW2]
                            + cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]
                              *(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1])
                            + DeltaR_EW1*DeltaR_EW1;
        double DeltaR_EW1QCD2 = - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD2];
        double DeltaR_EW2QCD1 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2QCD1]
                                + cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD1]
                                  *(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1])
                                + 2.0*DeltaR_EW1*DeltaR_EW1QCD1;
        double DeltaR_EW3 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3]
                            + cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]
                              *(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1])
                            + pow(DeltaR_EW1, 3.0)
                            + 2.0*DeltaR_EW1*(DeltaR_EW2 - DeltaR_EW1*DeltaR_EW1);

        std::cout << "  EW1             =  " << DeltaR_EW1 << std::endl;
        std::cout << "    DeltaAlphaL5q =  " << myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    -cW2/sW2*dRho1= " << - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]  << std::endl;
        std::cout << "    DeltaR1_rem   =  " << DeltaR_rem[EWSM::EW1] << std::endl;
        std::cout << "  EW1QCD1         =  " << DeltaR_EW1QCD1 << std::endl;
        std::cout << "  EW2(full)       =  " << DeltaR_EW2 << std::endl;
        std::cout << "    EW1*EW1       =  " << DeltaR_EW1*DeltaR_EW1 << std::endl;
        std::cout << "      dAle*dAle   =  " << myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "      others      = " << DeltaR_EW1*DeltaR_EW1 - myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    -cW2/sW2*dRho2=  " << - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2] << std::endl;
        std::cout << "    DeltaR2_rem   =  " << DeltaR_rem[EWSM::EW2] << std::endl;
        std::cout << "    others        =  " << cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]*(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1]) << std::endl;
        std::cout << "  EW1QCD2         =  " << DeltaR_EW1QCD2 << std::endl;
        std::cout << "  EW2QCD1         = " << DeltaR_EW2QCD1 << std::endl;
        std::cout << "  EW3             =  " << DeltaR_EW3 << std::endl;
        std::cout << "    -cW2/sW2*dRho3= " << - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3] << std::endl;
        std::cout << "    EW1^3         =  " << pow(DeltaR_EW1, 3.0) << std::endl;
        std::cout << "    2*EW1*(EW2-EW1^2)=" << 2.0*DeltaR_EW1*(DeltaR_EW2 - DeltaR_EW1*DeltaR_EW1) << std::endl;
        std::cout << "    others        = " << cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]*(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1]) << std::endl;

    } else if (myEWSM.schemeMw == myEWSM.OMSII) {

        // R = 1/(1 - Delta r)
        double R = 1.0/( (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)*(1.0 - myEWSM.DeltaAlphaL5q())
                         - (1.0 + cW2_TMP/sW2_TMP*DeltaRho_G)*DeltaR_rem[EWSM::EW1]
                         - DeltaR_rem[EWSM::EW1QCD1] - DeltaR_rem[EWSM::EW2] );

        // each contribution
        double DeltaR_EW1 = myEWSM.DeltaAlphaL5q() - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1] + DeltaR_rem[EWSM::EW1];
        double DeltaR_EW1QCD1 = - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD1] + DeltaR_rem[EWSM::EW1QCD1];
        double DeltaR_EW2 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]
                            + DeltaR_rem[EWSM::EW2]
                            + cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]
                              *(myEWSM.DeltaAlphaL5q() + DeltaR_rem[EWSM::EW1])
                            + DeltaR_EW1*DeltaR_EW1;
        double DeltaR_EW1QCD2 = - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD2];
        double DeltaR_EW2QCD1 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2QCD1]
                                + cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1QCD1]*myEWSM.DeltaAlphaL5q()
                                + 2.0*DeltaR_EW1*DeltaR_EW1QCD1;
        double DeltaR_EW3 = - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3]
                            + cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]*myEWSM.DeltaAlphaL5q()
                            + pow(DeltaR_EW1, 3.0)
                            + 2.0*DeltaR_EW1*(DeltaR_EW2 - DeltaR_EW1*DeltaR_EW1);

        std::cout << "1/(1-dr) - 1 (exact)                 = " << R - 1.0 << std::endl;
        std::cout << "     --> dr = " << 1.0 - 1.0/R << std::endl;
        std::cout << "1/(1-dr) - 1 (sum of expanded terms) = "
                  << DeltaR_EW1 + DeltaR_EW1QCD1 + DeltaR_EW2 + DeltaR_EW1QCD2
                     + DeltaR_EW2QCD1 + DeltaR_EW3 << std::endl;
        std::cout << "  EW1             =  " << DeltaR_EW1 << std::endl;
        std::cout << "    DeltaAlphaL5q =  " << myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    -cW2/sW2*dRho1= " << - cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]  << std::endl;
        std::cout << "    DeltaR1_rem   =  " << DeltaR_rem[EWSM::EW1] << std::endl;
        std::cout << "  EW1QCD1         =  " << DeltaR_EW1QCD1 << std::endl;
        std::cout << "  EW2(full)       =  " << DeltaR_EW2 << std::endl;
        std::cout << "    EW1*EW1       =  " << DeltaR_EW1*DeltaR_EW1 << std::endl;
        std::cout << "      dAle*dAle   =  " << myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "      others      = " << DeltaR_EW1*DeltaR_EW1 - myEWSM.DeltaAlphaL5q()*myEWSM.DeltaAlphaL5q() << std::endl;
        std::cout << "    -cW2/sW2*dRho2=  " << - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2] << std::endl;
        std::cout << "    DeltaR2_rem   =  " << DeltaR_rem[EWSM::EW2] << std::endl;
        std::cout << "    others        =  " << cW2_TMP/sW2_TMP*f_AlphaToGF*DeltaRho[EWSM::EW1]*(myEWSM.DeltaAlphaL5q()+DeltaR_rem[EWSM::EW1]) << std::endl;
        std::cout << "  EW1QCD2         =  " << DeltaR_EW1QCD2 << std::endl;
        std::cout << "  EW2QCD1         = " << DeltaR_EW2QCD1 << std::endl;
        std::cout << "  EW3             =  " << DeltaR_EW3 << std::endl;
        std::cout << "    -cW2/sW2*dRho3= " << - cW2_TMP/sW2_TMP*pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3] << std::endl;
        std::cout << "    EW1^3         =  " << pow(DeltaR_EW1, 3.0) << std::endl;
        std::cout << "    2*EW1*(EW2-EW1^2)=" << 2.0*DeltaR_EW1*(DeltaR_EW2 - DeltaR_EW1*DeltaR_EW1) << std::endl;
        std::cout << "    others        = " << cW2_TMP/sW2_TMP*pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2]*myEWSM.DeltaAlphaL5q() << std::endl;

    } else
        std::cout << "EWSM_Output::outputEachDeltaR(): Not implemented for schemeMw="
                  << myEWSM.schemeMw << std::endl;
}


void EWSM_Output::outputEachDeltaRhoZ_l(const StandardModel::lepton l, const double Mw_i) const
{
    std::cout << "================================================" << std::endl;
    std::cout << "rhoZ_l[(StandardModel::lepton)" << l << "]" << std::endl;
    std::cout << "Mw(input)   = " << Mw_i << std::endl;

    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;

    double DeltaRho[EWSM::orders_EW_size];
    DeltaRho[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaRho(Mw_i);

    /* compute delta rho_rem^f */
    complex deltaRho_rem_f[EWSM::orders_EW_size];
    deltaRho_rem_f[EWSM::EW1] = myEWSM.myOneLoopEW->deltaRho_rem_l(l,Mw_i);
    #ifdef WITHIMTWOLOOPQCD
    deltaRho_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaRho_rem_l(l,Mw_i).real(),
                                      myEWSM.myTwoLoopQCD->deltaRho_rem_l(l,Mw_i).imag(), false);
    #else
    deltaRho_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaRho_rem_l(l,Mw_i).real(), 0.0, false);
    #endif
    deltaRho_rem_f[EWSM::EW1QCD2] = complex(myEWSM.myThreeLoopQCD->deltaRho_rem_l(l,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW2] = complex(myEWSM.myTwoLoopEW->deltaRho_rem_l(l,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW2QCD1] = complex(myEWSM.myThreeLoopEW2QCD->deltaRho_rem_l(l,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW3] = complex(myEWSM.myThreeLoopEW->deltaRho_rem_l(l,Mw_i).real(), 0.0, false);

    /* compute Delta rbar_rem */
    double DeltaRbar_rem = myEWSM.myOneLoopEW->DeltaRbar_rem(Mw_i);

    double f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                         *sW2_TMP*cW2_TMP/M_PI/SM.getAle();

    /* Re[rho_Z^f] with or without resummation */
    double deltaRho_rem_f_real[EWSM::orders_EW_size];
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();

    double dummy[EWSM::orders_EW_size];
    outputEachDeltaRhoZ(f_AlphaToGF, DeltaRho, deltaRho_rem_f_real,
                        DeltaRbar_rem, false, dummy, 0.0);

    /* Im[rho_Z^f] without resummation */
    double ImRhoZf = 0.0;
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        ImRhoZf += deltaRho_rem_f[j].imag();
    std::cout << "ImRhoZf(with alpha)=" << ImRhoZf << std::endl;

    std::cout << "================================================" << std::endl;
}


void EWSM_Output::outputEachDeltaRhoZ_q(const StandardModel::quark q, const double Mw_i) const
{
    std::cout << "================================================" << std::endl;
    std::cout << "rhoZ_q[(StandardModel::quark)" << q << "]" << std::endl;
    std::cout << "Mw(input)   = " << Mw_i << std::endl;

    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;

    double DeltaRho[EWSM::orders_EW_size];
    DeltaRho[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaRho(Mw_i);

    /* compute delta rho_rem^f */
    complex deltaRho_rem_f[EWSM::orders_EW_size];
    deltaRho_rem_f[EWSM::EW1] = myEWSM.myOneLoopEW->deltaRho_rem_q(q,Mw_i);
    #ifdef WITHIMTWOLOOPQCD
    deltaRho_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaRho_rem_q(q,Mw_i).real(),
                                      myEWSM.myTwoLoopQCD->deltaRho_rem_q(q,Mw_i).imag(), false);
    #else
    deltaRho_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaRho_rem_q(q,Mw_i).real(), 0.0, false);
    #endif
    deltaRho_rem_f[EWSM::EW1QCD2] = complex(myEWSM.myThreeLoopQCD->deltaRho_rem_q(q,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW2] = complex(myEWSM.myTwoLoopEW->deltaRho_rem_q(q,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW2QCD1] = complex(myEWSM.myThreeLoopEW2QCD->deltaRho_rem_q(q,Mw_i).real(), 0.0, false);
    deltaRho_rem_f[EWSM::EW3] = complex(myEWSM.myThreeLoopEW->deltaRho_rem_q(q,Mw_i).real(), 0.0, false);

    /* compute Delta rbar_rem */
    double DeltaRbar_rem = myEWSM.myOneLoopEW->DeltaRbar_rem(Mw_i);

    /* conversion factor */
    double f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                         *sW2_TMP*cW2_TMP/M_PI/SM.getAle();

    /* Zbb */
    bool bool_Zbb = false;
    if (q==StandardModel::BOTTOM) bool_Zbb = true;
    double ZbbSubtract = 0.0;
    if (bool_Zbb)
        ZbbSubtract = - myEWSM.myCache->ale()/4.0/M_PI/sW2_TMP
                        *pow(myEWSM.myCache->Mt()/Mw_i, 2.0);
    double taub[EWSM::orders_EW_size];
    double Xt = myEWSM.myCache->Xt_alpha(Mw_i);
    if (bool_Zbb) {
        taub[EWSM::EW1] = -2.0*Xt;
        taub[EWSM::EW1QCD1] = 2.0/3.0*M_PI*Xt*myEWSM.myCache->alsMt();
        taub[EWSM::EW2]= -2.0*Xt*Xt*myEWSM.myTwoLoopEW->tau_2();
    }

    /* Re[rho_Z^f] with or without resummation */
    double deltaRho_rem_f_real[EWSM::orders_EW_size];
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();

    outputEachDeltaRhoZ(f_AlphaToGF, DeltaRho, deltaRho_rem_f_real,
                        DeltaRbar_rem, bool_Zbb, taub, ZbbSubtract);

    /* Im[rho_Z^f] without resummation */
    double ImRhoZf = 0.0;
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        ImRhoZf += deltaRho_rem_f[j].imag();
    std::cout << "ImRhoZf(with alpha)=" << ImRhoZf << std::endl;

    std::cout << "================================================" << std::endl;
}


void EWSM_Output::outputEachDeltaRhoZ(const double f_AlphaToGF,
                               const double DeltaRho[EWSM::orders_EW_size],
                               const double deltaRho_rem[EWSM::orders_EW_size],
                               const double DeltaRbar_rem,
                               const bool bool_Zbb,
                               const double taub[EWSM::orders_EW_size],
                               const double ZbbSubtract) const
{
    if (myEWSM.schemeRhoZ == myEWSM.APPROXIMATEFORMULA) {

    } else if (myEWSM.schemeRhoZ == myEWSM.NORESUM) {
        std::cout << "Leading contributions: alpha or Gmu" << std::endl;
        std::cout << "  DeltaRho[EW1]=" << DeltaRho[EWSM::EW1] << " "
                   << f_AlphaToGF*DeltaRho[EWSM::EW1] << std::endl;
        std::cout << "  DeltaRho[EW1QCD1]=" << DeltaRho[EWSM::EW1QCD1] << " "
                  << f_AlphaToGF*DeltaRho[EWSM::EW1QCD1] << std::endl;
        std::cout << "  DeltaRho[EW1QCD2]=" << DeltaRho[EWSM::EW1QCD2] << " "
                  << f_AlphaToGF*DeltaRho[EWSM::EW1QCD2] << std::endl;
        std::cout << "  DeltaRho[EW2]=" << DeltaRho[EWSM::EW2] << " "
                  << pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2] << std::endl;
        std::cout << "  DeltaRho[EW2QCD1]=" << DeltaRho[EWSM::EW2QCD1] << " "
                  << pow(f_AlphaToGF,2.0)*DeltaRho[EWSM::EW2QCD1] << std::endl;
        std::cout << "  DeltaRho[EW3]=" << DeltaRho[EWSM::EW3] << " "
                  << pow(f_AlphaToGF,3.0)*DeltaRho[EWSM::EW3] << std::endl;
        std::cout << "Remainder contributions: alpha or Gmu" << std::endl;
        std::cout << "  DeltaRbar_rem[EW1]=" << DeltaRbar_rem << " "
                  << f_AlphaToGF*DeltaRbar_rem << std::endl;
        std::cout << "  deltaRho_rem[EW1]=" << deltaRho_rem[EWSM::EW1] - ZbbSubtract << " "
                  << f_AlphaToGF*(deltaRho_rem[EWSM::EW1] - ZbbSubtract) << std::endl;
        std::cout << "  deltaRho_rem[EW1QCD1]=" << deltaRho_rem[EWSM::EW1QCD1] << " "
                  << f_AlphaToGF*deltaRho_rem[EWSM::EW1QCD1] << std::endl;
        std::cout << "  deltaRho_rem[EW2]=" << deltaRho_rem[EWSM::EW2] << " "
                  << pow(f_AlphaToGF,2.0)*deltaRho_rem[EWSM::EW2] << std::endl;
        if (bool_Zbb) {
            std::cout << "Taub: alpha or Gmu" << std::endl;
            std::cout << "  taub[EW1]=" << taub[EWSM::EW1] << " "
                      << f_AlphaToGF*taub[EWSM::EW1] << std::endl;
            std::cout << "  taub[EW1QCD1]=" << taub[EWSM::EW1QCD1] << " "
                      << f_AlphaToGF*taub[EWSM::EW1QCD1] << std::endl;
            std::cout << "  taub[EW2]=" << taub[EWSM::EW2] << " "
                      << pow(f_AlphaToGF,2.0)*taub[EWSM::EW2] << std::endl;
        }
        std::cout << "Each order: alpha or Gmu" << std::endl;
        double dRho_EW1 = DeltaRho[EWSM::EW1] + deltaRho_rem[EWSM::EW1] - ZbbSubtract;
        double dRho_EW1QCD1 = DeltaRho[EWSM::EW1QCD1] + deltaRho_rem[EWSM::EW1QCD1];
        double dRho_EW2 = DeltaRho[EWSM::EW1]*DeltaRho[EWSM::EW1]
                          + (deltaRho_rem[EWSM::EW1] - ZbbSubtract)*DeltaRho[EWSM::EW1]
                          - DeltaRho[EWSM::EW1]*DeltaRbar_rem
                          + DeltaRho[EWSM::EW2]
                          + deltaRho_rem[EWSM::EW2];
        double dRho_EW1QCD2 = DeltaRho[EWSM::EW1QCD2];
        double dRho_EW2QCD1 = DeltaRho[EWSM::EW2QCD1] + deltaRho_rem[EWSM::EW1QCD1]*DeltaRho[EWSM::EW1];
        double dRho_EW3 = DeltaRho[EWSM::EW3];
        //
        double dRho_EW1_TMP = dRho_EW1;
        double dRho_EW1QCD1_TMP = dRho_EW1QCD1;
        double dRho_EW2_TMP = dRho_EW2;
        double dRho_EW1QCD2_TMP = dRho_EW1QCD2;
        double dRho_EW2QCD1_TMP = dRho_EW2QCD1;
        double dRho_EW3_TMP = dRho_EW3;
        //
        if (bool_Zbb) {
            dRho_EW1 = dRho_EW1_TMP + 2.0*taub[EWSM::EW1];
            dRho_EW1QCD1 = dRho_EW1QCD1_TMP + 2.0*taub[EWSM::EW1QCD1];
            dRho_EW2 = dRho_EW2_TMP + 2.0*taub[EWSM::EW2] + taub[EWSM::EW1]*taub[EWSM::EW1]
                       + dRho_EW1_TMP*2.0*taub[EWSM::EW1];
            dRho_EW1QCD2 = dRho_EW1QCD2_TMP;
            dRho_EW2QCD1 = dRho_EW2QCD1_TMP + dRho_EW1_TMP*2.0*taub[EWSM::EW1QCD1]
                           + dRho_EW1QCD1_TMP*2.0*taub[EWSM::EW1] + 2.0*taub[EWSM::EW1]*taub[EWSM::EW1QCD1];
            dRho_EW3 = dRho_EW3_TMP + dRho_EW1_TMP*2.0*taub[EWSM::EW2]
                       + dRho_EW2_TMP*2.0*taub[EWSM::EW1] + 2.0*taub[EWSM::EW2]*taub[EWSM::EW1];
        }
        std::cout << "  EW1: " << dRho_EW1 << " " << f_AlphaToGF*dRho_EW1 << std::endl;
        std::cout << "  EW1QCD1: " << dRho_EW1QCD1 << " " << f_AlphaToGF*dRho_EW1QCD1 << std::endl;
        std::cout << "  EW2: " << dRho_EW2 << " " << pow(f_AlphaToGF,2.0)*dRho_EW2 << std::endl;
        std::cout << "  EW1QCD2: " << dRho_EW1QCD2 << " " << f_AlphaToGF*dRho_EW1QCD2 << std::endl;
        std::cout << "  EW2QCD1: " << dRho_EW2QCD1 << " " << pow(f_AlphaToGF,2.0)*dRho_EW2QCD1 << std::endl;
        std::cout << "  EW3: " << dRho_EW3 << " " << pow(f_AlphaToGF,3.0)*dRho_EW3 << std::endl;
        std::cout << "Total contribution: alpha or Gmu" << std::endl;
        std::cout << "  rhoZ="
                  << 1.0 + dRho_EW1 + dRho_EW1QCD1 + dRho_EW2
                     + dRho_EW1QCD2 + dRho_EW2QCD1 + dRho_EW3
                  << " "
                  << 1.0 + f_AlphaToGF*dRho_EW1 + f_AlphaToGF*dRho_EW1QCD1
                     + pow(f_AlphaToGF,2.0)*dRho_EW2
                     + f_AlphaToGF*dRho_EW1QCD2 + pow(f_AlphaToGF,2.0)*dRho_EW2QCD1
                     + pow(f_AlphaToGF,3.0)*dRho_EW3
                  << std::endl;
        if (bool_Zbb) {
            std::cout << "  rhoZ(taub resummed)="
                      << (1.0 + dRho_EW1_TMP + dRho_EW1QCD1_TMP + dRho_EW2_TMP
                          + dRho_EW1QCD2_TMP + dRho_EW2QCD1_TMP + dRho_EW3_TMP)
                         *pow(1.0 + taub[EWSM::EW1] + taub[EWSM::EW1QCD1] + taub[EWSM::EW2],2.0)
                      << " "
                      << (1.0 + f_AlphaToGF*dRho_EW1_TMP + f_AlphaToGF*dRho_EW1QCD1_TMP
                          + pow(f_AlphaToGF,2.0)*dRho_EW2_TMP
                          + f_AlphaToGF*dRho_EW1QCD2_TMP + pow(f_AlphaToGF,2.0)*dRho_EW2QCD1_TMP
                          + pow(f_AlphaToGF,3.0)*dRho_EW3_TMP)
                         *pow(1.0 + f_AlphaToGF*taub[EWSM::EW1] + f_AlphaToGF*taub[EWSM::EW1QCD1]
                              + pow(f_AlphaToGF,2.0)*taub[EWSM::EW2],2.0)
                      << std::endl;
        }
    } else
        std::cout << "EWSM_Output::outputEachDeltaRhoZ(): Not implemented for schemeRhoZ="
                  << myEWSM.schemeRhoZ << std::endl;
}


void EWSM_Output::outputEachDeltaKappaZ_l(const StandardModel::lepton l, const double Mw_i) const
{
    std::cout << "================================================" << std::endl;
    std::cout << "kappaZ_l[(StandardModel::lepton)" << l << "]" << std::endl;
    std::cout << "Mw(input)   = " << Mw_i << std::endl;

    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;

    double DeltaRho[EWSM::orders_EW_size];
    DeltaRho[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaRho(Mw_i);

    /* compute delta kappa_rem^f */
    complex deltaKappa_rem_f[EWSM::orders_EW_size];
    deltaKappa_rem_f[EWSM::EW1] = myEWSM.myOneLoopEW->deltaKappa_rem_l(l,Mw_i);
    #ifdef WITHIMTWOLOOPQCD
    deltaKappa_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaKappa_rem_l(l,Mw_i).real(),
                                        myEWSM.myTwoLoopQCD->deltaKappa_rem_l(l,Mw_i).imag(), false);
    #else
    deltaKappa_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaKappa_rem_l(l,Mw_i).real(), 0.0, false);
    #endif
    deltaKappa_rem_f[EWSM::EW1QCD2] = complex(myEWSM.myThreeLoopQCD->deltaKappa_rem_l(l,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW2] = complex(myEWSM.myTwoLoopEW->deltaKappa_rem_l(l,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW2QCD1] = complex(myEWSM.myThreeLoopEW2QCD->deltaKappa_rem_l(l,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW3] = complex(myEWSM.myThreeLoopEW->deltaKappa_rem_l(l,Mw_i).real(), 0.0, false);

    /* compute Delta rbar_rem */
    double DeltaRbar_rem = myEWSM.myOneLoopEW->DeltaRbar_rem(Mw_i);

    /* conversion factor */
    double f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                         *sW2_TMP*cW2_TMP/M_PI/SM.getAle();

    /* Re[Kappa_Z^f] with or without resummation */
    double deltaKappa_rem_f_real[EWSM::orders_EW_size];
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();

    /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
    double ReKappaZf = myEWSM.resumKappaZ(DeltaRho, deltaKappa_rem_f_real,
                                          DeltaRbar_rem, false);
    double Zgamma_EW2 = 35.0*myEWSM.alphaMz()*myEWSM.alphaMz()/18.0/sW2_TMP
                        *(1.0 - 8.0/3.0*ReKappaZf*sW2_TMP);

    double dummy[EWSM::orders_EW_size];
    outputEachDeltaKappaZ(f_AlphaToGF, cW2_TMP/sW2_TMP,
                          DeltaRho, deltaKappa_rem_f_real,
                          DeltaRbar_rem, false, dummy, 0.0, Zgamma_EW2);

    /* Im[kappa_Z^f] without resummation */
    double ImKappaZf = 0.0;
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        ImKappaZf += deltaKappa_rem_f[j].imag();
    std::cout << "ImKappaZf(with alpha)=" << ImKappaZf << std::endl;

    std::cout << "================================================" << std::endl;
}


void EWSM_Output::outputEachDeltaKappaZ_q(const StandardModel::quark q, const double Mw_i) const
{
    std::cout << "================================================" << std::endl;
    std::cout << "kappaZ_q[(StandardModel::quark)" << q << "]" << std::endl;
    std::cout << "Mw(input)   = " << Mw_i << std::endl;

    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;

    double DeltaRho[EWSM::orders_EW_size];
    DeltaRho[EWSM::EW1] = myEWSM.myOneLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD1] = myEWSM.myTwoLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW1QCD2] = myEWSM.myThreeLoopQCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2] = myEWSM.myTwoLoopEW->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW2QCD1] = myEWSM.myThreeLoopEW2QCD->DeltaRho(Mw_i);
    DeltaRho[EWSM::EW3] = myEWSM.myThreeLoopEW->DeltaRho(Mw_i);

    /* compute delta kappa_rem^f */
    complex deltaKappa_rem_f[EWSM::orders_EW_size];
    deltaKappa_rem_f[EWSM::EW1] = myEWSM.myOneLoopEW->deltaKappa_rem_q(q,Mw_i);
    #ifdef WITHIMTWOLOOPQCD
    deltaKappa_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaKappa_rem_q(q,Mw_i).real(),
                                        myEWSM.myTwoLoopQCD->deltaKappa_rem_q(q,Mw_i).imag(), false);
    #else
    deltaKappa_rem_f[EWSM::EW1QCD1] = complex(myEWSM.myTwoLoopQCD->deltaKappa_rem_q(q,Mw_i).real(), 0.0, false);
    #endif
    deltaKappa_rem_f[EWSM::EW1QCD2] = complex(myEWSM.myThreeLoopQCD->deltaKappa_rem_q(q,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW2] = complex(myEWSM.myTwoLoopEW->deltaKappa_rem_q(q,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW2QCD1] = complex(myEWSM.myThreeLoopEW2QCD->deltaKappa_rem_q(q,Mw_i).real(), 0.0, false);
    deltaKappa_rem_f[EWSM::EW3] = complex(myEWSM.myThreeLoopEW->deltaKappa_rem_q(q,Mw_i).real(), 0.0, false);

    /* compute Delta rbar_rem */
    double DeltaRbar_rem = myEWSM.myOneLoopEW->DeltaRbar_rem(Mw_i);

    /* conversion factor */
    double f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                         *sW2_TMP*cW2_TMP/M_PI/SM.getAle();

    /* Zbb */
    bool bool_Zbb = false;
    if (q==StandardModel::BOTTOM) bool_Zbb = true;
    double ZbbSubtract = 0.0;
    if (bool_Zbb)
        ZbbSubtract = myEWSM.myCache->ale()/8.0/M_PI/sW2_TMP
                      *pow(myEWSM.myCache->Mt()/Mw_i, 2.0);
    double taub[EWSM::orders_EW_size];
    double Xt = myEWSM.myCache->Xt_alpha(Mw_i);
    if (bool_Zbb) {
        taub[EWSM::EW1] = -2.0*Xt;
        taub[EWSM::EW1QCD1] = 2.0/3.0*M_PI*Xt*myEWSM.myCache->alsMt();
        taub[EWSM::EW2]= -2.0*Xt*Xt*myEWSM.myTwoLoopEW->tau_2();
    }

    /* Re[Kappa_Z^f] with or without resummation */
    double deltaKappa_rem_f_real[EWSM::orders_EW_size];
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();

    /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
    double ReKappaZf = myEWSM.resumKappaZ(DeltaRho, deltaKappa_rem_f_real,
                                          DeltaRbar_rem, bool_Zbb);
    double Zgamma_EW2 = 35.0*myEWSM.alphaMz()*myEWSM.alphaMz()/18.0/sW2_TMP
                        *(1.0 - 8.0/3.0*ReKappaZf*sW2_TMP);

    outputEachDeltaKappaZ(f_AlphaToGF, cW2_TMP/sW2_TMP,
                          DeltaRho, deltaKappa_rem_f_real,
                          DeltaRbar_rem, bool_Zbb, taub, ZbbSubtract, Zgamma_EW2);

    /* Im[kappa_Z^f] without resummation */
    double ImKappaZf = 0.0;
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        ImKappaZf += deltaKappa_rem_f[j].imag();
    std::cout << "ImKappaZf(with alpha)=" << ImKappaZf << std::endl;

    std::cout << "================================================" << std::endl;
}


void EWSM_Output::outputEachDeltaKappaZ(const double f_AlphaToGF,
                                 const double cW2overSW2,
                                 const double DeltaRho[EWSM::orders_EW_size],
                                 const double deltaKappa_rem[EWSM::orders_EW_size],
                                 const double DeltaRbar_rem,
                                 const bool bool_Zbb,
                                 const double taub[EWSM::orders_EW_size],
                                 const double ZbbSubtract,
                                 const double Zgamma_EW2) const
{
    /* rescale */
    double DeltaRho_new[EWSM::orders_EW_size];
    for (int j=0; j<EWSM::orders_EW_size; ++j)
        DeltaRho_new[j] = cW2overSW2*DeltaRho[j];

    if (myEWSM.schemeKappaZ == myEWSM.APPROXIMATEFORMULA) {
        std::cout << "Delta kappaZb (from the approximate formula of sin2thb) = "
                  << myEWSM.kappaZ_q_SM(StandardModel::BOTTOM)  - 1.0 << std::endl;
    } else if (myEWSM.schemeKappaZ == myEWSM.NORESUM) {
        std::cout << "Leading contributions: alpha or Gmu" << std::endl;
        std::cout << "  DeltaRho[EW1]=" << DeltaRho_new[EWSM::EW1] << " "
                   << f_AlphaToGF*DeltaRho_new[EWSM::EW1] << std::endl;
        std::cout << "  DeltaRho[EW1QCD1]=" << DeltaRho_new[EWSM::EW1QCD1] << " "
                  << f_AlphaToGF*DeltaRho_new[EWSM::EW1QCD1] << std::endl;
        std::cout << "  DeltaRho[EW1QCD2]=" << DeltaRho_new[EWSM::EW1QCD2] << " "
                  << f_AlphaToGF*DeltaRho_new[EWSM::EW1QCD2] << std::endl;
        std::cout << "  DeltaRho[EW2]=" << DeltaRho_new[EWSM::EW2] << " "
                  << pow(f_AlphaToGF,2.0)*DeltaRho_new[EWSM::EW2] << std::endl;
        std::cout << "  DeltaRho[EW2QCD1]=" << DeltaRho_new[EWSM::EW2QCD1] << " "
                  << pow(f_AlphaToGF,2.0)*DeltaRho_new[EWSM::EW2QCD1] << std::endl;
        std::cout << "  DeltaRho[EW3]=" << DeltaRho_new[EWSM::EW3] << " "
                  << pow(f_AlphaToGF,3.0)*DeltaRho_new[EWSM::EW3] << std::endl;
        std::cout << "EW2 from Z-gamma = " << Zgamma_EW2 << std::endl;
        std::cout << "Remainder contributions: alpha or Gmu" << std::endl;
        std::cout << "  DeltaRbar_rem[EW1]=" << DeltaRbar_rem << " "
                  << f_AlphaToGF*DeltaRbar_rem << std::endl;
        std::cout << "  deltaKappa_rem[EW1]=" << deltaKappa_rem[EWSM::EW1] - ZbbSubtract << " "
                  << f_AlphaToGF*(deltaKappa_rem[EWSM::EW1] - ZbbSubtract) << std::endl;
        std::cout << "  deltaKappa_rem[EW1QCD1]=" << deltaKappa_rem[EWSM::EW1QCD1] << " "
                  << f_AlphaToGF*deltaKappa_rem[EWSM::EW1QCD1] << std::endl;
        std::cout << "  deltaKappa_rem[EW1QCD2]=" << deltaKappa_rem[EWSM::EW1QCD2] << " "
                  << f_AlphaToGF*deltaKappa_rem[EWSM::EW1QCD2] << std::endl;
        std::cout << "  deltaKappa_rem[EW2]=" << deltaKappa_rem[EWSM::EW2] << " "
                  << pow(f_AlphaToGF,2.0)*deltaKappa_rem[EWSM::EW2] << std::endl;
        if (bool_Zbb) {
            std::cout << "Taub: alpha or Gmu" << std::endl;
            std::cout << "  taub[EW1]=" << taub[EWSM::EW1] << " "
                      << f_AlphaToGF*taub[EWSM::EW1] << std::endl;
            std::cout << "  taub[EW1QCD1]=" << taub[EWSM::EW1QCD1] << " "
                      << f_AlphaToGF*taub[EWSM::EW1QCD1] << std::endl;
            std::cout << "  taub[EW2]=" << taub[EWSM::EW2] << " "
                      << pow(f_AlphaToGF,2.0)*taub[EWSM::EW2] << std::endl;
        }
        std::cout << "Each order: alpha or Gmu" << std::endl;
        double dKappa_EW1 = DeltaRho_new[EWSM::EW1] + deltaKappa_rem[EWSM::EW1] - ZbbSubtract;
        double dKappa_EW1QCD1 = DeltaRho_new[EWSM::EW1QCD1] + deltaKappa_rem[EWSM::EW1QCD1];
        double dKappa_EW2 = (deltaKappa_rem[EWSM::EW1] - ZbbSubtract)*DeltaRho_new[EWSM::EW1]
                            - DeltaRho_new[EWSM::EW1]*DeltaRbar_rem
                            + DeltaRho_new[EWSM::EW2]
                            + deltaKappa_rem[EWSM::EW2];
        double dKappa_EW1QCD2 = DeltaRho_new[EWSM::EW1QCD2] + deltaKappa_rem[EWSM::EW1QCD2];
        double dKappa_EW2QCD1 = DeltaRho_new[EWSM::EW2QCD1] + deltaKappa_rem[EWSM::EW1QCD1]*DeltaRho_new[EWSM::EW1];
        double dKappa_EW3 = DeltaRho_new[EWSM::EW3];
        double dKappa_EW2QCD2 = deltaKappa_rem[EWSM::EW1QCD2]*DeltaRho_new[EWSM::EW1];
        //
        double dKappa_EW1_TMP = dKappa_EW1;
        double dKappa_EW1QCD1_TMP = dKappa_EW1QCD1;
        double dKappa_EW2_TMP = dKappa_EW2;
        double dKappa_EW1QCD2_TMP = dKappa_EW1QCD2;
        double dKappa_EW2QCD1_TMP = dKappa_EW2QCD1;
        double dKappa_EW3_TMP = dKappa_EW3;
        double dKappa_EW2QCD2_TMP = dKappa_EW2QCD2;
        //
        if (bool_Zbb) {
            dKappa_EW1 = dKappa_EW1_TMP - taub[EWSM::EW1];
            dKappa_EW1QCD1 = dKappa_EW1QCD1_TMP - taub[EWSM::EW1QCD1];
            dKappa_EW2 = dKappa_EW2_TMP - taub[EWSM::EW2] + taub[EWSM::EW1]*taub[EWSM::EW1]
                         - dKappa_EW1_TMP*taub[EWSM::EW1];
            dKappa_EW1QCD2 = dKappa_EW1QCD2_TMP;
            dKappa_EW2QCD1 = dKappa_EW2QCD1_TMP - dKappa_EW1_TMP*taub[EWSM::EW1QCD1]
                             - dKappa_EW1QCD1_TMP*taub[EWSM::EW1] + 2.0*taub[EWSM::EW1]*taub[EWSM::EW1QCD1];
            dKappa_EW3 = dKappa_EW3_TMP - dKappa_EW1_TMP*taub[EWSM::EW2]
                         - dKappa_EW2_TMP*taub[EWSM::EW1] + 2.0*taub[EWSM::EW2]*taub[EWSM::EW1];
            dKappa_EW2QCD2 = dKappa_EW2QCD2_TMP - dKappa_EW1QCD2_TMP*taub[EWSM::EW1]
                             - dKappa_EW1QCD1_TMP*taub[EWSM::EW1QCD1]
                             + taub[EWSM::EW1QCD1]*taub[EWSM::EW1QCD1];
        }
        std::cout << "  EW1: " << dKappa_EW1 << " " << f_AlphaToGF*dKappa_EW1 << std::endl;
        std::cout << "  EW1QCD1: " << dKappa_EW1QCD1 << " " << f_AlphaToGF*dKappa_EW1QCD1 << std::endl;
        std::cout << "  EW2: " << dKappa_EW2 + Zgamma_EW2
                  << " " << pow(f_AlphaToGF,2.0)*dKappa_EW2 + Zgamma_EW2 << std::endl;
        std::cout << "  EW1QCD2: " << dKappa_EW1QCD2 << " " << f_AlphaToGF*dKappa_EW1QCD2 << std::endl;
        std::cout << "  EW2QCD1: " << dKappa_EW2QCD1 << " " << pow(f_AlphaToGF,2.0)*dKappa_EW2QCD1 << std::endl;
        std::cout << "  EW3: " << dKappa_EW3 << " " << pow(f_AlphaToGF,3.0)*dKappa_EW3 << std::endl;
        std::cout << "  EW2QCD2: " << dKappa_EW2QCD2 << " " << pow(f_AlphaToGF,2.0)*dKappa_EW2QCD2 << std::endl;
        std::cout << "Total contribution: alpha or Gmu" << std::endl;
        std::cout << "  kappaZ="
                  << 1.0 + dKappa_EW1 + dKappa_EW1QCD1 + dKappa_EW2
                     + dKappa_EW1QCD2 + dKappa_EW2QCD1 + dKappa_EW3
                     + dKappa_EW2QCD2 + Zgamma_EW2
                  << " "
                  << 1.0 + f_AlphaToGF*dKappa_EW1 + f_AlphaToGF*dKappa_EW1QCD1
                     + pow(f_AlphaToGF,2.0)*dKappa_EW2 + f_AlphaToGF*dKappa_EW1QCD2
                     + pow(f_AlphaToGF,2.0)*dKappa_EW2QCD1
                     + pow(f_AlphaToGF,3.0)*dKappa_EW3
                     + pow(f_AlphaToGF,2.0)*dKappa_EW2QCD2
                     + Zgamma_EW2 << std::endl;
        if (bool_Zbb) {
            std::cout << "  kappaZ(taub resummed)="
                      << (1.0 + dKappa_EW1_TMP + dKappa_EW1QCD1_TMP
                          + dKappa_EW2_TMP + dKappa_EW1QCD2_TMP
                          + dKappa_EW2QCD1_TMP + dKappa_EW3_TMP
                          + dKappa_EW2QCD2_TMP)
                         /(1.0 + taub[EWSM::EW1] + taub[EWSM::EW1QCD1] + taub[EWSM::EW2])
                         + Zgamma_EW2
                      << " "
                      << (1.0 + f_AlphaToGF*dKappa_EW1_TMP
                          + f_AlphaToGF*dKappa_EW1QCD1_TMP
                          + pow(f_AlphaToGF,2.0)*dKappa_EW2_TMP
                          + f_AlphaToGF*dKappa_EW1QCD2_TMP
                          + pow(f_AlphaToGF,2.0)*dKappa_EW2QCD1_TMP
                          + pow(f_AlphaToGF,3.0)*dKappa_EW3_TMP
                          + pow(f_AlphaToGF,2.0)*dKappa_EW2QCD2_TMP)
                         /(1.0 + f_AlphaToGF*taub[EWSM::EW1]
                           + f_AlphaToGF*taub[EWSM::EW1QCD1]
                           + pow(f_AlphaToGF,2.0)*taub[EWSM::EW2])
                         + Zgamma_EW2
                      << std::endl;
        }
    } else
        std::cout << "EWSM_Output::outputEachDeltaKappaZ(): Not implemented for schemeKappaZ="
                  << myEWSM.schemeKappaZ << std::endl;
}



