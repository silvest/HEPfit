/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"
#include <EWSM.h>


double sigmaHadron::computeThValue() 
{ 
    double sigma_had;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        sigma_had = myEW.getMyEW_CHMN().sigma0_had();
    else if (myEWTYPE==EW::EWABC) 
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),false);
    else if (myEWTYPE==EW::EWABC2)
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb(),true)/GeVminus2_to_nb;
    else {   
        if (SM.IsFlagApproximateSigmaH())
            sigma_had = SM.getEWSM()->sigmaHadron_SM()/GeVminus2_to_nb;
        else
            sigma_had = myEW.sigma0_had();
        
        if(myEWTYPE==EW::EWBURGESS) {
            sigma_had = myEW.getMyEW_BURGESS().sigmaHadron(sigma_had,
                    myEW.Gamma_Z(), myEW.Gamma_had(), myEW.Gamma_l(SM.ELECTRON));
            return ( sigma_had*GeVminus2_to_nb );
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            bool nonZeroNP = false;

            double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
            for (int p=0; p<6; ++p) {
                delGVl[p] = SM.deltaGVl((StandardModel::lepton)p);
                delGAl[p] = SM.deltaGAl((StandardModel::lepton)p);
                delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
                delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
                if (delGVl[p]!=0.0 || delGAl[p]!=0.0
                        || delGVq[p]!=0.0 || delGAq[p]!=0.0)
                    nonZeroNP = true;
            }

            if (nonZeroNP) {
                double gVf, gAf;
                double Gl[6], deltaGl[6], Gq[6], deltaGq[6];
                double Gq_sum = 0.0, delGq_sum = 0.0;
                double Gf_sum = 0.0, delGf_sum = 0.0;
                for (int p=0; p<6; ++p) {
                    gVf = SM.StandardModel::gVl((StandardModel::lepton)p).real();
                    gAf = SM.StandardModel::gAl((StandardModel::lepton)p).real();
                    Gl[p] = gVf*gVf + gAf*gAf;
                    deltaGl[p] = 2.0*(gVf*delGVl[p] + gAf*delGAl[p]);

                    gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
                    gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
                    Gq[p] = gVf*gVf + gAf*gAf;
                    deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

                    Gq_sum += 3.0*Gq[p];
                    Gf_sum += Gl[p] + 3.0*Gq[p];
                    delGq_sum += 3.0*deltaGq[p];
                    delGf_sum += deltaGl[p] + 3.0*deltaGq[p];
                }

                sigma_had += 12.0*M_PI/SM.getMz()/SM.getMz()
                             *Gl[(int)SM.ELECTRON]*Gq_sum/Gf_sum/Gf_sum
                             *( deltaGl[(int)SM.ELECTRON]/Gl[(int)SM.ELECTRON]
                                + delGq_sum/Gq_sum - 2.0*delGf_sum/Gf_sum );
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("sigmaHadron::computeThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //sigma_had -= myEW.sigma0_had();
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


