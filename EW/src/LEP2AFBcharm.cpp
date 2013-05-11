/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBcharm.h"


double LEP2AFBcharm::getThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {        
        mq_cache = m_q(SM.CHARM, sqrt_s);

        double AFB_noBox, sigma;
        if (!flag[ISR])
            AFB_noBox = AFB_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2AFBcharm::Integrand_AFBnumeratorWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            double numerator = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << numerator << " +- " << ig.Error() << std::endl;
            // results: 7.5e-9 -- 3.0e-8
            
            // denominator
            myLEP2sigmaCharm.setFlags(flag);
            sigma = myLEP2sigmaCharm.getThValue()/GeVminus2_to_nb/1000.0;
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flag[WeakBox]) {
            // numerator
            ROOT::Math::Functor1D wf_F(this, &LEP2AFBcharm::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_F(wf_F, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_F.SetAbsTolerance(1.E-16); // desired absolute error
            ig_F.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_F = ig_F.Integral(0.0, 1.0); // interval
            //std::cout << sigma_box_F << " +- " << ig_F.Error()  << std::endl;
            // results: 8.6e-11 -- 3.1e-10
            //
            ROOT::Math::Functor1D wf_B(this, &LEP2AFBcharm::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_B(wf_B, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_B.SetAbsTolerance(1.E-16); // desired absolute error
            ig_B.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_B = ig_B.Integral(-1.0, 0.0); // interval
            //std::cout << sigma_box_B << " +- " << ig_B.Error()  << std::endl;
            // results: 5.0e-11 -- 8.5e-11
                        
            // denominator
            if (!flag[ISR]) {
                myLEP2sigmaCharm.setFlags(flag);
                sigma = myLEP2sigmaCharm.getThValue()/GeVminus2_to_nb/1000.0;
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }

        if ( myEW.checkSTUVWXY() && SM.IsFlagFixedAllSMparams() ){
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.AFB_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
            }
        }
    }
    double AFB_c = SMresult_cache;
    
    #ifdef LEP2TEST
    AFB_c = myTEST.AFBcharmTEST(sqrt_s);
    #endif
        
    if ( myEW.checkSTUVWXY() ){
        if ( SM.IsFlagFixedAllSMparams() ) {
            AFB_c += Coeff_cache[myLEP2oblique.Shat]*SM.obliqueShat()
                   + Coeff_cache[myLEP2oblique.That]*SM.obliqueThat()
                   + Coeff_cache[myLEP2oblique.Uhat]*SM.obliqueUhat()
                   + Coeff_cache[myLEP2oblique.V]*SM.obliqueV()
                   + Coeff_cache[myLEP2oblique.W]*SM.obliqueW()
                   + Coeff_cache[myLEP2oblique.X]*SM.obliqueX()
                   + Coeff_cache[myLEP2oblique.Y]*SM.obliqueY();
        } else {
            double ObParam[7] = {SM.obliqueShat(), SM.obliqueThat(), SM.obliqueUhat(),
                                 SM.obliqueV(), SM.obliqueW(), SM.obliqueX(), SM.obliqueY()};
            AFB_c += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
        }
    }
    
    return AFB_c;
}
        


