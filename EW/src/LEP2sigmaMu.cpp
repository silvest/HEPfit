/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaMu.h"


double LEP2sigmaMu::getThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        ml_cache = m_l(SM.MU);
        
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_l();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaMu::Integrand_sigmaWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            /* check
             MathMore library is required to use kADAPTIVESINGULAR. */
            //std::cout << ig.Options().Integrator() << std::endl;
            //std::cout << ig.Options().AbsTolerance() << std::endl;
            //std::cout << ig.Options().RelTolerance() << std::endl;
            //std::cout << SMresult_cache << " +- " << ig.Error() << std::endl;
            // results: 6.8e-9 -- 2.2e-8
        }
        
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaMu::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-17); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            SMresult_cache += sigma_box;
            //std::cout << sigma_box << " +- " << ig_box.Error() << std::endl;
            // results: 3.5e-12 -- +2.9e-10
        }

        if ( myEW.checkSTUVWXY() && !bSigmaForAFB && SM.IsFlagFixedAllSMparams() ) {
            //std::cout << "sqrt_s=" << sqrt_s << " in LEP2sigmaMu" << std::endl; // TEST
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.sigma_l_LEP2_NP(StandardModel::MU, s, ml_cache, ObParam);
            }
        }
    }
    double sigma_mu = SMresult_cache;
    
    #ifdef LEP2TEST
    sigma_mu = myTEST.sigmaMuTEST(sqrt_s)/GeVminus2_to_nb/1000.0;
    #endif 
    
    if ( myEW.checkSTUVWXY() && !bSigmaForAFB ) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            sigma_mu += Coeff_cache[myLEP2oblique.Shat]*myEW.Shat()
                      + Coeff_cache[myLEP2oblique.That]*myEW.That()
                      + Coeff_cache[myLEP2oblique.Uhat]*myEW.Uhat()
                      + Coeff_cache[myLEP2oblique.V]*myEW.V()
                      + Coeff_cache[myLEP2oblique.W]*myEW.W()
                      + Coeff_cache[myLEP2oblique.X]*myEW.X()
                      + Coeff_cache[myLEP2oblique.Y]*myEW.Y();
        } else {
            double ObParam[7] = {myEW.Shat(), myEW.That(), myEW.Uhat(),
                                 myEW.V(), myEW.W(), myEW.X(), myEW.Y()};
            sigma_mu += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::MU, s, ml_cache, ObParam);
        }
    }

    return ( sigma_mu*GeVminus2_to_nb*1000.0 );
}
        

