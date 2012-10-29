/* 
 * File:   LEP2sigmaMu.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaMu.h"


double LEP2sigmaMu::getThValue() { 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        ml_cache = m_l(SM.MU);
        
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_l();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaMu::Integrand_sigmaWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-15); // desired absolute error
            ig.SetRelTolerance(1.E-6); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << std::endl;
        }
        
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaMu::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-14); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            //std::cout << sigma_box << std::endl;
            SMresult_cache += sigma_box;
        }

        if ( myEW.checkModelForSTU() && !bSigmaForAFB && SM.FixedSMparams() ) {
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
    
    if ( myEW.checkModelForSTU() && !bSigmaForAFB ) {
        if ( SM.FixedSMparams() ) {
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
        

