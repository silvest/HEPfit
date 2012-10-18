/* 
 * File:   LEP2sigmaMu.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaMu.h"


double LEP2sigmaMu::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        if (!bRCs[LEP2TwoFermions::ISR])
            SMresult_cache = myTwoFermions.sigma_l(StandardModel::MU, s, Mw, GammaZ, bRCs);
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaMu::IntegrandISR_sigma_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-15); // desired absolute error
            ig.SetRelTolerance(1.E-6); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << std::endl;
        }
    }
    double sigma_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() )
        sigma_mu += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::MU, s);
    
    return ( sigma_mu*GeVminus2_to_nb*1000.0 );
}
        

