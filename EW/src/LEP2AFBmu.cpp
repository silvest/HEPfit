/* 
 * File:   LEP2AFBmu.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBmu.h"


double LEP2AFBmu::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        double AFB_noBox, sigma;
        if (!flag[ISR])
            AFB_noBox = AFB_NoISR_l();
        else {
            // numerator
            ROOT::Math::Functor1D wf(this, &LEP2AFBmu::Integrand_AFBnumeratorWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-5); // desired relative error
            double numerator = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "numerator = " << numerator << std::endl;

            // denominator
            myLEP2sigmaMu.setFlags(flag);
            sigma = myLEP2sigmaMu.getThValue()/GeVminus2_to_nb/1000.0;
            
            AFB_noBox = numerator/sigma;
        }    
        SMresult_cache = AFB_noBox;
        
        if (flag[WeakBox]) {
            // numerator
            ROOT::Math::Functor1D wf_F(this, &LEP2AFBmu::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_F(wf_F, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_F.SetAbsTolerance(1.E-16); // desired absolute error
            ig_F.SetRelTolerance(1.E-5); // desired relative error
            double sigma_box_F = ig_F.Integral(0.0, 1.0); // interval
            //std::cout << "sigma_box_F = " << sigma_box_F << std::endl;
            //
            ROOT::Math::Functor1D wf_B(this, &LEP2AFBmu::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_B(wf_B, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_B.SetAbsTolerance(1.E-16); // desired absolute error
            ig_B.SetRelTolerance(1.E-5); // desired relative error
            double sigma_box_B = ig_B.Integral(-1.0, 0.0); // interval
            //std::cout << "sigma_box_B = " << sigma_box_B << std::endl;
                        
            // denominator
            if (!flag[ISR]) {
                myLEP2sigmaMu.setFlags(flag);
                sigma = myLEP2sigmaMu.getThValue()/GeVminus2_to_nb/1000.0;
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    }
    double AFB_mu = SMresult_cache;
    
    if ( myEW.checkModelForSTU() )
        AFB_mu += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::MU, s);
    
    return AFB_mu;
}
        
