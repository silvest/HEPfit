/* 
 * File:   LEP2AFBcharm.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBcharm.h"


double LEP2AFBcharm::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        double AFB_noBox, sigma;
        if (!flag[ISR])
            AFB_noBox = AFB_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2AFBcharm::Integrand_AFBnumeratorWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-5); // desired relative error
            double numerator = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "numerator = " << numerator << std::endl;
            
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
            ig_F.SetAbsTolerance(1.E-14); // desired absolute error
            ig_F.SetRelTolerance(1.E-5); // desired relative error
            double sigma_box_F = ig_F.Integral(0.0, 1.0); // interval
            //std::cout << "sigma_box_F = " << sigma_box_F << std::endl;
            //
            ROOT::Math::Functor1D wf_B(this, &LEP2AFBcharm::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_B(wf_B, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_B.SetAbsTolerance(1.E-15); // desired absolute error
            ig_B.SetRelTolerance(1.E-5); // desired relative error
            double sigma_box_B = ig_B.Integral(-1.0, 0.0); // interval
            //std::cout << "sigma_box_B = " << sigma_box_B << std::endl;
                        
            // denominator
            if (!flag[ISR]) {
                myLEP2sigmaCharm.setFlags(flag);
                sigma = myLEP2sigmaCharm.getThValue()/GeVminus2_to_nb/1000.0;
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    }
    double AFB_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() )
        AFB_c += myLEP2oblique.AFB_q_LEP2_NP(StandardModel::CHARM, s);
    
    return AFB_c;
}
        


