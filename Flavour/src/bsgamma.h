/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMA_H
#define	BSGAMMA_H

#include <ThObservable.h>
#include "Flavour.h"
#include <StandardModel.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>

/**
 * @class Bsgamma
 * @ingroup flavour
 * @brief A class for the decay b -> s gamma. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Bsgamma : public ThObservable {
public: 
    /**
    * @brief constructor
    * @param[in] SM_i model
    * @param[in] obsFlag flag to choose which observable to compute
    */
    Bsgamma(const StandardModel& SM_i, int obsFlag);
    
    
    /**
    * @brief \f$ dY1 \f$
    * @param[in] E0 energy cutoff
    * @return return the function \delta Y^1 from arXiv:0805.3911
    */
    double dY1(double E0);
    
    
    /**
    * @brief \f$ Y1 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(1) from arXiv:0805.3911
    */
    double Y1(double E0, double Mu);
    
    
    /**
    * @brief \f$ Phi1 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_1 from arXiv:1005.5587
    */
    double Phi1(double rho);
    
    
    /**
    * @brief \f$ Phi2 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_2 from arXiv:1005.5587
    */
    double Phi2(double rho);
    
    
    /**
    * @brief \f$ Phi3 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_3 from arXiv:1005.5587
    */
    double Phi3(double rho);
    
    
    /**
    * @brief \f$ Phi4 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_4 from arXiv:1005.5587
    */
    double Phi4(double rho);
    
    
    /**
    * @brief \f$ Y2CF \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,CF) from arXiv:1005.5587
    */
    double Y2CF(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2CA \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,CA) from arXiv:1005.5587
    */
    double Y2CA(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NH \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NH) from arXiv:0805.3911
    */
    double Y2NH(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NV \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NV) from arXiv:0805.3911
    */
    double Y2NV(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NL \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NL) from arXiv:0805.3911
    */
    double Y2NL(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2) from arXiv:0805.3911 and arXiv:1005.5587
    */
    double Y2(double E0, double Mu);
    
    
    /**
    * @brief \f$ G78 \f$
    * @param[in] order perturbation theory order
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function G_{78} from arXiv:0805.3911 and arXiv:1005.5587
    */
    double G78(orders order, double E0, double Mu);
    
    
    /**
    * @brief \f$ computeCoeff \f$
    * @param[in] order perturbation theory order
    * @return compute the Wilson Coefficient part of the BR
    */
    void computeCoeff(orders order);
    
    
    /**
    * @brief \f$ computeGamma \f$
    * @param[in] order perturbation theory order
    * @return compute the decay rate \Gamma
    */
    void computeGamma(orders order);
    
    
    /**
    * @brief \f$ computeThValue \f$
    * @return compute the BR and A_{CP}
    */
    double computeThValue();
    
    
private:
    const StandardModel& mySM;
    
    double GF;
    double ale;
    double Mb;
    double Mc;
    double mu_b;
    double width;
    gslpp::complex lambda_t;
    double E0;
    
    double CF;
    double CA;
    double TR;
    double NL;
    double NH;
    double NV;
    
    double coeff;
    gslpp::vector<gslpp::complex> ** allcoeff;
    gslpp::complex C_2;
    gslpp::complex C_7;
    gslpp::complex C_8;
    int obs;
    
    double Gamma;
    double Gamma_conj;
};

#endif	/* BSGAMMA_H */