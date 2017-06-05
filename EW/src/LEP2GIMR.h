/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2GIMR_H
#define	LEP2GIMR_H

#include <stdexcept>
#include <gslpp.h>
#include "StandardModel.h"


using namespace gslpp;

/**
 * @class LEP2GIMR
 * @ingroup EW
 * @brief A class for NP analyses of LEP-II observables with the dimension 6 operators in the GIMR basis. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2GIMR {
public:
    
    enum Param {C_LL=0, C_LR, C_RL, C_RR, delta_GammaZ, delta_gLf, delta_gRf, delta_gLe, delta_gRe, delta_Mz2};
    
    /**
     * @brief LEP2GIMR constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    LEP2GIMR(const StandardModel& SM_i);
    
    double sigma_l_LEP2_GIMR(const QCD::lepton l, const double s,
                             const double GIMRParam_i[]) const;
    double sigma_q_LEP2_GIMR(const QCD::quark q, const double s,
                             const double GIMRParam_i[]) const;

    double sigmaFminusB_l_LEP2_GIMR(const QCD::lepton l, const double s,
                             const double GIMRParam_i[]) const;
    double sigmaFminusB_q_LEP2_GIMR(const QCD::quark q, const double s,
                             const double GIMRParam_i[]) const; 
    
    private:
    const StandardModel& SM;
    
    double gL_l(const QCD::lepton l) const;
    double gR_l(const QCD::lepton l) const;

    double gL_q(const QCD::quark q) const;
    double gR_q(const QCD::quark q) const;
    
    double deltaA1q(const QCD::quark q, const double GIMRParam_i[]) const;    
    double deltaA2q(const QCD::quark q, const double GIMRParam_i[]) const;    
    double deltaB1q(const QCD::quark q, const double GIMRParam_i[]) const;    
    double deltaB2q(const QCD::quark q, const double GIMRParam_i[]) const; 
  
    double deltaA1l(const QCD::lepton l, const double GIMRParam_i[]) const;    
    double deltaA2l(const QCD::lepton l, const double GIMRParam_i[]) const;    
    double deltaB1l(const QCD::lepton l, const double GIMRParam_i[]) const;    
    double deltaB2l(const QCD::lepton l, const double GIMRParam_i[]) const; 
    
};

#endif	/* LEP2GIMR_H */

