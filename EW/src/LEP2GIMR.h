/* 
 * File:   LEP2GIMR.h
 * Author: ggrillidc
 *
 * Created on September 28, 2015, 4:23 PM
 */

#ifndef LEP2GIMR_H
#define	LEP2GIMR_H

#include <stdexcept>
#include <gslpp.h>
#include <NPEffectiveGIMR.h>


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
    
    double sigma_l_LEP2_GIMR(const StandardModel::lepton l, const double s,
                             const double GIMRParam_i[]) const;
    double sigma_q_LEP2_GIMR(const QCD::quark q, const double s,
                             const double GIMRParam_i[]) const;

    double AFB_l_LEP2_GIMR(const StandardModel::lepton l, const double s, 
                           const double GIMRParam_i[]) const;
    double AFB_q_LEP2_GIMR(const QCD::quark q, const double s, 
                           const double GIMRParam_i[]) const;

    double sigmaF_l_LEP2_GIMR(const StandardModel::lepton l, const double s,
                             const double GIMRParam_i[]) const;
    double sigmaB_l_LEP2_GIMR(const StandardModel::lepton l, const double s,
                             const double GIMRParam_i[]) const;
    double sigmaF_q_LEP2_GIMR(const QCD::quark q, const double s,
                             const double GIMRParam_i[]) const;
    double sigmaB_q_LEP2_GIMR(const QCD::quark q, const double s,
                             const double GIMRParam_i[]) const;
    
    
    
    private:
    const StandardModel& SM;
    //
    
    
    
    double gL_l(const StandardModel::lepton l) const;
    double gR_l(const StandardModel::lepton l) const;

    double gL_q(const QCD::quark q) const;
    double gR_q(const QCD::quark q) const;

    
};

#endif	/* LEP2GIMR_H */

