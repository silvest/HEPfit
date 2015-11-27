/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2OBLIQUE_H
#define	LEP2OBLIQUE_H

#include <stdexcept>
#include "StandardModel.h"

/**
 * @class LEP2oblique
 * @ingroup EW
 * @brief A class for NP analyses of LEP-II observables with the extended oblique parameters. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2oblique {
public:

    enum Oblique {Shat=0, That, Uhat, V, W, X, Y};
    
    /**
     * @brief LEP2oblique constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    LEP2oblique(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////
    
    double sigma_l_LEP2_NP(const StandardModel::lepton l, const double s,
                           const double ml, const double ObParam_i[]) const;
    double sigma_q_LEP2_NP(const QCD::quark q, const double s,
                           const double mq, const double ObParam_i[]) const;

    double AFB_l_LEP2_NP(const StandardModel::lepton l, const double s, 
                         const double ml, const double ObParam_i[]) const;
    double AFB_q_LEP2_NP(const QCD::quark q, const double s, 
                         const double mq, const double ObParam_i[]) const;
    
    double R_q_LEP2_NP(const QCD::quark q, const double s, 
                       const double mq, const double ObParam_i[]) const;
        
    
    ////////////////////////////////////////////////////////////////////////    
private:
    const StandardModel& SM;
   
    double DeltaEpsilon_1(const double alpha0, const double ObParam_i[]) const;
    double DeltaEpsilon_2(const double alpha0, const double ObParam_i[]) const;    
    double DeltaEpsilon_3(const double alpha0, const double ObParam_i[]) const;
    
    double epsilonZZ(const double alpha0, const double ObParam_i[]) const;
    double epsilonGammaGamma(const double alpha0, const double ObParam_i[]) const;
    double epsilonGammaZ(const double alpha0, const double ObParam_i[]) const;
    
    double alpha_at_s(const double s) const 
    {
        //return ( SM.alphaMz() ); // TEST
        return ( SM.ale_OS(sqrt(s), FULLNLO) );
    }
    
    double Mw0(const double alpha0) const
    {
        double Mz = SM.getMz();
        return ( sqrt(c02(alpha0))*Mz );
    }
    
    double s02(const double alpha0) const 
    {
        
        throw std::runtime_error("Error in LEP2oblique::s02() needs to be modified with sW2 in the SM!");
        
        double GF = SM.getGF();
        double Mz = SM.getMz();
        return ( ( 1.0 - sqrt(1.0 - 4.0*M_PI*alpha0/sqrt(2.0)/GF/Mz/Mz) )/2.0 );
    }
    
    double c02(const double alpha0) const
    {
        return ( 1.0 - s02(alpha0) );
    }
    
    double vl(const StandardModel::lepton l, const double alpha0) const;
    double vq(const QCD::quark q, const double alpha0) const;
    double al(const StandardModel::lepton l, const double alpha0) const;
    double aq(const QCD::quark q, const double alpha0) const;

    double G1_l_NP(const StandardModel::lepton l, const double s, 
                   const double alpha0, const double ObParam_i[]) const;
    double G1_q_NP(const QCD::quark q, const double s, 
                   const double alpha0, const double ObParam_i[]) const;
    double G1_NP(const double s, const double alpha0, const double Qf, 
                 const double vf, const double af, const double ObParam_i[]) const;

    double G3_l_NP(const StandardModel::lepton l, const double s, 
                   const double alpha0, const double ObParam_i[]) const;
    double G3_q_NP(const QCD::quark q, const double s, 
                   const double alpha0, const double ObParam_i[]) const;
    double G3_NP(const double s, const double alpha0, const double Qf, 
                 const double vf, const double af, const double ObParam_i[]) const;
    
    double G1_l_SM0(const StandardModel::lepton l, const double s, 
                    const double alpha0) const;
    double G1_q_SM0(const QCD::quark q, const double s, 
                    const double alpha0) const;
    double G1_SM0(const double s, const double alpha0, const double Qf, 
                  const double vf, const double af) const;

    double G2_l_SM0(const StandardModel::lepton l, const double s, 
                    const double alpha0) const;
    double G2_q_SM0(const QCD::quark q, const double s, 
                    const double alpha0) const;
    double G2_SM0(const double s, const double alpha0, const double Qf, 
                  const double vf, const double af) const;

    double G3_l_SM0(const StandardModel::lepton l, const double s, 
                    const double alpha0) const;
    double G3_q_SM0(const QCD::quark q, const double s, 
                    const double alpha0) const;    
    double G3_SM0(const double s, const double alpha0, const double Qf, 
                  const double vf, const double af) const;
    
    double sigma_l_LEP2_SM0(const StandardModel::lepton l, const double s, 
                            const double alpha0, const double ml) const;
    double sigma_q_LEP2_SM0(const QCD::quark q, const double s, 
                            const double alpha0, const double mq) const;
    
};

#endif	/* LEP2OBLIQUE_H */

