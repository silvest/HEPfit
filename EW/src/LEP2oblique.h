/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2OBLIQUE_H
#define	LEP2OBLIQUE_H

#include <stdexcept>
#include <gslpp.h>
#include "EW.h"
using namespace gslpp;


/**
 * @class LEP2oblique
 * @brief a class for NP analyses of LEP-II observables with the extended oblique parameters
 */
class LEP2oblique {
public:
    
    /**
     * @brief LEP2oblique constructor
     * @param[in] EW_i an object of EW class
     */
    LEP2oblique(const EW& EW_i);

    
    ////////////////////////////////////////////////////////////////////////
    
    double sigma_l_LEP2_NP(const StandardModel::lepton l, const double s) const;
    double sigma_q_LEP2_NP(const StandardModel::quark q, const double s) const;

    double AFB_l_LEP2_NP(const StandardModel::lepton l, const double s) const;
    double AFB_q_LEP2_NP(const StandardModel::quark q, const double s) const;
    
    double R_q_LEP2_NP(const StandardModel::quark q, const double s) const;
        
    
    ////////////////////////////////////////////////////////////////////////    
private:
    const EW& myEW;
   
    double DeltaEpsilon_1() const;
    double DeltaEpsilon_2() const;    
    double DeltaEpsilon_3() const;
    
    double epsilonZZ() const;
    double epsilonGammaGamma() const;
    double epsilonGammaZ() const;
    
    double ml(const StandardModel::lepton l) const {
        return myEW.getSM().getLeptons(l).getMass();
    }

    double mq(const StandardModel::quark q, const double mu, const orders order=FULLNLO) const {
        switch(q) {
            case StandardModel::UP:
            case StandardModel::DOWN:
            case StandardModel::STRANGE:
                return myEW.getSM().Mrun(mu, myEW.getSM().getQuarks(q).getMass_scale(), 
                                         myEW.getSM().getQuarks(q).getMass(), order);
            case StandardModel::CHARM:
            case StandardModel::BOTTOM:
                return myEW.getSM().Mrun(mu, myEW.getSM().getQuarks(q).getMass(), order);
            case StandardModel::TOP:
                return myEW.getSM().getMtpole(); // the pole mass
            default:
                throw std::runtime_error("Error in LEP2oblique::mq()"); 
        }
    }
    
    double vl(const StandardModel::lepton l) const;
    double vq(const StandardModel::quark q) const;
    double al(const StandardModel::lepton l) const;
    double aq(const StandardModel::quark q) const;

    double G1_l_NP(const StandardModel::lepton l, const double s) const;
    double G1_q_NP(const StandardModel::quark q, const double s) const;
    double G1_NP(const double s, const double Qf, 
                 const double vf, const double af) const;

    double G3_l_NP(const StandardModel::lepton l, const double s) const;
    double G3_q_NP(const StandardModel::quark q, const double s) const;
    double G3_NP(const double s, const double Qf, 
                 const double vf, const double af) const;
    
    double G1_l_SM0(const StandardModel::lepton l, const double s) const;
    double G1_q_SM0(const StandardModel::quark q, const double s) const;
    double G1_SM0(const double s, const double Qf, 
                  const double vf, const double af) const;

    double G2_l_SM0(const StandardModel::lepton l, const double s) const;
    double G2_q_SM0(const StandardModel::quark q, const double s) const;
    double G2_SM0(const double s, const double Qf, 
                  const double vf, const double af) const;

    double G3_l_SM0(const StandardModel::lepton l, const double s) const;
    double G3_q_SM0(const StandardModel::quark q, const double s) const;    
    double G3_SM0(const double s, const double Qf, 
                  const double vf, const double af) const;
    
    double sigma_l_LEP2_SM0(const StandardModel::lepton l, const double s) const;
    double sigma_q_LEP2_SM0(const StandardModel::quark q, const double s) const;
    
};

#endif	/* LEP2OBLIQUE_H */

