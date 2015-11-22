/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LI_LJ_GAMMA_H
#define	LI_LJ_GAMMA_H

#include "gslpp.h"
#include "ThObservable.h"
#include "LeptonFlavour.h"

class li_lj_gamma : public ThObservable {
public:
    /**
     * constructor
     * @param LeptonFlavour
     */
    li_lj_gamma(const StandardModel& SM_i);

    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();

protected:

private:
    
};

class mu_e_gamma : public li_lj_gamma {
public:
    
    /**
     * @brief Constructor.
     */
    mu_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return mu_e_gamma
     */
    double computeThValue();
    
private:
    const StandardModel& mySM;

};

class tau_mu_gamma : public li_lj_gamma {
public:
    
    /**
     * @brief Constructor.
     */
    tau_mu_gamma(const StandardModel& SM_i);
    
    /**
     * @return tau_mu_gamma
     */
    double computeThValue();
    
private:
    const StandardModel& mySM;

};

class tau_e_gamma : public li_lj_gamma {
public:
    
    /**
     * @brief Constructor.
     */
    tau_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return tau_e_gamma
     */
    double computeThValue();
    
private:
    const StandardModel& mySM;

};


#endif	/* LI_LJ_GAMMA_H */
