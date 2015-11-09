/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LI_3LJ_H
#define	LI_3LJ_H

#include <gslpp.h>
#include <ThObservable.h>
#include "LeptonFlavour.h"

class li_3lj : public ThObservable {
public:
    /**
     * constructor
     * @param LeptonFlavour
     */
    li_3lj(const StandardModel& SM_i);

    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();

protected:

private:
    
};

class mu_3e : public li_3lj {
public:
    
    /**
     * @brief Constructor.
     */
    mu_3e(const StandardModel& SM_i);
    
    /**
     * @return mu_3e
     */
    double computeThValue ();
    
private:
    const StandardModel& mySM;

};

class tau_3mu : public li_3lj {
public:
    
    /**
     * @brief Constructor.
     */
    tau_3mu(const StandardModel& SM_i);
    
    /**
     * @return tau_3mu
     */
    double computeThValue ();
    
private:
    const StandardModel& mySM;

};

class tau_3e : public li_3lj {
public:
    
    /**
     * @brief Constructor.
     */
    tau_3e(const StandardModel& SM_i);
    
    /**
     * @return tau_3e
     */
    double computeThValue ();
    
private:
    const StandardModel& mySM;

};


#endif	/* LI_3LJ_H */
