/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LI_3LJ_H
#define	LI_3LJ_H

#include "gslpp.h"
#include "ThObservable.h"
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
     * @brief Calculates the value of the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @return value of the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$
     */
    double computeThValue();

protected:

private:
    
};

/**
 * @class mu_3e
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \mu \to eee \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The mu_3e class calculates the branching ratio of the process \f$ \mu \to eee \f$ in the model.
 */
class mu_3e : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \mu \to eee \f$.
     * @return value of the process \f$ \mu \to eee \f$
     */
    mu_3e(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \mu \to eee \f$
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class tau_3mu
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to \mu \mu \mu \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The tau_3mu class calculates the branching ratio of the process \f$ \tau \to \mu \mu \mu \f$ in the model.
 */
class tau_3mu : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to \mu \mu \mu \f$.
     * @return value of the process \f$ \tau \to \mu \mu \mu \f$
     */
    tau_3mu(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \tau \to \mu \mu \mu \f$
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class tau_3e
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to eee \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The tau_3e class calculates the branching ratio of the process \f$ \tau \to eee \f$ in the model.
 */
class tau_3e : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to eee \f$.
     * @return value of the process \f$ \tau \to eee \f$
     */
    tau_3e(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \tau \to eee \f$
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};


#endif	/* LI_3LJ_H */
