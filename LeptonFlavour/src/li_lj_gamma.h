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

//class li_lj_gamma : public ThObservable {
//public:
//    /**
//     * @brief Calculates the value of the process \f$ \ell_j \to \ell_i \gamma \f$.
//     * @return value of the process \f$ \ell_j \to \ell_i \gamma \f$.
//     */
//    li_lj_gamma(const StandardModel& SM_i);
//
//    /**
//     *
//     * @return value of the process \f$ \ell_j \to \ell_i \gamma \f$.
//     */
//    double computeThValue();
//
//protected:
//
//private:
//    
//};

/**
 * @class mu_e_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \mu \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The mu_e_gamma class calculates the branching ratio of the process \f$ \mu \to e \gamma \f$ in the model.
 */
class mu_e_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \mu \to e \gamma \f$.
     * @return value of the process \f$ \mu \to e \gamma \f$.
     */
    mu_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \mu \to e \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class tau_mu_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to \mu \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The tau_mu_gamma class calculates the branching ratio of the process \f$ \tau \to \mu \gamma \f$ in the model.
 */
class tau_mu_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to \mu \gamma \f$.
     * @return value of the process \f$ \tau \to \mu \gamma \f$.
     */
    tau_mu_gamma(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \tau \to \mu \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class tau_e_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The tau_e_gamma class calculates the branching ratio of the process \f$ \tau \to e \gamma \f$ in the model.
 */
class tau_e_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to e \gamma \f$.
     * @return value of the process \f$ \tau \to e \gamma \f$.
     */
    tau_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return value of the process \f$ \tau \to e \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};


#endif	/* LI_LJ_GAMMA_H */
