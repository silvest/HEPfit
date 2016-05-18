/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LI_LJ_GAMMA_H
#define	LI_LJ_GAMMA_H

#include "gslpp.h"
#include "ThObservable.h"
#include "LeptonFlavour.h"

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
 * @class log_meg
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \mu \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The log_meg class calculates the branching ratio of the process \f$ \mu \to e \gamma \f$ in the model.
 */
class log_meg : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \mu \to e \gamma \f$.
     * @return value of the process \f$ \mu \to e \gamma \f$.
     */
    log_meg(const StandardModel& SM_i);
    
    /**
     * @return decadic log of the process \f$ \mu \to e \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class Robs_mu_e_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the upper limit ratio for the process \f$ \mu \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class Robs_mu_e_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the upper limit ratio of the process \f$ \mu \to e \gamma \f$.
     */
    Robs_mu_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return @f$R=1+(BR_{\text{theo}}-BR_{\text{obs,90\%}})/BR_{\text{exp,90\%}}\cdot \sqrt{2}\text{InvErf}(0.90)@f$
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
 * @class log_tmg
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to \mu \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The log_tmg class calculates the branching ratio of the process \f$ \tau \to \mu \gamma \f$ in the model.
 */
class log_tmg : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to \mu \gamma \f$.
     * @return value of the process \f$ \tau \to \mu \gamma \f$.
     */
    log_tmg(const StandardModel& SM_i);
    
    /**
     * @return decadic log of the process \f$ \tau \to \mu \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class Robs_tau_mu_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the upper limit ratio for the process \f$ \tau \to \mu \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class Robs_tau_mu_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the upper limit ratio of the process \f$ \tau \to \mu \gamma \f$.
     */
    Robs_tau_mu_gamma(const StandardModel& SM_i);
    
    /**
     * @return @f$R=1+(BR_{\text{theo}}-BR_{\text{obs,90\%}})/BR_{\text{exp,90\%}}\cdot \sqrt{2}\text{InvErf}(0.90)@f$
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

/**
 * @class log_teg
 * @ingroup LeptonFlavour
 * @brief A class for calculating the branching ratio of the process \f$ \tau \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The log_teg class calculates the branching ratio of the process \f$ \tau \to e \gamma \f$ in the model.
 */
class log_teg : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \tau \to e \gamma \f$.
     * @return value of the process \f$ \tau \to e \gamma \f$.
     */
    log_teg(const StandardModel& SM_i);
    
    /**
     * @return decadic log of the process \f$ \tau \to e \gamma \f$.
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

/**
 * @class Robs_tau_e_gamma
 * @ingroup LeptonFlavour
 * @brief A class for calculating the upper limit ratio for the process \f$ \tau \to e \gamma \f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class Robs_tau_e_gamma : public ThObservable {
public:
    
    /**
     * @brief Calculates the upper limit ratio of the process \f$ \tau \to e \gamma \f$.
     */
    Robs_tau_e_gamma(const StandardModel& SM_i);
    
    /**
     * @return @f$R=1+(BR_{\text{theo}}-BR_{\text{obs,90\%}})/BR_{\text{exp,90\%}}\cdot \sqrt{2}\text{InvErf}(0.90)@f$
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};


#endif	/* LI_LJ_GAMMA_H */
