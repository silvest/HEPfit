/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMAELECTRON_H
#define	LEP2SIGMAELECTRON_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2dsigmadcosElectron
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to e^+e^-@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2dsigmadcosElectron : public LEP2ThDiffObservable {
public:

    /**
     * @brief LEP2dsigmadcosElectron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] cos_i the polar angle of the final state particle wrt e^-
     */
    LEP2dsigmadcosElectron(const StandardModel& SM_i, const double sqrt_s_i, const double cos_i) 
    : LEP2ThDiffObservable(SM_i, sqrt_s_i, cos_i) 
    {
    }

    /**
     * @return the cross section for e^+ e^- -> e^+ e^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffsigmaElectron
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to e^+e^-@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmaElectron : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaElectron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmaElectron(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> e^+e^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffRelectron
 * @ingroup EW
 * @brief A class for the ratio of cross sections of @f$R_e=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to e^+e^-)@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffRelectron : public eeffThObservable {
public:

    /**
     * @brief eeffRelectron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffRelectron(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the ratio @f$R_e=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to e^+e^-)@f$ at sqrt_s
     */
    double computeThValue();

private:
    
};


// t-channel subtracted version 

/**
 * @class eeffsigmatsubElectron
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to e^+e^-@f$ off the @f$Z@f$ pole, subtracting the SM t-channel.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmatsubElectron : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaElectron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmatsubElectron(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> e^+e^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffRelectrontsub
 * @ingroup EW
 * @brief A class for the ratio of cross sections of @f$R_e=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to e^+e^-)@f$ off the @f$Z@f$ pole, subtracting the SM t channel.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffRelectrontsub : public eeffThObservable {
public:

    /**
     * @brief eeffRelectron constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffRelectrontsub(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the ratio @f$R_e=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to e^+e^-)@f$ at sqrt_s
     */
    double computeThValue();

private:
    
};


#endif	/* LEP2SIGMAELECTRON_H */

