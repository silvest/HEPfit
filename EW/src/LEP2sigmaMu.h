/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaMu
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \mu^+\mu^-@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaMu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaMu(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i) 
    {
    }

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class LEP2dsigmadcosMu
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \mu^+\mu^-@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2dsigmadcosMu : public LEP2ThDiffObservableBin {
public:

    /**
     * @brief LEP2dsigmadcosMu constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] cos_i the polar angle of the final state particle wrt e^-
     * @param[in] cosmin_i the minimum polar angle of the final state particle wrt e^- in the bin
     * @param[in] cosmax_i the maximum polar angle of the final state particle wrt e^- in the bin
     */
    LEP2dsigmadcosMu(const StandardModel& SM_i, const double sqrt_s_i, const double cos_i, const double cosmin_i, const double cosmax_i) 
    : LEP2ThDiffObservableBin(SM_i, sqrt_s_i, cos_i, cosmin_i, cosmax_i) 
    {
    }

    /**
     * @return the differential cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};

/**
 * @class eeffsigmaMu
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \mu^+\mu^-@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmaMu : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaMu constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmaMu(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffRmuon
 * @ingroup EW
 * @brief A class for the ratio of cross sections of @f$R_mu=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to \mu^+\mu^-)@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffRmuon : public eeffThObservable {
public:

    /**
     * @brief eeffRmuon constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffRmuon(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the ratio @f$R_mu=\sigma(e^+e^-\to hadrons)/\sigma(e^+e^-\to \mu^+\mu^-)@f$ at sqrt_s
     */
    double computeThValue();

private:
    
};



#endif	/* LEP2SIGMAMU_H */

