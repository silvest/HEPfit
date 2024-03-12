/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LE_NuScatt_H
#define	LE_NuScatt_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */

/**
 * @class gLnuN2 
 * @ingroup EW
 * @brief An observable class for the effective neutrino nucleon LH coupling
 * @copyright GNU General Public License
 * @details This class is used to compute the effective neutrino nucleon LH coupling 
 * @f$g_L^2(\nu N)@f$
 *
 */
class gLnuN2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gLnuN2(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The effective neutrino nucleon LH coupling @f$g_L^2(\nu N)@f$
     * @return @f$g_L^2(\nu N)@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class gRnuN2 
 * @ingroup EW
 * @brief An observable class for the effective neutrino nucleon RH coupling
 * @copyright GNU General Public License
 * @details This class is used to compute the effective neutrino nucleon RH coupling 
 * @f$g_R^2(\nu N)@f$
 *
 */
class gRnuN2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gRnuN2(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The effective neutrino nucleon RH coupling @f$g_R^2(\nu N)@f$
     * @return @f$g_R^2(\nu N)@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class gVnue 
 * @ingroup EW
 * @brief An observable class for the effective (muon) neutrino-electron vector coupling
 * @copyright GNU General Public License
 * @details This class is used to compute the effective (muon) neutrino-electron vector coupling 
 * @f$g_V^{\nu_\mu e}@f$
 *
 */
class gVnue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gVnue(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The effective (muon) neutrino-electron vector coupling @f$g_V^{\nu_\mu e}@f$
     * @return @f$g_V^{\nu_\mu e}@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class gAnue 
 * @ingroup EW
 * @brief An observable class for the effective (muon) neutrino-electron axial-vector coupling
 * @copyright GNU General Public License
 * @details This class is used to compute the effective (muon) neutrino-electron axial-vector coupling 
 * @f$g_A^{\nu_\mu e}@f$
 *
 */
class gAnue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gAnue(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The effective (muon) neutrino-electron axial-vector coupling @f$g_A^{\nu_\mu e}@f$
     * @return @f$g_A^{\nu_\mu e}@f$
     */
    double computeThValue();

    
private:


};



/** 
 * @}
 */

#endif	/* LE_NuScatt_H */
