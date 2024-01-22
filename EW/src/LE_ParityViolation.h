/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LE_PV_H
#define	LE_PV_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */

/**
 * @class QWe 
 * @ingroup EW
 * @brief An observable class for the weak charge of the electron
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the weak charge of the electron 
 * @f$Q_W(e)@f$
 *
 */
class QWe : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    QWe(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The weak charge of the electron @f$Q_W(e)@f$
     * @return @f$Q_W(e)@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class QWp 
 * @ingroup EW
 * @brief An observable class for the weak charge of the proton
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the weak charge of the proton 
 * @f$Q_W(p)@f$
 *
 */
class QWp : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    QWp(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The weak charge of the proton @f$Q_W(p)@f$
     * @return @f$Q_W(p)@f$
     */
    double computeThValue();

    
private:


};


/**
 * @class QWn 
 * @ingroup EW
 * @brief An observable class for the weak charge of the neutron
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the weak charge of the neutron 
 * @f$Q_W(n)@f$
 *
 */
class QWn : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    QWn(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The weak charge of the neutron @f$Q_W(n)@f$
     * @return @f$Q_W(n)@f$
     */
    double computeThValue();

    
private:


};



/**
 * @class QWAPV 
 * @ingroup EW
 * @brief An observable class for the weak charge from atomic parity violation
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the weak charge from atomic parity violation
 * @f$Q_W(Z,N)@f$
 *
 */
class QWAPV : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    QWAPV(const StandardModel& SM_i, const double Z_i, const double N_i) 
    : ThObservable(SM_i), Z_at(Z_i), N_at(N_i)  
    {
    };

    /**
     * @brief The weak charge from atomic parity violation @f$Q_W(Z,N)@f$
     * @return @f$Q_W(Z,N)@f$
     */
    double computeThValue();
    
    const double Z_at, N_at;

    
private:


};






// The following should not be part of this class. Placed here temporarily for testing.

/**
 * @class Muon g-2 
 * @ingroup EW
 * @brief An observable class for the muon g-2
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the muon g-2 
 * @f$(g_\mu-2)/2@f$
 *
 */
class agminus2muon : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    agminus2muon(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The muon g-2 @f$(g_\mu-2)/2@f$
     * @return @f$(g_\mu-2)/2@f$
     */
    double computeThValue();

    
private:


};



/** 
 * @}
 */

#endif	/* LE_PV_H */
