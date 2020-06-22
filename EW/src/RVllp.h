/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RVLLP_H
#define	RVLLP_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class RWmue
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{W\mu/e}=\Gamma(W\to \mu\nu)/\Gamma(W\to e\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$W\to \mu\nu@f$ 
 * width to the @f$W\to e\nu@f$ width:
 * @f[
 * R_{W\mu/e}=\Gamma(W\to \mu\nu)/\Gamma(W\to e\nu)\,.
 * @f]
 *
 */
class RWmue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RWmue(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{W\mu/e}=\Gamma(W\to \mu\nu)/\Gamma(W\to e\nu)@f$.
     * @return @f$R_{W\mu/e}@f$
     */
    double computeThValue();
    
private:

};


/**
 * @class RWtaue
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{W\tau/e}=\Gamma(W\to \tau\nu)/\Gamma(W\to e\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$W\to \tau\nu@f$ 
 * width to the @f$W\to e\nu@f$ width:
 * @f[
 * R_{W\tau/e}=\Gamma(W\to \tau\nu)/\Gamma(W\to e\nu)\,.
 * @f]
 *
 */
class RWtaue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RWtaue(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{W\tau/e}=\Gamma(W\to \tau\nu)/\Gamma(W\to e\nu)@f$.
     * @return @f$R_{W\tau/e}@f$
     */
    double computeThValue();
    
private:

};


/**
 * @class RWtaumu
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{W\tau/\mu}=\Gamma(W\to \tau\nu)/\Gamma(W\to \mu\nu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$W\to \tau\nu@f$ 
 * width to the @f$W\to \mu\nu@f$ width:
 * @f[
 * R_{W\tau/\mu}=\Gamma(W\to \tau\nu)/\Gamma(W\to \mu\nu)\,.
 * @f]
 *
 */
class RWtaumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RWtaumu(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{W\tau/e}=\Gamma(W\to \tau\nu)/\Gamma(W\to \mu\nu)@f$.
     * @return @f$R_{W\tau/\mu}@f$
     */
    double computeThValue();
    
private:

};




/**
 * @class RZmue
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{Z\mu/e}=\Gamma(Z\to \mu\mu)/\Gamma(Z\to ee)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to \mu\mu@f$ 
 * width to the @f$Z\to ee@f$ width:
 * @f[
 * R_{Z\mu/e}=\Gamma(Z\to \mu\mu)/\Gamma(Z\to ee)\,.
 * @f]
 *
 */
class RZmue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RZmue(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{Z\mu/e}=\Gamma(Z\to \mu\mu)/\Gamma(Z\to ee)@f$.
     * @return @f$R_{Z\mu/e}@f$
     */
    double computeThValue();
    
private:

};


/**
 * @class RZtaue
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{Z\tau/e}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to ee)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to \tau\tau@f$ 
 * width to the @f$Z\to ee@f$ width:
 * @f[
 * R_{Z\tau/e}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to ee)\,.
 * @f]
 *
 */
class RZtaue : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RZtaue(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{Z\tau/e}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to ee)@f$.
     * @return @f$R_{Z\tau/e}@f$
     */
    double computeThValue();
    
private:

};


/**
 * @class RZtaumu
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{Z\tau/\mu}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to \mu\mu)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to \tau\tau@f$ 
 * width to the @f$Z\to \mu\mu@f$ width:
 * @f[
 * R_{Z\tau/\mu}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to \mu\mu)\,.
 * @f]
 *
 */
class RZtaumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RZtaumu(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{Z\tau/e}=\Gamma(Z\to \tau\tau)/\Gamma(Z\to \mu\mu)@f$.
     * @return @f$R_{Z\tau/\mu}@f$
     */
    double computeThValue();
    
private:

};


#endif	/* RVLLP_H */

