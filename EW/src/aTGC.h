/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ATGC_H
#define ATGC_H

#include <stdexcept>
#include <ThObservable.h>
#include "NPbase.h"

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */

/**
 * @class deltag1Z
 * @brief An observable class for the anomalous triple gauge coupling
 * @f$\delta g_{1,Z}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the anomalous triple gauge coupling
 * @f$\delta g_{1,Z}@f$.
 *
 */
class deltag1Z : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltag1Z(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltag1Z class.
     */
    virtual ~deltag1Z();

    /**
     * @brief The anomalous triple gauge coupling @f$\delta g_{1,Z}@f$.
     * @return @f$\delta g_{1,Z}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaKgamma
 * @brief An observable class for the anomalous triple gauge coupling
 * @f$\delta \kappa_{\gamma}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the anomalous triple gauge coupling
 * @f$\delta \kappa_{\gamma}@f$.
 *
 */
class deltaKgamma : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaKgamma(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaKgamma class.
     */
    virtual ~deltaKgamma();

    /**
     * @brief The anomalous triple gauge coupling @f$\delta \kappa_{\gamma}@f$.
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class lambdaZ
 * @brief An observable class for the anomalous triple gauge coupling
 * @f$\lambda_{Z}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the anomalous triple gauge coupling
 * @f$\lambda_{Z}@f$.
 *
 */
class lambdaZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    lambdaZ(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the lambdaZ class.
     */
    virtual ~lambdaZ();

    /**
     * @brief The anomalous triple gauge coupling @f$\lambda_{Z}@f$.
     * @return @f$\lambda_{Z}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/** 
 * @}
 */

#endif	/* ATGC_H */

