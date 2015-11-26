/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MUECONVERSION_H
#define	MUECONVERSION_H

#include "gslpp.h"
#include "ThObservable.h"
#include "LeptonFlavour.h"

//class mueconversion : public ThObservable {
//public:
//    /**
//     * constructor
//     * @param LeptonFlavour
//     */
//    mueconversion(const StandardModel& SM_i);
//
//    /**
//     *
//     * @brief 
//     * @return
//     */
//    double computeThValue();
//
//protected:
//
//private:
//    
//};

/**
 * @class mueconversion_Ti
 * @ingroup LeptonFlavour
 * @brief A class for calculating the decay rate of the process \f$ \mu \to e \f$ conversion in Titanium Nuclei.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The mueconversion_Ti class calculates the decay rate of the process \f$ \mu \to e \f$ conversion in Titanium Nuclei in the model.
 */
class mueconversion_Ti : public ThObservable {
public:
    
    /**
     * @brief Calculates the value of the process \f$ \mu \to e \f$ conversion in Titanium Nuclei.
     * @return returns the value of \f$ \mu \to e \f$ conversion rate in Titanium Nuclei
     */
    mueconversion_Ti(const StandardModel& SM_i);
    
    /**
     * @return returns the value of \f$ \mu \to e \f$ conversion rate in Titanium Nuclei
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

#endif	/* MUECONVERSION_H */
