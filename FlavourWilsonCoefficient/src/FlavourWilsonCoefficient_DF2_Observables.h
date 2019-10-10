/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_DF2_OBSERVABLES_H
#define FLAVOURWILSONCOEFFICIENT_DF2_OBSERVABLES_H

#include <ThObservable.h>


#include "FlavourWilsonCoefficient_DF2.h"

/**
 * @class FlavourWilsonCoefficient_DF2_CK
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the absolute value and phase of NP \f$\Delta S=2\f$ Wilson Coefficients
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */

class FlavourWilsonCoefficient_DF2_CK : public ThObservable {
public:
    FlavourWilsonCoefficient_DF2_CK(const StandardModel& SM_i, const int ind, const int absarg)
    : ThObservable(SM_i), index(ind), absa(absarg), myFlavourWilsonCoefficient_DF2(static_cast<const FlavourWilsonCoefficient_DF2*> (&SM_i))
    {
        if (myFlavourWilsonCoefficient_DF2->isModelFWC_DF2() == false)
            throw std::runtime_error("\nERROR: The Delta F=2 Wilson Coefficients are only defined in a FlavourWilsonCoefficient_DF2 model. Please check your observables list.\n");
    };
        
    double computeThValue()
    {
        if (absa == 0)
            return ((myFlavourWilsonCoefficient_DF2->getC_s()(index)).abs());
        else if (absa == 1)
            return ((myFlavourWilsonCoefficient_DF2->getC_s()(index)).arg()/M_PI*180.);
        else 
            throw std::runtime_error("\nERROR: the type of observable for FlavourWilsonCoefficient_DF2_CK must be either 0 (abs) or 1 (arg)\n");
    };


private:
    const int index, absa;
    const FlavourWilsonCoefficient_DF2 * myFlavourWilsonCoefficient_DF2;
};

class FlavourWilsonCoefficient_DF2_CD : public ThObservable {
public:
    FlavourWilsonCoefficient_DF2_CD(const StandardModel& SM_i, const int ind, const int absarg)
    : ThObservable(SM_i), index(ind), absa(absarg), myFlavourWilsonCoefficient_DF2(static_cast<const FlavourWilsonCoefficient_DF2*> (&SM_i))
    {
        if (myFlavourWilsonCoefficient_DF2->isModelFWC_DF2() == false)
            throw std::runtime_error("\nERROR: The Delta F=2 Wilson Coefficients are only defined in a FlavourWilsonCoefficient_DF2 model. Please check your observables list.\n");
    };
        
    double computeThValue()
    {
        if (absa == 0)
            return ((myFlavourWilsonCoefficient_DF2->getC_c()(index)).abs());
        else if (absa == 1)
            return ((myFlavourWilsonCoefficient_DF2->getC_c()(index)).arg()/M_PI*180.);
        else 
            throw std::runtime_error("\nERROR: the type of observable for FlavourWilsonCoefficient_DF2_CD must be either 0 (abs) or 1 (arg)\n");
    };


private:
    const int index, absa;
    const FlavourWilsonCoefficient_DF2 * myFlavourWilsonCoefficient_DF2;
};

class FlavourWilsonCoefficient_DF2_CBd : public ThObservable {
public:
    FlavourWilsonCoefficient_DF2_CBd(const StandardModel& SM_i, const int ind, const int absarg)
    : ThObservable(SM_i), index(ind), absa(absarg), myFlavourWilsonCoefficient_DF2(static_cast<const FlavourWilsonCoefficient_DF2*> (&SM_i))
    {
        if (myFlavourWilsonCoefficient_DF2->isModelFWC_DF2() == false)
            throw std::runtime_error("\nERROR: The Delta F=2 Wilson Coefficients are only defined in a FlavourWilsonCoefficient_DF2 model. Please check your observables list.\n");
    };
        
    double computeThValue()
    {
        if (absa == 0)
            return ((myFlavourWilsonCoefficient_DF2->getC_bd()(index)).abs());
        else if (absa == 1)
            return ((myFlavourWilsonCoefficient_DF2->getC_bd()(index)).arg()/M_PI*180.);
        else 
            throw std::runtime_error("\nERROR: the type of observable for FlavourWilsonCoefficient_DF2_CBd must be either 0 (abs) or 1 (arg)\n");
    };


private:
    const int index, absa;
    const FlavourWilsonCoefficient_DF2 * myFlavourWilsonCoefficient_DF2;
};

class FlavourWilsonCoefficient_DF2_CBs : public ThObservable {
public:
    FlavourWilsonCoefficient_DF2_CBs(const StandardModel& SM_i, const int ind, const int absarg)
    : ThObservable(SM_i), index(ind), absa(absarg), myFlavourWilsonCoefficient_DF2(static_cast<const FlavourWilsonCoefficient_DF2*> (&SM_i))
    {
        if (myFlavourWilsonCoefficient_DF2->isModelFWC_DF2() == false)
            throw std::runtime_error("\nERROR: The Delta F=2 Wilson Coefficients are only defined in a FlavourWilsonCoefficient_DF2 model. Please check your observables list.\n");
    };
        
    double computeThValue()
    {
        if (absa == 0)
            return ((myFlavourWilsonCoefficient_DF2->getC_bs()(index)).abs());
        else if (absa == 1)
            return ((myFlavourWilsonCoefficient_DF2->getC_bs()(index)).arg()/M_PI*180.);
        else 
            throw std::runtime_error("\nERROR: the type of observable for FlavourWilsonCoefficient_DF2_CBs must be either 0 (abs) or 1 (arg)\n");
    };


private:
    const int index, absa;
    const FlavourWilsonCoefficient_DF2 * myFlavourWilsonCoefficient_DF2;
};


#endif /* FLAVOURWILSONCOEFFICIENT_DF2_OBSERVABLES_H */

