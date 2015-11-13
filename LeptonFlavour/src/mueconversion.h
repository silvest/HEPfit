/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MUECONVERSION_H
#define	MUECONVERSION_H

#include <gslpp.h>
#include <complex>
#include <ThObservable.h>
#include "LeptonFlavour.h"

class mueconversion : public ThObservable {
public:
    /**
     * constructor
     * @param LeptonFlavour
     */
    mueconversion(const StandardModel& SM_i);

    /**
     *
     * @brief 
     * @return
     */
    double computeThValue();

protected:

private:
    
};

class mueconversion_Ti : public mueconversion {
public:
    
    /**
     * @brief Constructor.
     */
    mueconversion_Ti(const StandardModel& SM_i);
    
    /**
     * @return mueconversion_Ti
     */
    double computeThValue();
    
private:
    const StandardModel& mySM;

};

#endif	/* MUECONVERSION_H */
