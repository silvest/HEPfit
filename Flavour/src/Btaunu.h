/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BTAUNU_H
#define	BTAUNU_H

#include "ThObservable.h"
#include "QCD.h"

class StandardModel;

class Btaunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    Btaunu(const StandardModel& SM_i, QCD::meson meson_i);
    
    /**
     * 
     * @brief arXiv:1206.2634v2
     * @return theoretical value of |\f$ BR(B \rightarrow \tau \nu) \f$|
     */
    double computeThValue();
    
protected:
    
private:
    
    QCD::meson meson;
    
};

#endif	/* BTAUNU_H */