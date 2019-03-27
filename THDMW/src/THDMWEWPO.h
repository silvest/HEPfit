/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef THDMWEWPO_H
#define THDMWEWPO_H

#include "ThObservable.h"

class THDMW;

/**
 * @class Rb0
 * @ingroup THDMW
 * @brief 
 */
class  Rb0THDMW: public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Rb0THDMW(const StandardModel& SM_i);

    /**
     * @return Rb0GTHDM
     */
    double computeThValue ();
private:
    const THDMW& myTHDMW;
};




#endif /* THDMWEWPO_H */
