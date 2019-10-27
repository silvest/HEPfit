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
 * @class Rb0THDMW
 * @ingroup THDMW
 * @brief An observable class to calculate the Rb0THDMW observable in the %THDMW.
 */
class  Rb0THDMW: public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Rb0THDMW(const StandardModel& SM_i);

    /**
     * @return Rb0THDMW
     */
    double computeThValue ();
private:
    const THDMW& myTHDMW;
};




#endif /* THDMWEWPO_H */
