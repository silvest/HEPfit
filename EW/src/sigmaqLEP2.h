/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIGMAQLEP2_H
#define	SIGMAQLEP2_H

#include <ThObservable.h>
#include "EW.h"

class sigmaqLEP2 : public ThObservable {
public:
    
    /**
     * @brief sigmaqLEP2 constructor
     */
    sigmaqLEP2(const EW& EW_i);

    /**
     * @return the cross section of the quark channel for LEP2 energies
     */
    double getThValue();
    
private:
    
    double Sigmaq_LEP2;

};

#endif	/* SIGMAQLEP2_H */

