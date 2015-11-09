/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDS1_H
#define	AMPDS1_H

#include <gslpp_complex.h>
#include "Flavour.h"
#include <StandardModel.h>

class AmpDS1 {
public:
    /**
     * 
     * @brief comupte the amplitude for K_L decay in 2 pion
     * @param Flavour
     */
    AmpDS1(const StandardModel& SM_i);

protected:
    /**
     * 
     * @param order
     * @return the amplitude for K_L decay in 2 pion with 0 isospin change
     */
    gslpp::complex AmpDS1pp0(orders order);
    
    /**
     * 
     * @param order
     * @return the amplitude for K_L decay in 2 pion with double isospin change
     */
    gslpp::complex AmpDS1pp2(orders order);
    
private:
    
    const StandardModel& mySM;

};


#endif	/* AMPDS1_H */
