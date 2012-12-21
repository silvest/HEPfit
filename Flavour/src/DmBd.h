/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBD_H
#define	DMBD_H

#include <ThObservable.h> 
#include <Flavour.h>
#include <AmpDB2.h> 

class DmBd : public ThObservable, AmpDB2 {
public:
    DmBd(Flavour& ObsType) : ThObservable(ObsType), AmpDB2(ObsType) {};
    
    double getThValue();
//       { 
//        
//
//        return(2.*AmpBd(NLO).abs());
//    };
};

#endif	/* DMBD_H */
