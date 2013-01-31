/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAKAPPAZB_H
#define	DELTAKAPPAZB_H

class deltaKappaZb : public ThObservable {
public:

    deltaKappaZb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( SM.kappaZ_q(SM.BOTTOM).real() 
                 - SM.StandardModel::kappaZ_q(SM.BOTTOM).real() );
    };
    
private:

};

#endif	/* DELTAKAPPAZB_H */

