/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBS_H
#define	DMBS_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>


class DmBs : public ThObservable, AmpDB2 {
public:
    DmBs(Flavour& ObsType) : ThObservable(ObsType), AmpDB2(ObsType) {};
   
    double getThValue()
    {
        
        //std::cout << "Delta MB_s = " << 2.*AmpBs(NLO).abs() << std::endl;
        
        
        return(2.*AmpBs(NLO).abs());
    };
};



#endif	/* DMBS_H */

