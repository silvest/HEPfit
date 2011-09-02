/* 
 * File:   DmBs.h
 * Author: silvest
 *
 * Created on June 22, 2011, 5:55 PM
 */

#ifndef DMBS_H
#define	DMBS_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBs : public ThObservable, AmpDB2 {
public:
    DmBs(Flavour& ObsType) : ThObservable(ObsType), AmpDB2(ObsType) {};

    double getThValue() {
        return(2.*AmpBs(NLO).abs());
    };
};



#endif	/* DMBS_H */

