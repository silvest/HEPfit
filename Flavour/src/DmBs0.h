/* 
 * File:   DmBs0.h
 * Author: silvest
 *
 * Created on June 22, 2011, 5:57 PM
 */

#ifndef DMBS0_H
#define	DMBS0_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBs0 : public ThObservable, AmpDB2 {
public:
    DmBs0(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return(2.*AmpBs(LO).abs());
    };
};



#endif	/* DMBS0_H */

