/* 
 * File:   DmBd0.h
 * Author: marco
 *
 * Created on April 29, 2011, 12:52 PM
 */

#ifndef DMBD0_H
#define	DMBD0_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBd0 : public ThObservable, AmpDB2 {
public:
    DmBd0(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return(2.*AmpBd(LO).abs());
    };
};

#endif	/* DMBD0_H */
