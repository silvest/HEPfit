/* 
 * File:   DmBd.h
 * Author: marco
 *
 * Created on April 29, 2011, 12:45 PM
 */

#ifndef DMBD_H
#define	DMBD_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBd : public ThObservable, AmpDB2 {
public:
    DmBd(Flavour& ObsType) : ThObservable(ObsType), AmpDB2(ObsType) {};

    virtual double getThValue() {
        return(2.*Amp(NLO).abs());
    };
};

#endif	/* DMBD_H */
