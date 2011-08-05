/* 
 * File:   SJPsiK.h
 * Author: silvest
 *
 * Created on June 8, 2011, 5:23 PM
 */

#ifndef SJPSIK_H
#define	SJPSIK_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class SJPsiK : public ThObservable, AmpDB2 {
public:
    SJPsiK(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    virtual double getThValue() {
        return sin(-AmpBd(NLO).arg());
    };
};


#endif	/* SJPSIK_H */

