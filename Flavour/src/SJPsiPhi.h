/* 
 * File:   SJPsiPhi.h
 * Author: girardi_mac
 *
 * Created on 26 settembre 2012, 0.06
 */

#ifndef SJPSIPHI_H
#define	SJPSIPHI_H


#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class SJPsiPhi : public ThObservable, AmpDB2 {
public:
    SJPsiPhi(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return sin(AmpBs(NLO).arg());
    };
};


#endif	/* SJPSIPHI_H */

