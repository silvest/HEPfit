/* 
 * File:   Mneu1.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.30
 */

#ifndef MNEU1_H
#define	MNEU1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mneu1 : public ThObservable {
public:
    Mneu1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return ((static_cast<const SUSY*> (&SM))->getMneu()(0));
    };
private:

};

#endif	/* MNEU1_H */

