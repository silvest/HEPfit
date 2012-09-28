/* 
 * File:   Msq.h
 * Author: silvest
 *
 * Created on September 28, 2012, 1:54 PM
 */

#ifndef MSQ_H
#define	MSQ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msq : public ThObservable {
public:
    Msq(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        double res = (static_cast<const SUSY*> (&SM))->getMsd2()(0);
        return ((static_cast<const SUSY*> (&SM))->getMsd2()(0));
    };
private:

};

#endif	/* MSQ_H */

