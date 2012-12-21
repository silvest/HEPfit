/* 
 * File:   Mch1.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.25
 */

#ifndef MCH1_H
#define	MCH1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mch1 : public ThObservable {
public:
    Mch1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){        
 
        return ((static_cast<const SUSY*> (&SM))->getMch()(0));
    };
private:

};


#endif	/* MCH1_H */

