/* 
 * File:   Mch2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.26
 */

#ifndef MCH2_H
#define	MCH2_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mch2 : public ThObservable {
public:
    Mch2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
 
        return ((static_cast<const SUSY*> (&SM))->getMch()(1));
    };
private:

};

#endif	/* MCH2_H */

