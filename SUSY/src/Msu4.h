/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MU4_H
#define	MU4_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu4 : public ThObservable {
public:
    Msu4(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(3)));
    };
private:

};

#endif	/* MU4_H */

