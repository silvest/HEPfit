/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MU3_H
#define	MU3_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu3 : public ThObservable {
public:
    Msu3(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(2)));
    };
private:

};

#endif	/* MU3_H */

