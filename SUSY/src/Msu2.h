/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MU2_H
#define	MU2_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu2 : public ThObservable {
public:
    Msu2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(1)));
    };
private:

};

#endif	/* MU2_H */

