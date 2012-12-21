/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MU6_H
#define	MU6_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu6 : public ThObservable {
public:
    Msu6(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(5)));
    };
private:

};

#endif	/* MU6_H */

