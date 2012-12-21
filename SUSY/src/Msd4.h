/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MD4_H
#define	MD4_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd4 : public ThObservable {
public:
    Msd4(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(3)));
    };
private:

};

#endif	/* MD4_H */

