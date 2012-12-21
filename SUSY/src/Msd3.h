/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MD3_H
#define	MD3_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd3 : public ThObservable {
public:
    Msd3(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(2)));
    };
private:

};

#endif	/* MD3_H */

