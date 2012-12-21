/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MD6_H
#define	MD6_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd6 : public ThObservable {
public:
    Msd6(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(5)));
    };
private:

};

#endif	/* MD6_H */

