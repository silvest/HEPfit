/* 
 * File:   M2u2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.56
 */

#ifndef MD5_H
#define	MD5_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd5 : public ThObservable {
public:
    Msd5(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(4)));
    };
private:

};

#endif	/* MD5_H */

