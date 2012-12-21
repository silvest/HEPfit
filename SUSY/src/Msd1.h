/* 
 * File:   Mu1.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 15.50
 */

#ifndef MD1_H
#define	MD1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd1 : public ThObservable {
public:
    Msd1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(0)));
    };
private:

};


#endif	/* MD1_H */

