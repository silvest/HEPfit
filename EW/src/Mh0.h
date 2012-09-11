/* 
 * File:   Mh0.h
 * Author: girardi_mac
 *
 * Created on 16 agosto 2012, 14.24
 */

#ifndef MH0_H
#define	MH0_H

#include <ThObservable.h>
#include "EW.h"

class Mh0 : public ThObservable {
public: 
    Mh0(const EW& EW_i) : ThObservable(EW_i){ };
    double getThValue(){
        return SM.getMHl();
    }
private:

};

#endif	/* MH0_H */

