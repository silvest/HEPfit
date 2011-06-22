/* 
 * File:   Vcb.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef VCB_H
#define	VCB_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vcb : public ThObservable {
public:
    Vcb(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VCB_H */
