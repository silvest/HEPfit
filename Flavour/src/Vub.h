/* 
 * File:   Vub.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:41 PM
 */

#ifndef VUB_H
#define	VUB_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vub : public ThObservable {
public:
    Vub(const ThObsType& ObsType) : ThObservable(ObsType) {};
    Vub(const Vub& orig) : ThObservable(orig.ObsType) {};
    virtual ~Vub() {};

    double getThValue();
};

#endif	/* VUB_H */
