/* 
 * File:   alpha_2a.h
 * Author: silvest
 *
 * Created on April 4, 2011, 10:58 AM
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include <ThObservable.h>
#include <ThObsType.h>

class Alpha_2a : public ThObservable {
public:
    Alpha_2a(const ThObsType& ObsType) : ThObservable(ObsType) {};
    Alpha_2a(const Alpha_2a& orig) : ThObservable(orig.ObsType) {};
    virtual ~Alpha_2a() {};

    double getThValue();
};

#endif	/* ALPHA_2A_H */

