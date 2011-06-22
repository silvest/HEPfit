/* 
 * File:   Vus.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:17 PM
 */

#ifndef VUS_H
#define	VUS_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vus : public ThObservable {
public:
    Vus(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VUS_H */
