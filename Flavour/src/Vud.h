/* 
 * File:   Vud.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef VUD_H
#define	VUD_H

#include <ThObservable.h>
#include <ThObsType.h>

class Vud : public ThObservable {
public:
    Vud(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* VUD_H */
