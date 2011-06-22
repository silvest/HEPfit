/* 
 * File:   alpha.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include <ThObservable.h>
#include <ThObsType.h>

class Alpha : public ThObservable {
public:
    Alpha(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue();
};

#endif	/* ALPHA_H */
