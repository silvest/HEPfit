/* 
 * File:   DmBd.h
 * Author: marco
 *
 * Created on April 29, 2011, 12:45 PM
 */

#ifndef DMBD_H
#define	DMBD_H

#include <ThObservable.h>
#include <ThObsType.h>

class DmBd : public ThObservable {
public:
    DmBd(const ThObsType& ObsType) : ThObservable(ObsType) {};
    DmBd(const DmBd& orig) : ThObservable(orig.ObsType) {};
    virtual ~DmBd() {};

    virtual double getThValue();
protected:
    double getDmBd(int I);
};

#endif	/* DMBD_H */
