/* 
 * File:   DmBd0.h
 * Author: marco
 *
 * Created on April 29, 2011, 12:52 PM
 */

#ifndef DMBD0_H
#define	DMBD0_H

#include "DmBd.h"

class DmBd0 : public DmBd {
public:
    DmBd0(const ThObsType& ObsType) : DmBd(ObsType) {};
    DmBd0(const DmBd0& orig) : DmBd(orig.ObsType) {};
    virtual ~DmBd0() {};

    double getThValue();
};

#endif	/* DMBD0_H */
