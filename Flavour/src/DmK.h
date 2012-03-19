/* 
 * File:   DmK.h
 * Author: stefano
 *
 * Created on 10 gennaio 2012, 15.58
 */

#ifndef DMK_H
#define	DMK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

class DmK : public ThObservable, AmpDK2 {
public:
    /**
     * constructor
     * @param ObsType
     */
    DmK(Flavour& ObsType) : ThObservable(ObsType), AmpDK2(ObsType) {};
    
    /**
     * 
     * @return theoretical value of Delta M_K 
     */
    virtual double getThValue();
};

#endif	/* DMK_H */


