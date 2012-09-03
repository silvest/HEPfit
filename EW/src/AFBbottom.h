/* 
 * File:   AFBbottom.h
 * Author: mishima
 */

#ifndef AFBBOTTOM_H
#define	AFBBOTTOM_H

#include <ThObservable.h>
#include "EW.h"


class AFBbottom : public ThObservable {
public:

    /**
     * @brief AFBbottom constructor
     * @param[in] EW_i an object of EW class
     */
    AFBbottom(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* AFBBOTTOM_H */

