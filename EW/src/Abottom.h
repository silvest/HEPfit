/* 
 * File:   Abottom.h
 * Author: mishima
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <ThObservable.h>
#include "EW.h"


class Abottom : public ThObservable {
public:

    /**
     * @brief Abottom constructor
     * @param[in] EW_i an object of EW class
     */
    Abottom(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* ABOTTOM_H */

