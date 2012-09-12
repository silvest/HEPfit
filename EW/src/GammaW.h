/* 
 * File:   GammaW.h
 * Author: mishima
 */

#ifndef GAMMAW_H
#define	GAMMAW_H

#include <ThObservable.h>
#include "EW.h"


class GammaW : public ThObservable {
public:
    
    /**
     * @brief GammaW constructor
     * @param[in] EW_i an object of EW class
     */
    GammaW(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the total width of the W boson 
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* GAMMAW_H */

