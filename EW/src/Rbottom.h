/* 
 * File:   Rbottom.h
 * Author: mishima
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <ThObservable.h>
#include "EW.h"


class Rbottom : public ThObservable {
public:

    /**
     * @brief Rbottom constructor
     * @param[in] EW_i an object of EW class
     */
    Rbottom(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
};

#endif	/* RBOTTOM_H */

