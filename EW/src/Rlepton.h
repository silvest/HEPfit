/* 
 * File:   Rlepton.h
 * Author: mishima
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <ThObservable.h>
#include "EW.h"


class Rlepton : public ThObservable {
public:

    /**
     * @brief Rlepton constructor
     * @param[in] EW_i an object of EW class
     */
    Rlepton(const EW& EW_i) : ThObservable(EW_i), myEW(EW_i) {};

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double getThValue();

    
private:
    const EW& myEW; 
};

#endif	/* RLEPTON_H */

