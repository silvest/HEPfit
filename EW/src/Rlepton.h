/* 
 * File:   Rlepton.h
 * Author: mishima
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class Rlepton : public ThObservable {
public:

    /**
     * @brief Rlepton constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    Rlepton(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double getThValue();

    
private:
    const EW& myEW; 
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* RLEPTON_H */

