/* 
 * File:   Rbottom.h
 * Author: mishima
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class Rbottom : public ThObservable {
public:

    /**
     * @brief Rbottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     */
    Rbottom(const EW& EW_i, const bool bCHMN_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i) {};

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN;
};

#endif	/* RBOTTOM_H */

