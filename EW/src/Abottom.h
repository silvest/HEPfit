/* 
 * File:   Abottom.h
 * Author: mishima
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"
#include "EW_CHMN.h"


class Abottom : public ThObservable {
public:

    /**
     * @brief Abottom constructor
     * @param[in] EW_i an object of EW class
     * @param[in] bCHMN_i true if using EW_CHMN class 
     * @param[in] bBURGESS_i true if using the formula in hep-ph/9411257 by C.P. Burgess
     */
    Abottom(const EW& EW_i, const bool bCHMN_i=false, const bool bBURGESS_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myEW_CHMN(EW_i.getSM()), bCHMN(bCHMN_i), bBURGESS(bBURGESS_i) {
        if (bCHMN && bBURGESS)
            throw std::runtime_error("bCHMN and bBURGESS cannot be set to true simultaneously in Abottom()");
    };

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

    
private:
    const EW& myEW;
    const EW_CHMN myEW_CHMN;
    const bool bCHMN, bBURGESS;
};

#endif	/* ABOTTOM_H */

