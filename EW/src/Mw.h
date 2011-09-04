/* 
 * File:   Mw.h
 * Author: mishima
 */

#ifndef MW_H
#define	MW_H

#include <ThObservable.h>
#include "EW.h"


class Mw : public ThObservable {
public:

    /**
     * @brief Mw constructor
     * @param[in] EW_i an object of EW class
     */
    Mw(const EW& EW_i);

    /**
     * @return the W-boson mass
     */
    double getThValue();

    
private:
    double myMw;

};

#endif	/* MW_H */

