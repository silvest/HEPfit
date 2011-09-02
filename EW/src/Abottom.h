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
    Abottom(const EW& EW_i);

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

private:

};

#endif	/* ABOTTOM_H */

