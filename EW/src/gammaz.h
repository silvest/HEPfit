/* 
 * File:   GammaZ.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:41 PM
 */

#ifndef GAMMAZ_H
#define	GAMMAZ_H

#include <ThObservable.h>
#include "EW.h"

class GammaZ : public ThObservable {
public:

    /**
     * @brief GammaZ constructor
     * @param[in] myEW an object of EW class
     */
    GammaZ(const EW& myEW);

    /**
     * @return the total width of the Z boson
     */
    double getThValue();

private:

};

#endif	/* GAMMAZ_H */

