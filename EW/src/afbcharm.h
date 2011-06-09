/* 
 * File:   AFBcharm.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:42 PM
 */

#ifndef AFBCHARM_H
#define	AFBCHARM_H

#include <ThObservable.h>
#include "EW.h"

class AFBcharm : public ThObservable {
public:

    /**
     * @brief AFBcharm constructor
     * @param[in] myEW an object of EW class
     */
    AFBcharm(const EW& myEW);

    /**
     * @return the forward-backward asymmetry of the c-cbar channel
     */
    double getThValue();

private:

};

#endif	/* AFBCHARM_H */

