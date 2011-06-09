/* 
 * File:   AFBlepton.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:41 PM
 */

#ifndef AFBLEPTON_H
#define	AFBLEPTON_H

#include <ThObservable.h>
#include "EW.h"

class AFBlepton : public ThObservable {
public:
 
    /**
     * @brief AFBlepton constructor
     * @param[in] myEW an object of EW class
     */
    AFBlepton(const EW& myEW);

    /**
     * @return the forward-backward asymmetry of a leptonic channel
     */
    double getThValue();

private:

};

#endif	/* AFBLEPTON_H */

