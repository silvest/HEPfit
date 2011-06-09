/* 
 * File:   Acharm.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:43 PM
 */

#ifndef ACHARM_H
#define	ACHARM_H

#include <ThObservable.h>
#include "EW.h"

class Acharm : public ThObservable {
public:

    /**
     * @brief Acharm constructor
     * @param[in] myEW an object of EW class
     */
    Acharm(const EW& myEW);

    /**
     * @return the left-right asymmetry of the c-cbar channel
     */
    double getThValue();

private:

};

#endif	/* ACHARM_H */

