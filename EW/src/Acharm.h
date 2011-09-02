/* 
 * File:   Acharm.h
 * Author: mishima
 */

#ifndef ACHARM_H
#define	ACHARM_H

#include <ThObservable.h>
#include "EW.h"

class Acharm : public ThObservable {
public:

    /**
     * @brief Acharm constructor
     * @param[in] EW_i an object of EW class
     */
    Acharm(const EW& EW_i);

    /**
     * @return the left-right asymmetry of the c-cbar channel
     */
    double getThValue();

private:

};

#endif	/* ACHARM_H */

