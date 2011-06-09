/* 
 * File:   Rcharm.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:44 PM
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <ThObservable.h>
#include "EW.h"

class Rcharm : public ThObservable {
public:

    /**
     * @brief Rcharm constructor
     * @param[in] myEW an object of EW class
     */
    Rcharm(const EW& myEW);

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

private:

};

#endif	/* RCHARM_H */

