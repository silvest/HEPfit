/* 
 * File:   Rcharm.h
 * Author: mishima
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <ThObservable.h>
#include "EW.h"

class Rcharm : public ThObservable {
public:

    /**
     * @brief Rcharm constructor
     * @param[in] EW_i an object of EW class
     */
    Rcharm(const EW& EW_i);

    /**
     * @return the ratio of the c-cbar width to the hadronic width
     */
    double getThValue();

private:

};

#endif	/* RCHARM_H */

