/* 
 * File:   Rbottom.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:44 PM
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <ThObservable.h>
#include "EW.h"

class Rbottom : public ThObservable {
public:

    /**
     * @brief Rbottom constructor
     * @param[in] myEW an object of EW class
     */
    Rbottom(const EW& myEW);

    /**
     * @return the ratio of the b-bbar width to the hadronic width
     */
    double getThValue();

private:

};

#endif	/* RBOTTOM_H */

