/* 
 * File:   Abottom.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:43 PM
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <ThObservable.h>
#include "EW.h"

class Abottom : public ThObservable {
public:

    /**
     * @brief Abottom constructor
     * @param[in] myEW an object of EW class
     */
    Abottom(const EW& myEW);

    /**
     * @return the left-right asymmetry of the b-bbar channel
     */
    double getThValue();

private:

};

#endif	/* ABOTTOM_H */

