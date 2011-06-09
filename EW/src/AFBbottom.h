/* 
 * File:   AFBbottom.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:42 PM
 */

#ifndef AFBBOTTOM_H
#define	AFBBOTTOM_H

#include <ThObservable.h>
#include "EW.h"

class AFBbottom : public ThObservable {
public:

    /**
     * @brief AFBbottom constructor
     * @param[in] myEW an object of EW class
     */
    AFBbottom(const EW& myEW);

    /**
     * @return the forward-backward asymmetry of the b-bar channel
     */
    double getThValue();

private:

};

#endif	/* AFBBOTTOM_H */

