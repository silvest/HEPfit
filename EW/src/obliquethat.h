/* 
 * File:   obliqueThat.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:47 PM
 */

#ifndef OBLIQUETHAT_H
#define	OBLIQUETHAT_H

#include <ThObservable.h>
#include "EW.h"

class obliqueThat : public ThObservable {
public:

    /**
     * @brief obliqueThat constructor
     * @param[in] myEW an object of EW class
     */
    obliqueThat(const EW& myEW);

    /**
     * @return the oblique parameter T-hat
     */
    double getThValue();

private:

};

#endif	/* OBLIQUETHAT_H */

