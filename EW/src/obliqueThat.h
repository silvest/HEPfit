/* 
 * File:   obliqueThat.h
 * Author: mishima
 */

#ifndef OBLIQUETHAT_H
#define	OBLIQUETHAT_H

#include <ThObservable.h>
#include "EW.h"

class obliqueThat : public ThObservable {
public:

    /**
     * @brief obliqueThat constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueThat(const EW& EW_i);

    /**
     * @return the oblique parameter T-hat
     */
    double getThValue();

private:

};

#endif	/* OBLIQUETHAT_H */

