/* 
 * File:   obliqueShat.h
 * Author: mishima
 */

#ifndef OBLIQUESHAT_H
#define	OBLIQUESHAT_H

#include <ThObservable.h>
#include "EW.h"


class obliqueShat : public ThObservable {
public:

    /**
     * @brief obliqueShat constructor
     * @param[in] EW_i an object of EW class
     */
    obliqueShat(const EW& EW_i);

    /**
     * @return the oblique parameter S-hat
     */
    double getThValue();


private:

};

#endif	/* OBLIQUESHAT_H */

