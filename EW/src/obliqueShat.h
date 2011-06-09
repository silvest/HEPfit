/* 
 * File:   obliqueShat.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:46 PM
 */

#ifndef OBLIQUESHAT_H
#define	OBLIQUESHAT_H

#include <ThObservable.h>
#include "EW.h"

class obliqueShat : public ThObservable {
public:

    /**
     * @brief obliqueShat constructor
     * @param[in] myEW an object of EW class
     */
    obliqueShat(const EW& myEW);

    /**
     * @return the oblique parameter S-hat
     */
    double getThValue();

private:

};

#endif	/* OBLIQUESHAT_H */

