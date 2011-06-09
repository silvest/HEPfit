/* 
 * File:   GammaW.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:41 PM
 */

#ifndef GAMMAW_H
#define	GAMMAW_H

#include <ThObservable.h>
#include "EW.h"

class GammaW : public ThObservable {
public:
    
    /**
     * @brief GammaW constructor
     * @param[in] myEW an object of EW class
     */
    GammaW(const EW& myEW);

    /**
     * @return the total width of the W boson 
     */
    double getThValue();

private:

};

#endif	/* GAMMAW_H */

