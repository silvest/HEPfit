/* 
 * File:   mW.h
 * Author: mishima
 *
 * Created on June 9, 2011, 2:53 PM
 */

#ifndef MW_H
#define	MW_H

#include <ThObservable.h>
#include "EW.h"

class mW : public ThObservable {
public:

    /**
     * @brief mW constructor
     * @param[in] myEW an object of EW class
     */
    mW(const EW& myEW);

    /**
     * @return the W-boson mass
     */
    double getThValue();
    
private:

};

#endif	/* MW_H */

