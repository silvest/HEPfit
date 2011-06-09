/* 
 * File:   Rlepton.h
 * Author: mishima
 *
 * Created on June 9, 2011, 3:44 PM
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <ThObservable.h>
#include "EW.h"

class Rlepton : public ThObservable {
public:

    /**
     * @brief Rlepton constructor
     * @param[in] myEW an object of EW class
     */
    Rlepton(const EW& myEW);

    /**
     * @return the ratio of the hadronic width to the leptonic width
     */
    double getThValue();

private:

};

#endif	/* RLEPTON_H */

