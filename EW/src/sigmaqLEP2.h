/* 
 * File:   sigmaqLEP2.h
 * Author: giovannigrilli
 *
 * Created on 27 agosto 2012, 9.25
 */

#ifndef SIGMAQLEP2_H
#define	SIGMAQLEP2_H

#include <ThObservable.h>
#include "EW.h"

class sigmaqLEP2 : public ThObservable {
public:
    
    /**
     * @brief sigmaqLEP2 constructor
     */
    sigmaqLEP2(const EW& EW_i);

    /**
     * @return the cross section of the quark channel for LEP2 energies
     */
    double getThValue();
    
private:
    
    double Sigmaq_LEP2;

};

#endif	/* SIGMAQLEP2_H */

