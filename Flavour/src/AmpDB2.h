/* 
 * File:   AmpDB2.h
 * Author: marco
 *
 * Created on June 14, 2011, 2:40 PM
 */

#ifndef AMPDB2_H
#define	AMPDB2_H

#include <gslpp_complex.h>
#include <Flavour.h>

using namespace gslpp;

class AmpDB2 {
public:
    AmpDB2(Flavour& Flavour);

protected:
    complex Amp(orders order);

private:
    Flavour& myFlavour;

};

#endif	/* AMPDB2_H */

