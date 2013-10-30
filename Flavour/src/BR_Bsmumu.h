/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BSMUMU_H
#define	BR_BSMUMU_H

#include <ThObservable.h>
#include "Flavour.h"

class BR_Bsmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bsmumu(Flavour& Flavour, int obsFlag);
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    void setAmp(orders order);
    double getAmumu(orders order);
    double getSmumu(orders order);
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    void AmpSqBsmumu(orders order);
    
private:
    Flavour& myFlavour;
    double beta;
    double mBs;
    double mmu;
    double mb;
    double ms;
    double chiral;
    double absP;
    double argP;
    double absS;
    double argS;
    double ampSq;
    double Amumu;
    double Smumu;
    double phiNP;
    double timeInt;
    double ys;
    int obs;

};

#endif	/* BR_BSMUMU_H */