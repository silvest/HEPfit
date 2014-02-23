/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BDMUMU_H
#define	BDMUMU_H

#include <ThObservable.h>
#include "Flavour.h"

class Bdmumu : public ThObservable {
public:
    /**
     * constructor
     * @param Flavour
     */
    Bdmumu(Flavour& Flavour, int obsFlag);
    
    /**
     *
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_d \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    double computeAmumu(orders order);
    double computeSmumu(orders order);
    
    
protected:
    
    /**
     *
     * @param order
     * @param order_ew
     * @return the short distance contribution to the
     * |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    void computeAmpSq(orders order);
    void computeObs(orders order);
    
private:
    Flavour& myFlavour;
    double beta;
    double mBd;
    double mmu;
    double mb;
    double md;
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
    double yd;
    //double tF;
    int obs;
    
};

#endif	/* BDMUMU_H */