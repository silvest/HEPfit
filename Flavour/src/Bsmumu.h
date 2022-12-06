/* 
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSMUMU_H
#define	BSMUMU_H

class StandardModel;
class EvolBsmm;
#include <memory>
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class Bsmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    Bsmumu(const StandardModel& SM_i, int obsFlag, QCD::lepton lep_i=QCD::MU);
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    double computeAmumu(orders order);
    double computeSmumu(orders order);
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    void computeAmpSq(orders order, orders_qed order_qed, double mu);
    void computeObs(orders order, orders_qed order_qed);
    
private:
    
    QCD::lepton lep; /**< Final leptons type. */
    double beta;
    double mBs;
    double mW;
    double mlep;
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
    int obs;
    std::unique_ptr<EvolBsmm> evolbsmm;

};

#endif	/* BSMUMU_H */