/* 
 * Copyright (C) 2018 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MLL_H
#define MLL_H

class StandardModel;
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"

class Mll : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    Mll(const StandardModel& SM_i, int obsFlag, QCD::meson meson_i, QCD::lepton lep_i);
    
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
    QCD::meson meson;
    double ys;
    gslpp::complex CKM_factor;
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
    gslpp::complex C_10;
    gslpp::complex C_10p;
    gslpp::complex C_S;
    gslpp::complex C_Sp;
    gslpp::complex C_P;
    gslpp::complex C_Pp;
    bool FixedWCbtos;
    
    gslpp::vector<gslpp::complex> ** allcoeff;
    gslpp::vector<gslpp::complex> ** allcoeffprime;
    gslpp::vector<gslpp::complex> ** allcoeff_noSM;

};

#endif /* MLL_H */

