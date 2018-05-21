/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MVlnuObservables.h
 * Author: mauro_87
 *
 * Created on 19 maggio 2018, 20.51
 */

#ifndef MVLNUOBSERVABLES_H
#define MVLNUOBSERVABLES_H

#include "QCD.h"
#include "ThObservable.h"

class BR_MVlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    BR_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2);
    
    /**
    * @brief The binned observable @f$<BR>@f$ in @f$M \to V l^+l^-@f$.
    * @return @f$<BR>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
    
};


#endif /* MVLNUOBSERVABLES_H */

