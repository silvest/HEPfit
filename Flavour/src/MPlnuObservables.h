/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLNUOBSERVABLES_H
#define MPLNUOBSERVABLES_H

#include "QCD.h"
#include "ThObservable.h"

class RD_MPlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    RD_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$R_{D}@f$ in @f$M \to V l \nu@f$.
    * @return @f$R_{D}@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::lepton lep3; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final pseudoscalar meson type. */
    
};

#endif /* MPLNUOBSERVABLES_H */

