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

class Gammaw_MPlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_1 final leptons of the decay
     * @param[in] lep_2 final leptons of the decay
     */
    Gammaw_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2);
    
    /**
    * @brief The binned observable @f$<\Gamma><w>@f$ in @f$M \to P l \nu@f$.
    * @return @f$<\Gamma><w>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::lepton lep2; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final pseudoscalar meson type. */
    
};

class RD_MPlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_1 final leptons of the decay
     * @param[in] lep_2 final leptons of the decay
     * @param[in] lep_3 final leptons of the decay
     */
    RD_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3);
    
    /**
    * @brief The binned observable @f$R_{D}@f$ in @f$M \to P l \nu@f$.
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

class UnitarityV_MPlnu : public ThObservable{
public:
     
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    * @param[in] meson_i initial meson of the decay
    * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
    * @param[in] lep_i final leptons of the decay
    */
    UnitarityV_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
     
    /**
    * @brief Vectorial unitarity constraint for @f$M \to P l \nu@f$.
    * @return @f$ \Sum_i ag_i^2 @f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */
    
};
 
class UnitarityA_MPlnu : public ThObservable{
public:
     
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    * @param[in] meson_i initial meson of the decay
    * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
    * @param[in] lep_i final leptons of the decay
    */
    UnitarityA_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
    /**
    * @brief Axial unitarity constraint for @f$M \to P l \nu@f$.
    * @return @f$ \Sum_i (af_i^2 + aF1_i^2) @f$
    */
    double computeThValue ();
     
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */

};

class Unitarity_Strong_MPlnu : public ThObservable{
public:
     
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    * @param[in] meson_i initial meson of the decay
    * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
    * @param[in] lep_i final leptons of the decay
    */
    Unitarity_Strong_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
    /**
    * @brief Strong unitarity constraint for @f$M \to P l \nu@f$.
    * @return @f$ \Sum_i (af_i^2 + aF1_i^2) @f$
    */
    double computeThValue ();
     
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */

};

class FFplus_MPlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    FFplus_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$f_{+}@f$ in @f$M \to P l \nu@f$.
    * @return @f$<\Gamma><w>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */
    
};

class FF0_MPlnu : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    FF0_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$f_{0}@f$ in @f$M \to P l \nu@f$.
    * @return @f$<\Gamma><w>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */
    
};

class af0_0 : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    af0_0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
    /**
    * @brief The BGL parameter @f$ a_0^{f_0} @f$ in @f$M \to P l \nu@f$.
    * @return @f$ a_0^{f_0} @f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep1; /**< Final leptons type. */
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson pseudoscalarM; /**< Final vector meson type. */
    
};

#endif /* MPLNUOBSERVABLES_H */

