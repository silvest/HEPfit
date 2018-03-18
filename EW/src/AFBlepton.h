/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBLEPTON_H
#define	AFBLEPTON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class AFBlepton
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to \ell^+ \ell^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for 
 * an @f$e^+ e^-\to Z\to \ell^+ \ell^-@f$ process at the @f$Z@f$ pole in terms
 * of the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,\ell}=\frac 34 \mathcal{A}_e \mathcal{A}_\ell\,,
 * @f]
 * where @f$\ell@f$ denotes a charged lepton, and lepton-flavour universality
 * is assumed.
 *
 * @sa EW_NPZff::AFBlepton() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class AFBlepton : public ThObservable {
public:
 
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBlepton(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to \ell^+ \ell^-@f$
     * at the @f$Z@f$ pole, @f$A^{0,\ell}_{\mathrm{FB}}@f$
     * @return @f$A^{0,\ell}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

/**
 * @class AFBelectron
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to e^+ e^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for 
 * an @f$e^+ e^-\to Z\to e^+ e^-@f$ process at the @f$Z@f$ pole in terms
 * of the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,e}=\frac 34 \mathcal{A}_e \mathcal{A}_e\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class AFBelectron : public ThObservable {
public:
 
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBelectron(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to e^+ e^-@f$
     * at the @f$Z@f$ pole, @f$A^{0,e}_{\mathrm{FB}}@f$
     * @return @f$A^{0,e}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

/**
 * @class AFBmuon
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to \mu^+ \mu^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for 
 * an @f$e^+ e^-\to Z\to \mu^+ \mu^-@f$ process at the @f$Z@f$ pole in terms
 * of the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,\mu}=\frac 34 \mathcal{A}_e \mathcal{A}_\mu\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class AFBmuon : public ThObservable {
public:
 
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBmuon(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to \mu^+ \mu^-@f$
     * at the @f$Z@f$ pole, @f$A^{0,\mu}_{\mathrm{FB}}@f$
     * @return @f$A^{0,\mu}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

/**
 * @class AFBtau
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to \tau^+ \tau^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for 
 * an @f$e^+ e^-\to Z\to \tau^+ \tau^-@f$ process at the @f$Z@f$ pole in terms
 * of the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,\tau}=\frac 34 \mathcal{A}_e \mathcal{A}_\tau\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class AFBtau : public ThObservable {
public:
 
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBtau(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to \tau^+ \tau^-@f$
     * at the @f$Z@f$ pole, @f$A^{0,\tau}_{\mathrm{FB}}@f$
     * @return @f$A^{0,\tau}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* AFBLEPTON_H */

