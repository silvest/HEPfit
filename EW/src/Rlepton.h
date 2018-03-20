/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rlepton
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_\ell^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \ell^+ \ell^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * hadronic width to the @f$Z\to \ell^+ \ell^-@f$ width:
 * @f[
 * R_\ell = \frac{\Gamma_h}{\Gamma_\ell}\,,
 * @f]
 * where @f$\ell@f$ denotes a charged lepton, and lepton-flavour universality
 * is assumed.
 *
 * @sa EW_NPZff::Rlepton() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rlepton : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rlepton(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_\ell^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \ell^+ \ell^-)@f$.
     * @return @f$R_\ell^0@f$
     */
    double computeThValue();


private:


};

/**
 * @class Relectron
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_e^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to e^+ e^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * hadronic width to the @f$Z\to e^+ e^-@f$ width:
 * @f[
 * R_e = \frac{\Gamma_h}{\Gamma_e}\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class Relectron : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Relectron(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_e^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to e^+ e^-)@f$.
     * @return @f$R_e^0@f$
     */
    double computeThValue();


private:


};

/**
 * @class Rmuon
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_\mu^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \mu^+ \mu^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * hadronic width to the @f$Z\to \mu^+ \mu^-@f$ width:
 * @f[
 * R_\mu = \frac{\Gamma_h}{\Gamma_\mu}\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class Rmuon : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rmuon(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_\mu^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \mu^+ \mu^-)@f$.
     * @return @f$R_\mu^0@f$
     */
    double computeThValue();


private:


};

/**
 * @class Rtau
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_\tau^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \tau^+ \tau^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * hadronic width to the @f$Z\to \tau^+ \tau^-@f$ width:
 * @f[
 * R_\tau = \frac{\Gamma_h}{\Gamma_\tau}\,.
 * @f]
 * Lepton-flavour universality is not assumed.
 *
 *
 */
class Rtau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rtau(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_\tau^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \tau^+ \tau^-)@f$.
     * @return @f$R_\tau^0@f$
     */
    double computeThValue();


private:


};

#endif	/* RLEPTON_H */

