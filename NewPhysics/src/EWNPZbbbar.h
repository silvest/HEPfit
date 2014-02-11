/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWNPZBBBAR_H
#define	EWNPZBBBAR_H

#include <EWSM.h>

/**
 * @class EWNPZbbbar
 * @brief A class for the fermionic neutral-current couplings,
 * used together with NPZbbbar class.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the functions to compute the fermionic neutral-current
 * couplings:
 *
 * @li @f$\rho_Z^f@f$&nbsp;&nbsp; (with rhoZ_l() and rhoZ_q()),
 * @li @f$\kappa_Z^f@f$&nbsp;&nbsp; (with kappaZ_l() and kappaZ_q()),
 * @li @f$g_V^f@f$&nbsp;&nbsp; (with gVl() and gVq()),
 * @li @f$g_A^f@f$&nbsp;&nbsp; (with gAl() and gAq()),
 *
 * in presence of new physics (NP) contributing only to the bottom quark sector,
 * i.e. to \f$g_V^b\f$ and \f$g_A^b\f$ (or, equivalently, \f$\rho_Z^b\f$ and \f$\kappa_Z^b\f$).
 * The definitions of the above couplings are given in the description of EWSM
 * class.
 *
 * When the model flag @ref NPZbbbarFlags "NotLinearizedNP" defined in the model 
 * class NPZbbbar is set to TRUE (equivalently, NPZbbbar::FlagNotLinearizedNP=false),
 * the effective couplings for @f$f=b@f$ contain both SM and NP contributions.
 * On the other hand, when the flag is FALSE, all the couplings are identical
 * to the SM ones, and the NP contributions are included via NPZbbbar::deltaGVl(),
 * NPZbbbar::deltaGVq(), NPZbbbar::deltaGAl() and NPZbbbar::deltaGAq().
 *
 */
class EWNPZbbbar : public EWSM {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWNPZbbbar(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief EWSM::rhoZ_l()
     * @details SM contribution only. 
     * @copydetails EWSM::rhoZ_l()
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::rhoZ_q()
     * @details Non-SM contribution only for \f$q=b\f$.
     * @copydetails EWSM::rhoZ_q()
     */
    virtual complex rhoZ_q(const QCD::quark q) const;

    /**
     * @brief @copybrief EWSM::kappaZ_l()
     * @details SM contribution only.
     * @copydetails EWSM::kappaZ_l()
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::kappaZ_q()
     * @details Non-SM contribution only for \f$q=b\f$.
     * @copydetails EWSM::kappaZ_q()
     */
    virtual complex kappaZ_q(const QCD::quark q) const;

    /**
     * @brief @copybrief EWSM::gVl()
     * @details SM contribution only.
     * @copydetails EWSM::gVl()
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::gVq()
     * @details Non-SM contribution only for \f$q=b\f$.
     * @copydetails EWSM::gVq()
     */
    virtual complex gVq(const QCD::quark q) const;

    /**
     * @brief @copybrief EWSM::gAl()
     * @details SM contribution only.
     * @copydetails EWSM::gAl()
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief @copybrief EWSM::gAq()
     * @details Non-SM contribution only for \f$q=b\f$.
     * @copydetails EWSM::gAq()
     */
    virtual complex gAq(const QCD::quark q) const;


};

#endif	/* EWNPZBBBAR_H */

