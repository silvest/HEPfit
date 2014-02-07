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
 * @brief A class for the fermionic neutral-current couplings with new physics 
 * contributions to the \f$Zb\bar{b}\f$ vertices.
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
 * in presence of new physics contributing only to the bottom quark sector,
 * i.e. to \f$g_V^b\f$ and \f$g_A^b\f$ (or, equivalently, \f$\rho_Z^b\f$ and \f$\kappa_Z^b\f$).
 * The definitions of the above couplings are given in the description of EWSM
 * class.
 *
 * This class is called from NPZbbbar class.
 *
 * @callergraph
 *
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
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\rho_Z^l@f$ (SM contribution only)
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\rho_Z^q@f$ (Non-SM contribution only for \f$q=b\f$)
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\kappa_Z^l@f$ (SM contribution only)
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return@f$\kappa_Z^q@f$ (Non-SM contribution only for \f$q=b\f$)
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_V^l@f$ (SM contribution only)
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current vector coupling @f$g_V^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_V^q@f$ (Non-SM contribution only for \f$q=b\f$)
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_A^l@f$ (SM contribution only)
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current axial-vector coupling @f$g_A^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_A^q@f$ (Non-SM contribution only for \f$q=b\f$)
     */
    virtual complex gAq(const StandardModel::quark q) const;


};

#endif	/* EWNPZBBBAR_H */

