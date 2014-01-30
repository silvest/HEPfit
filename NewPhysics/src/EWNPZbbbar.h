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
 * @brief A class for new physics in the form of corrections to the \f$Zb\bar{b}\f$ vertices
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the functions to compute the fermionic neutral-current
 * couplings in presence of new physics contributing only to the bottom quark sector,
 * i.e. to \f$g_V^b\f$ and \f$g_A^b\f$ (or, equivalently, \f$\rho_Z^b\f$ and \f$\kappa_Z^b\f$).
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
     * @param[in] l name of a lepton
     * @return the (SM) effective neutral-current coupling @f$\rho_Z^l@f$
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$.
     * @param[in] q name of a quark
     * @return the effective neutral-current coupling @f$\rho_Z^q@f$ (Non-SM 
     * contributions only for \f$q=b\f$)
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$.
     * @param[in] l name of a lepton
     * @return the (SM) effective neutral-current coupling @f$\kappa_Z^l@f$
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$.
     * @param[in] q name of a quark
     * @return the effective neutral-current coupling @f$\kappa_Z^q@f$ (Non-SM 
     * contributions only for \f$q=b\f$)
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @brief The leptonic neutral-current vector coupling @f$g_V^l@f$.
     * @param[in] l name of a lepton
     * @return the (SM) neutral-current vector coupling @f$g_V^l@f$
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief The quark neutral-current vector coupling @f$g_V^q@f$.
     * @param[in] q name of a quark
     * @return the neutral-current vector coupling @f$g_V^q@f$ (Non-SM 
     * contributions only for \f$q=b\f$)
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @brief The leptonic neutral-current axial-vector coupling @f$g_A^l@f$.
     * @param[in] l name of a lepton
     * @return the (SM) neutral-current axial-vector coupling @f$g_A^l@f$
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief The quark neutral-current axial-vector coupling @f$g_A^q@f$.
     * @param[in] q name of a quark
     * @return the neutral-current axial-vector coupling @f$g_A^q@f$ (Non-SM 
     * contributions only for \f$q=b\f$)
     */
    virtual complex gAq(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
private:
    const StandardModel& SM;

};

#endif	/* EWNPZBBBAR_H */

