/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWNPZBBBAR_H
#define	EWNPZBBBAR_H

#include <EWSM.h>

/**
 * @class EWNPZbbbar
 * @brief A class for new physics in the form of corrections to the \f$Zb\bar{b}\f$ vertex
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class EWNPZbbbar : public EWSM {
public:

    EWNPZbbbar(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] l name of a lepton
     * @return the effective neutral-current coupling @f$\rho_Z^l@f$
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the effective neutral-current coupling @f$\rho_Z^q@f$
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the effective neutral-current coupling @f$\kappa_Z^l@f$
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the effective neutral-current coupling @f$\kappa_Z^q@f$
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the effective neutral-current vector coupling @f$g_V^l@f$
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the effective neutral-current vector coupling @f$g_V^q@f$
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the effective neutral-current axial-vector coupling @f$g_A^l@f$
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the effective neutral-current axial-vector coupling @f$g_A^q@f$
     */
    virtual complex gAq(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
private:
    const StandardModel& SM;

};

#endif	/* EWNPZBBBAR_H */

